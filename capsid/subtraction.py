#!/usr/bin/env python


# Copyright 2011(c) The Ontario Institute for Cancer Reserach. All rights reserved.
#
# This program and the accompanying materials are made available under the
# terms of the GNU Public License v3.0.
#
# You should have received a copy of the GNU General Public License along with
# this program.  If not, see <http://www.gnu.org/licenses/>.


from __future__ import division
from itertools import count, ifilter, imap
from collections import namedtuple
from functools import partial
import re

import pysam
from bx.intervals.intersection import Intersecter, Interval

from database import *

Counter = namedtuple('Counter', ['xeno_mapped', 'ref_mapped', 'ref_unmapped', 'maps_gene'])
Meta = namedtuple('Meta', ['sample', 'alignment'])
LookupGroup = namedtuple('LookupGroup', ['xeno', 'ref'])
Lookup = namedtuple('Lookup', ['file', 'column'])

db = None
logger = None
meta = None
mapq = None
intersecters = {}
genomes = {}
counter = Counter(count(), count(), count(), count())
regex = re.compile("gi\|(.+?)($|\|)|ref\|(.+?)(\.|$|\|)")


def get_meta(align_name):
    '''Gather the meta data for the alignment from the database'''

    alignment = db.Alignment.one({"name": align_name})
    try:
        sample = db.Sample.one({"name": alignment.sample})
    except AttributeError:
        exit(logger.error("Alignment {0} not found in Database.".format(align_name)))

    logger.debug('Subtraction for alignment: {0}'.format(align_name))

    return Meta(sample, alignment)


def update_isref(readId):
    ''' '''

    db.mapped.update({'readId': readId, 'alignment':meta.alignment.name},
                     {'$set': {'isRef': 1}}, False, False, False, True)


def extract_unmapped(align, fastq):
    '''Output unmapped alignments to Fastq file'''

    fastq.write('@'+str(align.qname)+'\n'+str(align.seq)+'\n+\n'+str(align.qqual)+'\n')


def maps_gene(mapped):
    '''Determine if the mapped alignment falls within a gene.'''
    global intersecters

    try:
        intersecter = intersecters[mapped['genome']]
    except KeyError:
        genes = db.feature.find({'genome': mapped['genome'], 'type': 'gene'})

        intersecter = Intersecter()

        # Interval end is exclusive, need to +1 to line up with actual position
        [intersecter.add_interval(Interval(gene['start'], gene['end'] + 1, gene['name']))
         for gene in genes]

        intersecters[mapped['genome']] = intersecter

    return intersecter.find(mapped['refStart'], mapped['refEnd'])


def build_mapped(align, genome, reference):
    '''Generates dict for mapped alignments'''

    scores = [ord(m)-33 for m in align.qqual]

    if align.is_proper_pair:
        align_length = align.isize
        ref_end = align.pos + align.isize + 1
    else:
        # If the cigar data is missing [(), ()] give it a length of 1
        align_length = align.alen
        ref_end = align.aend or align.pos + align_length + 1

    try:
        MD = align.opt('MD')
        mismatch = len(re.findall("\D", MD))
    except KeyError:
        MD = ''
        try: mismatch = int(align.rlen)
        except TypeError:
            logger.debug(align)
            logger.error('Aligned read with null length')
            exit()

    mapped = {
       "readId": align.qname
       , "refStrand": -1 if align.is_reverse else 1
       , "refStart": int(align.pos) + 1 # pysam is 0-based index
       , "refEnd": int(ref_end)
       , "alignLength": int(align_length)
       , "readLength": int(align.rlen)  # Total Length of the read
       , "mapq": int(align.mapq)
       , "minQual": min(scores)
       , "avgQual": sum(scores) / len(scores)
       , "miscalls": align.qqual.count('.')
       , "mismatch": mismatch
       , "pairEnd": 1 if align.is_proper_pair else 0
       , "genome": int(genome)
       , "project": meta.sample.project
       , "sample": meta.sample.name
       , "alignment": meta.alignment.name
       , "platform": meta.alignment.platform
       , "sequencingType": meta.alignment.type
       , "sequence": align.query
       , "cigar": align.cigar
       , "MD": MD
       }

    mapped_genes = maps_gene(mapped)
    if mapped_genes:
        counter.maps_gene.next()
        mapped['mapsGene'] = [gene.value for gene in mapped_genes]

    if reference:
        mapped['isRef'] = 1

    return mapped


def get_genome(gid):
    '''Query database for genome using both gi and accession'''

    if gid['gi']:
        genome = int(gid['gi'])
    elif gid['accession']:
        try: genome = db.genome.find_one({'accession':gid['accession']},
                                         {'_id':0, 'gi':1})['gi']
        except TypeError: genome = None

    return genome


def build_lookup(column, header, gid):
    ''' '''
    global genomes
    gids = {'gi': None, 'accession': None}
    gids[column] = gid

    try:
        genome = get_genome(gids)
        if genome:
            genomes[str(header)] = genome
    except ValueError: pass


def genome_lookup(lookup):
    ''' '''
    global genomes
    logger.info('Building Lookup Table...')
    logger.debug('Lookup file {} using {} as bridge'.format(lookup.file, lookup.column))

    genomes = {} # Remove any old values before re-populating
    with open(lookup.file, 'rU') as fh:
        [build_lookup(lookup.column, *line.strip().split(';')) for line in fh]


def determine_genome(genome_header):
    '''Use the genome header data in BAM file to determine the genome'''

    result = regex.findall(genome_header)

    gids = [{'gi':r[0], 'accession':r[2]} for r in result]
    try:
        genome = filter(None, [get_genome(gid) for gid in gids]).pop()
    except IndexError:
        genome = None

    return genome


def insert_mapped(mapped, process):
    '''Insert mapped alignments into Database'''

    # If it maps to a genome save in the database, otherwise
    # just return the intersecting read id
    if mapped['genome'] and process in ['both', 'mapped']:
        db.mapped.insert(mapped)

    return mapped['readId']


def valid_mapped(align):
    '''Returns true if single-end or pair-end with isize'''

    return not align.is_proper_pair or align.is_proper_pair and align.isize


def extract_mapped(align, bamfile, reference=False):
    '''Process mapped alignment and return dict'''
    global genomes

    if (0 <= align.mapq <= 3 or align.mapq >= mapq) and valid_mapped(align):
        try:
            genome = genomes[str(bamfile.getrname(align.tid))]
        except KeyError:
            genome = determine_genome(bamfile.getrname(align.tid))
            if genome: genomes[str(bamfile.getrname(align.tid))] = genome

        # If it doesn't map to a genome in the database, assume it is mapping to
        # a junction and just save as an intersecting readId
        if genome:
            counter.ref_mapped.next() if reference else counter.xeno_mapped.next()
            mapped = build_mapped(align, genome, reference)
        else:
            mapped = {'readId': align.qname, 'genome': None}

        return mapped


def maps_xeno(align, readids):
    '''Checks if alignment readId matches any in readids set'''

    return align.qname in readids


def extract_alignment(align, bamfile, readids, fastq):
    '''Process alignments '''

    in_xeno = maps_xeno(align, readids)

    if align.is_unmapped and not in_xeno and fastq:
        counter.ref_unmapped.next()
        extract_unmapped(align, fastq)
    elif not align.is_unmapped and in_xeno:
        return extract_mapped(align, bamfile, True)


def parse_ref(args, readids, process):
    '''Extract alignments from Reference BAM file'''
    global genomes

    if args.lookup.ref.file: genome_lookup(args.lookup.ref)

    logger.info('Finding alignments in Reference BAM file...')
    logger.debug('Reference BAM File: {0}'.format(args.ref))

    bamfile = pysam.Samfile(args.ref, 'rb')

    fastq = open(meta.alignment.name + '.unmapped.fastq', 'w') if process in ['both', 'unmapped'] else None

    return ifilter(None, (extract_alignment(align, bamfile, readids, fastq)
                          for align in bamfile.fetch(until_eof=True)))


def parse_xeno(args):
    '''Extract alignments from Xeno BAM file'''
    global genomes

    if args.lookup.xeno.file: genome_lookup(args.lookup.xeno)

    logger.info('Finding alignments in Xeno BAM file...')
    logger.debug('Xeno BAM File: {0}'.format(args.xeno))

    bamfile = pysam.Samfile(args.xeno, 'rb')
    return ifilter(None, (extract_mapped(align, bamfile)
                          for align in bamfile.fetch()))


def summary(xeno_mapped_readids, intersecting_mapped_readids, process):
    '''Logging summary of added records'''

    # Counter starts at 0, so >it = count(); >print it.next(); >0
    # Using counter.*.next here sets it to the correct value for printing.
    xeno_mapped = counter.xeno_mapped.next()
    ref_mapped = counter.ref_mapped.next()
    ref_unmapped = counter.ref_unmapped.next()
    maps_gene = counter.maps_gene.next()

    total_mapped = xeno_mapped + ref_mapped if process in ['mapped', 'both'] else 0

    logger.info('Total mapped alignments found in Xeno BAM file: {0}'.format(xeno_mapped))
    logger.info('Total mapped alignments found in Reference BAM file: {0}'.format(ref_mapped))
    logger.info('Total mapped alignments added to database: {0}'.format(total_mapped))
    logger.info('Xeno reads that hit a gene ("mapsGene":1): {0}'.format(maps_gene))
    logger.info('Reads that map to Xeno with unique Read IDs: {0}'.format(len(xeno_mapped_readids)))
    logger.info('Reads that map to both Xeno and Reference with unique Read IDs: {0}'.format(len(intersecting_mapped_readids)))
    logger.info('Unmapped alignments found in Reference BAM file: {0}'.format(ref_unmapped))


def main(args):
    ''' '''
    global db, logger, meta, mapq
    logger = args.logging.getLogger(__name__)

    if args.lookup:
        args.lookup = LookupGroup(Lookup._make(args.lookup),
                                  Lookup._make(args.lookup))
    else:
        args.lookup = LookupGroup(Lookup._make(args.xeno_lookup),
                                  Lookup._make(args.ref_lookup))

    db = connect(args)
    meta = get_meta(args.align)
    mapq = int(args.filter)
    process = args.process

    p_mapped = partial(insert_mapped, process=process)

    xeno_mapped = parse_xeno(args)
    if process in ['both', 'mapped']:
        logger.info('Inserting mapped alignments from Xeno BAM file...')
    xeno_mapped_readids = set(map(p_mapped, xeno_mapped))

    ref_mapped = parse_ref(args, xeno_mapped_readids, process)
    if process == 'mapped':
        logger.info('Inserting mapped alignments from Reference BAM file...')
    elif process == 'unmapped':
        logger.info('Outputting unmapped alignments from Reference BAM file...')
    else:
        logger.info('Inserting mapped and outputting unmapped from Reference BAM file...')
    intersecting_mapped_readids = set(map(p_mapped, ref_mapped))

    logger.info('Updating reads that map in both Xeno and Reference...')
    map(update_isref, intersecting_mapped_readids)

    summary(xeno_mapped_readids, intersecting_mapped_readids, process)


if __name__ == '__main__':
    print 'This program should be run as part of the capsid package:\n\t$ capsid subtraction -h\n\tor\n\t$ /path/to/capsid/bin/capsid subtraction -h'
