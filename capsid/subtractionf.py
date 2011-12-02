#!/usr/bin/env python


# Copyright 2011(c) The Ontario Institute for Cancer Reserach. All rights reserved.
#
# This program and the accompanying materials are made available under the
# terms of the GNU Public License v3.0.
#
# You should have received a copy of the GNU General Public License along with
# this program.  If not, see <http://www.gnu.org/licenses/>.


from __future__ import division
from itertools import count
from collections import namedtuple
import re

import pysam
from bx.intervals.intersection import Intersecter, Interval

import capsid


Counter = namedtuple('Counter', ['xeno_mapped', 'ref_mapped', 'ref_unmapped'])
GenomeIds = namedtuple('GenomeIds', ['gi', 'accession'])
Meta = namedtuple('Meta', ['sample', 'alignment'])

db = None
logger = None
meta = None
counter = Counter(count(), count(), count())
regex = re.compile("gi\|(.+?)($|\|)|ref\|(.+?)(\.|$|\|)")


def get_meta(align_name):
    '''Gather the meta data for the alignment from the database'''

    alignment = db.alignment.find_one({"name": align_name})
    try:
        sample = db.sample.find_one({"name": alignment['sample']})
    except TypeError:
        exit(logger.error("Alignment {0} not found in Database.".format(align_name)))
    except KeyError:
        exit(logger.error("Sample cannot be found. Please make sure the Alignment has been created properly."))

    logger.debug('Subtraction for alignment: {0}'.format(align_name))

    return Meta(sample, alignment)


def maps_gene(mapped):
    '''Determine if the mapped alignment falls within a gene.'''
    pass


def build_mapped(align, genome, reference=False):
    '''Generates dict for mapped alignments'''

    scores = [ord(m)-33 for m in align.qqual]

    if align.is_proper_pair:
        align_length = int(align.isize)
        ref_end = int(align.pos + align.isize) + 1
    else:
        # If the cigar data is missing [(), ()] give it a length of 1
        align_length = align.alen if align.alen else 1
        ref_end = align.aend or int(align.pos + 1) + align_length

    try: mismatch = len(re.findall("\D", align.opt('MD')))
    except KeyError: mismatch = align.alen or 0

    mapped = {
       "readId": align.qname
       , "refStrand": -1 if align.is_reverse else 1
       , "refStart": align.pos + 1 # pysam is 0-based index
       , "refEnd": int(ref_end)
       , "alignLength": int(align_length)
       , "readLength": int(align.rlen)  # Total Length of the read
       , "mapq": int(align.mapq)
       , "minQual": min(scores)
       , "avgQual": sum(scores) / len(scores)
       , "miscalls": align.qqual.count('.')
       , "mismatch": mismatch
       , "pairEnd": 1 if align.is_proper_pair else 0
       , "genome": genome
       , "project": meta.sample['project']
       , "sample": meta.sample['name']
       , "alignment": meta.alignment['name']
       , "platform": meta.alignment['platform']
       , "sequencingType": meta.alignment['type']
       }

    if maps_gene(mapped):
        mapped['mapsGene'] = 1

    if reference:
        mapped['isRef'] = 1

    return mapped


def query_genome(gid):
    '''Query database for genome using both gi and accession'''

    if gid.gi:
        genome = db.genome.find_one({'gi':int(gid.gi)}, {'_id':0, 'gi':1})
    elif gid.accession:
        genome = db.genome.find_one({'acession':gid.accession}, {'_id':0, 'gi':1})
    else:
        logger.error('Could not find Genome using gi or accession from genome header. Please make sure your header contains the standard formats for either GenInfo integrated database (gi|integer) or RefSeq (ref|accession|name).')

    return genome


def determine_genome(genome_header):
    '''Use the genome header data in BAM file to determine the genome'''

    result = regex.findall(genome_header)

    gids = [GenomeIds(r[0], r[2]) for r in result]

    return filter(None, [query_genome(gid) for gid in gids]).pop()['gi']


def insert_mapped(mapped):
    '''Insert mapped alignments into Database'''

    print 'insert_mapped: ', mapped


def extract_mapped(align, genome_header):
    '''Process alignment to determine if valid'''

    genome = determine_genome(genome_header)

    return build_mapped(align, genome)


def parse_ref(f):
    '''Extract alignments from Reference BAM file'''
    bamfile = pysam.Samfile(f, 'rb')
    pass


def parse_xeno(f):
    '''Extract alignments from Xeno BAM file'''

    bamfile = pysam.Samfile(f, 'rb')
    mapped = (extract_mapped(align, bamfile.getrname(align.tid)) for align in bamfile.fetch() if not align.is_proper_pair or align.is_proper_pair and align.isize > 0)
    ids = set(insert_mapped(align) for align in mapped)
    print 'done'

    return ids


def main(args):
    ''' '''
    global db, logger, meta

    logger = args.logging.getLogger(__name__)
    db = capsid.connect(args)
    meta = get_meta(args.align)

    xeno_mapped_readids = parse_xeno(args.xeno)
    #intersecting_mapped_readids = parse_ref(args.ref)


if __name__ == '__main__':
    print 'This program should be run as part of the capsid package:\n\t$ capsid subtraction -h\n\tor\n\t$ /path/to/capsid/bin/capsid subtraction -h'
