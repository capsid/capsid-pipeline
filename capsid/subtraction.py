#!/usr/bin/env python


# Copyright 2011(c) The Ontario Institute for Cancer Research. All rights reserved.
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
import os, sys
import subprocess
import csv

import pysam
from bx.intervals.intersection import Intersecter, Interval

from database import *
import alignment


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
temp = None


def check_align(args):

    finder = {"name": args.align, "sample": args.sample, "projectLabel": args.project}
    aln = db.alignment.find_one(finder)
    if not aln:
        logger.error("Alignment {0} not found in {1}/{2}.".format(args.align, args.project, args.sample))
        sys.exit(1)

    return aln

def get_meta(args):
    """Gather the meta data for the alignment from the database"""

    aln = check_align(args)
    sample = db.sample.find_one({"_id": aln['sampleId']})

    if not sample:
        logger.error("Alignment {0} not found in database. --project and --sample need to be passed to auto-create".format(args.align))
        exit(1)

    logger.debug('Subtraction for alignment: {0}'.format(args.align))

    return Meta(sample, aln)


def update_isref(readId):
    ''' '''

    db.mapped.update({'readId': readId, 'alignmentId':meta.alignment['_id']},
                     {'$set': {'isRef': 1}}, False, False, False, True)


def extract_unmapped(align, fastq):
    """Output unmapped alignments to Fastq file"""
    # here align.qname should be extracted in its original form i.e. with /1 and /2 if they exist 
    fastq.write("@{0:>s}\n{1:>s}\n+\n{2:>s}\n".format(str(align.qname), str(align.query), str(align.qqual)))


# new 
def tuple_to_CIG(readcig):
    #cig_ops=['M','I','D','N','S','H','P','=','X']
    cig_ops=['M','I','D','N','S','H','P','M','X']
    cig_string=str()
    for tple in readcig:
        cig_string = cig_string + str(tple[1]) + str(cig_ops[tple[0]])
    return cig_string    


# new 
def extract_mapped_sam(align, genome_sam, sam_file):
    """Output mapped xeno reads to a sam file"""
    sam_file.write("{0:>s}\t{1:>s}\t{2:>s}\t{3:>s}\t{4:>s}\t{5:>s}\t{6:>s}\t{7:>s}\t{8:>s}\t{9:>s}\t{10:>s}\t{11:>s}\n".format(str(align.qname.split("/")[0]),str(align.flag),genome_sam,str(align.pos + 1),str(align.mapq),str(tuple_to_CIG(align.cigar)),'*','0','0',str(align.query),'*','AS:i:' + str(align.opt('AS'))))

# new 
def extract_hg_readIds(align,hgids):
    """Output mapped hg readId to a file"""
    hgids.write("{0:>s}\n".format(align.qname.split("/")[0]))

def get_readscan_data(readscan_file):
    row_number = 0
    genome_identifier = re.compile("gi:(\d+)")
    with open(readscan_file, 'rU') as f:
        readscan_reader = csv.DictReader(f, delimiter='\t', quotechar='|')
        for row in readscan_reader:
            row['row'] = row_number
            row['index'] = int(row['index'])
            row['score'] = float(row['score'])
            row_number = row_number + 1
            if row.has_key(None):
                del row[None]
            match = genome_identifier.match(row['id'])
            if match:
                row['genome'] = int(match.group(1))
            yield row

# new
def get_only_xeno_reads(pathogen, human, args):
    """Obtain reads in sam format that only map to xeno"""    
    #print args.xeno.rsplit('/',1)[0]
    logger.info("Temporary directory: {0}".format(temp))
    logger.info("Sorting {0}".format(temp + pathogen))

    pathogen_input_file = temp + pathogen
    pathogen_sorted_file = temp + pathogen + '.sorted'
    human_input_file = temp + human
    human_sorted_file = temp + human + '.sorted'

    if subprocess.call(['sort', pathogen_input_file, '-t\t', '-d', '--key=1,1', '-T', temp, '-o', pathogen_sorted_file]) != 0:
        logger.error("Error sorting {0}".format(pathogen_input_file))
        sys.exit(1)
    logger.info("Sorting {0}".format(temp + human))
    if subprocess.call(['sort', human_input_file, '-t\t', '-d', '--key=1,1', '-T', temp, '-o', human_sorted_file]) != 0:
        logger.error("Error sorting {0}".format(human_input_file))
        sys.exit(1)
        
    os.remove(pathogen_input_file)
    os.remove(human_input_file)

    xeno_file = os.path.abspath(args.xeno)
    xeno_directory = os.path.dirname(xeno_file)

    gra = __file__.rsplit('/',1)[0] + '/gra.sh'
    logger.info("Output directory: {0}".format(xeno_directory))
    if subprocess.call([gra, pathogen_sorted_file, human_sorted_file, xeno_directory]) != 0:
        logger.error("Failed to calculate genome relative abundance successfully")
        sys.exit(1)

    readscan_file = xeno_directory + '/pathogen.gra.txt'
    result = list(get_readscan_data(readscan_file))

    selector = {"_id" : meta.alignment['_id']}
    updater = {"$set" : {"gra" : result}}
    db.alignment.update(selector, updater)
    
 
def maps_gene(mapped):
    """Determine if the mapped alignment falls within a gene."""
    global intersecters

    try:
        intersecter = intersecters[mapped['genome']]
    except KeyError:
        genes = db.feature.find({'genome': mapped['genome'], 'type': 'gene'})

        intersecter = Intersecter()

        # Interval end is exclusive, need to +1 to line up with actual position
        [intersecter.add_interval(Interval(gene['start'], gene['end'] + 1, gene['uid']))
         for gene in genes]

        intersecters[mapped['genome']] = intersecter

    return intersecter.find(mapped['refStart'], mapped['refEnd'])


def build_mapped(align, genome, reference):
    '''Generates dict for mapped alignments'''

    # Since align.qqual is in Sanger format. 
    scores = [ord(m)-33 for m in align.qqual]

    #if align.is_proper_pair:
        #align_length = align.isize
        #ref_end = align.pos + align.isize + 1
    #else:
        # If the cigar data is missing [(), ()] give it a length of 1
    #    align_length = align.alen
    #    ref_end = align.aend or align.pos + align_length + 1


    if align.alen is None:
        align_length = align.pos + align.qlen + 1        
    else:
        align_length = align.alen 
        

    if align.aend is None:
        ref_end = align.pos + align.qlen + 1
    else:
        ref_end = align.aend

 
    try:
        MD = align.opt('MD')
        mismatch = len(re.findall("\D", MD))
    except KeyError:
        MD = None
        try: mismatch = int(align.rlen)
        except TypeError:
            logger.debug(align)
            logger.error('Aligned read with null length')
            exit()

    try:
        AS = int(align.opt('AS'))
    except KeyError:
        AS = None

    try:
        PG = align.opt('PG')
    except KeyError:
        PG = None

    mapped = {
       # take only the read identifier exclude /1 and /2 
       "readId": align.qname.split("/")[0]
       , "refStrand": -1 if align.is_reverse else 1
       , "refStart": int(align.pos) + 1 # pysam is 0-based index
       , "refEnd": int(ref_end)
       , "alignLength": int(align_length)
       , "readLength": int(align.rlen)  # Total Length of the read
       , "mapq": int(align.mapq)
       , "minQual": min(scores)
       , "avgQual": sum(scores) / len(scores)
       , "qqual": align.qqual
       , "miscalls": align.qqual.count('.')
       , "mismatch": mismatch
       , "pairEnd": 1 if align.is_proper_pair else 0
       , "genome": int(genome)
       , "projectLabel": meta.sample['projectLabel']
       , "projectId": meta.sample['projectId']
       , "sample": meta.sample['name']
       , "sampleId": meta.sample['_id']
       , "alignment": meta.alignment['name']
       , "alignmentId": meta.alignment['_id']
       , "platform": meta.alignment['platform']
       , "sequencingType": meta.alignment['type']
       , "sequence": align.query
       , "cigar": align.cigar
    }

    if AS: mapped['alignScore'] = AS
    if MD: mapped['MD'] = MD
    if PG: mapped['PG'] = PG

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



def get_genome_gra(gid):
    '''Query database and return genome in proper format for gra'''
    try: 
        genome_db = db.genome.find_one({'gi':gid})
        genome_gra = 'gi' + '|' + str(genome_db['gi']) + '|ref|' + genome_db['accession'] + '.' + str(genome_db['version']) + '|' + genome_db['name'].split(' ')[0]  
    except TypeError: 
        genome_gra = False    
    return genome_gra



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
    # new do not insert into mongo right now
    if mapped['genome'] and process in ['both', 'mapped']:
        db.mapped.insert(mapped)
    
    return mapped['readId']


def valid_mapped(align):
    '''Returns true if valid single-end or pair-end alignment'''
    #return (not align.is_proper_pair or align.is_proper_pair and align.isize) and not align.is_unmapped
    if align.is_paired:
        return (align.is_proper_pair or not align.mate_is_unmapped or not align.is_unmapped)
    else:
        return (not align.is_unmapped)


def extract_mapped(align, bamfile, sam_file, reference=False):
    '''Process mapped alignment and return dict'''
    global genomes

    if (0 <= align.mapq <= 3 or align.mapq >= mapq) and valid_mapped(align):
        if gra and sam_file: 
          try:  
              #genome_sam = str(bamfile.getrname(align.tid))                                                                           
              #extract_mapped_sam(align,genome_sam,sam_file)   
              genome_sam = get_genome_gra(genomes[str(bamfile.getrname(align.tid))])   
          except KeyError:
              genome_sam = str(bamfile.getrname(align.tid))       
          except ValueError:
              pass
          else:
              if genome_sam:
                  extract_mapped_sam(align,genome_sam,sam_file)

        try:
            genome = genomes[str(bamfile.getrname(align.tid))]
        except KeyError:
            genome = determine_genome(bamfile.getrname(align.tid))
            if genome: genomes[str(bamfile.getrname(align.tid))] = genome
        except ValueError:     
            genome = False
            pass

        
        # If it doesn't map to a genome in the database, assume it is mapping to
        # a junction and just save as an intersecting readId
        if genome:
            counter.ref_mapped.next() if reference else counter.xeno_mapped.next()
            mapped = build_mapped(align, genome, reference)
        else:
            mapped = {'readId': align.qname.split("/")[0], 'genome': None}

        return mapped


def maps_xeno(align, readids):
    '''Checks if alignment readId matches any in readids set'''
    return align.qname.split("/")[0] in readids


def extract_alignment(align, bamfile, readids, fastq, hgids):
    '''Process alignments '''

    in_xeno = maps_xeno(align, readids)

    #if align.is_unmapped and not in_xeno and fastq:
    if not valid_mapped(align) and not in_xeno and fastq:
        counter.ref_unmapped.next()
        extract_unmapped(align, fastq)
    #elif not align.is_unmapped and in_xeno:
    elif valid_mapped(align) and in_xeno:
        if gra and hgids:
            # i.e map to both pathogen and human ref
            extract_hg_readIds(align,hgids)
        return extract_mapped(align, bamfile, False, True)


def parse_ref(args, readids, process):
    '''Extract alignments from Reference BAM file'''
    global genomes

    if args.lookup.ref.file: genome_lookup(args.lookup.ref)

    logger.info('Finding alignments in Reference BAM file...')
    logger.debug('Reference BAM File: {0}'.format(args.ref))

    bamfile = pysam.Samfile(args.ref, 'rb')

    fastq = open(meta.alignment['name'] + '.unmapped.fastq', 'w') if process in ['both', 'unmapped'] else None

    if args.gra:
        hgids = open(temp + meta.alignment['name'] + '.hg.mapped.txt', 'a') if process in ['both', 'mapped'] else None
        return ifilter(None, (extract_alignment(align, bamfile, readids, fastq, hgids)
                          for align in bamfile.fetch(until_eof=True)))
    else:
          return ifilter(None, (extract_alignment(align, bamfile, readids, fastq, hgids=False)
                          for align in bamfile.fetch(until_eof=True)))
      

def parse_xeno(args,process):
    '''Extract alignments from Xeno BAM file'''
    global genomes

    if args.lookup.xeno.file: genome_lookup(args.lookup.xeno)

    logger.info('Finding alignments in Xeno BAM file...')
    logger.debug('Xeno BAM File: {0}'.format(args.xeno))

    bamfile = pysam.Samfile(args.xeno, 'rb')

    # open the location for the sam file
    if args.gra:
        sam_file = open(temp + meta.alignment['name'] + '.vg.mapped.sam', 'a') if process in ['both', 'mapped'] else None
        return ifilter(None, (extract_mapped(align, bamfile, sam_file)
                          for align in bamfile.fetch(until_eof=True)))
    else:
        return ifilter(None, (extract_mapped(align, bamfile, sam_file=False)
                          for align in bamfile.fetch(until_eof=True)))
 


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
    logger.info('Total mapped alignments found in Reference (only for Refs stored in the db) and Xeno BAM files: {0}'.format(ref_mapped))
    logger.info('Total mapped alignments added to database: {0}'.format(total_mapped))
    logger.info('Xeno reads that hit a gene ("mapsGene":1): {0}'.format(maps_gene))
    logger.info('Reads that map to Xeno with unique Read IDs: {0}'.format(len(xeno_mapped_readids)))
    # Note: Reads that map to both Xeno and Reference _BUT_ also include those human ref that are not stored in the db with unique Read IDs
    logger.info('Reads that map to both Xeno and Reference with unique Read IDs : {0}'.format(len(intersecting_mapped_readids)))
    if process in ['both', 'unmapped']:
        logger.info('Unmapped alignments found in Reference BAM file: {0}'.format(ref_unmapped))


def main(args):
    ''' '''
    global db, logger, meta, mapq, temp, gra
    logger = args.logging.getLogger(__name__)

    if args.lookup:
        args.lookup = LookupGroup(Lookup._make(args.lookup),
                                  Lookup._make(args.lookup))
    else:
        args.lookup = LookupGroup(Lookup._make(args.xeno_lookup),
                                  Lookup._make(args.ref_lookup))

    db = connect(args)
    meta = get_meta(args)
    mapq = int(args.filter)
    process = args.process
    temp = args.temp
    gra = args.gra

    p_mapped = partial(insert_mapped, process=process)

    xeno_mapped = parse_xeno(args,process)
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

    if args.gra:
        logger.info('Outputting the set of reads mapping to Xeno only for GRA calculation')
        get_only_xeno_reads(meta.alignment['name'] + '.vg.mapped.sam', meta.alignment['name'] + '.hg.mapped.txt', args) if process in ['both', 'mapped'] else None

    logger.info('Updating reads that map in both Xeno and Reference...')
    map(update_isref, intersecting_mapped_readids)

    summary(xeno_mapped_readids, intersecting_mapped_readids, process)


if __name__ == '__main__':
    print 'This program should be run as part of the capsid package:\n\t$ capsid subtraction -h\n\tor\n\t$ /path/to/capsid/bin/capsid subtraction -h'
