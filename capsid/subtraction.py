#!/usr/bin/env python


# Copyright 2011(c) The Ontario Institute for Cancer Reserach. All rights reserved.
#
# This program and the accompanying materials are made available under the
# terms of the GNU Public License v3.0.
#
# You should have received a copy of the GNU General Public License along with
# this program.  If not, see <http://www.gnu.org/licenses/>.


from __future__ import division
import math
import re

from bson import ObjectId, errors
import pysam
from bx.intervals.intersection import Intersecter, Interval

import capsid


db, logger = None, None
summary = {'xeno_mapped': 0, 'human_mapped': 0, 'human_unmapped': 0}
intersecters = {}


def _file2dict(f, sep):
    '''
    Saves the file as a dictionary using the accession/gi as keys:
        {line[1]: line[0]}
    '''

    genomes = {}

    with open(f, 'rU') as fh:
        for line in fh.readlines():
            val, key = line.strip().split(sep)
            genomes[key] = val

    return genomes


def _db2dict(column):
    '''
    Returns a dictionary from the database to compare with _file2dict in the format:
        {'_id':genome['_id'], 'gi':genome['gi']}
    '''

    return db.genome.find({}, {'_id':1, column:1})


def _lookup_custom(f, sep, col):
    '''
    Will match the genome header info with genome Ids in the database using
    a shared column. For best results use either gi or accession
    '''
    logger.debug('Entering _lookup_custom')

    # parse lookup file
    # [{id:col_value},...]
    f_genomes = _file2dict(f, sep)
    # Query db
    # [{_id:_id,col:col_value},...]
    db_genomes = _db2dict(col)
    # merge lists to return _id:id
    genomes = dict((f_genomes[d[col]],d['_id']) for d in db_genomes if d[col] in f_genomes)

    logger.debug('Exiting _lookup_custom')
    return genomes


def _lookup(from_file, genomes, method):
    '''Take the input for the header of the BAM return to ObjectId of the Genome in the Database'''

    if method == 'default':
        # Already using ObjectIds
        g = ObjectId(from_file)
    elif method == 'fasta':
        # TODO parse out gi from header
        gi = re.search('^>gi:(.+?)\|', from_file)
        g = db.genome.find_one({'gi': gi}, {'_id': 1})['_id']
    elif method == 'custom':
        # If the genomeId in the bam file cannot be mapped to the db, assume it is from a junction and
        # assign gId = None. Later any reads that map to a junction and a virus will have isHuman:1 for the vg,
        # but the hg will not be saved
        # Replace id in file with id in MongoDB
        try:
            g = genomes[from_file]
        except KeyError:
            g = None

    return g


def _find_gene_hits(mapped):
    '''Find the mapped hits that fall within a gene'''

    global intersecters
    hits = {}

    # Loop through mapped to see if there is a hit
    for m in mapped:
            # Get tree
            try:
                intersecter = intersecters[m['genomeId']]
            except:
                genes = db.feature.find({"genomeId": m['genomeId'], "type": "gene"})

                # Initialize tree
                intersecter = Intersecter()

                for gene in genes:
                    # Interval end is exclusive, need to +1 to line up with actual position
                    try:
                        intersecter.add_interval(Interval(gene['start'], gene['end'] + 1))
                    except AssertionError:
                        logger.error('Gene: {0} - start is greater than stop. Skipping...'.format(gene['_id']))

                intersecters[m['genomeId']] = intersecter

            # Check if mapped hits any gene, save readId
            if intersecter.find(m['refStart'], m['refEnd']):
                hits[m['_id']] = 1

    return hits


def meta_data(align):
    '''Use the alignment name to find sample and project'''

    meta = {}

    meta['alignment'] = db.alignment.find_one({"name": align})
    try:
        meta['sample'] = db.sample.find_one({"name": meta['alignment']['sample']})
    except TypeError:
        exit(logger.error("Alignment {0} not found in Database.".format(align)))
    except KeyError:
        exit(logger.error("Sample cannot be found. Please make sure the Alignment has been created properly."))

    logger.debug('Subtraction for alignment: {0}'.format(align))

    return meta


def _mapped(align, genome, meta, human=False):
    '''Build the mapped hit document'''

    fsum = math.fsum
    scores = [ord(m)-33 for m in align.qqual]
    try:
        mismatch = len(re.findall("\D", align.opt('MD')))
    except KeyError:
        mismatch = align.alen

    mapped = {
        "_id": ObjectId()
      , "readId": align.qname
      , "refStrand": -1 if align.is_reverse else 1
      , "refStart": align.pos + 1 # pysam is 0-based index
      , "refEnd": int(align.aend)
      , "alignLength": align.alen  # Length of alignment, using cigar information
      , "readLength": align.rlen  # Total Length of the read
      , "mapq": int(align.mapq)
      , "minQual": min(scores)
      , "avgQual": fsum(scores) / len(scores)
      , "miscalls": align.qqual.count('.')
      , "mismatch": mismatch
      , "sequencingType": meta['alignment']['type']
      , "platform": meta['alignment']['platform']
      , "project": meta['sample']['project']
      , "genomeId": genome
      , "sample": meta['sample']['name']
      , "alignment": meta['alignment']['name']
    }

    if human:
        mapped['isHuman'] = 1

    return mapped


def _unmapped(align, meta):
    '''Build the unmapped hit document'''

    scores = [ord(m)-33 for m in align.qqual]
    unmapped = {
        "readId": align.qname,
        "sequencingType": meta['alignment']['type'],
        "readLength": align.rlen,
        "minQual": min(scores),
        "avgQual": sum(scores, 0.0) / len(scores),
        "miscalls": align.qqual.count('.'),
        "project": meta['sample']['project'],
        "sample": meta['sample']['name'],
        "alignment": meta['alignment']['name'],
    }

    return unmapped


def parse_human(f, mapped_ids, meta, lookup, fetch):
    '''Iterate through the Human BAM file'''
    logger.debug('Entering _parse_human')

    global summary

    mapped, unmapped, intersecting_mapped_ids = [], [], {}

    logger.info('Scanning Human BAM file...')
    logger.debug("{0}".format(f))
    bamfile = pysam.Samfile(f, 'rb')

    for align in bamfile.fetch(until_eof=True):
        # If the alignemnt is unmapped and the readId is not in the xeno mapped ids save
        if fetch['unmapped'] and align.is_unmapped and align.qname not in mapped_ids:
            summary['human_unmapped'] += 1
            unmapped.append(_unmapped(align, meta))
            # When it gets too big insert and clear
            if len(unmapped) >= 100:
                db.unmapped.insert(unmapped)
                unmapped = []
        elif fetch['mapped'] and not align.is_unmapped and align.qname in mapped_ids:
            # Returns the genome Id from the DB - will return None if a Junction
            genome = _lookup(bamfile.getrname(align.tid), lookup['human_genomes'], lookup['human_method'])
            intersecting_mapped_ids[align.qname] = 1
            # Will append intersecting human genome hits, unless from junction
            if genome:
                summary['human_mapped'] += 1
                mapped.append(_mapped(align, genome, meta, True))
                # When it gets too big insert and clear
                if len(mapped) >= 100:
                    db.mapped.insert(mapped)
                    mapped = []

    bamfile.close()

    try:
        db.mapped.insert(mapped)
    except:
        pass

    try:
        db.unmapped.insert(unmapped)
    except:
        pass

    logger.debug('Exiting _parse_human')
    return intersecting_mapped_ids


def parse_xeno(f, meta, lookup, fetch):
    '''Iterate through the Virus BAM file'''

    logger.debug('Entering _parse_xeno')

    global summary

    mapped, mapped_ids, hits_ids = [], {}, {}

    logger.info('Scanning Xeno BAM file...')
    logger.debug("{0}".format(f))
    bamfile = pysam.Samfile(f, 'rb')

    # Loop through the mapped alignments in Bamfile
    for align in bamfile.fetch():
        # Create mapped dict and add to list
        summary['xeno_mapped'] += 1
        mapped_ids[align.qname] = 1
        genome = _lookup(bamfile.getrname(align.tid), lookup['xeno_genomes'], lookup['xeno_method'])
        mapped.append(_mapped(align, genome, meta))
        # When it gets too big insert and clear
        if len(mapped) >= 100:
            #if fetch['mapped']:
            db.mapped.insert(mapped)
            hits_ids.update(_find_gene_hits(mapped))
            mapped = []

    bamfile.close()

    # Insert any left over
    try:
        hits_ids.update(_find_gene_hits(mapped))
        db.mapped.insert(mapped)
    except:
        pass

    logger.debug('Exiting _parse_xeno')
    return mapped_ids, hits_ids


def main(args):
    '''Processes bams files and populates the DB with mapped entries that
    contain relavent read information'''

    global db, logger

    logger = args.logging.getLogger(__name__)
    db = capsid.connect(args)

    meta = meta_data(args.align)
    lookup = dict(zip(['xeno_method', 'human_method',
                       'xeno_file', 'human_file',
                       'xeno_column', 'human_column',
                       'xeno_separator', 'human_separator',
                       'xeno_genomes', 'human_genomes'],
                      [args.lookup_method[0], args.lookup_method[1],
                       args.lookup_file[0], args.lookup_file[1],
                       args.lookup_column[0], args.lookup_column[1],
                       args.lookup_separator[0], args.lookup_separator[1],
                       None, None]))

    # Generate DB Ids from the Lookup file
    if lookup['xeno_method'] == 'custom':
        logger.debug('Xeno using lookup file {0}'.format(lookup['xeno_file']))
        lookup['xeno_genomes'] = _lookup_custom(lookup['xeno_file'], lookup['xeno_separator'], lookup['xeno_column'])

    if lookup['human_method'] == 'custom':
        logger.debug('Human using lookup file {0}'.format(lookup['human_file']))
        lookup['human_genomes'] = _lookup_custom(lookup['human_file'], lookup['human_separator'], lookup['human_column'])

    # Fetch Mapped or Unmapped
    fetch = {
            'mapped': not args.fetch_only_unmapped
        ,   'unmapped': args.fetch_unmapped or args.fetch_only_unmapped
            }

    # Create list of reads(mapped/unmapped) from xeno file
    logger.info('Finding mapped alignments in the Xeno BAM file...')
    xeno_mapped_ids, xeno_hits_ids = parse_xeno(args.xeno, meta, lookup, fetch)

    # Update mapped with mapped.mapsGene : 1 if it hits a gene
    # update({},{}, Upsert, Manipulate, Safe, Multi) http://api.mongodb.org/python/current/api/pymongo/collection.html
    logger.info('Updating database for reads that hit a gene...')
    for h in capsid.chunks(xeno_hits_ids.keys(), 1000):
        db.mapped.update({"_id": {"$in": h}, "alignment": meta['alignment']['name']}, {"$set": {"mapsGene": 1}}, False, False, False, True)

    #updated = db.system_js.mapsGene(args.align)

    # Create list of reads(mapped/unmapped) from xeno file
    logger.info('Finding mapped/unmapped alignments in the Human BAM file...')
    intersecting_mapped_ids = parse_human(args.human, xeno_mapped_ids, meta, lookup, fetch)

    # Updating mapped with mapped.isHuman : 1 for intersecting read Ids
    logger.info('Updating database for reads that map in both Xeno and Human...')
    for m in capsid.chunks(intersecting_mapped_ids.keys(), 1000):
        db.mapped.update({"readId": {"$in": m}, "alignment": meta['alignment']['name']}, {"$set": {"isHuman": 1}}, False, False, False, True)

    # Summary Stats
    logger.info('Total mapped alignments found in Xeno BAM file: {0}'.format(summary['xeno_mapped']))
    logger.info('Total mapped alignments found in Human BAM file: {0}'.format(summary['human_mapped']))
    logger.info('Total mapped alignments added to database: {0}'.format(summary['xeno_mapped'] + summary['human_mapped']))
    logger.info('Xeno reads that hit a gene ("mapsGene":1): {0}'.format(len(xeno_hits_ids)))
    logger.info('Reads that map in Xeno and Human: {0}'.format(len(intersecting_mapped_ids)))
    logger.info('Intersecting unmapped alignments added to database: {0}'.format(summary['human_unmapped']))


if __name__ == '__main__':
    print 'This program should be run as part of the capsid package:\n\t$ capsid subtraction -h\n\tor\n\t$ /path/to/capsid/bin/capsid subtraction -h'
