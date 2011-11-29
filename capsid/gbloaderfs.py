#!/usr/bin/env python
'''Functions for dealing with the records returned from SeqIO and creates lists of dictionary items that can be loaded into MongoDB'''


# Copyright 2011(c) The Ontario Institute for Cancer Reserach. All rights reserved.
#
# This program and the accompanying materials are made available under the
# terms of the GNU Public License v3.0.
#
# You should have received a copy of the GNU General Public License along with
# this program.  If not, see <http://www.gnu.org/licenses/>.

from itertools import imap, ifilter, ifilterfalse, repeat

from Bio import SeqIO
from bson import ObjectId
import gridfs
import multiprocessing

import capsid

logger = None


def extract_genome(record):
    '''Returns a dictionary representation of a genome that can be inserted into MongoDB'''

    genome = {
        "_id": ObjectId()
        , "gi": int(record.annotations['gi'])
        , "name": record.description
        , "accession": record.name
        , "version": record.annotations['sequence_version']
        , "length": record.features[0].location.nofuzzy_end
        , "strand":  record.features[0].strand
        , "taxonomy": record.annotations['taxonomy']
        , "organism": record.annotations['organism']
        }

    return genome

def exists(record, genomes):
    '''Checks if the record's GI exists in the list of GIs from the database'''
    return int(record.annotations['gi']) in genomes


def parse_record(record, genomes):
    '''Pull out genomic information from the GenBank file'''
    if exists(record, genomes):
        return None

    genome = extract_genome(record)

    if not is_human(record):
        genome['seqId'] = ObjectId()
        features, genome['features'] = extract_features(record)

    return record.annotations['gi']


def get_db_genomes():
    '''Get a list of genomes currently in the db, uses GI to prevent duplication.'''
    logger.info('Getting list of Genomes from the Database')
    return set(g['gi'] for g in db.genome.find({}, {"_id":0, "gi":1}))


def parse_gb_file(f):
    logger.info('Scanning GenBank File {0}'.format(f))

    try:
        logger.debug('Trying to parse {0} with SeqIO'.format(f))
        with open(f, 'rU') as fh:
            dbg = get_db_genomes()
            return filter(None, (parse_record(r, dbg) for r in SeqIO.parse(fh, 'gb')))
    except IOError as (errno, strerror):
        logger.warning("I/O Error({0}): {1} {2}. Skipping...".format(errno, strerror, f))


def main(args):
    '''
    Reads through GenBank files and loads data into MongoDB in the genome, feature and sequence collections.
    python gloader.py g1.gbff gb2.gbff',
    '''

    global logger

    logger = args.logging.getLogger(__name__)
    db = capsid.connect(args)
    gfs = gridfs.GridFS(db, 'sequence')

    records = imap(parse_gb_file, args.files)

    for recs in records:
        for r in recs:
            pass

if __name__ == '__main__':
    print 'This program should be run as part of the capsid package:\n\t$ capsid gbloader -h\n\tor\n\t$ /path/to/capsid/bin/capsid gbloader -h'
