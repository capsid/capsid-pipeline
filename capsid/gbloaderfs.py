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

import capsid

db, logger = None, None


def insert_records(record):
    '''Inserts Genome, Features and Sequence into MongoDB'''

    print record


def extract_features(record, genome):
    '''Returns a list of features and feature Ids belonging to the genome'''

    return ['a', 'b', 'c'], [1, 2, 3]


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


def is_human(record):
    '''Checks if the genome is human'''

    return record.annotations['organism'].capitalize() is 'Homo sapiens'


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
        features, genome['features'] = extract_features(record, genome)
        sequence = ({'_id': genome['seqId'], 'seq': record.seq.tostring()})

    return genome, features, sequence


def get_db_genomes():
    '''Get a list of genomes currently in the db, uses GI to prevent duplication.'''

    logger.info('Getting list of Genomes from the Database')
    return set(g['gi'] for g in db.genome.find({}, {"_id":0, "gi":1}))


def parse_gb_file(f):
    '''Use SeqIO to extract genome data from GenBank files'''
    logger.info('Scanning GenBank File {0}'.format(f))

    with open(f, 'rU') as fh:
        dbg = get_db_genomes()
        records = (parse_record(r, dbg) for r in SeqIO.parse(fh, 'gb'))
        # parse_record will return None if the genome is already in the Database,
        # and using filter(None, records) will put the entire thing in memory.
        # This way it only deals with 1 record at a time and skips the 'None's.
        [insert_records(r) for r in records if r]


def main(args):
    '''
    Reads through GenBank files and loads data into MongoDB in the genome, feature and sequence collections.
    python gloader.py g1.gbff gb2.gbff',
    '''

    global db, logger

    logger = args.logging.getLogger(__name__)
    db = capsid.connect(args)
    gfs = gridfs.GridFS(db, 'sequence')

    map(parse_gb_file, args.files)


if __name__ == '__main__':
    print 'This program should be run as part of the capsid package:\n\t$ capsid gbloader -h\n\tor\n\t$ /path/to/capsid/bin/capsid gbloader -h'
