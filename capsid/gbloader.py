#!/usr/bin/env python
'''Functions for dealing with the records returned from SeqIO and creates lists of dictionary items that can be loaded into MongoDB'''


# Copyright 2011(c) The Ontario Institute for Cancer Reserach. All rights reserved.
#
# This program and the accompanying materials are made available under the
# terms of the GNU Public License v3.0.
#
# You should have received a copy of the GNU General Public License along with
# this program.  If not, see <http://www.gnu.org/licenses/>.


from itertools import count
from collections import namedtuple
import re

from Bio import SeqIO

from database import *


Qualifiers = namedtuple('Qualifiers', ['name', 'geneId', 'locusTag'])
Counter = namedtuple('Counter', ['records', 'genomes', 'pending', 'features', 'sequences'])

db = None
logger = None
counter = Counter(count(), count(), count(), count(), count())


def valid_seq(record):
    '''Filters out unknown sequences that are all 'N' so they are not saved'''

    m = re.search('[AGCT]', record.seq.tostring())

    return bool(m)


def extract_sequence(record, genome, delete=False):
    '''Returns a dictionary of the genome sequence'''
    global counter

    if valid_seq(record):
        counter.sequences.next()
        if delete:
            seq_id = genome.fs.get_last_version(genome.gi)._id
            genome.fs.delete(seq_id)
        genome.fs.put(record.seq.tostring(), filename=str(genome.gi), chunkSize=80)


def get_qualifiers(qualifiers):
    '''Returns dictionary of useful qualifiers for the feature'''

    gene = qualifiers['gene'][0] if 'gene' in qualifiers else None

    locusTag = qualifiers['locus_tag'][0] if 'locus_tag' in qualifiers else None

    try:
        geneId = [int(refs[7:]) for refs in qualifiers["db_xref"] if 'GeneID' in refs][0]
    except (IndexError, KeyError):
        geneId = None

    name = gene or locusTag or geneId or 'NA'

    return Qualifiers(name, geneId or u'NA', locusTag or 'NA')


def build_subfeatures(feature, genome):
    '''Needed for features with locations that are 'join' or 'order'. Recreates the parent features multiple times using the subfeatures' location.'''

    [build_feature(feature, genome, sf.location) for sf in feature.sub_features]


def build_feature(feature, genome, sf_location = None):
    '''Returns a dictionary of the feature'''
    global counter
    counter.features.next()

    qualifiers = get_qualifiers(feature.qualifiers)
    feature.location = sf_location or feature.location

    db.Feature({
        "name": qualifiers.name
        , "genome": genome.gi
        , "geneId": qualifiers.geneId
        , "locusTag": qualifiers.locusTag
        , "start": feature.location.nofuzzy_start + 1
        , "end": feature.location.nofuzzy_end
        , "operator": feature.location_operator
        , "strand": feature.strand
        , "type": feature.type
        }).save()


def extract_feature(feature, genome):
    '''Determines the feature's location operator and calls the appropriate build function'''

    has_subs = feature.location_operator in ['join', 'order']

    build_subfeatures(feature, genome) if has_subs else build_feature(feature, genome)


def extract_features(record, genome, delete=False):
    '''Returns a list of features belonging to the genome'''

    if delete: db.feature.remove({'genome': genome.gi})

    [extract_feature(f, genome) for f in record.features[1:] if f.type in ['gene', 'CDS']]


def extract_genome(record, delete):
    '''Returns a dictionary of the genome'''
    global counter
    counter.genomes.next()

    if delete: db.genome.remove({'gi': int(record.annotations['gi'])})

    genome = db.Genome({
        "gi": int(record.annotations['gi'])
        , "name": record.description
        , "accession": record.name
        , "version": record.annotations['sequence_version']
        , "length": record.features[0].location.nofuzzy_end
        , "strand":  record.features[0].strand
        , "taxonomy": record.annotations['taxonomy']
        , "organism": record.annotations['organism']
        , "pending": "features"
        })

    genome.save()

    return genome


def get_genome(record):
    '''Returns genome from Database'''
    global counter
    counter.pending.next()

    return db.Genome.one({'gi': int(record.annotations['gi'])})


def exists(record, genomes):
    '''Checks if the record's GI exists in the list of GIs from the database'''

    return int(record.annotations['gi']) in genomes


def parse_record(record, saved_genomes, pending_genomes, repair):
    '''Pull out genomic information from the GenBank file'''
    global counter
    counter.records.next()

    if not repair and exists(record, saved_genomes):
        return None

    is_pending = exists(record, pending_genomes)
    delete = repair or is_pending

    genome = get_genome(record) if is_pending else extract_genome(record, delete)

    if genome.pending == 'features':
        extract_features(record, genome, delete)
        genome.pending = 'sequence'
        genome.save()

    extract_sequence(record, genome, delete)
    del genome['pending']

    genome.save()


def get_pending_genomes():
    '''Returns a set of genomes with pending transactions'''

    return set(genome['gi'] for genome in db.Genome.find({'pending': {'$exists': True}}))


def get_saved_genomes():
    '''Returns a set of genomes currently in the db, uses GI to prevent duplication.'''

    return set(genome['gi'] for genome in db.Genome.find({'pending': {'$exists': False}}).hint([('_id', 1)]))


def parse_gb_file(f, repair):
    '''Use SeqIO to extract genome data from GenBank files'''
    logger.info('Scanning GenBank File {0}'.format(f))

    with open(f, 'rU') as fh:
        pending_genomes = get_pending_genomes() if not repair else []
        saved_genomes = get_saved_genomes() if not repair else []
        [parse_record(r, saved_genomes, pending_genomes, repair) for r in SeqIO.parse(fh, 'gb')]
        summary()


def summary():
    '''Logging summary of added records'''
    global counter

    # Counter starts at 0, so >it = count(); >print it.next(); >0
    # Using counter.*.next here sets it to the correct value for printing.
    records = counter.records.next()
    genomes = counter.genomes.next()
    pending = counter.pending.next()
    features = counter.features.next()
    sequences = counter.sequences.next()

    if records:
        logger.info("{0} Genomes found, {1} new Genomes added.".format(records, genomes))
        if pending:
            logger.info('{0} Genome found with pending transactions'.format(pending))
        if features:
            logger.info('{0} Features added successfully!'.format(features))
        if sequences:
            logger.info('{0} Sequences added successfully!'.format(sequences))
    else:
        logger.info('No Genomes found, make sure this is a GenBank file.')

    # Reset the counter for the next file
    counter = Counter(count(), count(), count(), count(), count())


def main(args):
    '''
    Reads through GenBank files and loads data into MongoDB in the genome, feature and sequence collections.
    python gloader.py g1.gbff gb2.gbff',
    '''

    global db, logger

    logger = args.logging.getLogger(__name__)
    db = connect(args)

    [parse_gb_file(f, args.repair) for f in args.files]


if __name__ == '__main__':
    print 'This program should be run as part of the capsid package:\n\t$ capsid gbloader -h\n\tor\n\t$ /path/to/capsid/bin/capsid gbloader -h'
