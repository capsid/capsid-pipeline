#!/usr/bin/env python
'''Functions for dealing with the records returned from SeqIO and creates lists of dictionary items that can be loaded into MongoDB'''


# Copyright 2011(c) The Ontario Institute for Cancer Reserach. All rights reserved.
#
# This program and the accompanying materials are made available under the
# terms of the GNU Public License v3.0.
#
# You should have received a copy of the GNU General Public License along with
# this program.  If not, see <http://www.gnu.org/licenses/>.


from Bio import SeqIO
from bson import ObjectId

import capsid

db, logger = None, None


def upload(collection, data, chunk=100):
    '''Inserts data into MongoDB in chunks rather than all at once.'''

    for c in capsid.chunks(data, chunk):
        db[collection].insert(c)


def _build_qualifiers(qualifiers):
    '''Returns dictionary of useful qualifiers for the feature'''

    # Parse out the GeneID as an int
    try:
        geneId = int([refs[7:] for refs in qualifiers["db_xref"] if 'GeneID' in refs][0])
    except (IndexError, KeyError):
        geneId = 'NA'

    if 'locus_tag' in qualifiers:
        locusTag = qualifiers["locus_tag"][0]
    else:
        locusTag = 'NA'

    if 'gene' in qualifiers:
        name = qualifiers["gene"][0]
    elif locusTag != 'NA':
        name = locusTag
    elif geneId != 'NA':
        name = geneId
    else:
        name = 'NA'

    q = {
        "geneId": geneId
      , "name": name
      , "locusTag": locusTag
    }

    return q


def _build_feature_join(feature, genome):
    '''Required to deal with features with locations that are join() or order(). Inserts each of the subfeatures as features linked the Genome'''

    fs = []
    fIds = []

    # Creates features using subfeature locations
    for subfeature in feature.sub_features:
        featureId = ObjectId()
        feature.location = subfeature.location
        fs.append(_build_feature(feature, featureId, genome))
        fIds.append(featureId)

    return fs, fIds


def _build_feature(feature, featureId,  genome):
    '''Returns a dictionary representation of a feature that can be inserted into MongoDB'''

    qualifiers = _build_qualifiers(feature.qualifiers)

    f = {
        "_id": featureId
        , "name": qualifiers["name"]
        , "genomeId": genome["_id"]
        , "geneId": qualifiers["geneId"]
        , "locusTag": qualifiers["locusTag"]
        , "start": feature.location.nofuzzy_start + 1
        , "end": feature.location.nofuzzy_end
        , "operator": feature.location_operator
        , "strand": feature.strand
        , "type": feature.type
        , "seqId": genome["seqId"]
    }

    return f


def extract_features(record, genome):
    '''Returns a list of features and feature Ids belonging to the genome'''

    # List of features for DB import
    fs = []
    # List of feature Ids to push into genome
    fIds = []

    # Loop through the feature records
    for feature in record.features[1:]:
        if feature.type == "gene" or feature.type == "CDS":
            # Deal with positions that look like (1..100,150..200)
            if feature.location_operator == 'join' or feature.location_operator == 'order':
                feats, fIds = _build_feature_join(feature, genome)
                fs.extend(feats)
                fIds.extend(fIds)
            else:
                featureId = ObjectId()
                fs.append(_build_feature(feature, featureId, genome))
                fIds.append(featureId)

    return fs, fIds


def extract_genome(record):
    '''Returns a dictionary representation of a genome that can be inserted into MongoDB'''

    # Build genome data
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


def main(args):
    '''
    Reads through GenBank files and loads data into MongoDB in the genome, feature and sequence collections.
    python gloader.py g1.gbff gb2.gbff',
    '''

    global db, logger

    logger = args.logging.getLogger(__name__)
    db = capsid.connect(args)

    for f in args.files:
        records_found = 0
        dbgenomes = {}
        genomes = []
        features = []
        seqs = []

        logger.info('Scanning GenBank File {0}'.format(f))
        try:
            logger.debug('Trying to parse {0} with SeqIO'.format(f))
            with open(f, 'rU') as fh:
                # Get a list of genomes currently in the db, to prevent duplication use gi
                logger.info('Getting list of Genomes from the Database')
                dbgenomes = dict((g['gi'], 1) for g in db.genome.find({}, {"_id":0, "gi":1}))

                for record in SeqIO.parse(fh, 'gb'):
                    records_found += 1
                    if int(record.annotations['gi']) not in dbgenomes:
                        # Extract Genomic data from the GenBank file
                        genome = extract_genome(record)

                        # If there is seq data in the genbank file add it to the db
                        # Human seq are too long to bother with
                        if record.annotations['organism'].capitalize() != 'Homo sapiens':
                            genome['seqId'] = ObjectId()
                            seqs.append({"_id": genome['seqId'], "seq": record.seq.tostring()})

                            # Extract genome features (Gene, CDS)
                            fs, genome['features'] = extract_features(record, genome)
                            features.extend(fs)

                        genomes.append(genome)
                        logger.debug('Added: {0}'.format(record.description))
                    else:
                        logger.debug('Skipped: {0}'.format(record.description))
        except IOError as (errno, strerror):
            logger.warning("I/O Error({0}): {1} {2}. Skipping...".format(errno, strerror, f))
            continue

        # Features and Sequences within this if so they are not added if Genome insert fails
        if records_found:
            try:
                logger.info('Adding Genomes...')
                upload('genome', genomes)
                logger.info("{0} Genomes found, {1} new Genomes added.".format(records_found, len(genomes)))

                if features:
                    try:
                        logger.info('Adding Features...')
                        upload('feature', features)
                        logger.info('{0} Features added successfully!'.format(len(features)))
                    except:
                        raise
                else:
                    logger.info('No Features found.')

                if seqs:
                    try:
                        logger.info('Adding Sequences...')
                        upload('sequence', seqs)
                        logger.info('{0} Sequences added successfully!'.format(len(seqs)))
                    except:
                        raise
                else:
                    logger.info('No Sequences found.')

            except:
                raise
        else:
            logger.info('No Genomes found, make sure this is a GenBank file.')


if __name__ == '__main__':
    print 'This program should be run as part of the capsid package:\n\t$ capsid gbloader -h\n\tor\n\t$ /path/to/capsid/bin/capsid gbloader -h'

