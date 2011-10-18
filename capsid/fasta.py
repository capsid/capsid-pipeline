#!/usr/bin/env python
"""Output Genomes collection as fasta file"""


# Copyright 2011(c) The Ontario Institute for Cancer Reserach. All rights reserved.
#
# This program and the accompanying materials are made available under the
# terms of the GNU Public License v3.0.
#
# You should have received a copy of the GNU General Public License along with
# this program.  If not, see <http://www.gnu.org/licenses/>.


import re

import capsid


def main(args):
    '''Fasta Output of Genomes in the Database'''

    logger = args.logging.getLogger(__name__)
    db = capsid.connect(args)

    # By default the query will not include Human
    query = {"organism": {'$ne': "Homo sapiens"}}

    if args.organism:
        query['organism'] = {'$in': args.organism}
    if args.taxonomy:
        query['taxonomy'] = {'$in': args.taxonomy}

    genomes = db.genome.find(query)
    logger.info('Found {0} genomes'.format(genomes.count()))
    out = open(args.output, 'w')
    logger.info('Output genomes in Fasta format...')
    logger.debug('Writting to {0}...'.format(args.output))
    for genome in genomes:
        out.write('>{gi}|{accession} {name}\n'.format(**genome))
        seq = db.sequence.find_one({"_id":genome['seqId']}, {"_id":0, "seq":1})
        if seq:
            s = re.sub(r'(\w{80})', r'\1\n', seq['seq'])
            out.write('{0}\n'.format(s))

    logger.info('Done.')

if __name__ == '__main__':
    print 'This program should be run as part of the capsid package:\n\t$ capsid fasta -h\n\tor\n\t$ /path/to/capsid/bin/capsid fasta -h'
