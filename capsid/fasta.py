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

from database import *


def fasta_output(genome, out):
    '''Output Genome in Fasta format'''

    out.write('>gi|{gi}|ref|{accession}.{version}| {name}\n'.format(**genome))
    [out.write('{0}\n'.format(line)) for line in genome.fs.get_last_version(str(genome.gi))] 


def main(args):
    '''Fasta Output of Genomes in the Database'''

    logger = args.logging.getLogger(__name__)
    db = connect(args)

    # By default the query will not include Human
    query = {"organism": {'$ne': "Homo sapiens"}}

    if args.organism: query['organism'] = {'$in': args.organism}
    if args.taxonomy: query['taxonomy'] = {'$in': args.taxonomy}
    if args.ref: query['accession'] = args.ref
    if args.gi: query['gi'] = int(args.gi)

    genomes = db.Genome.find(query)
    logger.info('Found {0} genomes'.format(genomes.count()))

    logger.debug('Writing Fasta output to {0}...'.format(args.output))
    with open(args.output, 'w') as out:
        [fasta_output(genome, out) for genome in genomes]

    logger.info('{0} created.'.format(args.output))


if __name__ == '__main__':
    print 'This program should be run as part of the capsid package:\n\t$ capsid fasta -h\n\tor\n\t$ /path/to/capsid/bin/capsid fasta -h'
