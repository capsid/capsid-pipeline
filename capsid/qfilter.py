#!/usr/bin/env python


# Copyright 2011(c) The Ontario Institute for Cancer Reserach. All rights reserved.
#
# This program and the accompanying materials are made available under the
# terms of the GNU Public License v3.0.
#
# You should have received a copy of the GNU General Public License along with
# this program.  If not, see <http://www.gnu.org/licenses/>.

from itertools import ifilter, izip, imap
from collections import namedtuple

from Bio import SeqIO

PairEnd = namedtuple('PairEnd', ['first', 'second'])

logger = None
threshold = None
limit = None


def quality_check(record):
    '''Ensures that enough of the bases pass the quality threshold'''

    scores = record.letter_annotations['phred_quality']
    failed = len([q for q in scores if q < int(threshold)])

    return failed <= limit


def qual_eq_seq(record):
    '''Checks that the length of the phred_quality is equal to the sequence '''

    scores = record.letter_annotations['phred_quality']

    return len(scores) == len(record.seq)


def valid_record(record):
    '''Return True if the reocrds passes all tests'''

    return qual_eq_seq(record) and quality_check(record)


def filter_reads_pe(records):
    '''Filters out both records if one fails a test'''

    return valid_record(records.first) and valid_record(records.second)


def filter_reads_se(record):
    '''Filters out records that do not pass the quality thresholds'''

    return valid_record(record)


def pair_end(args):
    ''' '''

    fq1 = SeqIO.parse(open(args.single, 'rU'), 'fastq')
    fq2 = SeqIO.parse(open(args.pair, 'rU'), 'fastq')

    return ifilter(filter_reads_pe, imap(PairEnd._make, izip(fq1, fq2)))


def single_end(args):
    ''' '''

    fq = SeqIO.parse(open(args.single, 'rU'), 'fastq')

    return ifilter(filter_reads_se, fq)


def main(args):
    ''' '''
    global logger, threshold, limit
    logger = args.logging.getLogger(__name__)
    threshold = int(args.threshold)
    limit = int(args.limit)

    print args

    results = pair_end(args) if args.pair else single_end(args)
    print results.next()
    print results.next()

if __name__ == '__main__':
    print 'This program should be run as part of the capsid package:\n\t$ capsid qfilter -h\n\tor\n\t$ /path/to/capsid/bin/capsid qfilter -h'
