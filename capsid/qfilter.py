#!/usr/bin/env python

# Copyright 2011(c) The Ontario Institute for Cancer Reserach. All rights reserved.
#
# This program and the accompanying materials are made available under the
# terms of the GNU Public License v3.0.
#
# You should have received a copy of the GNU General Public License along with
# this program.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import division
from itertools import count, ifilter, izip, imap, repeat
from collections import namedtuple
import os
import subprocess

from Bio import SeqIO

Records = namedtuple('Records', ['single', 'pair'])
FileName = namedtuple('FileName', ['name', 'dot', 'ext'])
Counter = namedtuple('Counter', ['records', 'saved'])

logger = None
threshold = None
limit = None
temp = None
counter = Counter(count(), count())


def clean_ext(name):
    '''Pull off .fastq extension'''

    f = FileName._make(name.rpartition('.'))

    return f.name if f.ext == 'fastq' else name


def make_fastq_file(f):
    '''Turn 1-lined sorted temp file into fastq format'''

    fh = f + '.quality.fastq'
    logger.debug('Filtered Fastq: {0}'.format(fh))

    logger.info('Generating filtered Fastq file...')

    fq_sorted = open(f + '.f.c.temp')
    fastq = open(fh, 'w')

    [fastq.write('{0}\n{1}\n{2}\n{3}'.format(*line.split('\t')))
     for line in fq_sorted]


def collapse_file(f):
    ''' '''

    fc = f + '.f.c.temp'
    logger.debug('Collapsed File: {0}'.format(fc))
    logger.info('Collapsing {0}...'.format(fc))
    subprocess.call(['sort', '-t\t', '+1', '-2', '-u', '-T', temp, f + '.f.temp', '-o', fc])


def sortable_output(record, fq_single, fq_pair):
    '''Output fastq records as 1-line into temp file so it can be sorted'''
    counter.saved.next()

    fq_single.write('@{description}\t{seq}\t+{description}\t{quality}\n'.format(
            description = record.single.description,
            seq = record.single.seq,
            quality = SeqIO.QualityIO._get_sanger_quality_str(record.single))
            )

    if fq_pair:
        fq_pair.write('@{description}\t{seq}\t+{description}\t{quality}\n'.format(
                description = record.pair.description,
                seq = record.pair.seq,
                quality = SeqIO.QualityIO._get_sanger_quality_str(record.pair))
                )


def make_sortable_file(records, f_single, f_pair):
    ''' '''

    ft_single = f_single + '.f.temp'
    ft_pair = f_pair + '.f.temp' if f_pair else None

    logger.debug('Temp Filtered File: {0}'.format(ft_single))
    if f_pair: logger.debug('Temp Filtered File: {0}'.format(ft_pair))

    logger.info('Creating temporary files for collapsing...')

    fh_single = open(ft_single, 'w')
    fh_pair = open(ft_pair, 'w') if f_pair else None
    
    [sortable_output(record, fh_single, fh_pair) for record in records]
    
    fh_single.close()
    if fh_pair: fh_pair.close()


def sort_unique(records, args):
    '''Sort fastq file and clean up'''

    f_single = clean_ext(args.single)
    f_pair = clean_ext(args.pair) if args.pair else None

    make_sortable_file(records, f_single, f_pair)

    #collapse if not pair end
    if not f_pair:
        collapse_file(f_single)
    make_fastq_file(f_single)

    if f_pair:
        # Thinking about how best to do the collapsing for pair end
        #collapse_file(f_pair)
        make_fastq_file(f_pair)


def quality_check(scores):
    '''Ensures that enough of the bases pass the quality threshold'''

    failed = len([q for q in scores if q < int(threshold)])

    return failed <= limit


def qual_eq_seq(scores):
    '''Checks that the length of the phred_quality is equal to the sequence '''

    return len(scores) == len(record.seq)


def valid_record(record):
    '''Return True if the records passes all tests'''
    if not record:
        return

    scores = record.letter_annotations['phred_quality']

    return len(scores) == len(record.seq) and quality_check(scores)


def filter_reads(records):
    '''Filters out records if both fail test'''
    counter.records.next()
    
    return valid_record(records.single) or valid_record(records.pair)


def parse_fastq(args):
    ''' '''

    s = 'pair end' if args.pair else 'single end'
    logger.info('Reading FastQ files as {}...'.format(s))

    fq1 = SeqIO.parse(open(args.single, 'rU'), args.format)
    fq2 = SeqIO.parse(open(args.pair, 'rU'), args.format) if args.pair else repeat(None)

    return ifilter(filter_reads, imap(Records._make, izip(fq1, fq2)))


def clean_up(f):
    '''Remove Temp files'''

    if os.path.isfile(f + '.f.temp'):
        logger.debug('Deleting {0}'.format(f + '.f.temp'))
        os.remove(f + '.f.temp')
    if os.path.isfile(f + '.f.c.temp'):
        logger.debug('Deleting {0}'.format(f + '.f.c.temp'))
        os.remove(f + '.f.c.temp')


def summary(args):
    ''' '''

    # Counter starts at 0, so >it = count(); >print it.next(); >0
    # Using counter.*.next here sets it to the correct value for printing.
    records = counter.records.next()
    saved = counter.saved.next()


    if records:
        percent = (saved/records) * 100
        logger.info('{0} filtered and saved to {1}.quality.fastq'.format(args.single, clean_ext(args.single)))
        if args.pair:
            logger.info('{0} collapsed and saved to {1}.quality.fastq'.format(args.pair, clean_ext(args.pair)))
        logger.info('{0} of {1} ({2:.2f}%) records passed filter.'.format(saved, records, percent))
        if not args.pair:
            p = subprocess.Popen(["wc", "-l", clean_ext(args.single) + ".quality.fastq"], stdout=subprocess.PIPE)
            reads = int(p.communicate()[0].partition(' ')[0]) // 4
            percent_c = (reads/saved) * 100
            logger.info('{0} of {1} ({2:.2f}%) records left after collapsing.'.format(reads, saved, percent_c))
            logger.info('{0} of {1} ({2:.2f}%) records saved.'.format(reads, records, (reads/records)*100))
    else:
        logger.warning('No records found.')


def main(args):
    ''' '''
    global logger, temp, threshold, limit

    logger = args.logging.getLogger(__name__)
    temp = args.temp
    threshold = int(args.threshold)
    limit = int(args.limit)

    args.format = 'fastq-' + 'sanger'

    records = parse_fastq(args)
    sort_unique(records, args)

    logger.info('Removing temporary files...')
    clean_up(clean_ext(args.single))
    if args.pair: clean_up(clean_ext(args.pair))

    summary(args)


if __name__ == '__main__':
    print 'This program should be run as part of the capsid package:\n\t$ capsid qfilter -h\n\tor\n\t$ /path/to/capsid/bin/capsid qfilter -h'
