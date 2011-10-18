#!/usr/bin/env python


# Copyright 2011(c) The Ontario Institute for Cancer Reserach. All rights reserved.
#
# This program and the accompanying materials are made available under the
# terms of the GNU Public License v3.0.
#
# You should have received a copy of the GNU General Public License along with


import os

import capsid


def main(args):
    '''Converts Sam files into BAM files and collapses them using Picard'''

    logger = args.logging.getLogger(__name__)
    picard = args.picard or '.'
    files = args.files

    for f in files:
        # Convert Sam to Bam
        logger.info('Converting {0} to bam...'.format(f))
        os.system('samtools view -bS -o ' + f + '.bam ' + f)

    # Merge all Bam files
    logger.info('Merging bam files...')
    os.system('samtools merge merged.bam ' + " ".join([f+'.bam' for f in files]))
    # Sort bam file alignment by left most coordinate
    logger.info('Sorting...')
    os.system('samtools sort merged.bam merged.sorted')
    # Index sorted alignments for fast random access - creates .bai file
    logger.info('Indexing...')
    os.system('samtools index merged.sorted.bam')
    # Remove duplicates using Picard
    logger.info('Removing duplicates with Picard...')
    os.system('java -jar {0}/MarkDuplicates.jar INPUT=merged.sorted.bam OUTPUT=merged.sorted.marked.bam METRICS_FILE=duplicateMETRICS.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT'.format(picard))
    logger.info('Indexing final bam file...')
    os.system('samtools index merged.sorted.marked.bam')


if __name__ == "__main__":
    print 'This program should be run as part of the capsid package:\n\t$ capsid sam2bam -h\n\tor\n\t$ /path/to/capsid/bin/capsid sam2bam -h'
