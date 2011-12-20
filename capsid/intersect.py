#!/usr/bin/env python

# Copyright 2011(c) The Ontario Institute for Cancer Reserach. All rights reserved.
#
# This program and the accompanying materials are made available under the
# terms of the GNU Public License v3.0.
#
# You should have received a copy of the GNU General Public License along with
# this program.  If not, see <http://www.gnu.org/licenses/>.

import os
import subprocess

logger = None
temp = None

def clean_ext(name):
    '''Pull off extension'''

    return name.rpartition('.')[0]


def collapse_file(f, name):
    ''' '''
    logger.info('Collapsing {0}...'.format(f))
    logger.debug('Collapsed {0} saved as {1}'.format(f, name))

    p1 = subprocess.Popen(['awk', 'ORS=NR%4?"\\t":"\\n"', f], stdout=subprocess.PIPE)
    p2 = subprocess.Popen(['sort', '-u', '-T', temp, '-o', name], stdin=p1.stdout, stdout=None)
    p1.stdout.close()


def intersect_files(f):
    ''' '''

    collapse_file(f, 'temp.fq')

    logger.info('Intersecting with {0}...'.format(f))
    logger.debug('Intersecting reads output to out.fq')

    subprocess.call(['comm -12 saved.fq temp.fq > out.fq'], shell=True)

    logger.debug('out.fq moved to saved.fq')
    os.rename('out.fq', 'saved.fq')


def finalize(args):
    '''Clean up temp files'''
    logger.info('Cleaning up temporary files...')

    if os.path.isfile('temp.fq'):
        logger.debug('Deleting temp.fq...')
        os.remove('temp.fq')

    files = args.files
    f = '_'.join([clean_ext(f) for f in files]) + '.intersect.fastq'

    if os.path.isfile('saved.fq'):
        logger.info('Saving intersecting reads to {0}...'.format(f))
        subprocess.call(['tr "\t" "\n" < saved.fq > {0}'.format(f)], shell=True)
        logger.debug('Deleting saved.fq...')
        os.remove('saved.fq')


def main(args):
    ''' '''
    global logger, temp

    logger = args.logging.getLogger(__name__)
    temp = args.temp

    collapse_file(args.files[0], 'saved.fq')
    map(intersect_files, args.files[1:])

    p = subprocess.Popen(["wc", "-l", "saved.fq"], stdout=subprocess.PIPE)
    reads = p.communicate()[0].partition(' ')[0]

    finalize(args)

    logger.info('{0} intersecting unmapped reads found.'.format(reads))


if __name__ == '__main__':
    print 'This program should be run as part of the capsid package:\n\t$ capsid qfilter -h\n\tor\n\t$ /path/to/capsid/bin/capsid qfilter -h'
