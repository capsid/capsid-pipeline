#!/usr/bin/env python


# Copyright 2011(c) The Ontario Institute for Cancer Reserach. All rights reserved.
#
# This program and the accompanying materials are made available under the
# terms of the GNU Public License v3.0.
#
# You should have received a copy of the GNU General Public License along with
# this program.  If not, see <http://www.gnu.org/licenses/>.

logger = None

def pair_end(args)


def main(args):
    ''' '''
    global logger
    logger = args.logging.getLogger(__name__)
    print args

    results = pair_end(args) if args.pair else single_end(args)


if __name__ == '__main__':
    print 'This program should be run as part of the capsid package:\n\t$ capsid qfilter -h\n\tor\n\t$ /path/to/capsid/bin/capsid qfilter -h'
