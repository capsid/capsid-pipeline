#!/usr/bin/env python


# Copyright 2011(c) The Ontario Institute for Cancer Reserach. All rights reserved.
#
# This program and the accompanying materials are made available under the
# terms of the GNU Public License v3.0.
#
# You should have received a copy of the GNU General Public License along with
# this program.  If not, see <http://www.gnu.org/licenses/>.


from __future__ import division
from itertools import imap, count
import re

import pysam
from bx.intervals.intersection import Intersecter, Interval

import capsid


db, logger = None, None
x_it, r_it, u_it = count(), count(), count()


def parse_xeno(f):


def main(args):
    ''' '''
    global db, logger

    xeno_mapped_readids = parse_xeno(args.xeno)
    intersecting_mapped_readids = parse_ref(args.ref)
    logger = args.logging.getLogger(__name__)
    db = capsid.connect(args)


if __name__ == '__main__':
    print 'This program should be run as part of the capsid package:\n\t$ capsid subtraction -h\n\tor\n\t$ /path/to/capsid/bin/capsid subtraction -h'
