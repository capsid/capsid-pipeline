#!/usr/bin/env python


# Copyright 2011(c) The Ontario Institute for Cancer Reserach. All rights reserved.
#
# This program and the accompanying materials are made available under the
# terms of the GNU Public License v3.0.
#
# You should have received a copy of the GNU General Public License along with
# this program.  If not, see <http://www.gnu.org/licenses/>.

import logging, logging.handlers
import os, sys, errno, ConfigParser

import configure, qfilter, gbloader, subtraction, statistics, fasta, intersect


def chunks(l, n):
    '''Yield successive n-sized chunks from l.'''

    for i in xrange(0, len(l), n):
        yield l[i:i+n]


def log(args):

    import colorize

    try:
        os.makedirs(os.path.expanduser('~/.capsid'))
    except OSError, e:
        if e.errno != errno.EEXIST:
            raise

    root_logger = logging.getLogger('')
    root_logger.setLevel(logging.DEBUG)

    console_handler = colorize.ColorizingStreamHandler(sys.stdout)
    console_format = '%(asctime)s %(message)s'
    console_date = '%H:%M:%S'
    console_handler.setFormatter(logging.Formatter(console_format, console_date))
    console_handler.setLevel(getattr(logging, args.logging))
    root_logger.addHandler(console_handler)

    file_handler = logging.handlers.RotatingFileHandler(os.path.expanduser('~/.capsid/capsid.log'))
    file_format = '%(asctime)s %(name)-20s %(levelname)-8s %(message)s'
    file_handler.setFormatter(logging.Formatter(file_format))
    file_handler.setLevel(logging.DEBUG)
    root_logger.addHandler(file_handler)

    traceback_log = logging.getLogger('traceback')
    traceback_log.propagate = False
    traceback_log.setLevel(logging.ERROR)
    traceback_log.addHandler(file_handler)

    return logging


def log_exception(exc_type, exc_value, traceback):
    logging.getLogger(__name__).error('{0}: {1}'.format(exc_type.__name__,exc_value))
    logging.getLogger('traceback').error(
            '{0}: {1}'.format(exc_type.__name__, exc_value),
            exc_info=(exc_type, exc_value, traceback),
            )
sys.excepthook = log_exception


def get_version(v):
    version = "{0}.{1}".format(v[0], v[1])
    if v[2]:
        version = "{0}.{1}".format(version, v[2])
    if v[3] != "f":
        version = "{0}{1}{2}".format(version, v[3], v[4])
        if v[5]:
            version = "{0}.dev{1}".format(version, v[5])
    return version

__version__ = get_version((1, 2, 0, "f", 0, 0))
