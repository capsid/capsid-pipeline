#!/usr/bin/env python


# Copyright 2011(c) The Ontario Institute for Cancer Research. All rights reserved.
#
# This program and the accompanying materials are made available under the
# terms of the GNU Public License v3.0.
#
# You should have received a copy of the GNU General Public License along with
# this program.  If not, see <http://www.gnu.org/licenses/>.

import os, ConfigParser

from pymongo import Connection


def connect(args):
    '''Connects to MongoDB using settings from the config file. Exits if no connection can be made.'''

    logger = args.logging.getLogger(__name__)

    config = ConfigParser.ConfigParser()
    config.read(os.path.expanduser('~/.capsid/capsid.cfg'))

    try:
        address = config.get('MongoDB', 'host')
        port = int(config.get('MongoDB', 'port'))
        database = config.get('MongoDB', 'database')
        username = config.get('MongoDB', 'username')
        password = config.get('MongoDB', 'password')
    except:
        logger.error('There is an error in the CaPSID configuration file. Please run `capsid configure` again.`')
        exit()

    logger.debug('Connecting to {0}:{1} ({2})'.format(address, port, database))
    connection = Connection(address, port)
    admindb = connection.admin
    admindb.authenticate(username, password)

    return connection[database]
