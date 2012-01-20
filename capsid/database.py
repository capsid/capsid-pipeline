#!/usr/bin/env python


# Copyright 2011(c) The Ontario Institute for Cancer Reserach. All rights reserved.
#
# This program and the accompanying materials are made available under the
# terms of the GNU Public License v3.0.
#
# You should have received a copy of the GNU General Public License along with
# this program.  If not, see <http://www.gnu.org/licenses/>.

import os, ConfigParser

from mongokit import *


class Project(Document):
    __collection__ = 'project'
    structure = {
        "name": basestring
        , "label": basestring
        , "description": basestring
        , "wikiLink": basestring
        , "role": [basestring]
        }
    use_dot_notation = True

class Sample(Document):
    __collection__ = 'sample'
    structure = {
        "name": basestring
        , "description": basestring
        , "project": basestring
        , "cancer": basestring
        , "role": basestring
        , "source": basestring
        }
    use_dot_notation = True

class Alignment(Document):
    __collection__ = 'alignment'
    structure = {
        "name": basestring
        , "project": basestring
        , "sample": basestring
        , "aligner": basestring
        , "platform": basestring
        , "type": basestring
        , "infile": basestring
        , "outfile": basestring
        }
    use_dot_notation = True

class Genome(Document):
    __collection__ = 'genome'
    structure = {
        "gi": int
        , "name": basestring
        , "accession": basestring
        , "version": int
        , "length": int
        , "strand": int
        , "taxonomy": [basestring]
        , "organism": basestring
        , "pending": basestring
        }
    gridfs = {'files': ['sequence']}
    use_dot_notation = True

class Feature(Document):
    __collection__ = 'feature'
    structure = {
        "name": basestring
        , "genome": int
        , "geneId": OR(int, unicode)
        , "locusTag": basestring
        , "start": int
        , "end": int
        , "operator": basestring
        , "strand": int
        , "type": basestring
        }
    use_dot_notation = True

class Mapped(Document):
    __collection__ = 'mapped'
    use_schemaless = True
    structure = {
        "readId": basestring
        , "refStrand": int
        , "refStart": int
        , "refEnd": int
        , "alignLength": int
        , "readLength": int
        , "mapq": int
        , "minQual": int
        , "avgQual": float
        , "miscalls": int
        , "mismatch": int
        , "pairEnd": int
        , "genome": int
        , "project": basestring
        , "sample": basestring
        , "alignment": basestring
        , "platform": basestring
        , "sequencingType": basestring
        , "mapsGene": int
        , "isRef": int
        }
    required_fields = ['readId', 'refStrand', 'refStart', 'refEnd', 'alignLength',
                       'readLength', 'mapq', 'minQual', 'avgQual', 'miscalls',
                       'mismatch', 'pairEnd', 'genome', 'project', 'sample',
                       'alignment', 'platform', 'sequencingType']
    use_dot_notation = True


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

    connection.register([Project, Sample, Alignment, Genome, Feature, Mapped])

    return connection[database]
