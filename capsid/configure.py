#!/usr/bin/env python


# Copyright 2011(c) The Ontario Institute for Cancer Reserach. All rights reserved.
#
# This program and the accompanying materials are made available under the
# terms of the GNU Public License v3.0.
#
# You should have received a copy of the GNU General Public License along with
# this program.  If not, see <http://www.gnu.org/licenses/>.


import pymongo
import getpass
import os
import base64
import ConfigParser

from database import *

db, logger = None, None


def ensure_indexes():
    '''Creates the MongoDB Indices needed to effeciently run CaPSID'''

    # Project
    logger.debug('Adding Project Index')
    db.project.ensure_index('label', unique=True)

    # Sample
    logger.debug('Adding Sample Indices')
    db.sample.ensure_index('name', unique=True)
    db.sample.ensure_index('project')

    # Alignment
    logger.debug('Adding Alignment Index')
    db.alignment.ensure_index('name', unique=True)

    # Genome
    logger.debug('Adding Genome Indices')
    db.genome.ensure_index('gi', unique=True)
    db.genome.ensure_index('accession', unique=True)
    db.genome.ensure_index('pending', sparse=True)

    # Feature
    logger.debug('Adding Feature Indices')
    db.feature.ensure_index([('genome', pymongo.ASCENDING), ('type', pymongo.ASCENDING), ('strand', pymongo.ASCENDING)])
    db.feature.ensure_index('genome')
    db.feature.ensure_index('start')

    # Mapped
    logger.debug('Adding Mapped Indices')
    db.mapped.ensure_index('refStart')
    db.mapped.ensure_index([('genome', pymongo.ASCENDING), ('sample', pymongo.ASCENDING), ('mapsGene', pymongo.ASCENDING)], sparse=True)
    db.mapped.ensure_index([('genome', pymongo.ASCENDING), ('project', pymongo.ASCENDING), ('mapsGene', pymongo.ASCENDING)], sparse=True)
    db.mapped.ensure_index([('readId', pymongo.ASCENDING), ('_id', pymongo.ASCENDING)])

    # User
    logger.debug('Adding User Index')
    db.user.ensure_index('username', unique=True)

    # Role
    logger.debug('Adding Role Index')
    db.role.ensure_index('authority', unique=True)

    # UserRole
    logger.debug('Adding UserRole Index')
    db.userRole.ensure_index([('role', pymongo.ASCENDING), ('user', pymongo.ASCENDING)], unique=True)

    # GridFS
    #db.fs.chunks.ensure_Index([('files_id', pymongo.ASCENDING), ('n', pymongo.ASCENDING)], unique=True)

def genome_samples():
    '''Calculate how many samples hit the genomes'''

    logger.debug('Deleting old GenomeSamples function...')
    del db.system_js.gs
    logger.debug('Adding GenomeSamples function...')
    db.system_js.gs = """
function() {
  genomes=db.genome.find({}, {_id:0, gi:1});
  genomes.forEach(function(g){
    s=db.mapped.distinct('sample',{genome:g.gi});
    db.genome.update({'gi':g.gi}, {$set: {samples: s}});
  });
}
"""


def setup_config(args):
    '''Saves MongoDB settings to a configuration file'''

    config = ConfigParser.SafeConfigParser()

    config.add_section('MongoDB')

    print 'Please enter the settings for your MongoDB server:'
    config.set('MongoDB', 'host', args.host or raw_input('Host [localhost]: ') or 'localhost')
    config.set('MongoDB', 'port', args.port or raw_input('Port [27017]: ') or '27017')
    config.set('MongoDB', 'database', args.database or raw_input('Database [capsid]: ') or 'capsid')
    config.set('MongoDB', 'username', args.username or raw_input('Username [capsid]: ') or 'capsid')
    config.set('MongoDB', 'password', args.password or getpass.getpass('Password: '))

    # Writing our configuration file
    with open(os.path.expanduser('~/.capsid/capsid.cfg'), 'wb') as configfile:
        config.write(configfile)


def main(args):
    '''Setup MongoDB for use by CaPSID'''

    global db, logger

    logger = args.logging.getLogger(__name__)

    # Setup config files
    setup_config(args)

    db = connect(args)

    # Add all req indexes to MongoDB
    logger.info('Setting up MongoDB...')
    logger.info('Adding Indices...')
    ensure_indexes()
    logger.info('Adding JavaScript Functions...')
    genome_samples()

    logger.info('Done!')


if __name__ == '__main__':
    print 'This program should be run as part of the capsid package:\n\t$ capsid configure -h\n\tor\n\t$ /path/to/capsid/bin/capsid configure -h'
