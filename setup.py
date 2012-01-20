#!/usr/bin/env python


# Copyright 2011(c) The Ontario Institute for Cancer Reserach. All rights reserved.
#
# This program and the accompanying materials are made available under the
# terms of the GNU Public License v3.0.
#
# You should have received a copy of the GNU General Public License along with
# this program.  If not, see <http://www.gnu.org/licenses/>.i


from setuptools import setup
import sys, os


if sys.version_info < (2, 6):
    exec('raise Error, "Python 2.6 or later is required"')


def read(*path):
        return open(os.path.join(os.path.abspath(os.path.dirname(__file__)), *path)).read()


VERSION = '1.2'
README = read('README')
NEWS = read('NEWS.rst')
install_requires = ['cython', 'numpy', 'anyjson', 'pysam', 'pymongo', 'biopython', 'bx-python', 'mongokit']

if sys.version_info < (2, 7):
    install_requires.append('argparse')

# Get classifiers from http://pypi.python.org/pypi?%3Aaction=list_classifiers
classifiers = """
Development Status :: 5 - Production/Stable
License :: OSI Approved :: GNU General Public License (GPL)
Environment :: Console
Intended Audience :: Science/Research
Intended Audience :: Developers
Programming Language :: Python :: 2.6
Programming Language :: Python :: 2.7
Topic :: Scientific/Engineering :: Bio-Informatics
Operating System :: Unix
"""

config = {
    'name': 'capsid',
    'version': VERSION,
    'description': 'CaPSID: Computational Pathogen Sequence Identification',
    'long_description': README + '\n\n' + NEWS,
    'license': 'GNU General Public License, Version 3.0',
    'author': 'Shane Wilson',
    'author_email': 'shane.wilson@oicr.on.ca',
    'url': 'https://github.com/capsid/capsid',
    'download_url': 'https://github.com/capsid/capsid',
    'classifiers': filter(None, classifiers.split("\n")),
    'scripts': ['bin/capsid'],
    'packages': ['capsid'],
    'zip_safe': True,
   'install_requires': install_requires
}


def setup_package():
    """Setup Package"""

    setup(**config)


if __name__ == '__main__':
    setup_package()
