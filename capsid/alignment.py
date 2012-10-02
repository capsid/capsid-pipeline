#!/usr/bin/env python

# Copyright 2012(c) The Ontario Institute for Cancer Research. All rights reserved.
#
# This program and the accompanying materials are made available under the terms of the GNU Public License v3.0.
# You should have received a copy of the GNU General Public License along with 
# this program. If not, see <http://www.gnu.org/licenses/>.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY 
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
# OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT 
# SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, 
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED 
# TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
# OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER 
# IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN 
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

from pymongo.errors import DuplicateKeyError

from database import *
import sample

db, logger = None, None


def check_sample(args):
    if (not db.sample.find_one({"name": args.sample})):
        args.sdesc = None
        args.cancer = None
        args.role = None
        args.source = None
        sample.main(args)


def create_alignment(args):
    return {
        "aligner" : args.aligner,
        "infile" : args.infile,
        "name" : args.align,
        "outfile" : args.outfile,
        "platform" : args.platform,
        "project" : args.project,
        "sample" : args.sample,
        "type" : args.type,
        "version" : 0
        }


def main(args):
    '''Create Alignment from the command line'''

    global db, logger

    logger = args.logging.getLogger(__name__)
    db = connect(args)
    
    check_sample(args)
    
    try:
        db.alignment.insert(create_alignment(args), safe=True)
        logger.debug("alignment {0} inserted successfully".format(args.align))
        logger.info("Alignment {0} has been added to {1}/{2}".format(args.align, args.project, args.sample))
    except DuplicateKeyError:
        logger.info("Alignment {0} already exists".format(args.align))


if __name__ == '__main__':
    print 'This program should be run as part of the capsid package:\n\t$ capsid alignment -h\n\tor\n\t$ /path/to/capsid/bin/capsid alignment -h'