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
import sys

db, logger = None, None

def create_role(args):
    return {
        "authority": "ROLE_" + args.project.upper()
        }

def create_project(args):
    return {
        "description" : args.pdesc,
        "label" : args.project,
        "name" : args.pname or args.project,
        "roles" : ["ROLE_" + args.project.upper()],
        "version" : 0,
        "wikiLink" : args.link
        }

def main(args):
    '''Create Projects from the command line'''

    global db, logger

    logger = args.logging.getLogger(__name__)
    db = connect(args)
    
    try:
        db.project.insert(create_project(args), safe=True)
        logger.debug("project {0} inserted successfully".format(args.project))
        db.role.insert(create_role(args), safe=True)
        logger.debug("role {0} inserted successfully".format(args.project))
        logger.info("Project {0} added successfully to the database".format(args.project))
    except DuplicateKeyError:
        logger.error("Project {0} already exists".format(args.project))
        sys.exit(1)

if __name__ == '__main__':
    print 'This program should be run as part of the capsid package:\n\t$ capsid project -h\n\tor\n\t$ /path/to/capsid/bin/capsid project -h'