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
import project
import sys

db, logger = None, None


def check_project(args):
    project = db.project.find_one({"label": args.project})

    if (not project):
        logger.error("Project {0} does not exists".format(args.project))
        sys.exit(1)
    else:
        return project


def create_sample(project, args):
    return {
        "cancer" : args.cancer,
        "description" : args.sdesc,
        "name" : args.sample,
        "projectLabel" : args.project,
        "projectId" : project["_id"],
        "role" : args.role,
        "source" : args.source,
        "version" : 0
        }


def main(args):
    '''Create Projects from the command line'''

    global db, logger

    logger = args.logging.getLogger(__name__)
    db = connect(args)
    
    project = check_project(args)
    
    try:
        db.sample.insert(create_sample(project, args), safe=True)
        logger.debug("sample {0} inserted successfully".format(args.sample))
        logger.info("Sample {0} has been added to {1}".format(args.sample, args.project))
    except DuplicateKeyError as inst:
        logger.error("Error {0}".format(inst))
        logger.error("Sample {0} already exists in project: {1}".format(args.sample, project["label"]))
        sys.exit(1)


if __name__ == '__main__':
    print 'This program should be run as part of the capsid package:\n\t$ capsid sample -h\n\tor\n\t$ /path/to/capsid/bin/capsid sample -h'