#!/usr/bin/env python

# Copyright 2011(c) The Ontario Institute for Cancer Research. All rights reserved.
#
# This program and the accompanying materials are made available under the
# terms of the GNU Public License v3.0.
#
# You should have received a copy of the GNU General Public License along with
# this program.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import division
from collections import namedtuple
from functools import partial
import multiprocessing
import re

from bx.intervals.intersection import Intersecter, Interval

from database import *


db = None
logger = None

def merge(lst):
    ''' '''

    sorted_list = sorted([sorted(elem) for elem in lst])
    try:
        saved = sorted_list[0]
    except IndexError:
        return

    for start, end in sorted_list:
        if start <= saved[1]:
            saved[1] = max(saved[1], end)
        else:
            yield tuple(saved)
            saved[0] = start
            saved[1] = end

    yield tuple(saved)


def insert_stats(stats):
    '''Load sample statistics into the database'''

    db.statistics.save(stats)


def find_genes(genome):
    ''' '''

    return db.feature.find({"type": "gene", "genome": genome},
                           {"start":1, "end":1})


def init_intersecter(hits):
    ''' '''

    intersecter = Intersecter()

    for h in hits:
        intersecter.add_interval(Interval(h[0], h[1]))

    return intersecter


def calc_gene_coverage(intersecter, gene):
    ''' '''

    # Find hits that intersect with the gene
    overlap = intersecter.find(gene['start'], gene['end'])
    # Merge overlap to prevent counting the same bases more than once
    m_overlap = merge([(o.start, o.end) for o in overlap])
    # Get the total area over the gene that is covered
    length = sum([min(gene['end'], end) - max(gene['start'], start) + 1
                  for start, end in m_overlap])

    # Divide by the length of the gene to get the percent coverage
    return length / (gene['end'] - gene['start'] + 1)


def gene_coverage(query, genome):
    ''' '''

    # Merge the hits to prevent miscalculations due to overlapping alignments
    hits = merge([(hit['refStart'], hit['refEnd']) for hit in query])
    # Create Intersector Object and fill with merged hit locations
    intersecter = init_intersecter(hits)
    # Get a list of gene locations
    genes = find_genes(genome['gi'])

    coverage = [calc_gene_coverage(intersecter, gene) for gene in genes]

    try:
        maximum = max(coverage)
        mean = sum(coverage) / len(coverage)
    except ValueError, ZeroDivisionError:
        maximum = 0.0
        mean = 0.0

    return mean, maximum


def gene_hits(col, value, genome):
    ''' '''

    return db.mapped.find({col: value, "genome": genome, "mapsGene": {'$exists':  True}},
                          {"_id":0, "refStart":1, "refEnd":1})


def genome_coverage(query, genome):
    ''' '''

    # Merge the hits to prevent miscalculations due to overlapping alignments
    hits = merge([(hit['refStart'], hit['refEnd']) for hit in query])

    hit_lengths = [hit[1] - hit[0] + 1 for hit in hits]
    coverage = sum(hit_lengths) / genome['length']

    return coverage


def genome_hits(col, value, genome):
    '''Calculate the number of hits the sample has on the genome'''

    return db.mapped.find({col: value, "genome": genome},
                          {"_id":0, "refStart":1, "refEnd":1})


def get_common_stats(project, genome_hits_query, gene_hits_query, genome):
    '''Build dictionary containing common coverage statistiscs'''

    # Total number of hits on the genome
    genome_hit_count = genome_hits_query.count()

    # Don't bother continuing if there are no hits
    if not genome_hit_count:
        return None

    # Get the coverage of mapped alignments in the genome for that sample
    genome_coverage_percent = genome_coverage(genome_hits_query, genome)

    # Total number of hits on just the genes of the genome
    gene_hit_count = gene_hits_query.count()

    # The average/max of the mean coverage on each gene
    gene_coverage_avg, gene_coverage_max = gene_coverage(gene_hits_query, genome)

    stats = {
        "accession": genome['accession']
        ,  "genome": genome['name']
        ,  "gi": genome["gi"]
        ,  "projectLabel": project['label']
        ,  "project": project['name']
        ,  "projectId": project['_id']
        ,  "genomeHits": genome_hit_count
        ,  "geneHits": gene_hit_count
        ,  "genomeCoverage": genome_coverage_percent
        ,  "geneCoverageAvg": gene_coverage_avg
        ,  "geneCoverageMax": gene_coverage_max
        ,  "tags": []
        }

    if "left" in genome:
        stats["left"] = genome["left"]

    return stats


def filter_stats(stats):
    '''
    Adds filtering tags to the statistics object. These are in a filters array element,
    which can be used in indexes for performance
    '''

    logger.info('Calculating filter: {0}'.format(stats))
    phagePattern = re.compile('phage', re.IGNORECASE)

    logger.info('Calculating filter 1: {0}'.format(stats))
    if stats["geneCoverageMax"] < 0.5:
        stats["tags"].append("lowMaxCover")
    logger.info('Calculating filter 2: {0}'.format(stats))
    if phagePattern.match(stats["genome"]):
        stats["tags"].append("phage")

    return stats


def build_project_stats(project, genome):
    '''Build dictionary containing all project converage statistiscs'''

    genome_hits_query = genome_hits('projectId', project['_id'], genome['gi'])
    gene_hits_query = gene_hits('projectId', project['_id'], genome['gi'])

    stats = get_common_stats(project, genome_hits_query, gene_hits_query, genome)
    if not stats:
        return None

    # Replaces the calculated gene_coverage_max from above with the max coverage of all samples
    sample_coverage = db.statistics.find({'projectId': project['_id'], 'accession': genome['accession'], 'sampleId': {'$exists': 1}}, {'_id':0, 'geneCoverageMax':1})
    gene_coverage_max = max([x['geneCoverageMax'] for x in sample_coverage])

    stats["ownerType"] = "project"
    stats["ownerId"] = project['_id']
    stats["geneCoverageMax"] = gene_coverage_max
    return filter_stats(stats)


def build_sample_stats(project, sample, genome):
    '''Build dictionary containing all sample converage statistiscs'''

    genome_hits_query = genome_hits('sampleId', sample['_id'], genome['gi'])
    gene_hits_query = gene_hits('sampleId', sample['_id'], genome['gi'])

    stats = get_common_stats(project, genome_hits_query, gene_hits_query, genome)
    if not stats:
        return None

    stats["ownerType"] = "sample"
    stats["ownerId"] = sample['_id']
    stats["sample"] = sample['name']
    stats["sampleId"] = sample['_id']
    return filter_stats(stats)


def build_alignment_stats(project, alignment, genome):
    '''Build dictionary containing all alignment coverage statistiscs'''

    genome_hits_query = genome_hits('alignmentId', alignment['_id'], genome['gi'])
    gene_hits_query = gene_hits('alignmentId', alignment['_id'], genome['gi'])

    stats = get_common_stats(project, genome_hits_query, gene_hits_query, genome)
    if not stats:
        return None

    stats["ownerType"] = "alignment"
    stats["ownerId"] = alignment['_id']
    stats["sample"] = alignment['sample']
    stats["sampleId"] = alignment['sampleId']
    stats["alignment"] = alignment['name']
    stats["alignmentId"] = alignment['_id']
    return filter_stats(stats)


def project_statistics(project):
    '''Calculate the statistics for a project'''
    logger.debug('Calculating statistics for project: {0}'.format(project['label']))

    genomes = db.genome.find({}, timeout=False)

    project_stats = (build_project_stats(project, genome) for genome in genomes)

    #map(insert_stats, filter(None, project_stats))
    [insert_stats(stat) for stat in project_stats if stat]


def sample_statistics(sample, project):
    '''Calculate the statistics for a sample'''
    logger.debug('Calculating statistics for sample: {0}'.format(sample['name']))

    genomes = db.genome.find({}, timeout=False);
    sample_stats = (build_sample_stats(project, sample, genome) for genome in genomes)

    #map(insert_stats, filter(None, sample_stats))
    [insert_stats(stat) for stat in sample_stats if stat]



def alignment_statistics(alignment, project):
    '''Calculate the statistics for an alignment'''
    logger.debug('Calculating statistics for alignment: {0}'.format(alignment['name']))

    genomes = db.genome.find({}, timeout=False);
    alignment_stats = (build_alignment_stats(project, alignment, genome) for genome in genomes)

    #map(insert_stats, filter(None, sample_stats))
    [insert_stats(stat) for stat in alignment_stats if stat]



def generate_statistics(project):
    '''Generates the statistics for the project and all samples under it'''
    logger.info('Calculating statistics for project: {0}'.format(project['name']))

    logger.debug('Remove old statistics for the project: {0}'.format(project['label']))
    db.statistics.remove({'projectId': project['_id']})

    samples = db.sample.find({"projectId": project['_id']})
    logger.info("Found samples: {0}".format(samples))

    alignments = db.alignment.find({"projectId": project['_id']})
    logger.info("Found alignments: {0}".format(alignments))

    pool_size = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(pool_size)

    p_statistics = partial(sample_statistics, project=project)
    a_statistics = partial(alignment_statistics, project=project)

    pool.map(p_statistics, samples)
    pool.map(a_statistics, alignments)
    #map(p_statistics, samples)

    project_statistics(project)


def update_sample_count(genome):
    ''' '''
    s = db.mapped.find({'genome': genome['gi']}).distinct('sampleId')
    db.genome.update({'gi': genome['gi']}, {'$set': {'samples': s, 'sampleCount': len(s)}})


def main(args):
    '''Calculate Genome Coverage Statistics'''

    global db, logger

    logger = args.logging.getLogger(__name__)
    db = connect(args)

    projects = list(db.project.find({'label': {'$in': args.projects}}))
    logger.info("Found projects: {0}".format(projects))

    map(generate_statistics, projects)

    # Updating Genomes with the number of sample hits
    logger.info('Updating Genome collection to show which samples hit the genome...')
    genomes = db.genome.find({}, {'_id': 0, 'gi': 1})
    map(update_sample_count, genomes)
    logger.info('Done.')

if __name__ == '__main__':
    print 'This program should be run as part of the capsid package:\n\t$ capsid statistics -h\n\tor\n\t$ /path/to/capsid/bin/capsid statistics -h'
