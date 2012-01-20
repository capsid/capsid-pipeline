#!/usr/bin/env python

# Copyright 2011(c) The Ontario Institute for Cancer Reserach. All rights reserved.
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

def build_project_stats(project, genome):
    '''Build dictionary containing all project converage statistiscs'''

    genome_hits_query = genome_hits('project', project['label'], genome['gi'])
    # Total number of hits on the genome
    genome_hit_count = genome_hits_query.count()
 
    # Don't bother continuing if there are no hits
    if not genome_hit_count:
        return None   

    # Get the coverage of mapped alignments in the genome for that sample
    genome_coverage_percent = genome_coverage(genome_hits_query, genome)

    gene_hits_query = gene_hits('project', project['label'], genome['gi'])
    # Total number of hits on just the genes of the genome
    gene_hit_count = gene_hits_query.count()
    # The average/max of the mean coverage on each gene
    gene_coverage_avg, gene_coverage_max = gene_coverage(gene_hits_query, genome)

    stats = {
        "accession": genome['accession']
        ,  "genome": genome['name']
        ,  "label": project['label']
        ,  "project": project['name']
        ,  "genomeHits": genome_hit_count
        ,  "geneHits": gene_hit_count
        ,  "genomeCoverage": genome_coverage_percent
        ,  "geneCoverageAvg": gene_coverage_avg
        ,  "geneCoverageMax": gene_coverage_max
        }

    return stats


def build_sample_stats(project, sample, genome):
    '''Build dictionary containing all sample converage statistiscs'''

    genome_hits_query = genome_hits('sample', sample['name'], genome['gi'])
    # Total number of hits on the genome
    genome_hit_count = genome_hits_query.count()

    # Don't bother continuing if there are no hits
    if not genome_hit_count:
        return None

    # Get the coverage of mapped alignments in the genome for that sample
    genome_coverage_percent = genome_coverage(genome_hits_query, genome)

    gene_hits_query = gene_hits('sample', sample['name'], genome['gi'])
    # Total number of hits on just the genes of the genome
    gene_hit_count = gene_hits_query.count()
    # The average/max of the mean coverage on each gene
    gene_coverage_avg, gene_coverage_max = gene_coverage(gene_hits_query, genome)

    stats = {
        "accession": genome['accession']
        ,  "genome": genome['name']
        ,  "label": project['label']
        ,  "project": project['name']
        ,  "sample": sample['name']
        ,  "genomeHits": genome_hit_count
        ,  "geneHits": gene_hit_count
        ,  "genomeCoverage": genome_coverage_percent
        ,  "geneCoverageAvg": gene_coverage_avg
        ,  "geneCoverageMax": gene_coverage_max
        }

    return stats


def project_statistics(project):
    '''Calculate the statistics for a project'''
    logger.info('Calculating statistics for project: {0}'.format(project['name']))

    genomes = db.genome.find()
    
    project_stats = (build_project_stats(project, genome) for genome in genomes)
    map(insert_stats, filter(None, project_stats))


def sample_statistics(sample, project):
    '''Calculate the statistics for a sample'''
    logger.debug('Calculating statistics for sample: {0}'.format(sample['name']))

    genomes = db.genome.find()
    sample_stats = (build_sample_stats(project, sample, genome) for genome in genomes)
    map(insert_stats, filter(None, sample_stats))


def generate_statistics(project):
    '''Generates the statistics for the project and all samples under it'''

    samples = db.sample.find({"project": project['label']})

    project_statistics(project)

    pool_size = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(pool_size)
    
    p_statistics = partial(sample_statistics, project=project)
    
    pool.map(p_statistics, samples)
    

def main(args):
    '''Calculate Genome Coverage Statistics'''

    global db, logger

    logger = args.logging.getLogger(__name__)
    db = connect(args)

    projects = db.project.find({'label': {'$in': args.projects}})

    map(generate_statistics, projects)

    # Updating Genomes with the number of sample hits
    logger.info('Updating Genome collection to show which samples hit the genome...')
    db.system_js.gs()
    logger.info('Done.')

if __name__ == '__main__':
    print 'This program should be run as part of the capsid package:\n\t$ capsid statistics -h\n\tor\n\t$ /path/to/capsid/bin/capsid statistics -h'
