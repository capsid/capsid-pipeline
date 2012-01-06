#!/usr/bin/env Python


# Copyright 2011(c) The Ontario Institute for Cancer Reserach. All rights reserved.
#
# This program and the accompanying materials are made available under the
# terms of the GNU Public License v3.0.
#
# You should have received a copy of the GNU General Public License along with
# this program.  If not, see <http://www.gnu.org/licenses/>.


#TODO Threshold
from __future__ import division
import math

from bson import ObjectId, errors
from bx.intervals.intersection import Intersecter, Interval

import capsid


db, logger = None, None


def merge(array):
    '''
    Accepts a list containing lists or tuples with a size of two(2) and merges the overlapping ranges. Returns a generator object.

        >>>L = [(1,15), (10,20), (25,30)]
        >>>merge(L)
        <generator object merge at 0x26e7230>
        >>>list(merge(L))
        [(1, 20), (25, 30)]
    '''
    try:
        saved = [min(list(array))[0], min(list(array))[0]]
    except ValueError:
        return

    for start, end in sorted([sorted(elem) for elem in array]):
        if start <= saved[1]:
            saved[1] = max(saved[1], end)
        else:
            yield tuple(saved)
            saved[0] = start
            saved[1] = end

    yield tuple(saved)


def _get_intersector(gene_hits):
    '''Build a intersector of hits. Use gene start, end to find hits on the gene'''

    intersector = Intersecter()

    for h in gene_hits:
        intersector.add_interval( Interval(h[0], h[1]) )

    return intersector


def find_genes(genome):
    '''Returns the genes of the genome.'''

    genes = db.feature.find({"type": "gene", "genomeId": genome['_id']}, {"start":1, "end":1, "name": 1})

    return list(genes)


def _hits(genome, col, value):
    #logger.debug('Entering _hits')

    # Get mapped alignments from the db
    hits_query = db.mapped.find({col: value, "genomeId": genome['_id']}, {"refStart":1, "refEnd":1})

    # Hit Counts
    hit_count = hits_query.count()

    # Merge the hits to prevent miscalculations due to overlapping alignments
    hits = merge([(h['refStart'], h['refEnd']) for h in hits_query])

    #logger.debug('Exiting _hits')
    return hits, hit_count


def find_hits_by_sample(genome, name):
    return _hits(genome, "sample", name)


def find_hits_by_project(genome, label):
    return _hits(genome, "project", label)


def _gene_hits(genome, col, value):
    #logger.debug('Entering _gene_hits')

    # Get mapped alignments that hit a gene from DB
    hits_query = db.mapped.find({col: value, "genomeId": genome['_id'], "mapsGene": 1}, {"refStart":1, "refEnd":1})

    # Hit Counts
    gene_hits_count = hits_query.count()

    # Merge the hits to prevent miscalculations due to overlapping alignments
    gene_hits = merge([(h['refStart'], h['refEnd']) for h in hits_query])

    #logger.debug('Exiting _gene_hits')
    return gene_hits, gene_hits_count


def find_gene_hits_by_sample(genome, name):
    return _gene_hits(genome, "sample", name)


def find_gene_hits_by_project(genome, label):
    return _gene_hits(genome, "project", label)


def calc_coverage(genome, hits):
    '''Calculate the percentage of bases that alignments map using the merged mapped alignments'''
    #logger.debug('Entering coverage')

    #map(lambda h: h[1] - h[0], hits['plus'])
    # This is ~50-100% faster than (i)map/lambda
    hit_lengths = [hit[1] - hit[0] + 1 for hit in hits]
    coverage = sum(hit_lengths)/genome['length']

    #logger.debug('Exiting coverage')
    return coverage


def calc_gene_coverage(genes, gene_hits):
    ''' Get the average/maximum coverage of genes '''
    #logger.debug('Entering gene_coverage')

    intersector = _get_intersector(gene_hits)
    coverage = []

    for gene in genes:
        # Finds the hits that intersect with the gene
        overlap = intersector.find(gene['start'], gene['end'])
        # Merge them to remove multiple count over the same bases
        m_overlap = merge([(o.start, o.end) for o in overlap])
        # Get the total area over the gene that is covered
        length = sum([min(gene['end'], end) - max(gene['start'], start) + 1 for start, end in m_overlap])
        # Divide by the length of the gene to get the percent coverage
        p_cov = length / (gene['end'] - gene['start'] + 1)
        coverage.append(p_cov)

    try:
        maximum = max(coverage)
    except ValueError:
        maximum = 0.0

    try:
        mean = math.fsum(coverage)/len(coverage)
    except ZeroDivisionError:
        mean = 0.0

    #logger.debug('Exiting gene_coverage')
    return mean, maximum


def get_stats(genome, samples, project):
    '''Runs the statistics of each sample and project on every Genome'''

    p_hits, p_ghits, p_avg, p_max, p_cov = 0, 0, [], [], []

    logger.debug('Calculating the statistics for {0}'.format(genome['name']))
    genes = find_genes(genome)

    # For each sample perform calculations
    for sample in samples:
        logger.debug('Using the hits from {0}'.format(sample['name']))
        # Total number of hits on the genome
        hits, hit_count = find_hits_by_sample(genome, sample['name'])
        # Total number of hits on just the genes of the genome
        gene_hits, gene_hit_count = find_gene_hits_by_sample(genome, sample['name'])
        # Get the coverage of mapped alignments in the genome for that sample
        coverage = calc_coverage(genome, hits)
        # The average of the mean coverage on each gene
        avg_gene_coverage, max_gene_coverage = calc_gene_coverage(genes, gene_hits)

        logger.debug("Stats for Sample:{0}, Genome:{1}".format(sample['name'], genome['name']))
        logger.debug("Hits: {0}".format(hit_count))
        logger.debug("Gene Hits: {0}".format(gene_hit_count))
        logger.debug("Coverage: {0:.2f}%".format(coverage*100))
        logger.debug("Average Gene Coverage: {0:.2f}%".format(avg_gene_coverage*100))
        logger.debug("Max Gene Coverage: {0:.2f}%".format(max_gene_coverage*100))

        if hit_count:
            db.statistics.insert({
                    "genomeId": genome['_id']
                ,   "genomeAccession": genome['accession']
                ,   "genomeName": genome['name']
                ,   "projectId": project['_id']
                ,   "projectLabel": project['label']
                ,   "projectName": project['name']
                ,   "sampleId": sample['_id']
                ,   "sampleName": sample['name']
                ,   "hits": hit_count
                ,   "geneHits": gene_hit_count
                ,   "totalCoverage": coverage
                ,   "geneCoverage": avg_gene_coverage
                ,   "maxCoverage": max_gene_coverage
                ,   "threshold": 20
                })

    # Total number of hits on the genome
    p_hits, p_hit_count = find_hits_by_project(genome, project['label'])
    # Total number of hits on just the genes of the genome
    p_gene_hits, p_gene_hit_count = find_gene_hits_by_project(genome, project['label'])
    # Get the coverage of mapped alignments in the genome for that sample
    p_coverage = calc_coverage(genome, p_hits)
    # The average of the mean coverage on each gene
    p_avg_gene_coverage, p_max_gene_coverage = calc_gene_coverage(genes, p_gene_hits)

    logger.debug("Stats for Project:{0}, Genome:{1}".format(project['name'], genome['name']))
    logger.debug("Hits: {0}".format(p_hit_count))
    logger.debug("Gene Hits: {0}".format(p_gene_hit_count))
    logger.debug("Coverage: {0:.2f}%".format(p_coverage*100))
    logger.debug("Average Gene Coverage: {0:.2f}%".format(p_avg_gene_coverage*100))
    logger.debug("Max Gene Coverage: {0:.2f}%".format(p_max_gene_coverage*100))

    if p_hit_count:
        db.statistics.insert({
                "genomeId": genome['_id']
            ,   "genomeAccession": genome['accession']
            ,   "genomeName": genome['name']
            ,   "projectId": project['_id']
            ,   "projectLabel": project['label']
            ,   "projectName": project['name']
            ,   "sampleId": 0
            ,   "hits": p_hit_count
            ,   "geneHits": p_gene_hit_count
            ,   "totalCoverage": p_coverage
            ,   "geneCoverage": p_avg_gene_coverage
            ,   "maxCoverage": p_max_gene_coverage
            ,   "threshold": 20
            })

def main(args):
    '''Calculates the mapped alignment statistics for each Genome'''

    global db, logger

    logger = args.logging.getLogger(__name__)
    db = capsid.connect(args)

    projects = db.project.find({'label': {'$in': args.projects}})

    if projects.count():
        genomes = list(db.genome.find())

        for project in projects:
            logger.info('Calculating the statistics for {0}'.format(project['name']))
            samples = list(db.sample.find({"project": project['label']}))

            try:
                logger.debug('Looking for Multiprocessing...')
                import multiprocessing
                logger.debug('Multiprocessing found.')

                multiprocessing.log_to_stderr(args.logging.WARNING)

                pool_size = multiprocessing.cpu_count()
                pool = multiprocessing.Pool(pool_size)
                for genome in genomes:
                        pool.apply_async(get_stats, [genome, samples, project])
                pool.close()
                pool.join()
            except ImportError:
                logger.debug('Multiprocessing not found.')
                for genome in genomes:
                    get_stats(genome, samples, project)

        # Updating Genomes with the number of sample hits
        logger.info('Updating Genome collection to show which samples hit the genome...')
        db.system_js.gs()
    else:
        logger.info('No projects found using the list: {0}'.format(args.projects))

if __name__ == '__main__':
    print 'This program should be run as part of the capsid package:\n\t$ capsid statistics -h\n\tor\n\t$ /path/to/capsid/bin/capsid statistics -h'
