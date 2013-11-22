#!/usr/bin/env python
'''Functions for dealing with the NCBI taxonomy information and loading them into MongoDB'''


# Copyright 2011(c) The Ontario Institute for Cancer Research. All rights reserved.
#
# This program and the accompanying materials are made available under the
# terms of the GNU Public License v3.0.
#
# You should have received a copy of the GNU General Public License along with
# this program.  If not, see <http://www.gnu.org/licenses/>.

import pymongo 
import os.path
import subprocess

from database import *

db = None
fs = None
logger = None

def main(args):
    '''
    Reads through NCBI taxonomy data files and loads data into MongoDB in the taxa collection.
    '''

    global db, fs, logger

    logger = args.logging.getLogger(__name__)
    db = connect(args)

    load_taxonomy(args.directory, args.repair)


def load_taxonomy_names(name_file):  
    '''
    Loads NCBI taxonomy names into a couple of dicts that are returned, 
    as a scientific name collection and a common name collection.
    '''

    logger.info("Loading taxonomy names: {0}".format(name_file))

    NamesFile = open(name_file)

    SciNameDict  = {}
    ComNameDict  = {}
    #AllNamesDict = {}
    for line in NamesFile:
        taxID, Name, Uname, NameClass = line.split('\t|\t')
        taxID = int(taxID)
        if NameClass[0:3] == 'sci': # Process scientific name
            if Name == 'environmental samples':
                # recode as ENV-xxx with xxx being phylogrouping
                # example: environmental samples <Bacillariophyta> becomes ENV-Bacillariophyta
                groups = Uname.split('<')
                if len(groups) == 2:
                    Name = 'ENV-' + groups[1][:-1]
                else:
                    Name = 'ENV-Undefined'
            SciNameDict[taxID] = Name
        elif NameClass[0:11] == 'genbank com': # Process Genbank common name
            ComNameDict[taxID] = Name  
            #AllNamesDict[taxID] = AllNamesDict.get(taxID, []) + [Name]

    return SciNameDict, ComNameDict


def load_taxonomy_data(name_file, node_file):
    '''
    Loads NCBI taxonomy node data into MongoDB in the taxa collection.
    '''

    SciNameDict, ComNameDict = load_taxonomy_names(name_file)

    logger.info("Loading taxonomy nodes: {0}".format(node_file))

    NodesFile = open(node_file)

    for line in NodesFile:
        fields  = line.split('\t|\t')
        taxID   = int(fields[0])
        if int(fields[1]) == 1 and fields[2] == 'no rank':
            first = 0
            #first = None
        else:
            first = int(fields[1])
        SciName = SciNameDict.get(taxID, "<%d>" % taxID)
        NCBIdata = {
            "_id": taxID,  
            "parent"  : first,  # Parent taxID
            "rank"    : fields[2], # Taxonomic rank name
            "sciName" : SciName,                # Scientific name
            "comName" : ComNameDict.get(taxID, SciName) # Genbank common name
            #"allNames": AllNamesDict.get(taxID, SciName)   # List of all names
            #"GenCode" : int(fields[6]),            # Genomic genetic code
            #"MitCode" : int(fields[8])         # Mitochondrial code
            }
        db.taxa.update({'_id': taxID},NCBIdata,True)

    logger.info("Adding parent index")
    db.taxa.ensure_index('parent')

    logger.info("Calculating nested set values")
    update_tree(db.taxa)

    logger.info("Adding left index")
    db.taxa.ensure_index('left')


def load_taxonomy_genomes(genome_file):

    taxon_viral_full = db.taxa.find({ 'rank': 'superkingdom', 'sciName': 'Viruses'})

    viral = taxon_viral_full[0]
    r = taxon_viral_full[0]['right']
    l = taxon_viral_full[0]['left']

    logger.info('Retrieving all the viral taxa')

    viral_taxa_all = []
    viral_taxa_dict = {}

    taxon_ids = db.taxa.find( {'left': { '$gt': l }, 'right': { '$lt': r }}, { '_id': 1} )
    tax_array = []
    for taxid in taxon_ids:
        tax_array.append(taxid['_id'])
        viral_taxa_dict[taxid['_id']] = True
    viral_taxa_all = viral_taxa_all + tax_array

    logger.info("Total number of viral taxa: {0}".format(len(viral_taxa_all)))

    logger.info("Loading viral genome identifiers")

    # We could sort the file, but it's actually a bad idea. That's because the sort 
    # works on the whole file, which is huge. We are only interested in a subset of
    # genomes. Also, we switched to use a dict to test for the match, which is 
    # much faster in this loop. 
    gi_file = open(genome_file)

    db.gitaxid.remove()

    list_gi = []
    current = None

    for line in gi_file:
        elements = line.split('\t')
        gi  = int(elements[0])
        taxid = int(elements[1])

        if taxid != 0 and (taxid in viral_taxa_dict):
            if current == None:
                current = taxid

            if taxid == current:
                list_gi.append(gi)   
            else:
                db.gitaxid.update({'_id': current}, {'$addToSet': { 'gi': {'$each': list_gi }}}, True)
                current = taxid
                list_gi = [gi] 

    db.gitaxid.update({'_id': current}, {'$addToSet': { 'gi': {'$each': list_gi }}}, True)

    logger.info("Finished Loading viral genome identifiers")


def update_genomes():
    '''
    Adds the preorder tree traversal values to the genomes
    '''

    logger.info("Calculating nested set values for genomes")
    genomes = db.gitaxid.find({}, { '_id': 1, 'gi': 1})
    for genome in genomes:
        taxon = db.taxa.find_one({'_id': genome['_id']}, {'left': 1})
        db.genome.update({'gi': {'$in': genome['gi']}}, {'$set': {'left': taxon['left']}}, multi=True)


def update_tree(collection): 
    '''
    Calculates the modified preorder tree traversal for the taxa collection.
    '''
    update_node(collection, 0, 1)


def update_node(collection, node_id, left):
    '''
    Calculates the modified preorder tree traversal for a node in the taxa collection.
    '''
    right = left + 1
    cursor = collection.find({ "parent": node_id }, { '_id': 1})
    for item in cursor:
        right = update_node(collection, item['_id'], right) + 1
    collection.update({'_id': node_id}, {'$set' : { "left" : left, "right" : right}}) 
    return right


def load_taxonomy(directory, repair):
    '''
    Loads NCBI taxonomy data files into MongoDB in the taxa collection.
    '''

    name_file = os.path.join(directory, 'names.dmp')
    node_file = os.path.join(directory, 'nodes.dmp')

    load_taxonomy_data(name_file, node_file)

    genome_file = os.path.join(directory, 'gi_taxid_nucl.dmp')
    load_taxonomy_genomes(genome_file)

    update_genomes()


if __name__ == '__main__':
    print 'This program should be run as part of the capsid package:\n\t$ capsid taxonomy -h\n\tor\n\t$ /path/to/capsid/bin/capsid taxonomy -h'
