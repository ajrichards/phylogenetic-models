#!/usr/bin/env python
"""
get data into a format that LDA understands
https://pythonhosted.org/ete2/tutorial/tutorial_trees.html#reading-newick-trees

the data are organized into splits so to keep things more managable 
we first move individual files into corresponding splits

extract the tar file into 'dataDir' and run this script

NOTES: 

    (1) in this example the tree is contstrained
    (2) the total positions are broken into N 'splits'
    (3) each split contains 30,000 files
    (4) so that means that there are 30,000 * N total positions

"""

import csv, os,re,sys,shutil,cPickle,getopt,gc,time
from multiprocessing import Pool, cpu_count
from vertebratesLib import *

def get_counts(args):
    """
    function called by multiprocessing
    """
    
    position = args    
    treeList = get_trees(split,position,dataDir)
    treeMatrix = np.zeros((len(treeList),TRANSITIONS.size),dtype=int)
    speciesMatrix = {}
    fixedTree,treeSummary = fix_tree(treeList[0])
    leaves = np.array([leaf for  leaf in fixedTree.get_leaves()])[[10,20,30,60,70,80]]
             
    for leaf in leaves:
        actualLeaf = re.sub("_+[A-Z]$","",leaf.name)
        speciesMatrix[actualLeaf] = np.zeros((len(treeList),TRANSITIONS.size),dtype=int)

    ## loop through all samples of the tree
    for t,pbTree in enumerate(treeList):
        fixedTree,treeSummary = fix_tree(pbTree)
        treePairs = []
        for item in treeSummary.itervalues():
            treePairs.extend(item['pairs'])
        treeCounts = transitions_to_counts(treePairs)
        treeMatrix[t,:] = treeCounts
        
        for leaf in leaves:
            actualLeaf = re.sub("_+[A-Z]$","",leaf.name)
            leafPairs = get_species_transitions(treeSummary,leaf)
            leafCounts = transitions_to_counts(leafPairs)
            speciesMatrix[actualLeaf][t,:] = leafCounts 

    ## update the matrix
    treeCounts = np.array([int(round(i)) for i in treeMatrix.mean(axis=0)])
    speciesCounts = {}
    for leaf in leaves:
        actualLeaf = re.sub("_+[A-Z]$","",leaf.name)
        speciesCounts[actualLeaf] = np.array([int(round(i)) for i in speciesMatrix[actualLeaf].mean(axis=0)])
    return speciesCounts,treeCounts

if __name__ == "__main__":
    ## read in input file
    if len(sys.argv) < 3:
        print sys.argv[0] + " -s split -f first -l last"
        sys.exit()
    try:
        optlist, args = getopt.getopt(sys.argv[1:], 'f:l:s:d:')
    except getopt.GetoptError:
        raise Exception(sys.argv[0] + "-f first -l last -s split -d dataDir")

    first,last,split = None,None,None
    for o,a in optlist:
        if o == '-f':
            first = int(a)
        if o == '-l':
            last  = int(a)
        if o == '-s':
            split = a
        if o == '-d':
            dataDir = a
            
    ## load the matrices
    print("extrating transitions for %s"%split)
    outputDir = os.path.join("..","data","hv-compressed")
    outputFile = os.path.join(outputDir,"%s.npz"%(split))
    npz = np.load(outputFile)
    summaryTree = npz['tree']
    summarySpecies = {}
    for key,item in npz.iteritems():
        summarySpecies[key] = item

    ## identify positions that are in scope and run
    positions = get_positions(split,dataDir)
    positionsToRun = positions[first:last]

    """
    toRun = range(first,last)
    for p,position in enumerate(positions):
        if p not in toRun:
            continue
        if p%100 == 0:
            print("%s/%s"%(p,len(positions)))
        
        speciesCounts,treeCounts = get_counts(p)
        summaryTree[p,:] = treeCounts
        for key,item in speciesCounts.iteritems():
            summarySpecies[key][p,:] = item
    """

    print("...running via multiprocessing")
    po = Pool(processes=cpu_count()-1)
    _results = po.map_async(get_counts,positionsToRun)
    print("...getting multiprocessing results")
    results =  _results.get()

    print("...updating results")
    ## save counts to loaded matrices
    toRun = np.arange(first,last)
    for indx, result in enumerate(results):
        p = toRun[indx]
        speciesCounts,treeCounts = result
        summaryTree[p,:] = treeCounts
        for key,item in speciesCounts.iteritems():
            summarySpecies[key][p,:] = item

    ##
    #toRun = np.arange(first,last)
    #for indx, result in enumerate(results):
    #    p = toRun[indx]
    #    speciesCounts,treeCounts = result
    #    npz['tree'][p,:] = treeCounts
    #    for key,item in speciesCounts.iteritems():
    #        npz[key][p,:] = item

    ## save it
    args = summarySpecies
    args['tree'] = summaryTree
    
    print("...saving")
    np.savez_compressed(outputFile,**args)
    print("done.")
