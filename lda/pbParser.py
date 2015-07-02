#!/usr/bin/env python
"""
to be called via multiprocessing

"""

import csv, os,re,sys,shutil,cPickle,getopt,gc,time
from multiprocessing import Pool, cpu_count
from vertebratesLib import *


## define function for multiprocessing
def get_counts_branches(args):
    """
    function called by multiprocessing
    """
    
    split,position = args
    treeList = get_trees(split,position,dataDir)
    branchSamples = {}
    for branch in branches:
        branchSamples[branch] = np.zeros((len(treeList),TRANSITIONS.size),dtype=int)
    
    ## loop through all samples of the tree
    for t,pbTree in enumerate(treeList):
        fixedTree,treeSummary = fix_tree(pbTree)
        
        for branch in branches:
            branchPairs = treeSummary[branch.split("#")[0]]['pairs']
            branchCounts = transitions_to_counts(branchPairs)
            branchSamples[branch][t,:] = branchCounts

    ## get the average for each sample matrix
    branchResults = {}
    for branch in branches:
        branchResults[branch] = np.array([int(round(i)) for i in branchSamples[branch].mean(axis=0)])
                    
    return branchResults


def get_counts_positions(args):
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
        print sys.argv[0] + " -s split -d dataDir"
        sys.exit()
    try:
        optlist, args = getopt.getopt(sys.argv[1:], 's:f:l:m:d:')
    except getopt.GetoptError:
        raise Exception(sys.argv[0] + "-s split -f first -l last -m matType -d dataDir")

    split,first,last,dataDir,matType = None,None,None,None,None
    for o,a in optlist:
        if o == '-f':
            first = int(a)
        if o == '-l':
            last  = int(a)
        if o == '-s':
            split = a
        if o == '-d':
            dataDir = a
        if o == '-m':
            matType = a
            
    if matType not in ['branches','positions']:
        raise Exception("Invalid matrix type specified")

    ## identify positions to run
    positions = get_positions(split,dataDir)
    positionsToRun = positions[first:last]
    print("...running via multiprocessing")

    if matType == 'branches':
        ## load the matrices
        print("extrating transitions for %s"%split)
        outputDir = os.path.join("..","data","hv-compressed")
        outputFile = os.path.join(outputDir,"branches.npz")
        print 'loading', outputFile
        npz = np.load(outputFile)
        mat = npz['matrix']
        branches = npz['rows']

        ## identify positions that are in scope and run
        po = Pool(processes=cpu_count()-1)
        _results = po.map_async(get_counts_branches,[(split,p) for p in positionsToRun])
        print("...getting multiprocessing results")
        results =  _results.get()

        print("...appending results to matrix")
        for p,pos in enumerate(positionsToRun):
            branchResults = results[p]
            if p % 200 == 0:
                print("%s/%s"%(p,len(positionsToRun)))
            for branch,branchMeans in branchResults.iteritems():
                branchInd = np.where(branches==branch)[0][0]
                mat[branchInd,:] += branchMeans
                
        ## save it
        print("...saving")
        args = {}
        args['matrix'] = mat
        args['rows'] = branches
        args['columns'] = TRANSITIONS
        np.savez_compressed(outputFile,**args)
        print("done.")

    elif matType == 'positions':
        ## load the matrices
        print("extrating transitions for %s"%split)
        outputDir = os.path.join("..","data","hv-compressed")
        outputFile = os.path.join(outputDir,"%s-positions.npz"%(split))
        
        npz = np.load(outputFile)
        summaryTree = npz['tree']
        summarySpecies = {}
        for key,item in npz.iteritems():
            summarySpecies[key] = item

        ## multiprocess it
        po = Pool(processes=cpu_count()-1)
        _results = po.map_async(get_counts_positions,positionsToRun)
        print("...getting multiprocessing results")
        results =  _results.get()

        ## save counts to loaded matrices
        print("...updating results")
        toRun = np.arange(first,last)
        for indx, result in enumerate(results):
            p = toRun[indx]
            speciesCounts,treeCounts = result
            summaryTree[p,:] = treeCounts
            for key,item in speciesCounts.iteritems():
                summarySpecies[key][p,:] = item

        ## save it
        args = summarySpecies
        args['tree'] = summaryTree
    
        print("...saving")
        np.savez_compressed(outputFile,**args)
        print("done.")
