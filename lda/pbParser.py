#!/usr/bin/env python
"""
to be called via multiprocessing

"""

import csv, os,re,sys,shutil,cPickle,getopt,gc,time
from multiprocessing import Pool, cpu_count
from vertebratesLib import *


## define function for multiprocessing
def get_counts(args):
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


if __name__ == "__main__":
    ## read in input file
    if len(sys.argv) < 3:
        print sys.argv[0] + " -s split -d dataDir"
        sys.exit()
    try:
        optlist, args = getopt.getopt(sys.argv[1:], 's:f:l:d:')
    except getopt.GetoptError:
        raise Exception(sys.argv[0] + "-s split -f first -l last -d dataDir")

    split,first,last,dataDir = None,None,None,None
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
    outputFile = os.path.join(outputDir,"branches.npz")
    print 'loading', outputFile
    npz = np.load(outputFile)
    mat = npz['matrix']
    branches = npz['rows']

    ## identify positions that are in scope and run
    positions = get_positions(split,dataDir)
    positionsToRun = positions[first:last]
    print("...running via multiprocessing")
    po = Pool(processes=cpu_count()-1)
    _results = po.map_async(get_counts,[(split,p) for p in positionsToRun])
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
