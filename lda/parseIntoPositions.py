#!/usr/bin/env python
"""
Parse the phylobayes outputted newick trees into a matrix of branch specific transitions

"""

import csv, os,re,sys,shutil,getopt,time
import subprocess
from vertebratesLib import *
        
## check that the files are extracted
dataDir = None
for ddir in [os.path.join("..","data","herve-vertebrates"),\
             os.path.join("/","media","ganda","mojo","phylogenetic-models","herve-vertebrates")]:
    if os.path.isdir(ddir):
        dataDir = ddir
    print ddir

if not dataDir:
    raise Exception("data directory needs to be set")
print("setting data dir as \n...%s"%(dataDir))

## if necessary move the files into subdirectories accoriding to split
ensure_files_properly_located(dataDir)

outputDir = os.path.join("..","data","hv-compressed")
if not os.path.isdir(outputDir):
    os.mkdir(outputDir)

def create_empty_matrix(split):
    ## get the leaves 
    treeList = get_trees(split,'1000',dataDir)
    fixedTree,treeSummary = fix_tree(treeList[0])
    leaves = np.array([leaf for  leaf in fixedTree.get_leaves()])[[10,20,30,60,70,80]]
   
    ## loop through the files
    outputFile1 = os.path.join(outputDir,"%s.npz"%(split))
    positions = get_positions(split,dataDir)

    summaryTree = np.zeros((positions.size,TRANSITIONS.size),)
    summarySpecies = {}
    for leaf in leaves:
        actualLeaf = re.sub("_+[A-Z]$","",leaf.name)
        summarySpecies[actualLeaf] = np.zeros((positions.size,TRANSITIONS.size),)

    print("...saving empty matrix")
    args = summarySpecies
    args['tree'] = summaryTree
    np.savez_compressed(outputFile1,**args)
    print("done.")

## run via subprocessing -- use tail -f file.log
timeStart = time.time()
for split in ["SPLIT%s"%i for i in range(1,11)]: 
    ## save the split identifiers
    splitPositions = get_positions(split,dataDir)
    np.savez(os.path.join(outputDir,"%s-rows.npz"%(split)),rows=splitPositions)
    
    #create_empty_matrix(split)
    #logFile = open("%s.log"%split,'w')
    #blocks = ((0,5000),(5000,10000),(10000,15000),(15000,20000),(20000,25000),(25000,30000))
    #for first,last in blocks:
    #    cmd = "python dataMungeVertebrates.py -s %s -f %s -l %s -d %s"%(split,first,last,dataDir)
    #    print 'running...\n%s'%cmd
    #    proc = subprocess.Popen(cmd,shell=True,stdout=logFile,stdin=subprocess.PIPE)
    #    proc.communicate()

print("end: %s"%time.strftime('%H:%M:%S',time.gmtime(time.time()-timeStart)))
