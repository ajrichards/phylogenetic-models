#!/usr/bin/env python
"""
Parse the phylobayes outputted newick trees into a matrices of position specific transitions

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
    leaves = np.array([leaf for  leaf in fixedTree.get_leaves()])[[10,20,30,40,50,60,70,80,90]]
   
    ## loop through the files
    outputFile = os.path.join(outputDir,"%s-positions.npz"%(split))
    positions = get_positions(split,dataDir)

    summaryTree = np.zeros((positions.size,TRANSITIONS.size),)
    summarySpecies = {}
    for leaf in leaves:
        actualLeaf = re.sub("_+[A-Z]$","",leaf.name)
        summarySpecies[actualLeaf] = np.zeros((positions.size,TRANSITIONS.size),)

    print("...saving empty matrix")
    args = summarySpecies
    args['tree'] = summaryTree
    args['rows'] = positions
    args['columns'] = TRANSITIONS
    np.savez_compressed(outputFile,**args)
    print("done.")

## run via subprocessing -- use tail -f file.log
timeStart = time.time()
for split in SPLITS:
    print(split)
    positions = get_positions(split,dataDir)
    create_empty_matrix(split)
    chunks = 2
    n = len(positions)
    chunkSize = int(round(float(n)/float(chunks)))
    startPoints = np.arange(0,n,chunkSize)
    blocks = [(sp,sp+chunkSize) for sp in startPoints]

    for first,last in blocks:  
        print("...running via multiprocessing")
        cmd = "python pbParser.py -s %s -f %s -l %s -m positions -d %s"%(split,first,last,dataDir)
        print '...%s'%cmd
        p = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE, bufsize=1)
        with p.stdout:
            for line in iter(p.stdout.readline, b''):
                print line,
        p.wait() 

print("end: %s"%time.strftime('%H:%M:%S',time.gmtime(time.time()-timeStart)))
