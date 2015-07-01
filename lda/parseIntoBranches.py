#!/usr/bin/env python
"""
Parse the phylobayes outputted newick trees into a matrix of branch specific transitions

"""

import csv, os,re,sys,shutil,getopt,time
import numpy as np
from multiprocessing import Pool, cpu_count
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

## create empty matrix
timeStart = time.time()
split = "SPLIT1"
treeList = get_trees(split,'1000',dataDir)
fixedTree,treeSummary = fix_tree(treeList[0])
branches = np.sort(np.array(["%s#%s"%(b1,b2['parent']) for b1,b2 in treeSummary.iteritems()]))

## loop through the files
outputFile = os.path.join(outputDir,"branches.npz")
mat = np.zeros((branches.size,TRANSITIONS.size),)

## save the empty matrix to file
args = {}
args['matrix'] = mat
args['rows'] = branches
args['columns'] = TRANSITIONS
np.savez_compressed(outputFile,**args)
logFile = open("branches.log",'w')

for split in SPLITS:
    print split

    positions = get_positions(split,dataDir)
    chunks = 5
    n = len(positions)
    chunkSize = int(round(float(n)/float(chunks)))
    startPoints = np.arange(0,n,chunkSize)
    blocks = [(sp,sp+chunkSize) for sp in startPoints]

    for first,last in blocks:  
        logFile = open("branches.log",'a')
        print("...running via multiprocessing")
        cmd = "python pbParser.py -s %s -f %s -l %s -d %s"%(split,first,last,dataDir)
        print '...%s'%cmd
        #proc = subprocess.Popen(cmd,shell=True,stdout=logFile,stdin=subprocess.PIPE)
        p = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE, bufsize=1)
        with p.stdout:
            for line in iter(p.stdout.readline, b''):
                print line,
        p.wait() # wait for the subprocess to exit
        #proc.communicate()
        #logFile.close()

print 'saving as a csv...'
npz = np.load(outputFile)
mat = npz['matrix']

outfid = open(os.path.join(outputDir,'branches.csv'),'w')
writer = csv.writer(outfid)
writer.writerow(["branch"] + TRANSITIONS.tolist())
for b,branch in enumerate(branches):
    writer.writerow([branch]+ [int(i) for i in mat[b,:].tolist()])
outfid.close()

print("end: %s"%time.strftime('%H:%M:%S',time.gmtime(time.time()-timeStart)))
