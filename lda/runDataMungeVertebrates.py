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

import csv, os,re,sys,shutil,cPickle,getopt,time
import subprocess
from htsint import run_subprocess
from vertebratesLib import *
        
## check that the files are extracted
splits = ["SPLIT%s"%i for i in range(1,11)] 
dataDir = None
for ddir in [os.path.join("..","data","herve-vertebrates"),\
             os.path.join("/","media","ganda","mojo","phylogenetic-models","herve-vertebrates")]:
    if os.path.isdir(ddir):
        dataDir = ddir
    print ddir

if not dataDir:
    raise Exception("data directory needs to be set")
print("setting data dir as \n...%s"%(dataDir))

## move the files into subdirectories accoriding to split
for _split in splits:
    splitPath = os.path.join(dataDir,_split)
    if not os.path.isdir(splitPath):
        os.mkdir(splitPath)
fileCount = 0
for fileName in os.listdir(dataDir):
    filePath = os.path.join(dataDir,fileName)
    if not re.search("\.map$",fileName):
        continue

    _split = re.findall("SPLIT\d+",fileName)[0]
    splitPath = os.path.join(dataDir,_split)
    newFilePath = os.path.join(splitPath,fileName)
    shutil.move(filePath,newFilePath)
    fileCount += 1

if fileCount != 0:
    print('Moved %s files...'%fileCount)

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
#split = "SPLIT2"
for split in ["SPLIT%s"%i for i in range(3,11)]: 
    create_empty_matrix(split)
    logFile = open("%s.log"%split,'w')
    blocks = ((0,5000),(5000,10000),(10000,15000),(15000,20000),(20000,25000),(25000,30000))
    for first,last in blocks:
        cmd = "python dataMungeVertebrates.py -s %s -f %s -l %s -d %s"%(split,first,last,dataDir)
        print 'running...\n%s'%cmd
        proc = subprocess.Popen(cmd,shell=True,stdout=logFile,stdin=subprocess.PIPE)
        proc.communicate()

print("end: %s"%time.strftime('%H:%M:%S',time.gmtime(time.time()-timeStart)))
