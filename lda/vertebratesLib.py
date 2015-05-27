#!/usr/bin/env python
"""
get data into a format that LDA understands
https://pythonhosted.org/ete2/tutorial/tutorial_trees.html#reading-newick-trees

the data are organized into splits so to keep things more managable 
we first move individual files into corresponding splits

extract the tar file into 'dataDir' and run this script

NOTES: 

    (1) in this example the tree is contstrained


"""

import csv, os,re,sys,shutil
from ete2 import Tree

## check that the files are extracted
dataDir = os.path.join("..","data","herve-vertebrates")
splits = ["SPLIT%s"%i for i in range(1,11)] 

## convenience functions
def get_tree_file_path(split,position):
    """
    return the map file path for a split and a position
    """

    if int(position) < 0 or int(position) > 29999:
        raise Exception("invalid numId arg %s"%position)
    if split not in splits:
        raise Exception("invalid split arg %s"%split)

    filePath = os.path.join(dataDir,split,"Amphibia-noo-100x1285039-%s-CATG_%s.map"%(split,position))
    if not os.path.exists(filePath):
        raise Exception("cannot find '%s'"%(filePath))
    return filePath

def get_trees(split,position):
    """
    return the trees as a list for a given split and position
    """

    filePath = get_tree_file_path(split,position)
    fid = open(filePath,'r')
    reader = csv.reader(fid,delimiter=';')
    treeList = []
    for linja in reader:
        if len(linja) == 0:
            continue
        treeList.append(linja[0]+";")
    return treeList
