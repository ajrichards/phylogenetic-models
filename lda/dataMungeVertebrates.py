#!/usr/bin/env python
"""
get data into a format that LDA understands
https://pythonhosted.org/ete2/tutorial/tutorial_trees.html#reading-newick-trees

the data are organized into splits so to keep things more managable 
we first move individual files into corresponding splits

extract the tar file into 'dataDir' and run this script

NOTES: 

    (1) in this example the tree is contstrained
    (2) the total position are broken into N 'splits'
    (3) each split contains 30,000 files
    (4) so that means that there are 30,000 * N total positions

"""

import csv, os,re,sys,shutil
from ete2 import Tree,TreeStyle,PhyloTree
from ete2 import faces, AttrFace
from vertebratesLib import *

## check that the files are extracted
dataDir = os.path.join("..","data","herve-vertebrates")

## move the files into subdirectories accoriding to split
splits = ["SPLIT%s"%i for i in range(1,11)] 
for split in splits:
    splitPath = os.path.join(dataDir,split)
    if not os.path.isdir(splitPath):
        os.mkdir(splitPath)
fileCount = 0
for fileName in os.listdir(dataDir):
    filePath = os.path.join(dataDir,fileName)
    if not re.search("\.map$",fileName):
        continue

    split = re.findall("SPLIT\d+",fileName)[0]
    splitPath = os.path.join(dataDir,split)
    newFilePath = os.path.join(splitPath,fileName)
    shutil.move(filePath,newFilePath)
    fileCount += 1

if fileCount != 0:
    print('Moved %s files...'%fileCount)

## playing
treeList = get_trees("SPLIT1","25000")
nwTree = treeList[0]
t = Tree(nwTree,format=0)

## setup a tree
nwTree = treeList[0]
ptree = PhyloTree(nwTree)
rootNode = "Acipenser__R"
ptree.set_outgroup(rootNode)

## iterate through tree an give each branch a identifier
branches = {}
level = 0
branchName = 'B'+str(level)

for node in ptree.traverse("levelorder"): #levelorder | postorder
    if node.name == rootNode:
        continue
    if node.name == 'NoName':
        level += 1
        branchName = 'B'+str(level)
        node.name = branchName

    if not branches.has_key(branchName):
        branches[branchName] = []
    branches[branchName].append(node.name)
    distance = node.get_distance("Acipenser__R")

#def my_layout(node):
#    if node.is_leaf():
#        name_face = AttrFace("name",fsize=25)
#        faces.add_face_to_node(name_face, node, column=0, position="branch-right")
#    else:
#        name_face = AttrFace("name", fsize=20, fgcolor="red")
#        faces.add_face_to_node(name_face, node, column=0, position="branch-right")

#ts = TreeStyle()
#ts.show_leaf_name = False
#ts.show_branch_length = True
#ts.show_branch_support = True
#ts.scale =  180
#ts.layout_fn = my_layout
#out = ptree.render("foo.png", units="mm", tree_style=ts,dpi=400)

## figure out how to get children nodes for a given parent (i.e. branch)
for node in ptree.traverse("levelorder"):
    if node.name != 'B3':
        continue

    descendants = [n.name for n in node.get_descendants()]
    print node.name,descendants
    
## from the current parent and a list of children we identify the transitions
print aas
trans = get_transitions("Acipenser__R",descendants)
print trans
