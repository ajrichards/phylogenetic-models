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

#def my_layout(node):
#    if node.is_leaf():
#         name_face = AttrFace("name",fsize=25)
#         faces.add_face_to_node(name_face, node, column=0, position="branch-right")
#
#ts = TreeStyle()
#ts.show_leaf_name = False
#ts.show_branch_length = True
#ts.show_branch_support = True
#ts.scale =  180
#ts.layout_fn = my_layout
#out = t.render("foo.png", units="mm", tree_style=ts,dpi=400)

## use the number of distances from root to figure out how many levels there are
branches = {}
debug = 0
rootNode = "Acipenser__R"
branches["0"] = [rootNode]
currentParent = rootNode
level = 1
ptree = PhyloTree(nwTree)

## root the tree
ptree.set_outgroup(rootNode)

#for x in  dir(pTree):
#    if re.search("root",x):
#        print x
#print pTree.get_tree_root()
#sys.exit()


for node in ptree.traverse("postorder"): # levelorder
    if node.name == rootNode:
        continue
    if node.name == 'NoName':
        level += 1
    
    if not branches.has_key(str(level)):
        branches[str(level)] = []
    branches[str(level)].append(node.name)
    debug += 1
    distance = node.get_distance("Acipenser__R")
    descendants = node.get_descendants()
    ancestors = node.get_ancestors()
    actualDescendants = []
    for d in descendants:
        if d.name != 'NoName':
            actualDescendants.append(d)

    #evolEvents = [event.etype for event in node.get_my_evol_events()]
    print "..."
    print node.name, level, distance,len(ancestors)
    event = node.get_my_evol_events()[0]
    print event.branch_supports, event.dup_score, event.e_newick, event.etype, event.famSize, event.fam_size 
    #print event.in_seqs, event.inparalogs, event.node, event.orthologs, event.out_seqs, event.outgroup
    print event.outgroup_spcs, event.root_age, event.seed, event.sos
    #sys.exit()
    #print node.name, level, distance,len(ancestors),evolEvents #ptree.get_descendant_evol_events(node.name)
    #print "...", node.get_my_evol_events()
    
    if debug > 20:
        break

#for key in range(len(branches.keys())):
#    print key, branches[str(key)]



print get_tree_file_path("SPLIT1","5000")
print nwTree
print dir(node)
#print dir(pTree)
sys.exit()

print(dir(node))
#print dir(t)
#print dir(ptree)

#t.show()
#print tree
#print filePath
#print 'total lines', totalLines
#sys.exit()

## create transition matrix
t = Tree( '((H:1,I:1):0.5, A:1, (B:1,(C:1,D:1):0.5):0.5);' )
print t
#                    /-H
#          /--------|
#         |          \-I
#         |
#---------|--A
#         |
#         |          /-B
#          \--------|
#                   |          /-C
#                    \--------|
#                              \-D

# I get D
D = t.search_nodes(name="D")[0]

# I get all nodes with distance=0.5
nodes = t.search_nodes(dist=0.5)
print len(nodes), "nodes have distance=0.5"

# We can limit the search to leaves and node names (faster method).
D = t.get_leaves_by_name(name="D")
print D




sys.exit()
