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

import csv, os,re,sys,shutil,cPickle
from ete2 import Tree,TreeStyle,PhyloTree
from ete2 import faces, AttrFace
from vertebratesLib import *
import matplotlib.pyplot as plt

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

## prep to save the transition matrices
outputDir = os.path.join("..","data","hv-compressed")
if not os.path.isdir(outputDir):
    os.mkdir(outputDir)

## loop through the files
#for split in splits:
split = "SPLIT1"
position = "25000"
outputFile = os.path.join(outputDir,"%s-%s.npz"%(split,position))
treeList = get_trees(split,position,dataDir)

countMatrix = np.zeros((len(treeList),len(TRANSITIONS)),)
t = 0
for t,pbTree in enumerate(treeList):
    fixedTree,treeSummary = fix_tree(pbTree)
    tlist = []
    for item in treeSummary.itervalues():
        tlist.extend(item['pairs'])
    counts = transitions_to_counts(tlist)
    countMatrix[t,:] = counts

print countMatrix.shape

# summarizing plot
cmap = plt.cm.PuBuGn
fig = plt.figure(figsize=(10,4))
ax = fig.add_subplot(1,1,1)
hmap = ax.imshow(countMatrix, interpolation='nearest',aspect='auto',cmap=cmap)
_ = ax.set_title("%s-%s"%(split,position))
#_ = ax.set_xticks(range(countMatrix.shape[1]))
#_ = ax.set_yticks(range(countMatrix.shape[0]))
#_ = ax.set_xticklabels(TRANSITIONS,fontsize=7)
#_ = ax.set_yticklabels(["s%s"%i for i in range(countMatrix.shape[0])],fontsize=7)
xlabs = ax.get_xticklabels()
plt.setp(xlabs, rotation=90)
#ax.set_aspect(1./ax.get_data_ratio())
plt.show()
#print counts
#print counts.shape

sys.exit()


"""
treeId = 0
for pbTree in treeList:
    print "tree%s"%treeId 
    t,ts = fix_tree(nwTree)
    tlist = []
    for item in ts.itervalues():
        tlist.extend(item['pairs'])

    counts = transitions_to_counts(tlist)
    indices = np.where(counts != 0)[0]
    print transitions[indices]
    print counts[indices]
"""

    
#fixedTree,treeSummary = fix_tree(nwTree)
#print treeSummary
#t = Tree(fixedTree,format=0)
    
## plot tree
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
#out = fixedTree.render("foo.png",units="mm",tree_style=ts,dpi=400)


    

## loop through the trees in a position


#a = np.arange(10)
#b = np.arange(10)
#np.savez_compressed('file.npz', a=a, b=b)
#print len(transitions)
#print len(list(set(transitions)))
sys.exit()




## write tree
#outfid = open("%s-%s.csv"%(split,position), 'wb')
#writer = csv.writer(outfid)
#writer.writerow(["parent","child","aa","dist","transitions","pairs"])
#for node in t.traverse("levelorder"):
#    key = node.name
#    if key == 'N1':
#        continue
#    item = ts[key]
#    writer.writerow([item['parent'],key,item['aa'],item['dist'],";".join(item['transitions']),";".join(item['pairs'])])
#outfid.close()

#for nwTree in treeList:
#    t = fix_tree(nwTree)

sys.exit()

#outfid = open("%s-%s.pkl"%(split,position), 'wb')
#cPickle.dump(t, output, -1)
#outfid.close()

#print node.attributes
#print node.aa

#node.add_feature('something',0.5)

#print node.something
#print dir(node)
#t.show()
#sys.exit()









#sys.exit()
## plot it 


sys.exit()


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
