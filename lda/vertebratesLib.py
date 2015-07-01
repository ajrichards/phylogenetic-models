#!/usr/bin/env python
"""
get data into a format that LDA understands
https://pythonhosted.org/ete2/tutorial/tutorial_trees.html#reading-newick-trees

the data are organized into splits so to keep things more managable 
we first move individual files into corresponding splits

extract the tar file into 'dataDir' and run this script

there are basically 30,000 positions in each split
for a total of 30000 * X splits


NOTES: 

    (1) in this example the tree is contstrained


"""

import csv, os,re,sys,shutil
import numpy as np
import matplotlib.pyplot as plt
from ete2 import Tree,TreeStyle
from ete2 import faces, AttrFace

## variables
aa2codon = {'C': ['TGT','TGC'],\
            'D': ['GAT','GAC'],\
            'S': ['TCT', 'TCG', 'TCA', 'TCC', 'AGC', 'AGT'],\
            'Q': ['CAA', 'CAG'],\
            'M': ['ATG'],\
            'N': ['AAC', 'AAT'],\
            'P': ['CCT', 'CCG', 'CCA', 'CCC'],\
            'K': ['AAG', 'AAA'],\
            'STOP': ['TAG', 'TGA', 'TAA'],\
            'T': ['ACC', 'ACA', 'ACG', 'ACT'],\
            'F': ['TTT', 'TTC'],\
            'A': ['GCA', 'GCC', 'GCG', 'GCT'],\
            'G': ['GGT', 'GGG', 'GGA', 'GGC'],\
            'I': ['ATC', 'ATA', 'ATT'],\
            'L': ['TTA', 'TTG', 'CTC', 'CTT', 'CTG', 'CTA'],\
            'H': ['CAT', 'CAC'],\
            'R': ['CGA', 'CGC', 'CGG', 'CGT', 'AGG', 'AGA'],\
            'W': ['TGG'],\
            'V': ['GTA', 'GTC', 'GTG', 'GTT'],\
            'E': ['GAG', 'GAA'],\
            'Y': ['TAT', 'TAC']}

aas = aa2codon.keys()
aas.sort()

## transition vocab
base = ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"]
TRANSITIONS = []
for t1 in base:
    for t2 in base:
        TRANSITIONS.append("%s%s"%(t1,t2))
TRANSITIONS.sort()
TRANSITIONS = np.array(TRANSITIONS)

## check that the files are extracted
#dataDir = os.path.join("..","data","herve-vertebrates")
SPLITS = ["SPLIT%s"%i for i in range(1,11)] 

## move the files into subdirectories accoriding to split
def ensure_files_properly_located(dataDir):
    for _split in SPLITS:
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

## convenience functions
def get_positions(split,dataDir):
    positions = []
    debug = 0
    pattern = "_\d+\."
    prog = re.compile(pattern)
    
    
    for fileName in os.listdir(os.path.join(dataDir,split)):
        position = prog.findall(fileName)[0][1:-1]
        positions.append(position)
        
    return np.sort(np.array(positions))

def get_tree_file_path(split,position,dataDir):
    """
    return the map file path for a split and a position
    """

    ## error checking
    if split not in SPLITS:
        raise Exception("invalid split arg %s"%split)    
    positions = get_positions(split,dataDir)
    
    if position not in positions:
        raise Exception("invalid position arg %s"%position)

    if split == "SPLIT10":
        filePath = os.path.join(dataDir,split,"Amphibia-noo-100x1285039-%s-CATG_%s.map"%(split,position))
    else:
        filePath = os.path.join(dataDir,split,"Amphibia-noo-100x1285039-%s-CATG_%s.map"%(split,position))
    if not os.path.exists(filePath):
        raise Exception("cannot find '%s'"%(filePath))
    return filePath

def get_trees(split,position,dataDir):
    """
    return the trees as a list for a given split and position
    """

    filePath = get_tree_file_path(split,position,dataDir)
    fid = open(filePath,'r')
    reader = csv.reader(fid,delimiter=';')
    treeList = []
    for linja in reader:
        if len(linja) == 0:
            continue
        treeList.append(linja[0]+";")
    return treeList

def transitions_to_counts(transitionList):
    """
    given a list of transition pairs translate it into a row in the transition count matrix
    """
    
    if type(transitionList) == type([]):
        transitionList = np.array(transitionList)

    uniqueList = np.unique(transitionList)
    counts = np.zeros(TRANSITIONS.size,)
    for ut in uniqueList:
        indx = np.where(TRANSITIONS==ut)[0]
        count = np.where(transitionList == ut)[0].size

        if indx.size == 0:
            raise Exception("Invalid value in transitionList '%s'"%(ut))
        indx = np.where(TRANSITIONS==ut)[0]
        counts[indx] = count

    return counts

def standardize_tree(nwTree):
    pattern1 = "[\w|:]+.+?[,|\)|\(]"
    parsedTree = nwTree
    
    for r in re.finditer(pattern1,nwTree):
        matched = nwTree[r.start(0):r.end(0)]
        parsed = re.sub(":[A-Z]","",matched)
        multiples = re.findall(":\d\.[\de-]+",parsed)
        trans = "".join(multiples)
        
        if len(multiples) > 1:
            parsed = parsed.replace(trans,multiples[0])
        parsedTree = parsedTree.replace(matched,parsed)

    return re.sub("\)_[A-Z]:",")",parsedTree)[:-3]+";"
                            
def write_tree_to_file(outFileName,treeSummary):
    ## write a treeSummary to file
    outfid = open(outFileName, 'wb')
    writer = csv.writer(outfid)
    writer.writerow(["parent","child","aa","dist","transitions","pairs"])
    for node in t.traverse("levelorder"):
        key = node.name     
        if key == 'N1':
            continue
        item = ts[key]
        writer.writerow([item['parent'],key,item['aa'],item['dist'],";".join(item['transitions']),";".join(item['pairs'])])
    outfid.close()   


def fix_tree(pbTree):
    """
    takes a phylobayes formatted tree
    """
    nwTree = standardize_tree(pbTree)
    t = Tree(nwTree,format=0)
    rootNode = "Acipenser__R"
    #t.set_outgroup(rootNode)

    ## iterate through tree an give each node an identifier
    internalNodes = {}
    level = 0
    iNodeName = 'N'+str(level)

    ## give the nodes names
    for node in t.traverse("levelorder"): #levelorder | postorder
        if node.name == rootNode:
            continue
        if node.name == '':
            level += 1
            iNodeName = 'N'+str(level)
            node.name = iNodeName

        if not internalNodes.has_key(iNodeName):
            internalNodes[iNodeName] = []
        internalNodes[iNodeName].append(node.name)

    pattern1 = "[\w|:]+.+?[,|\)|\(]"
    prog1 = re.compile(pattern1)
    result = prog1.findall(pbTree)

    pattern2 =  "\_[A-Z]"
    prog2 = re.compile(pattern2)

    ## traverse the tree to fix distances and extract information
    n =0
    treeSummary = {}
    for node in t.traverse("postorder"):
        actualNode = re.sub("_+[A-Z]$","",node.name)
        if node.name == 'N1':
            continue

        res = re.sub("[\(|\)|,]","",result[n])
        ## convert scientific notation
        if re.search("\d+\.\d+e-\d+",res):
            for match in re.findall("\d+\.\d+e-\d+",res):
                res = res.replace(match,"%f"%(float(match)))
                
        transitions = re.findall(":\d\.\d+\:[A-Z]",res)
        transitions = [tr[1:] for tr in transitions]
        current = transitions[0][-1]
        
        ## leaf node (sainity check)
        if prog2.search(node.name):
            if not re.search(node.name,res):
                raise Exception("did not pass sanity check")
        ## internal node
        else:
            node._set_dist(float(transitions[0].split(":")[0]))
        
        treeSummary[actualNode] = {'aa':current,\
                                   'dist':node.dist,\
                                   'transitions':transitions}
        n += 1

    ## again traverse this time identifing transition pairs
    transitionPairs = {}
    debug = 0
    for node in t.traverse("postorder"):
        actualNode = re.sub("_+[A-Z]$","",node.name)
        if node.name == 'N1':
            continue
        parent = [n.name for n in node.get_ancestors()][0]
        treeSummary[actualNode]['parent'] = parent

        if parent == 'N1':
            parentAA = treeSummary[actualNode]['aa']
        else:
            parentTransitions = treeSummary[parent]['transitions']
            parentAA = parentTransitions[0].split(":")[-1]
        
        childTransitions = treeSummary[actualNode]['transitions']
        childAAs = [child.split(":")[-1] for child in childTransitions]
        childAAs.reverse()

        transitionPairs = []
        p = parentAA
        for c in childAAs:
            transitionPairs.append(p+c)
            p = c

        treeSummary[actualNode]['pairs'] = transitionPairs
        debug += 1
        
    return t,treeSummary

def get_species_transitions(treeSummary,leaf):
    """
    get all species transitions for a givn leaf
    """

    actualLeaf = re.sub("_+[A-Z]$","",leaf.name)
    leafPairs = treeSummary[actualLeaf]['pairs']
    for ancestor in leaf.iter_ancestors():
        if ancestor.name == "N1":
            continue
        leafPairs.extend(treeSummary[ancestor.name]['pairs'])

    return leafPairs
                         
def profile_heatmap_plot(countMatrix,figName,figTitle=None,cmap=plt.cm.bone,rowIds=None,fontSize=(7,7)):
    """
    create a transition profile plot
    """

    fig = plt.figure(figsize=(9,4))
    ax = fig.add_subplot(1,1,1)
    nonZeroInds = np.where(countMatrix.sum(axis=0) != 0)[0]
    hmap = ax.imshow(countMatrix[:,nonZeroInds], interpolation='nearest',aspect='auto',cmap=cmap)
    if figTitle:
        _ = ax.set_title(figTitle)
    _ = ax.set_xticks(range(nonZeroInds.size))
    _ = ax.set_yticks(range(countMatrix.shape[0]))
    _ = ax.set_xticklabels(TRANSITIONS[nonZeroInds],fontsize=fontSize[0])

    if rowIds:
        ax.set_yticklabels(rowIds,fontsize=fontSize[1])
    else:
        ax.set_yticklabels(["%s"%(i) for i in range(countMatrix.shape[0])],fontsize=fontSize[1])

    xlabs = ax.get_xticklabels()
    plt.setp(xlabs, rotation=90)
    cbar = fig.colorbar(hmap,orientation='vertical')

    plt.savefig(figName,dpi=500)
    
def profile_box_plot(countMatrix,figName,figTitle=None,fontSize=(7,7)):
    """
    create a transition profile plot
    """

    fig = plt.figure(figsize=(9,4))
    ax = fig.add_subplot(1,1,1)
    nonZeroInds = np.where(countMatrix.sum(axis=0) != 0)[0]

    bp = ax.boxplot(countMatrix[:,nonZeroInds], notch=0, sym='+', vert=1, whis=1.5)
    plt.setp(bp['boxes'], color='black')
    plt.setp(bp['whiskers'], color='black')
    plt.setp(bp['fliers'], color='red', marker='+')

    if figTitle:
        _ = ax.set_title(figTitle)
    _ = ax.set_xticks(range(1,nonZeroInds.size+1))
    _ = ax.set_xticklabels(TRANSITIONS[nonZeroInds],fontsize=fontSize[0])

    buff = 0.05
    yrange = countMatrix.max() - countMatrix.min()
    _ = ax.set_ylim(0.0 - (0.05*yrange), countMatrix.max() + (0.05 * yrange))
    
    xlabs = ax.get_xticklabels()
    plt.setp(xlabs, rotation=90,fontsize=fontSize[0])
    ylabs = ax.get_yticklabels()
    plt.setp(ylabs,fontsize=fontSize[1])

    plt.savefig(figName,dpi=500)

def get_split_data(split,dataDir=os.path.join("..","data","hv-compressed")):
    outputFile1 = os.path.join(dataDir,"%s.npz"%(split))
    npz = np.load(outputFile1)
    outputFile2 = os.path.join(dataDir,"%s-rows.npz"%(split))
    splitPositions = np.load(outputFile2)['rows']
    summaryTree = npz['tree']
    summarySpecies = {}
    for key,item in npz.iteritems():
        if key == 'tree':
            continue
        summarySpecies[key] = item
    return summaryTree,summarySpecies,splitPositions

def plot_tree(fixedTree,fileName):

    ## plot the tree
    def my_layout(node):
        if node.is_leaf():
            name_face = AttrFace("name",fsize=25)
            faces.add_face_to_node(name_face, node, column=0, position="branch-right")
        else:
            name_face = AttrFace("name", fsize=20, fgcolor="red")
            faces.add_face_to_node(name_face, node, column=0, position="branch-right")

    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.show_branch_length = True
    ts.show_branch_support = True
    ts.scale =  180
    ts.layout_fn = my_layout
    out = fixedTree.render(fileName, units="mm",tree_style=ts,dpi=400)
