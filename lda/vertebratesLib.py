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
import numpy as np
from ete2 import Tree

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
transitions = []
for t1 in base:
    for t2 in base:
        transitions.append("%s%s"%(t1,t2))
transitions.sort()
transitions = np.array(transitions)

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

def transitions_to_counts(transitionList):
    """
    given a list of transition pairs translate it into a row in the transition count matrix
    """
    
    if type(transitionList) == type([]):
        transitionList = np.array(transitionList)

    uniqueList = np.unique(transitionList)
    counts = np.zeros(transitions.size,)
    for ut in uniqueList:
        indx = np.where(transitions==ut)[0]
        count = np.where(transitionList == ut)[0].size

        if indx.size == 0:
            raise Exception("Invalid value in transitionList '%s'"%(ut))
        indx = np.where(transitions==ut)[0]
        counts[indx] = count

    return counts


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

def fix_tree(nwTree):
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
        if node.name == 'NoName':
            level += 1
            iNodeName = 'N'+str(level)
            node.name = iNodeName

        if not internalNodes.has_key(iNodeName):
            internalNodes[iNodeName] = []
        internalNodes[iNodeName].append(node.name)

    pattern1 = "[\w|:]+.+?[,|\)|\(]"
    prog1 = re.compile(pattern1)
    result = prog1.findall(nwTree)

    pattern2 =  "\d\.\d+:[A-Z]"
    prog2 = re.compile(pattern2)
    
    pattern3 =  "\_[A-Z]"
    prog3 = re.compile(pattern3)

    ## traverse the tree to fix distances and extract information
    n =0
    treeSummary = {}
    for node in t.traverse("postorder"):

        if node.name == 'N1':
            continue
    
        res = re.sub("[\(|\)|,]","",result[n])
        res = re.sub(node.name+":","",res)
        transitions = prog2.findall(res)

        ## leaf node
        if prog3.search(node.name):
            current = prog3.findall(node.name)[0][-1]
        ## internal node
        else:
            current = prog3.findall(res)[0][-1]
            node._set_dist(float(transitions[0].split(":")[0]))
        
        treeSummary[node.name] = {'aa':current,\
                                  'dist':node.dist,\
                                  'transitions':transitions}
        n += 1

    ## again traverse this time identifing transition pairs
    transitionPairs = {}
    debug = 0
    for node in t.traverse("postorder"):
        if node.name == 'N1':
            continue
        parent = [n.name for n in node.get_ancestors()][0]
        treeSummary[node.name]['parent'] = parent

        if parent == 'N1':
            parentAA = treeSummary[node.name]['aa']
        else:
            parentTransitions = treeSummary[parent]['transitions']
            parentAA = parentTransitions[0].split(":")[-1]
        
        childTransitions = treeSummary[node.name]['transitions']
        childAAs = [child.split(":")[-1] for child in childTransitions]
        childAAs.reverse()

        transitionPairs = []
        p = parentAA
        for c in childAAs:
            transitionPairs.append(p+c)
            p = c

        treeSummary[node.name]['pairs'] = transitionPairs
        debug += 1
        
    return t,treeSummary
