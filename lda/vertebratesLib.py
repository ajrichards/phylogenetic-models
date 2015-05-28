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

def get_transitions(parent,children):
    """
    given a parent node get the transitions to children
    parent - (str) can be a valid amino acid (upper-case) or a 'somename_AA'
    children - (list) if 'somename_AA' then AA is extracted other wise assumed the same as parent
    """
    
    def extract(name):
        if re.search("\_[ACDEFGHIKLMNPQRSTVWY]$",name):
            return re.findall("\_[ACDEFGHIKLMNPQRSTVWY]$",name)[0][-1]
        else:
            return

    if parent not in aas:
        paa = extract(parent)
        if not paa:
            raise Exception("invalid parent passed")
    else:
        paa = parent

    caas = []
    for child in children:
        echild = extract(child)
        if child in aas:
            caas.append(child)
        elif echild:
            caas.append(echild)
        else:
            caas.append(paa)

    return [paa + caa for caa in caas]
