#!/usr/bin/env python
"""
TODO

change aa2pp to aa2ppc 
and use http://www.geneinfinity.org/sp/sp_aaprops.html

"""

__author__ = "Adam Richards"

class Codon(object):
    """
    A class to handle codon related functions
    Amino acids are handled using the single character id, with the exception of STOP.
    """

    def __init__(self):
        """
        initialize dictionaries

        """

        ## dictionaries to convert between amino acids and codons
        self.aa2codon = {'C': ['TGT','TGC'],\
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

        self.codon2aa = {}
        for key,val in self.aa2codon.iteritems():
            for c in val:
                self.codon2aa[c] = key

        ## dictionaries to convert between amino acids physical property class
        self.pp2aa = {"neg":["D","E"],\
                      "pos":["K","R","H"],\
                      "pnc":["S","T","C","M","N","Q"],\
                      "aro":["F","Y","W"],\
                      "npa":["G","A","V","L","I","P"]}
        self.aa2pp = {}
        for key,val in self.pp2aa.iteritems():
            for c in val:
                self.aa2pp[c] = key

        ## dictionaries to convert between short and long versions of amino acids
        self.long2short = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D",\
                           "CYS":"C","GLU":"E","GLN":"Q","GLY":"G",\
                           "HIS":"H","ILE":"I","LEU":"L","LYS":"K",\
                           "MET":"M","PHE":"F","PRO":"P","SER":"S",\
                           "THR":"T","TRP":"W","TYR":"Y","VAL":"V"}

        self.short2long = dict([(value,key) for key,value in self.long2short.iteritems()])


if __name__ == "__main__":
    print "Running..."

    cd = Codon()

    if cd.aa2pp["F"] != "aro":
        raise Exception("Failed aa to pp test")
