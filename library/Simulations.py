#!/usr/bin/env python
"""
TODO
"""

import sys,os,csv,time
import numpy as np
from Codon import Codon

__author__ = "Adam Richards"


class Simulator(Codon):
    """
    constructor

    """

    def __init__(self,logFile="simulation.log"):
        Codon.__init__(self)

        ## set a non-zero probability of all codon transitions
        self.baseMean = 0.001
        self.baseSd   = 0.0001

        ## set the within class distribution
        self.withinClsMean = 30
        self.withinClsSd = 3

        ## non-synomymous to synonymous rate (non-zero and <= 1.0)
        self.nsRate = 0.3

        ## distribution for no transitions
        self.ntMean = 80.0
        self.ntSd = 1.0

        ## initialize logfile
        self.log = csv.writer(open(logFile,'w'))
        self.log.writerow([sys.argv[0]])
        self.log.writerow([time.asctime()])

    def get_base_matrix(self):
        #self.codons = self.cd.aa2pp.keys()
        self.codons = []
        for allaa in self.pp2aa.itervalues():
            for aa in allaa:
                self.codons.extend(self.aa2codon[aa])
        self.codons = np.array(self.codons)
       
        mat = np.random.normal(self.baseMean,self.baseSd,(self.codons.size,self.codons.size))
        return mat

    def generate_class_matrix(self,mat,classCodons):
        """
        generate a class specific rate matrix
        takes a base matrix as a class
        """

        for cd in classCodons:
            if cd not in self.codons:
                raise Exception("invalid class codons specified: %s"%cd)

        x = np.random.normal(self.withinClsMean,self.withinClsSd) # like omega
        cmat = mat.copy()

        ## loop through all codons
        for i,codonI in enumerate(self.codons):
            codonI = self.codons[i]
            for j,codonJ in enumerate(self.codons):
                codonJ = self.codons[j]

                if i > j:
                    continue
                
                ## no transition
                if i == j:
                    cmat[i,i] = np.random.normal(self.ntMean,self.ntSd,1)
                ## only add probabilities to codons that are in class
                elif codonI not in classCodons or codonJ not in classCodons:
                    continue
                ## synonymous mutation
                elif self.codon2aa[codonI] == self.codon2aa[codonJ]:
                    prob = np.random.normal(x,1,1)
                    cmat[i,j] = prob
                    cmat[j,i] = prob
                ## non-synonymous mutation
                else:
                    prob = np.random.normal(self.nsRate*x,self.nsRate*1,1)
                    cmat[i,j] = prob
                    cmat[j,i] = prob

        cmat = cmat/cmat.sum(axis=1).mean()
        return cmat

    def generate_codon_sequence(self,n,cMatrices,class2codon,verbose=False):
        """
        cMatrices - a dictionary of codon matrices with keys serving as identifiers
        we are ignoring start and stop codons
        class2codon - is a dictionary {classId:[codon1, codon2]}

        """

        ## generate mean and sd to draw class probabilites
        classProbs = {}
        for className in cMatrices.keys():
            classProbs[className] = (np.random.uniform(.5,1),0.1)
        self.log.writerow(["\nclasses"])
        for key,item in class2codon.iteritems():
            self.log.writerow([key,item])

        if verbose:
            print("\nclasses")
            for key,item in class2codon.iteritems():
                print key,item

            print("\nclass probabilities")
            for key,item in classProbs.iteritems():
                print key,item
        self.log.writerow(["\nclass probabilites"])
        for key,item in classProbs.iteritems():
            self.log.writerow([key,item])
        
        ## generate a probability for each class
        codonProbs = {}
        for className,classCodons in class2codon.iteritems():
            for codon in classCodons:
                cpMean,cpSd = classProbs[className]
                cp = np.random.normal(cpMean,cpSd)
                if cp > 1.0:
                    cp = 0.99999
                if cp < 0.0:
                    cp = 0.00001

                codonProbs[codon] = cp

        ## save and print values
        self.log.writerow(["\ncodon probabilites"])
        for key,item in codonProbs.iteritems():
            self.log.writerow([key,item])

        ## get first codon
        classNames = np.array(class2codon.keys())
        seqClasses = np.array([None]*n)
        seqCodons = np.array([None]*n)

        self.log.writerow(["\nseed sequence"])

        ## generate a sequence
        for indx in range(n):
            ## use class probabilities to set the underlying class
            if indx == 0:
                seqClasses[indx] = classNames[np.random.randint(0,classNames.size)]
                prevCodon = None
            else:
                prevClass = seqClasses[indx-1]
                prevCodon = seqCodons[indx-1]

                if np.random.uniform(0,1) <= codonProbs[prevCodon]:
                    seqClasses[indx] = prevClass
                else:
                    seqClasses[indx] = classNames[np.random.randint(0,classNames.size)]

            ## choose based on prob density from transition matrix                
            cmat = cMatrices[seqClasses[indx]]
            normalizedProbs = cmat.sum(axis=0) / cmat.sum(axis=0).sum()
            pairs = [(normalizedProbs[i],self.codons[i]) for i in range(len(normalizedProbs))]
            probs = np.random.multinomial(1, zip(*pairs)[0])
            seqCodons[indx] = self.codons[np.where(probs==1)[0][0]]

        return seqCodons,seqClasses

    def generate_new_sequence(self,seedSeq,seedSeqClasses,cMatrices,class2codon):
        """
        given a seed sequence and it's known classes generate a new sequence
        """
        
        n = seedSeq.size
        seqClasses = np.array([None]*n)
        seqCodons = np.array([None]*n)

        for seqIndx,codon in enumerate(seedSeq):     
            cmat = cMatrices[seedSeqClasses[seqIndx]]
            codonIndx = np.where(self.codons == codon)[0][0]
            normalizedProbs = cmat[:,codonIndx] / cmat[:,codonIndx].sum()
            pairs = [(normalizedProbs[i],self.codons[i]) for i in range(len(normalizedProbs))]
            probs = np.random.multinomial(1, zip(*pairs)[0])
            seqCodons[seqIndx] = self.codons[np.where(probs==1)[0][0]]

        return seqCodons

if __name__ == "__main__":
    print "Running..."

    import matplotlib.pyplot as plt
    
    cmap = plt.cm.PuBuGn
    sim = Simulator()
    mat = sim.get_base_matrix()

    ## generate the matrices
    cMatrices = {}
    classCodons = {}
    for c,classAa in sim.pp2aa.iteritems():
        _classCodons = []
        for aa in classAa:
            for codon in sim.aa2codon[aa]:
                _classCodons.append(codon)
        cMatrices[c] = sim.generate_class_matrix(mat,_classCodons)
        classCodons[c] = _classCodons

    ## generate the seed sequence
    seedSeq,seedSeqClasses = sim.generate_codon_sequence(10,cMatrices,classCodons,verbose=True)
    nextSeq = sim.generate_new_sequence(seedSeq,seedSeqClasses,cMatrices,classCodons)

    for i in range(seedSeq.size):
        print seedSeqClasses[i], seedSeq[i], nextSeq[i]

    c = 'pos'
    cmat = cMatrices[c]
    fig = plt.figure()        
    ax = fig.add_subplot(111)
    hmap = ax.imshow(cmat, interpolation='nearest',aspect='auto',cmap=cmap)
    _ = ax.set_title("class=%s"%c)
    _ = ax.set_xticks(range(sim.codons.size))
    _ = ax.set_yticks(range(sim.codons.size))
    _ = ax.set_xticklabels(sim.codons,fontsize=7)
    _ = ax.set_yticklabels(sim.codons,fontsize=7)
    xlabs = ax.get_xticklabels()
    plt.setp(xlabs, rotation=90)
    ax.set_aspect(1./ax.get_data_ratio())
    cbar = fig.colorbar(hmap, orientation='vertical')
    plt.show()
