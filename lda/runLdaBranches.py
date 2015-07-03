#!/usr/bin/env python

import os,csv
import lda
import numpy as np

from vertebratesLib import *

#vocab = TRANSITIONS.tolist()


outputFile = os.path.join("..","data","hv-compressed","branches.npz")
npz = np.load(outputFile)
mat = npz['matrix']
rows = npz['rows']
vocab = npz['columns']

#print("loading data...")
#split = "SPLIT1"
#summaryTree, summarySpecies, splitPositions = get_split_data(split)

X = mat.astype(int)
model = lda.LDA(n_topics=20, n_iter=500, random_state=1)
model.fit(X) 
topic_word = model.topic_word_  # model.components_ also works
n_top_words = 8
outfile1 = open("lda-topics-branches.csv",'w')
writer1 = csv.writer(outfile1)
writer1.writerow(["topic"] + ["w"+str(i) for i in range(n_top_words-1)])
for i, topic_dist in enumerate(topic_word):
    topic_words = np.array(vocab)[np.argsort(topic_dist)][:-n_top_words:-1]
    print('Topic {}: {}'.format(i, ' '.join(topic_words)))
    writer1.writerow(["topic"+str(i)] + topic_words.tolist())
                     
## print top 10
doc_topic = model.doc_topic_
print doc_topic.shape
for i in range(10):
    sortedInds = np.argsort(doc_topic)
    #print sortedInds.shape
    print("{} (top topic: {})".format("position - " + splitPositions[i], doc_topic[i].argmax()))

#print doc_topic[i]
sys.exit()





"""
vocab = lda.datasets.load_reuters_vocab()
titles = lda.datasets.load_reuters_titles()
print X.shape
print X.sum()
model = lda.LDA(n_topics=20, n_iter=1500, random_state=1)
model.fit(X)  # model.fit_transform(X) is also available
topic_word = model.topic_word_  # model.components_ also works
n_top_words = 8
for i, topic_dist in enumerate(topic_word):
    topic_words = np.array(vocab)[np.argsort(topic_dist)][:-n_top_words:-1]
    print('Topic {}: {}'.format(i, ' '.join(topic_words)))
"""        
