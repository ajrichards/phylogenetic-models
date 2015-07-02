#!/usr/bin/env python

import os,csv
import gensim
import numpy as np

from vertebratesLib import *

#vocab = TRANSITIONS.tolist()


outputFile = os.path.join("..","data","hv-compressed","branches.npz")
npz = np.load(outputFile)
mat = npz['matrix'].astype(int)
rows = npz['rows']
vocab = npz['columns']

#id2word = {}
#for i in range(len(vocab)):
#    id2word[i] = vocab[i]
#    print i, vocab[i]

## create the documents from the matrix
texts = []
usedInds = set([])
for r in range(mat.shape[0]):
    hitInds = np.where(mat[r,:] > 1)[0]
    text = []
    for hi in hitInds:
        word = vocab[hi]
        
        if word[0] == word[1]:
            continue

        usedInds.update([hi])
        text.extend([vocab[hi]] * mat[r,hi])
    texts.append(text)
usedInds = list(usedInds)

print usedInds
print len(usedInds)
sys.exit()
    
dictionary = gensim.corpora.Dictionary(texts)
dictionary.save('/tmp/branchs.dict')
#print(dictionary)
#print(dictionary.token2id)

## create a corpus from the documents
corpus = [dictionary.doc2bow(text) for text in texts]
gensim.corpora.MmCorpus.serialize('/tmp/branches.mm', corpus)
mm = gensim.corpora.MmCorpus('/tmp/branches.mm')
print mm

#corpus = [dictionary.doc2bow(text) for text in texts]    
#id2word = gensim.corpora.Dictionary.load_from_text('wiki_en_wordids.txt')
# id2word = id2word
#corpus = gensim.matutils.Dense2Corpus(mat)

chunksize = 10  
passes = 10

lda = gensim.models.ldamodel.LdaModel(corpus=corpus, num_topics=100, update_every=1, chunksize=chunksize, passes=passes)
print lda.print_topics(20)
print('done')
sys.exit()
#lda = gensim.models.ldamodel.LdaModel(corpus=corpus, id2word=id2word, num_topics=100, update_every=1, chunksize=10000, passes=1)
#print lda.print_topics(20)


#print dir(corpus)
#print corpus.dense
#print dir(corpus.dense)
print 'done'

sys.exit()

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
