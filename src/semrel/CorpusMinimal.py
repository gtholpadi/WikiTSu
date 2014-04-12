#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import os
import cPickle as pickle
import gzip
from Common import *
import math
# ========================

class Corpus:

    def __init__(self, lang):  # language of the corpus
        self.lang = lang
        self.index = {} # for each term, for each doc, the number of occurrences
        self.docsz = {} # for each doc, the doc size
        self.vocab = {} # for each term, the number of occurrences
        self.total = 0.0 # total number of tokens
        self.V = 0 # number of unique terms
        self.wdf = {} #for each word, its document frequency
        self.ndocs = 0.0 #number of documents
        self.tii = {} #tf-idf index
    def save(self, corpfile):
        f = gzip.open(corpfile, 'wb')
        pickle.dump(self, f)
        f.close()
    @staticmethod
    def load(corpfile):
        f = gzip.open(corpfile)
        corp = pickle.load(f)
        f.close()
        return corp
    def readdir(self, corpdir, regex=None):
        import re
        if regex: pat = re.compile(regex)
        for fname in os.listdir(corpdir):
            if regex and not pat.search(fname): continue
            fil = open(os.path.join(corpdir,fname))
            for line in fil:
                line = line.decode('utf-8').rstrip()
                for tok in line.split():
                    self.update((fname, tok))
            fil.close()
    def update(self, dt):  # (doc, token)
        (doc, word) = dt
        if word not in self.index:
            self.index[word] = {}
            self.vocab[word] = 0.0
        if doc not in self.index[word]:
            self.index[word][doc] = 0.0
        if doc not in self.docsz:
            self.docsz[doc] = 0.0
            self.ndocs += 1.0

        self.index[word][doc] += 1.0
        self.docsz[doc] += 1.0
        self.vocab[word] += 1.0
        self.total += 1.0
    def setwdf(self):
        self.wdf = {} #for each word, its document frequency
        vocab, index, wdf = self.vocab, self.index, self.wdf
        for w in vocab:
            wdf[w] = len(index[w])
    def settfidf(self):
        if not self.wdf:
            self.setwdf()
        self.tii = {} #tf-idf index
        index, tii, ndocs, wdf = self.index, self.tii, self.ndocs, self.wdf
        for w in index:
            tii[w]={}
            idf = math.log(ndocs/wdf[w])
            for d in index[w]:
                # tfidf(w,d) = tf(w,d) log ndocs/df(w)
                # tf(w,d) = 1 + log count(w,d) if count(w,d)>0 else 0
                tii[w][d] = (1+math.log(index[w][d])) * idf
    def printstats(self, outstr = sys.stdout):
        header = '#NUM_DOCS=%d\tNUM_WORDS=%d\tNUM_TOKENS=%.0f' % \
                 ( len(self.docsz), len(self.vocab), self.total )
        print >> outstr, header.encode('utf-8')
        numdocs = len(self.docsz)
        print >> outstr, '#word\tcoll.freq\tdoc.freq\tavg.doc.ct\tidf\tcfidf\ttcfidf'.encode('utf-8')
        for w in self.vocab:
            cf = self.vocab[w]
            df = len(self.index[w])
            adc = float(cf)/df
            idf =  math.log(numdocs / df)
            cfidf = cf * idf
            tcfidf = math.log(cf) * idf #tapered cfidf
            line = '%s\t%.0f\t%.0f\t%.4f\t%.4f\t%.4f\t%.4f' % (w, cf, df, adc, idf, cfidf, tcfidf)
            print >> outstr, line.encode('utf-8')
    @staticmethod
    def loadwstat(wstatfile):
        wstat = {}
        for line in open(wstatfile):
            line = line.rstrip().decode('utf-8')
            if line.startswith(u'#'): continue
            rec = line.split(u'\t')
            wd = rec.pop(0)
            wstat[wd] = rec
        return wstat

class MLCorpus:  # multilingual corpus

    def __init__(self, langs):
        self.langs = langs
        self.corp = dict( [(l, Corpus(l)) for l in langs] )
#        self.corp = dict( map(lambda x: (x, Corpus(x)), langs) )

    def loadcorpus(self, malletstatefile, langs):
        '''
        Input: unzipped state file from mallet train-topics
               a dict mapping 0 to en, 1 to hi, etc.
        '''
        f = open(malletstatefile)
        f.readline() # skip first line
        for line in f:
            (
                doc,
                l,
                pos,
                t,
                word,
                topic,
                ) = line.decode('utf-8').split()
            self.corp[langs[l]].update((doc, word))
        f.close()
        for l in self.corp:
            self.corp[l].postprocess()

    def pw(self, ws):  # list of (l, w) pairs

#        # Try 1: use independence assumption
#        pr = 1.0
#        for (l, w) in ws:
#            pr *= self.corp[l].pw(w)
#        return pr
        # Try 2: without independence assumption
        pr = 1.0
        for i, w in enumerate(ws):
            pr *= self.pwgws(w, ws[i+1:])
        return pr

    def pwgws(self, w1, ws):
        # p(w1 | w2,...,wk)
        # Computed as follows:
        # Intersect the document sets of w2,...,wk.
        # Compute p(w) in this restricted corpus.
        corp = self.corp
        l1, w1 = w1
        # no conditioning
        if not ws:
            return self.corp[l1].pw(w1)
        # get restricted corpus
        l, w = ws[0]
        rescorp = set(corp[l].index[w].keys()) # restricted document set
        for l, w in ws[1:]:
            if len(rescorp) == 0: break
            rescorp = rescorp & set(corp[l].index[w].keys())
        # fallback to full corpus, if restricted corpus empty
        if len(rescorp) == 0:
            return self.corp[l1].pw(w1)
        # p(w1|restricted) with add-1 smoothing
        ntok = sum([corp[l1].docsz[doc] for doc in rescorp])
        nw1 = sum(corp[l1].index[w1].values())
        return (nw1+1)/(ntok+corp[l1].V)

    def vocab(self, l):
        return self.corp[l].vocab.keys()

def main():
    (stfile, corpfile) = sys.argv[1:3]
    id2lang = sys.argv[3:]  # e.g. 0:en 1:hi 2:bn ...
    id2lang = dict( [idla.split(':') for idla in id2lang] )
#    id2lang = dict(map(lambda x: x.split(':'), id2lang))

    mlc = MLCorpus(id2lang.values())
    mlc.loadcorpus(stfile, id2lang)
    f = open(corpfile, 'w')
    pickle.dump(mlc, f)
    f.close()


if __name__ == '__main__':
    main()
