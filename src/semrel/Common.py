#!/usr/bin/python
# -*- coding: utf-8 -*-

#========================
import gzip
class Lang:
    '''Language ids.'''
    en = 'en'
    hi = 'hi'
    bn = 'bn'
    kn = 'kn'

class Sep:
    '''Separators in data files.'''
    w1w2 = '\t' # two words
    w2wt = ' ' # word and weight (probability)
    lw = ' ' # language id and word
    nid = ":" # in node ids (in Graph.py)
    text = ' ' # while parsing free text

class Lexicon:
    def __init__(self, l1, l2, dfile=None): # dfile = dictionary
        self.lex = {l1:{}, l2:{}}
        if dfile:
            self.loaddfile(l1, l2, dfile)

    def add(self, l1, w1, l2, w2):
        lex = self.lex
        if w1 not in lex[l1]:
            lex[l1][w1] = {}
        if w2 not in lex[l2]:
            lex[l2][w2] = {}
        lex[l1][w1][w2] = 1.0
        lex[l2][w2][w1] = 1.0

    def loaddfile(self, l1, l2, dfile):
        # Dictionary file format:
        # w11 <tab> w21
        # w11 <tab> w22
        # w12 <tab> w23
        # ...
        i = 1
        for line in open(dfile):
            line = line.decode('utf-8').rstrip()
            line = line.lower()
            w1, w2 = line.split(Sep.w1w2)
            self.add(l1, w1, l2, w2)
            i += 1

    def haslangpair(self, l1, l2):
        if l1 in self.lex and l2 in self.lex:
            return True
        else: return False

    def printlex(self, l1):
        lex = self.lex
        for w1 in lex[l1]:
            for w2 in lex[l1][w1]:
                print (w1 + Sep.w1w2 + w2).encode('utf-8')

class L1L2PMF: # for each word in L1, a pmf over words in L2
    def __init__(self, l1, l2, l1l2prfile):
        # l1l2prfile has format:
        # w11 <tab> w21 pr <tab> w22 pr ...
        # w12 <tab> w23 pr <tab> w21 pr ...

        self.l1, self.l2 = l1, l2
        self.pmf = {}

        f = gzip.open(l1l2prfile, 'rb')
        for line in f:
            line = line.decode('utf-8').rstrip()
            rec = line.split('\t')
            w1 = rec.pop(0)
            self.pmf[w1] = {}
            for w2pr in rec:
                w2, pr = w2pr.split(' ')
                self.pmf[w1][w2] = float(pr)
        f.close()
    @staticmethod
    def truncate_pmf(cands, truncprob, renormalize=True, maxcands=99999):
        '''Truncate pmf and normalize.
        Input: a pmf cands = [(w1 pr1), (w2 pr2), ...] sorted in descending order of pr.
        Output: truncated to truncprob and renormalized to be a pmf.'''
        tot = 0
        for i in range(len(cands)):
            tot += cands[i][1]
            if tot >= truncprob:
                break
        cands = cands[:min(i+1,maxcands)]
        #re-normalize
        if renormalize:
            tot = sum([sc for w, sc in cands])
            if tot != 0:
                cands = [(w, sc/tot) for w, sc in cands]
        return cands