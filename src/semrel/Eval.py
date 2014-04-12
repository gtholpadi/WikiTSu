#!//usr/bin/python
# -*- coding: utf-8 -*-
import sys
from Common import L1L2PMF, Sep, Lexicon
#import rpy
from scipy.stats import mstats
#from scipy.stats import stats
#import math
from MyUtils import MyUtils
import gzip
# Evaluation
# Given an output file of the form:
#   w11 <tab> w21 prob <tab> w22 prob <tab> ...
#   w12 <tab> ...
#   ...
# compute the MRR as follows:
# RR for w11 = 1/min(rank(w21), rank(w22))
# MRR(lang1) = MEAN(RR for each word)

class Eval:
    def __init__(self, lex, V2):
        # lex = Lexicon object, V2 = l2 vocab size
        self.lex = lex
        self.V2 = V2

    def evaluate(self, l1, l2, ofile, qfile, colnum):
        # ofile = output file with ranked l2 words for each l1 word
        # qfile = query file with l1 words which should be evaluated
        qws = {}
        for line in open(qfile):
            line = line.decode('utf-8').rstrip()
            qws[line.split()[colnum-1]] = True
        lex = self.lex
        sco = {}
        for line in gzip.open(ofile):
            line = line.decode('utf-8').rstrip()
            rec = line.split(Sep.w1w2)
            w1 = rec.pop(0)
            rec = [(cand.split(Sep.w2wt)) for cand in rec]
            #if w1 in qws and w1 in lex.lex[l1]:
                #sco[w1] = []
                #sco[w1].append(self.rr(l1, l2, w1, rec))
                ##sco[w1].append(self.precatk(l1, l2, w1, rec, 1))
                ##sco[w1].append(self.precatk(l1, l2, w1, rec, 5))
                ##sco[w1].append(self.rprec(l1, l2, w1, rec))
                #sco[w1].append(self.presatk(l1, l2, w1, rec, 1))
                #sco[w1].append(self.presatk(l1, l2, w1, rec, 5))
            if w1 in qws:
                sc = self.evalword(l1, l2, w1, rec)
                if sc:
                    sco[w1] = sc
        return sco

    def evalword(self, l1, l2, w1, cands):
        lex = self.lex
        if w1 in lex.lex[l1]:
            sco = []
            sco.append(self.rr(l1, l2, w1, cands))
            #sco[w1].append(self.precatk(l1, l2, w1, rec, 1))
            #sco[w1].append(self.precatk(l1, l2, w1, rec, 5))
            #sco[w1].append(self.rprec(l1, l2, w1, rec))
            sco.append(self.presatk(l1, l2, w1, cands, 1))
            sco.append(self.presatk(l1, l2, w1, cands, 5))
        else:
            sco = None
        return sco
#    def score(self, l1, l2, w1, cands):
#        return self.rr(l1, l2, w1, cands)

    def rr(self, l1, l2, w1, cands):
        # Reciprocal rank
        found = False
        for i, cand in enumerate(cands):
            w2, pr = cand
            if w2 in self.lex.lex[l1][w1]:
                found = True
                break
        if found:
            return 1.0/(i+1)
        else:
            return 1.0/self.V2

    def precatk(self, l1, l2, w1, cands, k):
        # Precision @ k
        return [w2 in self.lex.lex[l1][w1] for w2, pr in cands[:k]].count(True)/float(k)

    def rprec(self, l1, l2, w1, cands):
        # R-Precision: same as prec@k, with k = #relevant docs for current query
        return self.precatk(l1, l2, w1, cands, len(self.lex.lex[l1][w1]))

    def presatk(self, l1, l2, w1, cands, k):
        # Presence @ k : returns 1 if present within first k
        present = [w2 in self.lex.lex[l1][w1] for w2, pr in cands[:k]].count(True)
        if present > 0: return 1
        else: return 0

    @staticmethod
    def spear_corr(X, Y):
        return mstats.spearmanr(X, Y, use_ties=True)
#        return stats.spearmanr(X, Y)
#        return rpy.r.corr(X, Y, method="spearman")
    @staticmethod
    def pears_corr(X, Y):
        return mstats.pearsonr(X,Y)
#        return stats.pearsonr(X,Y)
#        return rpy.r.corr(X, Y, method="pearson")
    @staticmethod
    def ktau_corr(X, Y):
        return mstats.kendalltau(X, Y, use_ties=True, use_missing=False)

def evalble(dfile, dl1, dl2, ofile, ol1, ol2, ol1qfile, ol1colnum, ol2vocfile, evtype):
    '''Evaluate bilingual lexicon extraction.
    dfile = dictionary file (see Lexicon for format)
    dl1, dl2 = languages l1 and l2 for dictionary file
    ofile = ranking of (language) ol2 words, for each ol1 word (gzipped)
    ol1qfile = ol1 words which should be evaluated
    ol1colnum = column in ol1qfile that contain the ol1 words
    ol2vocfile = ol2 vocabulary list, to determine V2
    evtype = evaluation type (SUM/DET)
    '''
    V2 = len( open(ol2vocfile).readlines() )
    lex = Lexicon(dl1, dl2, dfile)
    ev = Eval(lex, V2)
    sco = ev.evaluate(ol1, ol2, ofile, ol1qfile, ol1colnum)
    if evtype == 'SUM':
        vals = zip(*sco.values())
        for vallist in vals:
            print sum(vallist) / float(len(vallist)),
        print len(sco)
#        print sum([rec[0] for w, rec in sco.iteritems()]) / len(sco),
    else:
        for w, rec in sco.iteritems():
            line = w+'\t'+'\t'.join([unicode(sc) for sc in rec])
            print line.encode('utf-8')
#        print ('\n'.join([w+'\t'+'\t'.join([unicode(sc) for sc in rec]) for w, rec in sco.iteritems()])).encode('utf-8')

def evalsemrel(l1, l2, pairsfile, l1colnum, l2colnum, l2gpmffile, l1l2methpmffile):
    '''Evaluate semantic relatedness.
    l1 and l2 = languages s and t in p(t|s)
    pairsfile = translation pairs in l1 and l2, l1 words in column <colnum>
    l2gpmffile = gold pmf over l2 words, for l2 words (including those in pairsfile)
    l1l2methpmffile = method-induced pmf over l2 words, for l1 words (including those
                      in pairsfile)
    '''
    l2gpmf = L1L2PMF(l2, l2, l2gpmffile)
    l1l2pmf = L1L2PMF(l1, l2, l1l2methpmffile)
#    print "#JSDiv Spearmanr"
    jsds, rhos = [], []
    for line in open(pairsfile):
        line = line.decode('utf-8').rstrip()
        pair = line.split()
        w1, w2 = pair[l1colnum-1], pair[l2colnum-1]
        gpmf = l2gpmf.pmf[w2] #gold pmf
        mpmf = l1l2pmf.pmf[w1] #method pmf
        vecs = [ (gpmf[x2], mpmf[x2]) for x2 in gpmf if x2 in mpmf]
        vecs.extend( [ (gpmf[x2], 0.0) for x2 in gpmf if x2 not in mpmf] )
        vecs.extend( [ (0.0, mpmf[x2]) for x2 in mpmf if x2 not in gpmf] )
        gvec, mvec = zip(*vecs)
        jsd = MyUtils.jsd(gvec, mvec, base=2)
        rho, pval = mstats.spearmanr(gvec, mvec, use_ties=True)
        jsds.append(jsd)
        rhos.append(rho)
        print "%f\t%f\t%f" % (jsd, rho, pval)
    print "\t\t\t%f\t%f" % ( sum(jsds)/len(jsds), sum(rhos)/len(rhos) )
#    print "%f" % ( sum(jsds)/len(jsds) )

if __name__ == "__main__":
    cmd, args = sys.argv[1], sys.argv[2:]
    if cmd == "evalble":
        dfile, dl1, dl2, ofile, ol1, ol2, ol1qfile, ol1colnum, ol2vocfile, evtype  = args
        evalble(dfile, dl1, dl2, ofile, ol1, ol2, ol1qfile, int(ol1colnum), ol2vocfile, evtype)
    elif cmd == "evalsemrel":
        l1, l2, pairsfile, l1colnum, l2colnum, l2gpmffile, l1l2methpmffile = args
        evalsemrel(l1, l2, pairsfile, int(l1colnum), int(l2colnum),
                   l2gpmffile, l1l2methpmffile)
    else:
        print "invalid command:", cmd
