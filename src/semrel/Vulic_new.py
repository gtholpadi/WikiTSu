'''
Created on 17-Jan-2013

@author: gtholpadi
'''

from PTM_new import PTM
from Common import L1L2PMF
import math
import sys
import gzip
import time

class VulicMethod:
    def sim(self, l1, w1, l2, w2):
        '''dummy'''
        return 0
    def get_similar(self, l1, w1, l2):
        '''Order words in l2 by similarity to w1 and return
        '''
        ptm = self.ptm
        cands = [(w2, self.sim(l1, w1, l2, w2)) for w2 in ptm.vocabs[l2]]
        cands.sort(key=lambda x: x[1], reverse=True)
        return cands

class Cue(VulicMethod):
    '''Vulic's Cue method
    '''
    def __init__(self, ptm):
        '''Input: PTM learned earlier
        '''
        self.ptm = ptm
    def sim(self, l1, w1, l2, w2):
        '''Same as p(w2|w1) in the paper
        '''
        ptm = self.ptm
        norm = sum( [ptm.phi[l1][t][w1] for t in ptm.topics] )
        sc = sum( [ptm.phi[l2][t][w2] * ptm.phi[l1][t][w1] for t in ptm.topics] ) / norm
        return sc

class TI(VulicMethod):
    '''Vulic's TI method.
    '''
    def __init__(self, ptm):
        '''Input: PTM learned earlier
        '''
        self.ptm = ptm
        self.ti = {} # tfitf values; indexed by language, word and topic
        # set tfitf
        ti = self.ti
        T = float( len(ptm.topics) )
        tfd = {} # TF formula denominators
        for l in ptm.langs:
            tfd[l] = {}
            for t in ptm.topics:
                tfd[l][t] = sum( ptm.tw[l][t].itervalues() )
        for l in ptm.langs:
            ti[l] = {}
            for w in ptm.vocabs[l]:
                ti[l][w] = {}
                nt4w = sum( [1 for t in ptm.topics if w in ptm.tw[l][t]] )
                if nt4w == 0: nt4w = 1
                elif nt4w == T: nt4w = T - 1
                itf = math.log10(T/nt4w)
                for t in ptm.topics:
                    if w in ptm.tw[l][t]: ti[l][w][t] = ptm.tw[l][t][w] * itf / tfd[l][t]
                    else: ti[l][w][t] = 0
                # normalize ti vector for w, so that dot product = cosine similarity
                norm = math.sqrt( sum( [ti[l][w][t]**2 for t in ptm.topics] ) )
                for t in ptm.topics: ti[l][w][t] /= norm
    def sim(self, l1, w1, l2, w2):
        '''Same as cos(w1,w2) in the paper
        '''
        ptm, ti = self.ptm, self.ti
        return sum( [ti[l1][w1][t] * ti[l2][w2][t] for t in ptm.topics] )

class TICue(VulicMethod):
    '''Vulic's TI+Cue method.
    '''
    def __init__(self, ptm, Lambda=10):
        self.ptm = ptm
        self.Lambda = Lambda
        self.ti = TI(ptm)
        self.cue = Cue(ptm)
    def sim(self, l1, w1, l2, w2):
        '''same as in paper
        '''
        return self.ti.sim(l1, w1, l2, w2) + self.Lambda * self.cue.sim(l1, w1, l2, w2)

def savecands(l1, qfile, l2, meth, truncprob, candfile):
    f = gzip.open(candfile, 'wb')
    for w1 in open(qfile).read().decode('utf-8').split():
        cands = meth.get_similar(l1, w1, l2)
        # convert to pmf
        tot = sum([sc for w, sc in cands])
        cands = [(w, sc/tot) for w, sc in cands]
        # truncate and renormalize pmf
        cands = L1L2PMF.truncate_pmf(cands, truncprob, renormalize=True)
        line = w1 + '\t' + '\t'.join( [w2+' '+str(sc) for w2, sc in cands] )
        print >> f, line.encode('utf-8')
    f.close()
    
def runticue(l1, qfile, l2, ptmfile, Lambda, truncprob, candfile):
    tic = TICue(PTM.load(ptmfile), Lambda)
    savecands(l1, qfile, l2, tic, truncprob, candfile)

def runcue(l1, qfile, l2, ptmfile, truncprob, candfile):
    cue = Cue(PTM.load(ptmfile))
    savecands(l1, qfile, l2, cue, truncprob, candfile)

if __name__ == "__main__":
    cmd, args = sys.argv[1], sys.argv[2:]
    if cmd == "runticue":
        l1, qfile, l2, ptmfile, Lambda, truncprob, candfile = args
        runticue(l1, qfile, l2, ptmfile, float(Lambda), float(truncprob), candfile)
    elif cmd == "runcue":
        l1, qfile, l2, ptmfile, truncprob, candfile = args
        runcue(l1, qfile, l2, ptmfile, float(truncprob), candfile)
    else:
        print "invalid command: ", cmd