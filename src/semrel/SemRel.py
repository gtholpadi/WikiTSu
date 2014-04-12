'''
Created on 17-Feb-2013
Collects algorithms for semantic relatedness score computation.
@author: gtholpadi
'''
from Common import L1L2PMF
from MyUtils import MyUtils
import sys
from Eval import Eval

class SemRelPTM:
    @staticmethod 
    def SRptm(l1, l2, l1l2pairsfile, l1l2candfile, l2l1candfile):
        '''
        Input: p(t|s), p(s|t), {(s,t)} pairs from human annotation task
        Output: {(s,t,score)} for all the pairs
        Algo:
            For each s,t pair
                Get distributions p(t'|s) and p(s'|t).
                l-INF normalize the distributions to get r(t'|s) and r(s'|t)
                Take the max of r(t|s) and r(s|t); this is the score.
        '''
        l1l2, l2l1 = L1L2PMF(l1, l2, l1l2candfile), L1L2PMF(l2, l1, l2l1candfile)
        for line in open(l1l2pairsfile):
            line = line.decode('utf-8').rstrip()
            w1, w2 = line.split()
            if w2 in l1l2.pmf[w1]:
#                rl1l2 = l1l2.pmf[w1][w2] / max( l1l2.pmf[w1].itervalues() )
                rl1l2 = l1l2.pmf[w1][w2]
            else: rl1l2 = 0
            if w1 in l2l1.pmf[w2]:
#                rl2l1 = l2l1.pmf[w2][w1] / max( l2l1.pmf[w2].itervalues() )
                rl2l1 = l2l1.pmf[w2][w1]
            else: rl2l1 = 0
            score = max(rl1l2, rl2l1)
            print (w1+'\t'+w2+'\t'+unicode(score)).encode('utf-8')

def evalscores(goldfile, methfile):
    goldsc = {}
    for line in open(goldfile):
        line = line.decode('utf-8').rstrip()
        w1, w2, sc = line.split()
        goldsc[w1+'#'+w2] = float(sc)
    methsc = {}
    for line in open(methfile):
        line = line.decode('utf-8').rstrip()
        w1, w2, sc = line.split()
        methsc[w1+'#'+w2] = float(sc)
    scores = [(goldsc[w1w2], sc) for w1w2, sc in methsc.iteritems()]
    scorr, spval = Eval.spear_corr([s1 for s1, s2 in scores], [s2 for s1, s2 in scores])
    pcorr, ppval = Eval.pears_corr([s1 for s1, s2 in scores], [s2 for s1, s2 in scores])
    kcorr, kpval = Eval.ktau_corr([s1 for s1, s2 in scores], [s2 for s1, s2 in scores])
    print "%f\t%f\t%f\t%f\t%f\t%f" % (scorr, spval, pcorr, ppval, kcorr, kpval)
        
if __name__ == "__main__":
    cmd, args = sys.argv[1], sys.argv[2:]
    if cmd == "ptm":
        l1, l2, l1l2pairsfile, l1l2candfile, l2l1candfile = args
        SemRelPTM.SRptm(l1, l2, l1l2pairsfile, l1l2candfile, l2l1candfile)
    elif cmd == "eval":
        goldfile, methfile = args
        evalscores(goldfile, methfile)
    else:
        print "invalid command:", cmd
