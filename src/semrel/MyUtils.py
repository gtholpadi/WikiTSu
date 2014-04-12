#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
import math
import random
from Common import *
class MyUtils:
    @staticmethod
    def pnorm(v, p=2):
        '''p-Norm of a vector (implemented as a sequence'''
        if p == 'INF':
            norm = max( [math.fabs(x) for x in v] )
        elif p == 1:
            norm = sum( [math.fabs(x) for x in v] )
        elif p == 2:
            norm = math.sqrt( sum( [x*x for x in v] ) )
        else:
            norm = sum( [math.fabs(x)**p for x in v] ) ** (1.0/p)
        return norm
    @staticmethod
    def normalize(v, p=1):
        '''Lp-normalize the vector v.'''
        pnorm = MyUtils.pnorm(v, p)
        return [x/pnorm for x in v]
    @staticmethod
    def cumu_cutoff(a, t):
        '''Get the shortest prefix of the list a whose sum exceeds the cutoff threshold t.'''
        tot = 0
        for i in range(len(a)):
            tot += a[i]
            if tot > t:
                break
        return a[:i+1]
    @staticmethod
    def abs_disc(q):
        '''Return distribution q after smoothing using absolute discounting.
        '''
        eps = min([qi for qi in q if qi>0])/100.
        q = [qi+eps for qi in q]
        return MyUtils.normalize(q,p=1)
    @staticmethod
    def kld(p, q, base=math.e, smooth='ABS_DISC'):
        '''Return KL Divergence of pmf p from q.
        Ignore very low probabilities (below 1e-6).
        '''
        if 0 in q: q = MyUtils.abs_disc(q)
        return sum( [ p[i] * (math.log(p[i], base) - math.log(q[i], base) ) \
                     for i in range(len(p)) if p[i] > 1e-7] )
    @staticmethod
    def symkld(p, q, base=math.e):
        '''Return symmetrized KL Divergence of pmfs p and q.
        Ignore very low probabilities (below 1e-6).
        '''
        return sum( [(p[i]-q[i]) * math.log(p[i]/q[i], base) \
                     for i in range(len(p))] )
    @staticmethod
    def jsd(p, q, base=2):
        '''Return Jensen-Shannon Divergence between pmf's p and q'''
        m = [(p[i]+q[i])/2 for i in range(len(p))]
        return ( MyUtils.kld(p, m, base=base) + MyUtils.kld(q, m, base=base) ) / 2
    @staticmethod
    def dir_sample(K, alpha=None):
        '''Get a sample from a Dirichlet distribution'''
        if not alpha: alpha = [1.0 for j in range(K)] # uniform symmetric dirichlet
        sample = [random.gammavariate(a,1) for a in alpha]
        sample = [v/sum(sample) for v in sample]
        return sample
    @staticmethod
    def spear_corr(X, Y):
        pass
#        return mstats.spearmanr(X, Y, use_ties=True)
#        return stats.spearmanr(X, Y)
#        return rpy.r.corr(X, Y, method="spearman")
#        N = len(X)
#        x, y = Eval.val2rank(X), Eval.val2rank(Y)
#        xm, ym = sum(x)/float(N), sum(y)/float(N)
#        numr = sum( [ (x[i]-xm)*(y[i]-ym) for i in range(N)] )
#        denrx = sum( [ (x[i]-xm)**2 for i in range(N)] )
#        denry = sum( [ (y[i]-ym)**2 for i in range(N)] )
#        r = numr / math.sqrt(denrx*denry)
#        return (r, 1)
#    @staticmethod
#    def val2rank(X):
#        '''Convert values to ranks, with tie-breaking '''
#        N = len(X)
#        temp = [ [X[i], i, None] for i in range(N)] #val, id, rank
#        temp.sort(key=lambda a: a[0]) # sort on val
#        for i in range(N): temp[i][2] = i+1 # assign rank
#        val = [v for v, i, r in temp]
#        ids = [i for v, i, r in temp]
#        rnk = [r for v, i, r in temp]
#
#        i = 0
#        while i < N:
#            j = i+1
#            while j < N and val[j] == val[i]:
#                j += 1
#            tot = sum(rnk[i:j])
#            n = j-i
#            rnk[i:j] = [tot/n]*n
#            i = j
#        x = [None]*N
#        for i in range(N):
#            x[ ids[i] ] = rnk[i]
#        return x

if __name__ == '__main__':
    MyUtils.kld([1,0,0],[0.99, 0.01, 0])
