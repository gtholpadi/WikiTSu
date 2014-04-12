# -*- coding: utf-8 -*-
'''
Created on 30-Nov-2012

@author: gtholpadi
'''
import sys
import gzip
import time
import math
import numpy as np

from MyUtils import MyUtils
from Common import L1L2PMF, Lexicon
from PTM_new import PTM
from Eval import Eval

class SourceCatDist():
    '''Category-specific distributions over source words.'''
    def __init__(self):
        # phi = a dict over categs; each categ has a dict over words
        #       with probabilities
        self.phi = {}
    def set_from_PTM(self, ptm, lang):
        '''Set phi using a PTM object'''
        phi = ptm.phi[lang]
        self.T = len(phi) # number of categs
        self.cat = [None for k in range(self.T)] # categs
        self.phi = [None for k in range(self.T)] # categ-specific word distributions
        for k, phic in enumerate(phi.iteritems()):
            self.cat[k] = phic[0]
            self.phi[k] = phic[1]

class MixModelFine():
    '''Class for computing p*(t|s) using p(t|s) and p_a(t|s) for different 'a's.
    Also includes fine-grained choice of experts based on categories.'''
    def __init__(self, sl, tl, st_pmffiles, scd):
        '''Input: source lang, target lang, p(t|s) pmf file names,
        source categ dist object'''
        self.sl, self.tl = sl, tl
        self.beta = [] # p(t|s) pmfs for different experts
        self.tvocs = [] # target vocabularies for different experts
        self.OneByVt = [] # 1/V_t for different experts
        for fil in st_pmffiles:
            self.beta.append(L1L2PMF(sl, tl, fil).pmf)
            self.tvocs.append( set([ t for s in self.beta[-1].iterkeys()
                                    for t in self.beta[-1][s].iterkeys() ]) )
            self.OneByVt.append( 1.0/len(self.tvocs[-1]) )
        self.E = len(self.beta) # number of experts
        self.scd = scd
        self.T = len(scd.phi) # number of categories
#        self.seteps([1.0/self.E for i in range(self.E)]) # default expert mixture
    def setlatent(self, eps, theta):
        self.eps = eps
        self.theta = theta
    def savelatent(self, latentfile):
        eps, theta, scd = self.eps, self.theta, self.scd
        lf = open(latentfile, 'w')
        print >> lf, "E=%d\tT=%d" % (self.E, self.T)
        for epsk in eps:
            print >> lf, ' '.join( [str(pr) for pr in epsk] )
        print >> lf, '\n'.join( [scd.cat[k]+'\t'+str(pr) for k, pr in enumerate(theta)] )
        lf.close()
#    def readeps(self, epsfile):
#        self.eps = [float(pr) for pr in open(epsfile).read().split()]
    def settrdata(self, trfile, l1, l2):
        self.trlex = Lexicon(l1, l2, trfile)
    def learnEM(self, eps0=None, theta0=None):
        '''Learn latent variable values using EM'''
        sl, beta, E, phi, T, trl = self.sl, self.beta, self.E, self.scd.phi, self.T, self.trlex.lex
        trdata = [ (s,t) for s in trl[sl].iterkeys() for t in trl[sl][s].iterkeys() ]
        N = len(trdata)
        # starting point for optimization
#        # uniform init
        if not eps0: eps0 = [ [1.0/E for j in range(E)] for k in range(T) ]
        if not theta0: theta0 = [1.0/T for k in range(T)]
#        # random init
#        if not eps0: eps0 = [ MyUtils.dir_sample(E) for k in range(T) ]
#        if not theta0: theta0 = MyUtils.dir_sample(T)

        b = [ [ [0 for j in range(E)] for k in range(T)] for n in range(N) ]
        MINCHANGE = 1e-3
        epst1 = eps0
        thetat1 = theta0
        while True:
            # set current value of eps
            epst = epst1
            thetat = thetat1
            # compute the b's
            for n in range(N):
                s, t = trdata[n]
                tot = sum( [ beta[j][s][t] * phi[k][s] * epst[k][j] * thetat[k]
                            for j in range(E) for k in range(T) if t in beta[j][s] ] )
                for k in range(T):
                    for j in range(E):
                        b[n][k][j] = (beta[j][s][t] * phi[k][s] * epst[k][j] * thetat[k])/tot if t in beta[j][s] else 0
            # compute new value of eps, theta
            epst1 = [None for j in range(T)]
            for i in range(T):
                tot = sum( [b[n][i][j] for n in range(N) for j in range(E)] )
                epst1[i] = [ 0 if tot == 0 else sum([b[n][i][l] for n in range(N)])/tot for l in range(E) ]
            tot = sum( [sum(bnk) for bn in b for bnk in bn] )
            thetat1 = [ sum([b[n][i][j] for n in range(N) for j in range(E)]) / tot for i in range(T) ]
            # check for convergence
            if np.linalg.norm(np.subtract(epst, epst1)) < MINCHANGE: # Frobenius norm
                break
        self.setlatent(epst1, thetat1)
#        return epst1
    def pmm(self, t, s):
        '''p*(t|s) = \sum_kj eps_kj p_j(t|s) theta_k --- mixture of experts'''
        beta, tvocs, OneByVt, E, T, eps, theta = \
            self.beta, self.tvocs, self.OneByVt, self.E, self.T, self.eps, self.theta
        # for experts where both s and t are in the vocab
        pr = sum( [beta[j][s][t] * eps[k][j] * theta[k] \
                   for k in range(T) for j in range(E) \
                   if s in beta[j] and t in beta[j][s] ] )
        # for experts where t is there, but s is not there
        pr = pr + sum( [OneByVt[j] * eps[k][j] * theta[k] \
                        for k in range(T) for j in range(E)\
                        if s not in beta[j] and t in tvocs[j] ] )
        # for experts where t is not there, beta_jst = 0, hence ignore those experts.
        return pr

class MixModel():
    '''Class for computing p*(t|s) using p(t|s) and p_a(t|s) for different 'a's.'''
    def __init__(self, sl, tl, st_pmffiles):
        '''Input: source lang, target lang, p(t|s) pmf file names, training lexicon'''
        self.sl, self.tl = sl, tl
        self.beta = [] # p(t|s) pmfs for different experts
        self.tvocs = [] # target vocabularies for different experts
        self.OneByVt = [] # 1/V_t for different experts
        for fil in st_pmffiles:
            self.beta.append(L1L2PMF(sl, tl, fil).pmf)
            self.tvocs.append( set([ t for s in self.beta[-1].iterkeys()
                                    for t in self.beta[-1][s].iterkeys() ]) )
            self.OneByVt.append( 1.0/len(self.tvocs[-1]) )
        self.E = len(self.beta) # number of experts
        self.eps = [1.0/self.E for i in range(self.E)] # default expert mixture
    def saveeps(self, epsfile):
        print >> open(epsfile, 'w'), ' '.join([str(pr) for pr in self.eps])
    def readeps(self, epsfile):
        self.eps = [float(pr) for pr in open(epsfile).read().split()]
    def settrdata(self, trfile, l1, l2):
        self.trlex = Lexicon(l1, l2, trfile)
        sl, trl = self.sl, self.trlex.lex
        self.trdata = [ (s,t) for s in trl[sl].iterkeys() for t in trl[sl][s].iterkeys() ]
    def learnepsgrid(self, tlvoc):
        # Choose a grid of param values to try
        # For each param value
        #   Set eps to that point
        #   Compute cands for tr data instances
        #   Compute MRR for cand set
        #   If MRR best so far, save eps
        self.tlvoc = tlvoc
        self.ev = Eval(self.trlex, len(tlvoc))
        self.wtrange = [0., .1, .2, .3, .4, .5, .6, .7, .8, .9, 1.]
        self.bestmrr = 0
        self.besteps = None
        self.gridsearch([], 0, self.E)
        self.eps = self.besteps
    def gridsearch(self, curreps, tot, coordsleft):
        if coordsleft == 1:
            curreps.append(1.-tot)
            mrr = self.computemrr(curreps)
            if mrr > self.bestmrr:
                self.bestmrr = mrr
                self.besteps = list(curreps)
            curreps.pop()
        else:
            for wt in self.wtrange:
                if tot + wt <= 1.:
                    curreps.append(wt)
                    self.gridsearch(curreps, tot+wt, coordsleft-1)
                    curreps.pop()
                else:
                    break
    def computemrr(self, eps):
        rrs = []
        for s, tgt in self.trdata:
            # compute cands
            cands = [(t, self.pmm(t,s,eps)) for t in self.tlvoc]
            cands.sort(key=lambda x: x[1], reverse=True)
            # compute rr for s,t
            rr, presat1, presat5 = self.ev.evalword(self.sl, self.tl, s, cands)
            rrs.append(rr)
        # compute mrr
        return sum(rrs)/float(len(rrs))
    def learnepsEM(self, eps0=None):
        sl, beta, E, trdata = self.sl, self.beta, self.E, self.trdata
        #trdata = [ (s,t) for s in trl[sl].iterkeys() for t in trl[sl][s].iterkeys() ]
        N = len(trdata)
        # starting point for optimization
        # uniform init
        #if not eps0: eps0 = [1.0/E for i in range(E)]
        # random init
        #if not eps0: eps0 = MyUtils.dir_sample(E)
        bestlhood = float("-inf")
        besteps = None
        for ntries in range(10):
    #        # random init
            eps0 = MyUtils.dir_sample(E)

            b = [ [0 for j in range(E)] for n in range(N) ]
            MINCHANGE = 1e-3
            epst1 = eps0
            while True:
                # set current value of eps
                epst = epst1
                # compute the b's
                for n in range(N):
                    s, t = trdata[n]
                    tot = sum( [ beta[j][s][t] * epst[j] for j in range(E) if t in beta[j][s] ] )
                    for j in range(E):
                        b[n][j] = (beta[j][s][t] * epst[j])/tot if t in beta[j][s] else 0
    #                b[n] = [ (beta[j][s][t] * epst[j])/tot if t in beta[j][s] else 0 for j in range(E) ]
                # compute new value of eps
                tot = sum( [sum(bn) for bn in b] )
                epst1 = [ sum([bn[k] for bn in b])/tot for k in range(E) ]
                # check for convergence
                if MyUtils.pnorm( [t-t1 for t, t1 in zip(epst,epst1)], 2 ) < MINCHANGE:
                    break
            # keep this epsilon if it is the best so far
            currlhood = self.loglhood(trdata, epst1)
            if currlhood > bestlhood:
                bestlhood = currlhood
                besteps = epst1
        #self.eps = epst1
        self.eps = besteps
        #return epst1
        return besteps
    def loglhood(self, trdata, eps):
        E, beta = self.E, self.beta
        #tot = 0
        #for n in range(N):
            #tot1 = 0
            #for j in range(E):
                #if t in beta[j][s]:
                    #tot1 += eps[j] * beta[j][s][t]
            #tot += math.log(tot1)
        #return tot
        try:
            return sum( [ math.log( sum( [ beta[j][s][t] * eps[j]
                                       for j in range(E) if t in beta[j][s] ] ) ) \
                      for s,t, in trdata ] )
        except ValueError:
            return -99999;
    def pmm(self, t, s, eps=None):
        '''p*(t|s) = \sum_j eps_j p_j(t|s) --- mixture of experts'''
        beta, tvocs, OneByVt, E = self.beta, self.tvocs, self.OneByVt, self.E
        if not eps: eps = self.eps
        # for experts where both s and t are in the vocab
        pr = sum( [ eps[j] * beta[j][s][t] for j in range(E) \
                     if s in beta[j] and t in beta[j][s] ] )
        # for experts where t is there, but s is not there
        pr = pr + sum( [ eps[j] * OneByVt[j] for j in range(E) \
                     if s not in beta[j] and t in tvocs[j] ] )
        return pr

class Paux():
    '''Class for computing p_a(t|s) using third/auxiliary language'''
    def __init__(self, sl, tl, al, at_pmffile, sa_pmffile, at_tvocfile, as_avocfile):
        '''Input: source lang, target lang, auxiliary lang, p(t|a), p(a|s), |Vat_t|, |Vas_a|.
        '''
        self.sl, self.tl, self.al = sl, tl, al
        self.PtGa = L1L2PMF(al, tl, at_pmffile) #p(t|a)
        self.PaGs = L1L2PMF(sl, al, sa_pmffile) #p(a|s)
        # vocab size of target in at_pmf
        self.Vat_t = set( open(at_tvocfile).read().decode('utf-8').split() )
        self.OneByVat_t = 1/float(len(self.Vat_t))
        self.Vas_a = set( open(as_avocfile).read().decode('utf-8').split() )
        self.OneByVas_a = 1/float(len(self.Vas_a))
    def p_a(self, t, s):
        '''p_a(t|s) = \sum_a p(t|a) * p(a|s)'''
        if t not in self.Vat_t: return 0

        PtGa, PaGs, OneByVat_t = self.PtGa.pmf, self.PaGs.pmf, self.OneByVat_t
        if s in PaGs:
            pr = sum( [ PtGa[a][t] * PaGs[s][a] for a in PaGs[s] \
                                if a in PtGa and t in PtGa[a] ] )
            pr = pr + OneByVat_t * sum( [ PaGs[s][a] for a in PaGs[s] \
                                if a not in PtGa ] )
        else: # dont do this, ignore such s (done below in getpaux())
            OneByVas_a, Vas_a = self.OneByVas_a, self.Vas_a
            pr = sum( [ PtGa[a][t] for a in Vas_a if a in PtGa and t in PtGa[a] ] )
            pr =  pr + OneByVat_t * sum( [ 1 for a in Vas_a if a not in PtGa ] )
            pr = OneByVas_a * pr
        return pr

def getpaux(sl, qfile, tl, tlvocfile, al, at_pmffile, sa_pmffile, at_tvocfile, as_avocfile,
            truncprob, candfile):
    pa = Paux(sl, tl, al, at_pmffile, sa_pmffile, at_tvocfile, as_avocfile)
    tlvoc = open(tlvocfile).read().decode('utf-8').split()
    f = gzip.open(candfile, 'wb')
    for s in open(qfile).read().decode('utf-8').split():
        if s in pa.PaGs.pmf: #ignore source words without data
#            start = time.clock()
            cands = [(t, pa.p_a(t,s)) for t in tlvoc]
            cands.sort(key=lambda x: x[1], reverse=True)
            cands = L1L2PMF.truncate_pmf(cands, truncprob, renormalize=True)
            line = s + '\t' + '\t'.join( [t+' '+str(pr) for t, pr in cands] )
            print >> f, line.encode('utf-8')
#            print >> sys.stderr, "%d secs, %d words" % (time.clock() - start, len(cands))
    f.close()

#def getpmm(sl, goldfile, colnum, tl, tlvocfile, st_pmffile, sta_pmffile,
#           trlexfile, l1, l2, truncprob, candfile, epsfile):

#def getpmm(sl, tl, tlvocfile, st_pmffile, sta_pmffile,
#           tr_ts_cand_lat_filelist, truncprob):

def getpmm(sl, tl, tlvocfile, st_pmffiles,
           tr_ts_cand_lat_filelist, truncprob, learnmeth):
#    mm = MixModel(sl, tl, [st_pmffile, sta_pmffile])
    # st_pmffiles = comma separated list of pmf files, starting with
    #               the base dist, followed by aux lang dist
    mm = MixModel(sl, tl, st_pmffiles.split(','))

    tlvoc = open(tlvocfile).read().decode('utf-8').split()
    for line in open(tr_ts_cand_lat_filelist):
        line = line.rstrip()
        tsfile, colnum, trfile, l1, l2, candfile, epsfile = line.split()
        # trfile = (s,t) pairs such that they belong to vocab of st, at, and as corpora
        #          (so that impact of bringing in a can be measured accurately)
        # tsfile = should not contain words present in trfile
        colnum = int(colnum)

        mm.settrdata(trfile, l1, l2)
        if learnmeth == "EM":
            mm.learnepsEM()
        elif learnmeth == "GRID":
            mm.learnepsgrid(tlvoc)
        else:
            mm.learnepsgrid(tlvoc)
        mm.saveeps(epsfile)

        f = gzip.open(candfile, 'wb')

        for s in [line.decode('utf-8').strip().split()[colnum-1]
                  for line in open(tsfile).readlines()]:
    #        start = time.clock()
            cands = [(t, mm.pmm(t,s)) for t in tlvoc]
            cands.sort(key=lambda x: x[1], reverse=True)
            cands = L1L2PMF.truncate_pmf(cands, truncprob, renormalize=True)
            line = s + u'\t' + u'\t'.join( [t+u' '+unicode(pr) for t, pr in cands] )
            print >> f, line.encode('utf-8')
    #        print >> sys.stderr, "%d secs, %d words" % (time.clock() - start, len(cands))
        f.close()

def getpmmfine(sl, tl, tlvocfile, st_pmffile, sta_pmffile, tr_ts_cand_lat_filelist,
               st_ptmfile, truncprob):
    stptm = PTM.load(st_ptmfile)
    scd = SourceCatDist()
    scd.set_from_PTM(stptm, sl)
    mmf = MixModelFine(sl, tl, [st_pmffile, sta_pmffile], scd)
    tlvoc = open(tlvocfile).read().decode('utf-8').split()
    for line in open(tr_ts_cand_lat_filelist):
        line = line.rstrip()
        tsfile, colnum, trfile, l1, l2, candfile, latentfile = line.split()
        # trfile = (s,t) pairs such that they belong to vocab of st, at, and as corpora
        #          (so that impact of bringing in a can be measured accurately)
        # tsfile = should not contain words present in trfile
        colnum = int(colnum)

        mmf.settrdata(trfile, l1, l2)
        mmf.learnEM()
        mmf.savelatent(latentfile)

        f = gzip.open(candfile, 'wb')

        for s in [line.decode('utf-8').strip().split()[colnum-1]
                  for line in open(tsfile).readlines()]:
    #        start = time.clock()
            cands = [(t, mmf.pmm(t,s)) for t in tlvoc]
            cands.sort(key=lambda x: x[1], reverse=True)
            cands = L1L2PMF.truncate_pmf(cands, truncprob, renormalize=True)
            line = s + u'\t' + u'\t'.join( [t+u' '+unicode(pr) for t, pr in cands] )
            print >> f, line.encode('utf-8')
    #        print >> sys.stderr, "%d secs, %d words" % (time.clock() - start, len(cands))
        f.close()

def getepsEM(sl, tl, st_pmffile, sta_pmffile, trfile, l1, l2, epsfile):
    mm = MixModel(sl, tl, [st_pmffile, sta_pmffile], trfile, l1, l2)
    mm.learnepsEM()
    mm.saveeps(epsfile)

def getpmm4eps(sl, tl, tlvocfile, st_pmffiles, epsfile,
           tr_ts_cand_lat_filelist, tr_or_ts, truncprob):
#    mm = MixModel(sl, tl, [st_pmffile, sta_pmffile])
    # st_pmffiles = comma separated list of pmf files, starting with
    #               the base dist, followed by aux lang dist
    mm = MixModel(sl, tl, st_pmffiles.split(','))
    mm.readeps(epsfile)

    tlvoc = open(tlvocfile).read().decode('utf-8').split()
    for line in open(tr_ts_cand_lat_filelist):
        line = line.rstrip()
        tsfile, colnum, trfile, l1, l2, candfile, epsfile1 = line.split()
        # trfile = (s,t) pairs such that they belong to vocab of st, at, and as corpora
        #          (so that impact of bringing in a can be measured accurately)
        # tsfile = should not contain words present in trfile
        colnum = int(colnum)

        f = gzip.open(candfile, 'wb')
        qfile = tsfile if tr_or_ts == "TEST" else trfile
        for s in [line.decode('utf-8').strip().split()[colnum-1]
                  for line in open(qfile).readlines()]:
            cands = [(t, mm.pmm(t,s)) for t in tlvoc]
            cands.sort(key=lambda x: x[1], reverse=True)
            cands = L1L2PMF.truncate_pmf(cands, truncprob, renormalize=True)
            line = s + u'\t' + u'\t'.join( [t+u' '+unicode(pr) for t, pr in cands] )
            print >> f, line.encode('utf-8')
        f.close()


def main():
    cmd, args = sys.argv[1], sys.argv[2:]
    if cmd == "getpaux":
        # qfile = sa_svoc (or smaller), tlvocfile = ta_tvoc
        sl, qfile, tl, tlvocfile, al, at_pmffile, sa_pmffile, \
            at_tvocfile, as_avocfile, truncprob, candfile = args
        getpaux(sl, qfile, tl, tlvocfile, al, at_pmffile, sa_pmffile,
                at_tvocfile, as_avocfile, float(truncprob), candfile)
    elif cmd == "getpmm":
#        sl, tl, tlvocfile, st_pmffile, sta_pmffile, \
#            tr_ts_cand_lat_filelist, truncprob = args
#        getpmm(sl, tl, tlvocfile, st_pmffile, sta_pmffile,
#               tr_ts_cand_lat_filelist, float(truncprob))
        sl, tl, tlvocfile, st_pmffiles, \
            tr_ts_cand_lat_filelist, truncprob, learnmeth = args
        getpmm(sl, tl, tlvocfile, st_pmffiles,
               tr_ts_cand_lat_filelist, float(truncprob), learnmeth)
    elif cmd == "getpmmfine":
        sl, tl, tlvocfile, st_pmffile, sta_pmffile, \
            tr_ts_cand_lat_filelist, st_ptmfile, truncprob = args
        getpmmfine(sl, tl, tlvocfile, st_pmffile, sta_pmffile,
                   tr_ts_cand_lat_filelist, st_ptmfile, float(truncprob))
    elif cmd == "getepsEM":
        sl, tl, st_pmffile, sta_pmffile, \
            trfile, l1, l2, epsfile = args
        getepsEM(sl, tl, st_pmffile, sta_pmffile, trfile, l1, l2, epsfile)
    elif cmd == "getpmm4eps":
        sl, tl, tlvocfile, st_pmffiles, epsfile, \
           tr_ts_cand_lat_filelist, tr_or_ts, truncprob = args
        getpmm4eps(sl, tl, tlvocfile, st_pmffiles, epsfile,
           tr_ts_cand_lat_filelist, tr_or_ts, float(truncprob))
    else:
        print "invalid command: ", cmd

if __name__ == "__main__":
    main()
