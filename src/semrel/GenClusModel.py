import sys
from PTM_new import PTM
from CorpusMinimal import Corpus
from Common import L1L2PMF, Lexicon
import numpy as np
np.seterr(all='raise')
import math
import gzip
class MixModelClus:
    def __init__(self, sl, tl, slvecfile, st_pmffiles):
        '''Input: source lang, target lang, source word feature vectors,
         p(t|s) pmf file names, and number of clusters K.'''
        self.sl, self.tl = sl, tl
        self.slfeat = WordFeat.loadfeat(slvecfile)
        self.slfeat = self.scalefeat(self.slfeat)
        self.beta = [] # p(t|s) pmfs for different experts
        self.tvocs = [] # target vocabularies for different experts
        self.OneByVt = [] # 1/V_t for different experts
        for fil in st_pmffiles:
            self.beta.append(L1L2PMF(sl, tl, fil).pmf)
            self.tvocs.append( set([ t for s in self.beta[-1].iterkeys()
                                    for t in self.beta[-1][s].iterkeys() ]) )
            self.OneByVt.append( 1.0/len(self.tvocs[-1]) )
        self.E = len(self.beta) # number of experts
    def scalefeat(self, slfeat):
        wds = []
        feat = []
        for w, f in slfeat.iteritems():
            wds.append(w)
            feat.append(f)
        feat = np.array(feat)
        #normalize first 6 features
        feat[:,0:6] = feat[:,0:6]  / np.amax(feat[:,0:6], axis=0)
        #feat[:,0:6] = feat[:,0:6]  * np.array([1, 0, 1, 1, 0, 1])
        #for f in feat:
            #f[0] = 0 #math.log(f[0])
            #f[1] = 0
            #f[3] = 0
            #f[4] = 0
            #f[6:] = [0 for i in range(len(f)-6)]
        for i, w in enumerate(wds):
            slfeat[w] = feat[i]
            #slfeat[w] = np.array(feat[i])
        return slfeat
    def settrdata(self, trfile, l1, l2):
        self.trlex = Lexicon(l1, l2, trfile)
    @staticmethod
    def init_mu(D, K):
        '''Initialize K mu's, each of size D.'''
        mus = []
        for k in range(K):
            mu = [0 for d in range(D)]
            mu[k] = 2
            mus.append(mu)
        return np.array(mus)
    def train(self, K, eps0=None, pi0=None, mu0=None, sigma0=None):
        sl, slfeat, beta, E, self.K, trl = self.sl, self.slfeat, self.beta, self.E, K, self.trlex.lex
        trdata = [ (s,t, np.array(slfeat[s]) ) for s in trl[sl].iterkeys() for t in trl[sl][s].iterkeys() ]
        X = np.array( [x for s, t, x in trdata] )
        N = len(trdata)

        ##########initialize parameters
        if not eps0: eps0 = np.array( [ [1.0/E for j in range(E)] for k in range(K) ] )
        if not pi0: pi0 = [1.0/K for k in range(K)]
        NbyK = N/K #rounded
        # take a few data points and take their mean as the mean for each Gaussian
        #if not mu0: mu0 = np.array( [ np.mean( [ x for s,t,x in trdata[NbyK*k:NbyK*(k+1)] ], axis=0 ) \
                                        #for k in range(K)] )
        if not mu0: mu0 = MixModelClus.init_mu(len(X[0]), K)
        # take a variance of 1 for now; TODO: use sample variance of a few data points instead
        if not sigma0: sigma0 = [5000000.0 for k in range(K)]

        ##########EM iterations
        T = np.array( [ [ [0.0 for j in range(E)] for k in range(K)]for n in range(N)] )
        MINCHANGE = 1e-3
        epst1, pit1, mut1, sigmat1 = eps0, pi0, mu0, sigma0
        while True:
            epst, pit, mut, sigmat = epst1, pit1, mut1, sigmat1
            #########compute p(e_n=j,c_n=k|s_n,t_n,x_n,eps0,pi0,mu0) = T_nkj
            for n in range(N):
                s,t,x = trdata[n]
                for k in range(K):
                    diff = x-mut[k]
                    for j in range(E):
                        if t not in beta[j][s]:
                            T[n][k][j] = 0
                        else:
                            T[n][k][j] = beta[j][s][t] * epst[k][j] * pit[k] \
                                * math.exp( - np.dot(diff,diff) / (2*(sigmat[k]**2)) ) \
                                / (sigmat[k])
                                #* ( - np.dot(diff,diff) / (2*(sigmat[k]**2)) - math.log(sigmat[k]) )
                            #T[n][k][j] = math.exp( math.log(beta[j][s][t]) + math.log(epst[k][j]) + math.log(pit[k]) \
                                #+ (- np.dot(diff,diff) / (2*(sigmat[k]**2)) ) \
                                #- math.log( abs(sigmat[k]) ) )
                tot = np.sum(T[n])
                if tot != 0.: T[n] = T[n]/tot
            #########max E [ log p(s,t,x,c,e|eps,pi,mu) ]
            Ttot = np.sum(T)
            Tk = np.sum(np.sum(T, axis=0), axis=1) # \sum_n,j T_nkj
            epst1 = np.sum(T, axis=0) /  Tk[:,np.newaxis]
            pit1 = Tk/N
            Tnk = np.sum(T, axis=2)
            for k in range(K):
                mut1[k] = np.sum(X * Tnk[:,k][:,np.newaxis], axis=0) / Tk[k]
                diff = X - mut1[k]
                sigmat1[k] = math.sqrt( np.sum(np.diag(np.dot(diff, diff.T)) * Tnk[:,k]) / Tk[k] )
            #########check for convergence
            if np.linalg.norm(np.subtract(epst, epst1)) < MINCHANGE: # Frobenius norm
                break
        self.eps, self.pi, self.mu, self.sigma, self.N = epst1, pit1, mut1, sigmat1, N
    def save(self, paramfile):
        N, K, E, eps, pi, mu, sigma = self.N, self.K, self.E, self.eps, self.pi, self.mu, self.sigma
        pf = open(paramfile, 'w')
        # Format:
        # N=.. K=.. E=..
        # EPS
        # eps_1.
        # eps_2.
        # ...
        # PI
        # pi_1 pi_2 ...
        # MU
        # mu_1.
        # mu_2.
        # ...
        # SIGMA
        # sigma_1 sigma_2 ...
        print >> pf, "N=%d\tK=%d\tE=%d" % (N, K, E)
        print >> pf, "EPS"
        for epsk in eps:
            print >> pf, ' '.join( [str(pr) for pr in epsk] )
        print >> pf, "PI"
        print >> pf, ' '.join( [str(pr) for pr in pi] )
        print >> pf, "MU"
        for muk in mu:
            print >> pf, ' '.join( [str(mukd) for mukd in muk] )
        print >> pf, "SIGMA"
        print >> pf, ' '.join( [str(sk) for sk in sigma] )
        pf.close()
    def pmm(self, t, s):
        '''p(t|s,x)'''
        beta, eps, pi, mu, sigma, K, E, tvocs, OneByVt = self.beta, self.eps, \
            self.pi, self.mu, self.sigma, self.K, self.E, self.tvocs, self.OneByVt
        x = self.slfeat[s]
        tot = 0.
        for k in range(K):
            diff = x - mu[k]
            tempk = pi[k] * math.exp( - np.dot(diff,diff) / (2*(sigma[k]**2)) ) \
                    / sigma[k]
            #tempk = pi[k] * ( - np.dot(diff,diff) / (2*(sigma[k]**2)) - math.log(sigma[k]))
            #tempk = math.exp( math.log(pi[k]) + ( - np.dot(diff,diff) / (2*(sigma[k]**2)) ) \
                    #-  math.log(sigma[k]) )
            tot1 = 0.
            for j in range(E):
                if s in beta[j]:
                    # for experts where both s and t are in the vocab
                    if t in beta[j][s]: pr = beta[j][s][t]
                    else: pr = 0.
                else:
                    # for experts where t is there, but s is not there
                    if t in tvocs[j]: pr = OneByVt[j]
                    else: pr = 0.
                tot1 += pr * eps[k][j]
            try:
                tot += tempk * tot1
            except FloatingPointError:
                pass
        return tot

def getpmmc(sl, tl, slvecfile, tlvocfile, st_pmffiles, K,
           tr_ts_cand_par_filelist, truncprob):
    mmc = MixModelClus(sl, tl, slvecfile, st_pmffiles.split(','))
    tlvoc = open(tlvocfile).read().decode('utf-8').split()
    for line in open(tr_ts_cand_par_filelist):
        line = line.rstrip()
        tsfile, colnum, trfile, l1, l2, candfile, paramfile = line.split()
        # trfile = (s,t) pairs such that they belong to vocab of st, at, and as corpora
        #          (so that impact of bringing in `a' can be measured accurately)
        # tsfile = should not contain words present in trfile
        colnum = int(colnum) # column in tsfile, from which to read test words

        mmc.settrdata(trfile, l1, l2)
        mmc.train(K)
        mmc.save(paramfile)

        f = gzip.open(candfile, 'wb')
        for s in [line.decode('utf-8').strip().split()[colnum-1]
                  for line in open(tsfile).readlines()]:
    #        start = time.clock()
            cands = [(t, mmc.pmm(t,s)) for t in tlvoc]
            cands.sort(key=lambda x: x[1], reverse=True)
            cands = L1L2PMF.truncate_pmf(cands, truncprob, renormalize=True)
            line = s + u'\t' + u'\t'.join( [t+u' '+unicode(pr) for t, pr in cands] )
            print >> f, line.encode('utf-8')
    #        print >> sys.stderr, "%d secs, %d words" % (time.clock() - start, len(cands))
        f.close()

def trainmmc(sl, tl, slvecfile, st_pmffiles, K, trfile, l1, l2, paramfile):
    mmc = MixModelClus(sl, tl, slvecfile, st_pmffiles.split(','))
    # trfile = (s,t) pairs such that they belong to vocab of st, at, and as corpora
    #          (so that impact of bringing in `a' can be measured accurately)
    mmc.settrdata(trfile, l1, l2)
    mmc.train(K)
    mmc.save(paramfile)


class WordFeat:
    def __init__(self, l, ptmfile, wstatfile):
        self.l = l
        self.ptm = PTM.load(ptmfile)
        self.wstat = Corpus.loadwstat(wstatfile)
    def getfeat(self, w):
        l, ptm, wstat = self.l, self.ptm, self.wstat
        vec = []
        # stats for w
        if w in wstat:
            vec.extend( self.wstat[w] )
        else:
            return None
        # theta for w
        norm = sum( [ptm.phi[l][t][w] for t in ptm.topics] )
        vec.extend( [unicode(ptm.phi[l][t][w]/norm) for t in ptm.topics] )
        return vec
    @staticmethod
    def printfeat(w, vec):
        print ( w+u'\t'+u' '.join(vec) ).encode('utf-8')
    @staticmethod
    def loadfeat(featfile):
        feat = {}
        for line in open(featfile):
            line = line.rstrip().decode('utf-8')
            w, x = line.split(u'\t')
            x = [ float(xi) for xi in x.split(u' ') ]
            feat[w] = x
        return feat

def getwordfeat(l, vocfile, ptmfile, wstatfile):
    wf = WordFeat(l,ptmfile, wstatfile)
    for w in open(vocfile).read().decode('utf-8').split():
        vec = wf.getfeat(w)
        if vec:
            WordFeat.printfeat(w, vec)

def main():
    cmd, args = sys.argv[1], sys.argv[2:]
    if cmd == "getwordfeat":
        l, vocfile, ptmfile, wstatfile = args
        getwordfeat(l, vocfile, ptmfile, wstatfile)
    elif cmd == "getpmmc":
        sl, tl, slvecfile, tlvocfile, st_pmffiles, K, tr_ts_cand_par_filelist, truncprob = args
        getpmmc(sl, tl, slvecfile, tlvocfile, st_pmffiles, int(K), tr_ts_cand_par_filelist, float(truncprob))
    elif cmd == "trainmmc":
        sl, tl, slvecfile, st_pmffiles, K, trfile, l1, l2, paramfile = args
        trainmmc(sl, tl, slvecfile, st_pmffiles, int(K), trfile, l1, l2, paramfile)
    else:
        print "invalid command: ", cmd

if __name__ == "__main__":
    main()
