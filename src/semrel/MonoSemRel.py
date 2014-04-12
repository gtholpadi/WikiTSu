'''
Created on 07-Mar-2013
Collects algorithms for monolingual semantic relatedness.
@author: gtholpadi
'''

'''
Input : corpus, query list (words must exist in corpus)
Output: for each word in query list, compute pmi with entire
vocab, rank by pmi, normalize and return cands.
'''
from CorpusMinimal import Corpus
from Adhoc import readwordlistaslist
from Common import L1L2PMF
import math
import gzip
import sys
class ESA:
    WINDOWSIZE = 100 #window size for postprocessing
    MINDIFF = .05 #minimum difference in window as fraction of top score 
    def __init__(self, corp):
        self.corp = corp
        corp.settfidf()
        self.tii = corp.tii
        self.postprocess()
    def postprocess(self):
        tii = self.tii
        for w in tii:
            # sort decreasing by score, move a window of size WINDOWSIZE 
            # till diff between first and last in window falls below MINDIFF
            # fraction of top score.
            cvec = sorted(tii[w].iteritems(), key=lambda x:x[1], reverse=True)
#            cvec.sort(key=lambda x:x[1], reverse=True)
            maxoffset = len(cvec)-ESA.WINDOWSIZE
            maxoffset = maxoffset if maxoffset > 0 else 0
            i = 0
            for i in range(1, maxoffset):
                diff = cvec[i][1]-cvec[i+ESA.WINDOWSIZE][1]
                if diff/cvec[0][1] < ESA.MINDIFF: break
            if i:
                cvec = cvec[:i+ESA.WINDOWSIZE-1] #leave out the last one in window
            #l2-normalize and assign to tfidf index
            tot = math.sqrt( sum([sc*sc for d, sc in cvec]) )
            tii[w] = dict( [(d, sc/tot) for d, sc in cvec] )
    def sim(self, w1, w2):
        # use cosine similarity
        tii = self.tii
        tot = 0
        for d in tii[w1]:
            if d in tii[w2]:
                tot += tii[w1][d] * tii[w2][d]
        return tot
    def getrel(self, w):
        c = self.corp
        cands = [(x, self.sim(w, x)) for x in c.vocab]# if x != w]
        cands.sort(key=lambda x: x[1], reverse=True)
        return cands
class PMI:
    def __init__(self, corp):
        self.corp = corp #corpus object
        self.ndocs = corp.ndocs
        corp.setwdf()
        self.wdf = corp.wdf #for each word, its document frequency
    def pmi(self, w1, w2):
        # FORMULAS USED:
        # pmi(w1,w2) = log [ p(w1,w2) / p(w1)p(w2) ]
        # p(w1) = #docs with w1 / #docs
        # p(w2,w2) = #docs with w1 and w2 / #docs
        # pmi(w1,w2) = log [ (#docs with w1 and w2 * #docs) / (#docs with w1 * #docs with w2) ]
        c, wdf, ndocs = self.corp, self.wdf, self.ndocs

#        # add 1 for smoothing
#        nw1 = wdf[w1] + 1
#        nw2 = wdf[w2] + 1
#        nw1w2 = len(set(c.index[w1].keys()).intersection(set(c.index[w2].keys()))) + 1
#        ndocs += 1
#        frac = (nw1w2 * ndocs) / (nw1 * nw2)
##        return math.log(frac, 2)
#        return math.log(frac, 2) * nw1w2 # to mitigate rare words problem

        # without smoothing
        # npmi(w1,w2) = pmi(w1,w2) / -log p(w1,w2)
        nw1 = wdf[w1]
        nw2 = wdf[w2]
        nw1w2 = len(set(c.index[w1].keys()).intersection(set(c.index[w2].keys()))) + 1e-4
        frac = (nw1w2 * ndocs) / (nw1 * nw2)
        npmi = math.log(frac) #usual
#        npmi = math.log(frac)/-math.log(nw1w2/ndocs) #normalized
#        return npmi
        return npmi * nw1w2 # to mitigate rare words problem
    def getrel(self, w):
        c = self.corp
        cands = [(x, self.pmi(w, x)) for x in c.vocab]# if x != w]
        cands.sort(key=lambda x: x[1], reverse=True)
        return cands

def getpmi(lang, qfile, colnum, corpdir, candfile, truncprob):
    '''Given a list of words (column <colnum> in qfile), get a list of 
    related words in the same language, using a background corpus stored
    in <corpdir>. The list of words is sorted by the pmi score, and 
    l1-normalized (to get a probability distribution) and truncated, and
    then stored in <candfile>.  
    '''
    corp = Corpus(lang)
    corp.readdir(corpdir)
    pmi = PMI(corp)
    getcands4meth(qfile, colnum, candfile, truncprob, pmi.getrel)
#    f = gzip.open(candfile, 'wb')
#    for line in open(qfile):
#        line = line.decode('utf-8').rstrip()
#        w1 = line.split()[colnum-1]
#        cands = pmi.getrel(w1)
#        # convert to pmf
#        minsc = min([sc for w, sc in cands])
#        cands = [(w, sc-minsc) for w, sc in cands] # shift all scores to make them positive
#        tot = sum([sc for w, sc in cands])
#        cands = [(w, sc/tot) for w, sc in cands]
#        # truncate and renormalize pmf
#        cands = L1L2PMF.truncate_pmf(cands, truncprob, renormalize=True)
#        line = w1 + u'\t' + u'\t'.join( [w+u' '+unicode(sc) for w, sc in cands] )
#        print >> f, line.encode('utf-8')
#    f.close()
def getesa(lang, qfile, colnum, corpdir, candfile, truncprob):
    '''Given a list of words (column <colnum> in qfile), get a list of 
    related words in the same language, using a background corpus stored
    in <corpdir>. The list of words is sorted by the esa score, and 
    l1-normalized (to get a probability distribution) and truncated, and
    then stored in <candfile>.  
    '''
    corp = Corpus(lang)
    corp.readdir(corpdir)
    esa = ESA(corp)
    getcands4meth(qfile, colnum, candfile, truncprob, esa.getrel)
#    f = gzip.open(candfile, 'wb')
#    for line in open(qfile):
#        line = line.decode('utf-8').rstrip()
#        w1 = line.split()[colnum-1]
#        cands = esa.getrel(w1)
#        # convert to pmf
#        minsc = min([sc for w, sc in cands])
#        cands = [(w, sc-minsc) for w, sc in cands] # shift all scores to make them positive
#        tot = sum([sc for w, sc in cands])
#        cands = [(w, sc/tot) for w, sc in cands]
#        # truncate and renormalize pmf
#        cands = L1L2PMF.truncate_pmf(cands, truncprob, renormalize=True)
#        line = w1 + u'\t' + u'\t'.join( [w+u' '+unicode(sc) for w, sc in cands] )
#        print >> f, line.encode('utf-8')
#    f.close()
def getcands4meth(qfile, colnum, candfile, truncprob, getrelmeth):
    f = gzip.open(candfile, 'wb')
    for line in open(qfile):
        line = line.decode('utf-8').rstrip()
        w1 = line.split()[colnum-1]
        cands = getrelmeth(w1)
        # convert to pmf
        minsc = min([sc for w, sc in cands])
        cands = [(w, sc-minsc) for w, sc in cands] # shift all scores to make them positive
        tot = sum([sc for w, sc in cands])
        cands = [(w, sc/tot) for w, sc in cands]
        # truncate and renormalize pmf
        cands = L1L2PMF.truncate_pmf(cands, truncprob, renormalize=True)
        line = w1 + u'\t' + u'\t'.join( [w+u' '+unicode(sc) for w, sc in cands] )
        print >> f, line.encode('utf-8')
    f.close()
def main():
    cmd, args = sys.argv[1], sys.argv[2:]
    if cmd == "getpmi":
        lang, qfile, colnum, corpdir, candfile, truncprob = args
        getpmi(lang, qfile, int(colnum), corpdir, candfile, float(truncprob))
    elif cmd == "getesa":
        lang, qfile, colnum, corpdir, candfile, truncprob = args
        getesa(lang, qfile, int(colnum), corpdir, candfile, float(truncprob))
    else:
        print "invalid command:", cmd
if __name__ == "__main__":
    main()