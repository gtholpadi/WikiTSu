'''
Created on 17-Jan-2013

@author: gtholpadi
'''
import cPickle as pickle
import sys
import gzip

class PTM:
    '''Polylingual topic model
    '''
    def __init__(self):
        '''
        Input: State file given by mallet
        Output: Pickled PTM object
        '''
        self.beta = 0
        self.langs = set([])
        self.vocabs = {} #indexed by language
        self.topics = set([])
        self.tw = {} # topic word counts; indexed by language, topic and word
        self.phi = {} # p(w|t); indexed by language, topic and word
    def save(self, ptmfile):
        f = gzip.open(ptmfile, 'wb')
        pickle.dump(self, f)
        f.close()
    @staticmethod
    def load(ptmfile):
        f = gzip.open(ptmfile)
        ptm = pickle.load(f)
        f.close()
        return ptm
    def loadstatefile(self, stfile, id2lang, beta=0.01):
        self.beta = beta
        langs, topics, vocabs, tw, phi = self.langs, self.topics, self.vocabs, self.tw, self.phi
        f = gzip.open(stfile, 'rb')
        f.readline() # skip first line
        for line in f:
            doc, lid, pos, wid, w, t = line.decode('utf-8').split()
            l = id2lang[lid]
            # metadata
            if l not in langs:
                langs.add(l)
                vocabs[l] = set([])
            if w not in vocabs[l]: vocabs[l].add(w)
            if t not in topics: topics.add(t)
            # data
            if l not in tw: tw[l] = {}
            if t not in tw[l]: tw[l][t] = {}
            if w not in tw[l][t]: tw[l][t][w] = 0
            tw[l][t][w] += 1
        f.close()
        # set phi
        for l in langs:
            phi[l] = {}
            for t in topics:
                phi[l][t] = {}
                denom = sum(tw[l][t].itervalues()) + len(vocabs[l]) * beta # denom = \sum_w n(l,t,w) + V_l * beta
                for w in vocabs[l]:
                    if w in tw[l][t]: phi[l][t][w] = (tw[l][t][w] + beta) / denom
                    else : phi[l][t][w] = beta / denom

def main():
    stfile, beta, ptmfile = sys.argv[1:4]
    beta = float(beta)
    id2lang = sys.argv[4:]  # e.g. 0:en 1:hi 2:bn ...
    id2lang = dict( [x.split(':') for x in id2lang] )

    ptm = PTM()
    ptm.loadstatefile(stfile, id2lang, beta)
    f = gzip.open(ptmfile, 'wb')
    pickle.dump(ptm, f)
    f.close()


if __name__ == '__main__':
    main()
