#!/usr/bin/python
# -*- coding: utf-8 -*-
from PTM_new import PTM
import sys
import os
import random
from Common import *
import time
import gzip
from lib2to3.pgen2.token import STAR
from CorpusMinimal import Corpus
import numpy

def createdict(d1file, l11, l12, d2file, l21, l22):
    '''Create a new dictionary given two existing dictionaries,
    by pivoting on a common language.
    '''
    ls1 = set([l11, l12])
    ls2 = set([l21, l22])
    l = (ls1 & ls2).pop()
    l1 = (ls1-ls2).pop()
    l2 = (ls2-ls1).pop()

    lex1 = Lexicon(l11, l12, d1file)
    lex2 = Lexicon(l21, l22, d2file)
    lex12 = Lexicon(l1, l2)

    for w1 in lex1.lex[l1]:
        for w in lex1.lex[l1][w1]:
            if w not in lex2.lex[l]: continue
            for w2 in lex2.lex[l][w]:
                lex12.add(l1, w1, l2, w2)
    lex12.printlex(l1)

def readwordlist(wlfile):
    '''Read a list of words (one per line) from a file.
    '''
    wl = {}
    for line in open(wlfile):
        w = line.decode('utf-8').rstrip()
        wl[w] = True
    return wl

def readwordlistaslist(wlfile):
    '''
    Input: a file containing a list of words, one per line (ignore
    empty and '#'-prefixed lines)
    Output: a list containing the words
    '''
    wl = []
    for line in open(wlfile):
        line = line.decode('utf-8').strip()
        if line and not line.startswith('#'): wl.append(line)
    return wl

def readnamevaluelist(nvlfile, NVSEP='\t'):
    '''Read name-value pairs, one per line
    '''
    nv = {}
    f = open(nvlfile)
    for line in f:
        line = line.decode('utf-8').rstrip()
        if line.startswith('#'): continue
        rec = line.split(NVSEP)
        if len(rec) == 2:
            name, value = rec
            nv[name] = value
    f.close()
    return nv

def dictstats(dfile, l1, l2, l1vocfile, l2vocfile, l1qfile, l2qfile):
    '''Compute dictionary statistics; use this to see how many 1-1 maps are there.
    '''
    lex = Lexicon(l1, l2, dfile)
    voc1 = readwordlist(l1vocfile)
    voc2 = readwordlist(l2vocfile)
    q1 = readwordlist(l1qfile)
    q2 = readwordlist(l2qfile)
    print l1
    for w1 in q1:
        ct = 0
        for w2 in lex.lex[l1][w1]:
            if w2 in voc2:
                ct += 1
        line = "%s\t%d" % (w1, ct)
        print line.encode('utf-8')
    print l2
    for w2 in q2:
        ct = 0
        for w1 in lex.lex[l2][w2]:
            if w1 in voc1:
                ct += 1
        line = "%s\t%d" % (w2, ct)
        print line.encode('utf-8')

def trunccands(candfile, truncprob, maxperline, newcandfile):
    '''Read a candidate file, and truncate it to within a specified probability mass,
    and/or to within a max number of candidates per word.
    '''
    ncf = gzip.open(newcandfile, 'wb')
    for line in gzip.open(candfile):
        line = line.decode('utf-8').rstrip()
        cands = line.split(Sep.w1w2)
        w1 = cands.pop(0)
        cands = [ (w2, float(pr)) for w2pr in cands for w2, pr in [tuple(w2pr.split())] ]
        cands = L1L2PMF.truncate_pmf(cands, truncprob, renormalize=True, maxcands=maxperline) # renormalize
        line = w1 + u'\t' + u'\t'.join( [w2+u' '+unicode(pr) for w2, pr in cands] )
#        tot = 0
#        for i, w2pr in enumerate(cands):
#            w2, pr = w2pr.split(Sep.w2wt)
#            pr = float(pr)
#            if pr < 1e-6:
#                i -= 1
#                break
#            tot += pr
#            if tot >= truncprob: break
#        cands = cands[:min(i+1, maxperline)]
#        line = w1 + Sep.w1w2 + Sep.w1w2.join(cands)
        print >> ncf, line.encode('utf-8')
    ncf.close()

def filteronquery(qfile, candfile):
    '''Keep only those lines in the candidates file that correspond
    to the words in the query list.
    '''
    f = open(qfile)
    qws = dict([(qw, True) for qw in f.read().decode('utf-8').split()])

    f.close()
    for line in open(candfile):
        line = line.decode('utf-8').rstrip()
        pos = line.find(Sep.w1w2)
        if line[:pos] in qws:
            print line.encode('utf-8')

def docstatfromstate(stfile, l1, l2, id2lang):
    '''Compute document statistics from the state file given by PTM (Mallet).
    '''
    corp = dict( [(lang, {})  for lang in id2lang.values()] )
    f = gzip.open(stfile)
    f.readline()
    for line in f:
        line = line.decode('utf-8').rstrip()
        d, lid, pos, wid, w, t = line.split()
        l = id2lang[lid]
        if d not in corp[l]:
            corp[l][d] = 0
        corp[l][d] += 1
    f.close()
#    for l in corp:
#        print (l + ' corpus').encode('utf-8')
#        print ( u'\n'.join([d+u'\t'+unicode(ct) for d, ct in sorted(corp[l].items(), key=lambda x: x[0])]) ).encode('utf-8')
    ndocs = len(corp[l1])
    assert ndocs == len(corp[l2])
    ntoksl1 = [ ct for d, ct in sorted(corp[l1].iteritems(), key=lambda x:x[0]) ]
    ntoksl2 = [ ct for d, ct in sorted(corp[l2].iteritems(), key=lambda x:x[0]) ]
    ntokspfl1 = [ numpy.mean(ntoksl1), numpy.std(ntoksl1) ]
    ntokspfl2 = [ numpy.mean(ntoksl2), numpy.std(ntoksl2) ]
    ntokratios = [ float(max(ntoksl1[i],ntoksl2[i]))/min(ntoksl1[i],ntoksl2[i])
                      for i in range(len(ntoksl1))]
    ntokratiopfl = [ numpy.mean(ntokratios), numpy.std(ntokratios) ]
#    print "%s%s\t%d\t%.2f$\pm$%.2f\t%.2f$\pm$%.2f\t%.2f$\pm$%.2f" % \
#        (l1, l2, ndocs, ntokspfl1[0], ntokspfl1[1], ntokspfl2[0], ntokspfl2[1],
#         ntokratiopfl[0], ntokratiopfl[1])
    print "%s%s\t%d\t%d$\pm$%d\t%d$\pm$%d\t%.1f$\pm$%.1f" % \
        (l1, l2, ndocs, ntokspfl1[0], ntokspfl1[1], ntokspfl2[0], ntokspfl2[1],
         ntokratiopfl[0], ntokratiopfl[1])

def getcorpfromids(docidfile, col, srcdir, destdir):
    '''Copy only the articles specified in a docidfile to a new directory,
    at the same time renaming them by their line number in the docidfile.
    '''
    col -= 1
    lineno = 1
    for line in open(docidfile):
        rec = line.split()
        docids = rec[col].split(',')
        for docid in docids:
#            cmd = 'cat %s/%s >> %s/%s' % (srcdir, docid, destdir, lineno)#this works
            cmd = 'cp %s/%s %s/%s' % (srcdir, docid, destdir, lineno) #just trying
            os.system(cmd)
        lineno += 1

def filtcorponvocab(vocfile, srcdir, destdir):
    '''Create a new corpus from the old one after filtering out all words
    except those given in vocfile. destdir must already exist.
    '''
    from TokenFilter import KeepWord
    kpw = KeepWord(kpfiles=[vocfile])
    for fname in os.listdir(srcdir):
        f = open(os.path.join(srcdir, fname))
        text = f.read().decode('utf-8')
        f.close()
        f = open(os.path.join(destdir, fname), 'w')
        # assume text cleaned already
        print >> f, ( ' '.join( kpw.filter( text.split() ) ) ).encode('utf-8')
        f.close()

def printcorpwordstat(corpdir, lang):
    '''Compute word statistics for a corpus.
    '''
    corp = Corpus(lang)
    corp.readdir(corpdir)
    corp.printstats()

def printcorpdocstat(corpdir):
    '''Compute document statistics for a corpus.
    '''
    print "#file\tnchars\tntokens\tntypes"
    for fname in os.listdir(corpdir):
        f = open(os.path.join(corpdir, fname))
        text = f.read().decode('utf-8')
        f.close()
        nchar = len(text)
        tokens = text.split()
        ntokens = len(tokens)
        ntypes = len(set(tokens))
        print "%s\t%d\t%d\t%d" % (fname, nchar, ntokens, ntypes)

def getwnsim(envocfile, qfile):
    '''Compute English WordNet similarity.
    Input: en vocab, query words (en).
    Output: for each query word, similarity with each word in en vocab.
    '''
    from nltk.corpus import wordnet as wn
    f = open(envocfile)
    voc = dict([(w, wn.synsets(w)) for w in f.read().decode('utf-8').split()])
    f.close()
    f = open(qfile)
    qws = f.read().decode('utf-8').split()
    f.close()
    for w1 in qws:
        w1sim = []
        if len(voc[w1]) == 0: # word not in wordnet
            w1sim.append( (w1, 1) )
        else:
            for w2 in voc:
                # choose max similarity from any of their sense pairs
                maxsim = 0
                for s1 in voc[w1]:
                    for s2 in voc[w2]:
                        sim = s1.path_similarity(s2)
                        if sim > maxsim: maxsim = sim
                if maxsim > 0:
                    w1sim.append( (w2, maxsim) )
            #normalize
            tot = sum([sim for w2, sim in w1sim])
            w1sim = [(w2, sim/tot) for w2, sim in w1sim]
            w1sim = sorted(w1sim, key=lambda x: x[1], reverse=True)
        print (w1 + Sep.w1w2 + Sep.w1w2.join([w2+Sep.w2wt+str(pr) for w2, pr in w1sim])).encode('utf-8')

def preprocess(docfile, stopfile, minsize, lang):
    '''Preprocessing for clean wiki articles:
    remove stopwords, remove words smaller than minsize characters, remove words
    with foreign scripts
    Input: file to preprocess, stopfile, minsize of allowed words, lang
    '''
    # create filters
    from TokenFilter import StopWord, SmallWord, KeepScript
    stw = StopWord(stfiles=[stopfile]) if stopfile != "NONE" else StopWord()
    smw = SmallWord(minsize=3)
    ksc = KeepScript(lang=lang)
    for line in open(docfile):
        line = line.decode('utf-8').rstrip()
        tokens = line.split()
        tokens = stw.filter(tokens)
        tokens = smw.filter(tokens)
        tokens = ksc.filter(tokens)
        print ' '.join(tokens).encode('utf-8')

def get_t1t2_map(i1i2mapfile, itmap1file, itmap2file, IDTITSEP='\t'):
    '''Read an id1-id2 map and print the corresponding title1-title2 map.
    '''
    itmap1 = readnamevaluelist(itmap1file, IDTITSEP)
    itmap2 = readnamevaluelist(itmap2file, IDTITSEP)
    i1i2map = readnamevaluelist(i1i2mapfile, IDTITSEP)
    for i1, i2 in i1i2map.iteritems():
        if not ( i1 in itmap1 and i2 in itmap2 ): continue
        # restrict to single word titles
        tit1 = itmap1[i1].split()
        tit2 = itmap2[i2].split()
        if not ( len(tit1) == 1 and len(tit2) == 1 ): continue
        print (tit1[0]+'\t'+tit2[0]).encode('utf-8')

def print_useful_words(wordstatsfile, mincf=20, maxdffrac=0.1, maxvocsize=5000, validwordsfile=None):
    '''
    Input: output of printcorpwordstat() above. (a corpus word statistics file)
    Output: a list of words (size <= maxvocsize) where each word satisfies-- cf >= mincf, df <= maxdf.
    '''
    f = open(wordstatsfile)
    numdocs, = [pv.split('=')[1] for pv in f.readline()[1:].split() if pv.startswith('NUM_DOCS')]
    maxdf = int(float(numdocs) * maxdffrac)
    wl = []
    for line in f:
        line = line.decode('utf-8').rstrip()
        if line.startswith('#'): continue
        w, cf, df, adc, idf, cfidf, tcfidf = line.split()
        if int(cf) >= mincf and int(df) <= maxdf:
            wl.append( (w, float(adc)) )
    if len(wl) > maxvocsize:
        wl.sort(key=lambda x: x[1], reverse=True)

    if validwordsfile:
        vwl = set(readwordlistaslist(validwordsfile))
        wl = [(w, adc) for w, adc in wl if w in vwl]

    print (u'\n'.join( [w for w, adc in wl[:maxvocsize]] )).encode('utf-8')

def print_useful_docs(docstatfile, mintypes=10, mintokens=30):
    '''
    Input: output of printcorpdocstat() above. (a corpus document statistics file)
    Output: a list of files where each file satisfies-- ntypes >= mintypes, ntokens >= mintokens.
    '''
    for line in open(docstatfile):
        line = line.decode('utf-8').rstrip()
        if line.startswith('#'): continue
        fname, nchar, ntokens, ntypes = line.split()
        if int(ntypes) >= mintypes and int(ntokens) >= mintokens:
            print fname.encode('utf-8')

def intersect_useful_map(l1docfile, l2docfile, l1l2idmapfile):
    '''Restrict l1-l2 id mapping to ids present in l1.doc and l2.doc (useful docs).
    '''
    l1ids = set( open(l1docfile).read().split() )
    l2ids = set( open(l2docfile).read().split() )
    for line in open(l1l2idmapfile):
        l1id, l2id = line.split()
        if l1id in l1ids and l2id in l2ids:
            print l1id + "\t" + l2id

def get_eval_pairs(l1, l2, l1l2dictfile, l1vocfile, l2vocfile):
    '''Get evaluation word pairs for a language pair.
    Given the dictionary and the vocabularies, get pairs from the
    dictionary that exist in the vocabularies.'''
    l1voc = set(readwordlistaslist(l1vocfile))
    l2voc = set(readwordlistaslist(l2vocfile))
    #for line in open(l1l2dictfile):
        #w1, w2 = line.decode('utf-8').split()
        #if w1 in l1voc and w2 in l2voc:
            #print (w1+"\t"+w2).encode('utf-8')
    lex = Lexicon(l1,l2, l1l2dictfile).lex
    for w1 in lex[l1]:
        if w1 in l1voc:
            for w2 in lex[l1][w1]:
                if w2 in l2voc:
                    print (w1+"\t"+w2).encode('utf-8')

def get_vocab_from_ptm(ptmfile, lang):
    '''Get vocab for a language from the ptm file (pickle-zipped).'''
    print ( '\n'.join( PTM.load(ptmfile).vocabs[lang] ) ).encode('utf-8')

def gen_trdata_splits(datafile, splitsdir, trdatasize, numsamples, prefix='NEW'):
    '''
    Input: datafile, containing one data instance per line
    Output: create (keep if existing) numsamples subdirectories under splitsdir, and
            store 2 files in each---<prefix>.tr.gold and <prefix>.ts.gold.
            The tr.gold file contains trdatasize data instances randomly sampled
            from datafile. The ts.gold file contains the remaining instances.
    '''
    data = [line.rstrip() for line in open(datafile).readlines()]
    N = len(data)
    if trdatasize > N: trdatasize = N/2
    SEED = 2013
    random.seed(SEED)
    for i in range(numsamples):
        #random sample
        trdataidx = random.sample(xrange(N), trdatasize)
        #save tr and ts
        dirname = os.path.join(splitsdir, str(i))
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        f = open(os.path.join(dirname, prefix+'.tr.gold'), 'w')
        print >> f, '\n'.join( [data[j] for j in trdataidx] )
        f.close()
        f = open(os.path.join(dirname, prefix+'.ts.gold'), 'w')
        print >> f, '\n'.join( [data[j] for j in range(N) if j not in trdataidx] )
        f.close()

def gen_semrel_pairs(l1l2candfile, l2l1candfile, numpairsperfile, pairsfile):
    '''Scan candidates file. For each word, take its highest ranked candidate
    and form a pair. Do this for the reverse candidates file too. Only take
    numpairsperfile pairs from each file, though.
    '''
    pairs = set([])
    i = 0
    for line in gzip.open(l1l2candfile):
        line = line.decode('utf-8').rstrip()
        rec = line.split('\t')
        w1, w2 = rec[0], rec[1].split(' ')[0]
        pairs.add(w1+"#"+w2)
        i += 1
        if i == numpairsperfile: break
    i = 0
    for line in gzip.open(l2l1candfile):
        line = line.decode('utf-8').rstrip()
        rec = line.split('\t')
        w2, w1 = rec[0], rec[1].split(' ')[0]
        pair = w1+"#"+w2
        if pair not in pairs:
            pairs.add(w1+"#"+w2)
            i += 1
        if i == numpairsperfile: break
    f = open(pairsfile, 'w')
    for pair in pairs:
        w1, w2 = pair.split("#")
        print >> f, (w1+'\t'+w2).encode('utf-8')
    f.close()

def print_wid_map(vocfile, startid=1):
    '''Given a vocabulary list, create a word-id mapping.
    '''
    i = startid
    for w in readwordlistaslist(vocfile):
        print (unicode(i)+'\t'+w).encode('utf-8')
        i += 1

def print_did_map(corpdir, startid=1):
    '''Given a directory, list its files and create a file-id mapping.
    '''
    i = startid
    for fname in os.listdir(corpdir):
        print unicode(i)+'\t'+fname
        i += 1

def corp2mat_sparse(lang, corpdir, widfile, matfile):
    '''Convert a corpus to a term-doc matrix (input for CCA).
    Save in sparse matrix format of matlab.
    '''
    widmap = {}
    for line in open(widfile):
        line = line.decode('utf-8').rstrip()
        wid, w = line.split('\t')
        widmap[w] = wid
    corp = Corpus(lang)
    corp.readdir(corpdir)
    f = open(matfile, 'w')
    for w in corp.index:
        for d in corp.index[w]:
            print >> f, widmap[w], d, int(corp.index[w][d])
    f.close()

def formattrszeval(trszevalfile):
    # input file format:
    # trsz1
    # <mrr> <pres@1> <pres@5> <#ts wds>
    # <mrr> <pres@1> <pres@5> <#ts wds>
    # ...
    # trsz2
    # <mrr> ...
    # ...

    mrrs, presat1s, presat5s, numtssamps = None, None, None, None
    trsz = None
    for line in open(trszevalfile):
        if line.startswith('#'): continue
        line = line.rstrip()
        rec = line.split()
        if len(rec) == 1:
            if mrrs: # for first trsz, this is skipped
                print "%s\t%.4f\t%.4f\t%.4f" % (trsz,
                                                numpy.mean(mrrs), #numpy.std(mrrs)
                                                numpy.mean(presat1s), #numpy.std(presat1s)
                                                numpy.mean(presat5s), #numpy.std(presat5s)]
                                                #numpy.mean(numtssamps), numpy.std(numtssamps)
                                                )
            trsz = rec[0]
            mrrs, presat1s, presat5s, numtssamps = [], [], [], []
        else:
            mrr, presat1, presat5, numtssamp = rec
            mrrs.append(float(mrr))
            presat1s.append(float(presat1))
            presat5s.append(float(presat5))
            numtssamps.append(int(numtssamp))
    if mrrs: # for last trsz
        print "%s\t%.4f\t%.4f\t%.4f" % (trsz,
                                        numpy.mean(mrrs), #numpy.std(mrrs)
                                        numpy.mean(presat1s), #numpy.std(presat1s)
                                        numpy.mean(presat5s), #numpy.std(presat5s)]
                                        #numpy.mean(numtssamps), numpy.std(numtssamps)
                                        )

def formatevalblesum(evalblesumfile, fmt='TABLE'):
    # input file format:
    # l1l2
    # <mrr> <pres@1> <pres@5> <#ts wds>
    # <mrr> <pres@1> <pres@5> <#ts wds>
    # ...
    # l3l4
    # <mrr> ...
    # ...
    metrics = ['MRR', 'Pres@1', 'Pres@5']
    MRR, PRESat1, PRESat5 = metrics
    NUMTSSAMP = 'numtssamp'

    lp = {} # language pairs
    ls = set() # languages
    l1l2 = None
    mrrs, presat1s, presat5s, numtssamps = None, None, None, None
    for line in open(evalblesumfile):
        if line.startswith('#'): continue
        line = line.rstrip()
        rec = line.split()
        if len(rec) == 1:
            if mrrs: # for first l1l2, this is skipped
                lp[l1l2] = {}
                lp[l1l2][MRR] = [numpy.mean(mrrs), numpy.std(mrrs)]
                lp[l1l2][PRESat1] = [numpy.mean(presat1s), numpy.std(presat1s)]
                lp[l1l2][PRESat5] = [numpy.mean(presat5s), numpy.std(presat5s)]
                lp[l1l2][NUMTSSAMP] = [numpy.mean(numtssamps), numpy.std(numtssamps)]
            l1l2 = rec[0]
            l1, l2 = l1l2[:-2], l1l2[-2:]
            ls = ls.union( [l1,l2] )
            mrrs, presat1s, presat5s, numtssamps = [], [], [], []
        else:
            mrr, presat1, presat5, numtssamp = rec
            mrrs.append(float(mrr))
            presat1s.append(float(presat1))
            presat5s.append(float(presat5))
            numtssamps.append(int(numtssamp))
    if mrrs: # for last l1l2
        lp[l1l2] = {}
        lp[l1l2][MRR] = [numpy.mean(mrrs), numpy.std(mrrs)]
        lp[l1l2][PRESat1] = [numpy.mean(presat1s), numpy.std(presat1s)]
        lp[l1l2][PRESat5] = [numpy.mean(presat5s), numpy.std(presat5s)]
        lp[l1l2][NUMTSSAMP] = [numpy.mean(numtssamps), numpy.std(numtssamps)]
    if fmt=='TABLE':
        llist = sorted(list(ls))
        for met in metrics:
            print met+"\t"+"\t".join(llist)
            for l1 in llist:
                row = []
                for l2 in llist:
                    l1l2 = l1+l2
                    if l1l2 in lp:
    #                    cell = "%0.4f $\\pm$ %.4f (%d)" % (lp[l1l2][met][0], lp[l1l2][met][1], lp[l1l2][NUMTSSAMP][0])
    #                    cell = "%0.4f (%d)" % (lp[l1l2][met][0], lp[l1l2][NUMTSSAMP][0])
    #                    cell = "%0.4f $\\pm$ %.4f" % (lp[l1l2][met][0], lp[l1l2][met][1])
                        cell = "%0.4f" % (lp[l1l2][met][1])
    #                    cell = "%0.4f" % (lp[l1l2][met][0])
                    else:
                        cell = "--"
                    row.append(cell)
                print l1+"\t"+"\t".join(row)
    elif fmt=='LIST':
        for l1l2 in sorted(lp.iterkeys()):
            print l1l2 + "\t" + " ".join( [str(sco[0]) for sco in lp[l1l2].itervalues()] )
    else:
        pass
def formatevalbledet(evalbledetfile, wstatfile):
    # evalbledetfile file format:
    # w1 <mrr> <pres@1> <pres@5>
    # w2 <mrr> <pres@1> <pres@5>
    # ...
    # output (after averaging repeats):
    # cf(w1) df(w1) adc(w1) <mrr> <pres@1> <pres@5>
    # cf(w1) df(w2) adc(w3) <mrr> <pres@1> <pres@5>
    pm4w = {} # perf measure for word
    for line in open(evalbledetfile):
        line = line.decode('utf-8').rstrip()
        rec = line.split()
        w = rec.pop(0)
        if w not in pm4w:
            pm4w[w] = []
            for i in range(len(rec)):
                pm4w[w].append( [] )
                pm4w[w][i].append( float(rec[i]) )
        else:
            for i in range(len(rec)):
                pm4w[w][i].append( float(rec[i]) )
    for w in pm4w:
        for i in range(len(pm4w[w])): # for each perf measure
            pm4w[w][i] = numpy.mean(pm4w[w][i])
    for line in open(wstatfile):
        line = line.decode('utf-8').rstrip()
        if line.startswith('#'): continue
        w, cf, df, adc, idf, cfidf, tcfidf = line.split()
        if w in pm4w:
            line = "%s\t%s\t%s\t%s\t%s" % (w, cf, df, adc,
                                            "\t".join([unicode(pm) for pm in pm4w[w]]) )
            print line.encode('utf-8')

def goldbyexpert(goldfile, srccol, tgtcol, src2tgtpmffiles, l1, l2, srcvecfile):
    '''pmmfiles comma separated'''
    from GenClusModel import WordFeat
    vecs = WordFeat.loadfeat(srcvecfile)
    pmffiles = src2tgtpmffiles.split(',')
    pmfs = {}
    for pmffile in pmffiles:
        pmfs[pmffile] = L1L2PMF(l1, l2, pmffile)
    for line in open(goldfile):
        line = line.decode('utf-8').rstrip()
        rec = line.split()
        s, t = rec[srccol-1], rec[tgtcol-1]
        pfile, pr = max( [(pmffile, pmf.pmf[s][t] if t in pmf.pmf[s] else 0) for pmffile, pmf in pmfs.iteritems()] , key=lambda x: x[1] )
        beg = pfile.find('ptm/') + 4
        end = pfile.find('.ticue')
        print ("%s\t%s\t%s\t%s" % (s, t, unicode(pfile[beg:end]), \
            u'\t'.join( [unicode(f) for f in vecs[s][:6]] ) )).encode('utf-8')



if __name__ == "__main__":
    cmd, args = sys.argv[1], sys.argv[2:]
    if cmd == "createdict":
        d1file, l11, l12, d2file, l21, l22 = args
        createdict(d1file, l11, l12, d2file, l21, l22)
    elif cmd == "trunccands":
        candfile, probmass, maxperline, newcandfile = args
        trunccands(candfile, float(probmass), int(maxperline), newcandfile)
    elif cmd == "filteronquery":
        qfile, candfile = args
        filteronquery(qfile, candfile)
    elif cmd == "docstatfromstate":
        stfile, l1, l2, id2lang = args
        id2lang = id2lang.split(',')  # e.g. 0:en,1:hi,2:bn ...
        id2lang = dict( [x.split(':') for x in id2lang] )
        docstatfromstate(stfile, l1, l2, id2lang)
    elif cmd == "getcorpfromids":
        docidfile, col, srcdir, destdir = args
        getcorpfromids(docidfile, int(col), srcdir, destdir)
    elif cmd == "filtcorponvocab":
        vocfile, srcdir, destdir = args
        filtcorponvocab(vocfile, srcdir, destdir)
    elif cmd == "getwnsim":
        envocfile, qfile = args
        getwnsim(envocfile, qfile)
    elif cmd == "printcorpwordstat":
        corpdir, lang = args
        printcorpwordstat(corpdir,lang)
    elif cmd == "printcorpdocstat":
        corpdir, = args
        printcorpdocstat(corpdir)
    elif cmd == "dictstats":
        dfile, l1, l2, l1vocfile, l2vocfile, l1qfile, l2qfile = args
        dictstats(dfile, l1, l2, l1vocfile, l2vocfile, l1qfile, l2qfile)
    elif cmd == "preprocess":
        docfile, stopfile, minsize, lang = args
        if stopfile == 'NONE': stopfile = None
        preprocess(docfile, stopfile, int(minsize), lang)
    elif cmd == "get_t1t2_map":
        i1i2mapfile, itmap1file, itmap2file = args
        get_t1t2_map(i1i2mapfile, itmap1file, itmap2file, IDTITSEP='\t')
    elif cmd == "print_useful_words":
        wordstatsfile, mincf, maxdffrac, maxvocsize, validwordsfile = args
        print_useful_words(wordstatsfile, int(mincf), float(maxdffrac), int(maxvocsize),
                           None if validwordsfile == "NONE" else validwordsfile)
    elif cmd == "print_useful_docs":
        docstatfile, mintypes, mintokens = args
        print_useful_docs(docstatfile, int(mintypes), int(mintokens))
    elif cmd == "intersect_useful_map":
        l1docfile, l2docfile, l1l2idmapfile = args
        intersect_useful_map(l1docfile, l2docfile, l1l2idmapfile)
    elif cmd == "get_eval_pairs":
        l1, l2, l1l2dictfile, l1vocfile, l2vocfile = args
        get_eval_pairs(l1, l2, l1l2dictfile, l1vocfile, l2vocfile)
    elif cmd == "get_vocab_from_ptm":
        ptmfile, lang = args
        get_vocab_from_ptm(ptmfile, lang)
    elif cmd == "gen_trdata_splits":
        datafile, splitsdir, trdatasize, numsamples, prefix = args
        gen_trdata_splits(datafile, splitsdir, int(trdatasize), int(numsamples), prefix)
    elif cmd == "gen_semrel_pairs":
        l1l2candfile, l2l1candfile, numpairsperfile, pairsfile = args
        gen_semrel_pairs(l1l2candfile, l2l1candfile, int(numpairsperfile), pairsfile)
    elif cmd == "print_wid_map":
        vocfile, startid = args
        print_wid_map(vocfile, int(startid))
    elif cmd == "print_did_map":
        corpdir, startid = args
        print_did_map(corpdir, int(startid))
    elif cmd == "corp2mat_sparse":
        lang, corpdir, widfile, matfile = args
        corp2mat_sparse(lang, corpdir, widfile, matfile)
    elif cmd  == "formatevalblesum":
        evalblesumfile, fmt = args
        formatevalblesum(evalblesumfile, fmt)
    elif cmd == "formattrszeval":
        trszevalfile, = args
        formattrszeval(trszevalfile)
    elif cmd == "formatevalbledet":
        evalbledetfile, wstatfile = args
        formatevalbledet(evalbledetfile, wstatfile)
    elif cmd == "goldbyexpert":
        goldfile, srccol, tgtcol, src2tgtpmffiles, l1, l2, srcvecfile = args
        goldbyexpert(goldfile, int(srccol), int(tgtcol), src2tgtpmffiles, l1, l2, srcvecfile)
    else:
        print "invalid command:", cmd
