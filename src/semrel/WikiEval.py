import sys
import gzip

def gettitsugg(l1l2candfile, l1titfile, l2titfile, numsugg):
    l1tit = set( [ line.decode('utf-8').rstrip().split(u'\t')[1] for line in open(l1titfile).readlines() ] )
    l2tit = set( [ line.decode('utf-8').rstrip().split(u'\t')[1] for line in open(l2titfile).readlines() ] )
    for line in gzip.open(l1l2candfile):
        rec = line.decode('utf-8').rstrip().split(u'\t')
        w1 = rec.pop(0)
        if w1 in l1tit: continue # skip l1 words that are titles
        suggtit = []
        for w2pr in rec:
            w2, pr = w2pr.split()
            if w2 in l2tit: suggtit.append(w2)
            if len(suggtit) >= numsugg: break
        if suggtit: print (w1+u'\t'+u' '.join(suggtit)).encode('utf-8')

def getnolltitsugg(l1l2candfile, l1titfile, llidmapfile, colnum, l2titfile, numsugg):
    l1idtit = dict( [ line.decode('utf-8').rstrip().split(u'\t') for line in open(l1titfile).readlines() ] )
    l1llid = set( [ line.decode('utf-8').rstrip().split(u'\t')[colnum-1] for line in open(llidmapfile).readlines() ] )
    l1tit = set( [ tit for aid, tit in l1idtit.iteritems() if aid not in l1llid] )
    l2tit = set( [ line.decode('utf-8').rstrip().split(u'\t')[1] for line in open(l2titfile).readlines() ] )
    for line in gzip.open(l1l2candfile):
        rec = line.decode('utf-8').rstrip().split(u'\t')
        w1 = rec.pop(0)
        if w1 not in l1tit: continue # skip l1 words that are not titles
        suggtit = []
        for w2pr in rec:
            w2, pr = w2pr.split()
            if w2 in l2tit: suggtit.append(w2)
            if len(suggtit) >= numsugg: break
        if suggtit: print (w1+u'\t'+u' '.join(suggtit)).encode('utf-8')

def main():
    cmd, args = sys.argv[1], sys.argv[2:]
    if cmd == "gettitsugg":
        l1l2candfile, l1titfile, l2titfile, numsugg = args
        gettitsugg(l1l2candfile, l1titfile, l2titfile, int(numsugg))
    elif cmd == "getnolltitsugg":
        l1l2candfile, l1titfile, llidmapfile, colnum, l2titfile, numsugg = args
        getnolltitsugg(l1l2candfile, l1titfile, llidmapfile, int(colnum), l2titfile, int(numsugg))
    else:
        print "invalid command:", cmd

if __name__ == "__main__":
    main()
