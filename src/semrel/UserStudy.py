import gzip
import sys
import random
from Common import L1L2PMF

# Cluster labeling user study:
# Given a set of cluster labels in one language, 
# use TCI to get a set of labels in another language

def getTCIlabels(labfile, l1, l2, pmffile):
    # read labels.csv
    # load l1l2 pmf
    # for each l1 label set,
    #   for each l1 label
    #       get top prob l2 word from pmf; add to l2 label set
    #       if l1 label not present in pmf, use a top prob l2 word ...
    #           ... from orig. l2 label set instead

    l1l2pmf = {}
    for line in gzip.open(pmffile, 'rb'):
        line = line.rstrip().decode('utf-8')
        rec = line.split(u'\t')
        l1w = rec.pop(0)
        l1l2pmf[l1w] = [l2wpr.split()[0] for l2wpr in rec]

    labf = open(labfile)
    while True:
        line = labf.readline()
        if not line: break
        line = line.rstrip().decode('utf-8')
        if line == u'':
            continue
        #skip en labels
        l1labs = labf.readline().rstrip().decode('utf-8').split()
        l2labs = labf.readline().rstrip().decode('utf-8').split()
        K = len(l2labs)
        l2labsn = set([])
        while len(l2labsn) < K:
            lab1 = random.choice(l1labs)
            if lab1 in l1l2pmf:
                for lab2 in l1l2pmf[lab1]:
                    if lab2 not in l2labsn:
                        l2labsn.add(lab2)
                        break
            else:
                l1labs.remove(lab1)
        print (u' '.join(l2labsn)).encode('utf-8')


def main():
    cmd, args = sys.argv[1], sys.argv[2:]
    if cmd == "getTCIlabels":
        labfile, l1, l2, pmffile = args
        getTCIlabels(labfile, l1, l2, pmffile)
    else:
        print "invalid command: ", cmd

if __name__ == "__main__":
    main()
