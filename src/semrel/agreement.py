# Author: Tom Lippincott <tom@cs.columbia.edu>

"""calculate agreement coefficients for labeling results using various distance metrics

HOW TO USE agreement.py 
(Thomas M. Lippincott) 

The script can calculate alpha, kappa, pi and S using either simple
binary equality or MASI as the distance metric between labelings.

The file example.txt contains an example data file with the
"integrated coding example" from Artstein and Poesio's paper.

It creates a DataSet class that resembles the "common notation" from
Artstein & Poesio, so any metric definable with the i/k/c/etc notation
is easy to add.

Running "python agreement.py -h" will print usage info.

e.g. "python agreement.py -f example.txt -a alpha -d binary" prints
the alpha value from the paper. 

Right now the script only imports data in the format of the example
file, but this process is separate from the class itself so if there
is a particular data source/file format you can send me an example or
description and I'll add it to the import methods. 

"""

import sys, sets
import math

class DataSet:

    def __init__(self,distance):
        self.distance = distance
        self.I = sets.Set()
        self.K = sets.Set()
        self.C = sets.Set()
        self.data = []

    def __str__(self):
        return "\r\n".join(map(lambda x:"%s\t%s\t%s"%(x['coder'],x['item'].replace('_',"\t"),",".join(x['labels'])),self.data))

    def get(self,key):
        if(key=="pi"):
            return self.pi()
        elif(key=="kappa"):
            return self.kappa_avg()
        elif(key=="alpha"):
            return self.alpha()
        elif(key=="S"):
            return self.S()

    # and array of values of form "coder_item_label1:label2"
    def load_array(self,array,sep=':'):
        for i in array:
            coder,item,temp = i.split('_')
            labels = sets.Set(temp.split(sep))
            self.C.add(coder)
            self.K.add(labels)
            self.I.add(item)
            self.data.append({'coder':coder,'labels':labels,'item':item})

    # remove items that aren't fully annotated
    def sanitize(self):
        temp_I = self.I.copy()
        temp_data = self.data
        for i in temp_I:
            if(len(filter(lambda x:x['item']==i,temp_data))!=len(self.C)):
                self.I = self.I.difference(sets.Set([i]))
                self.data = filter(lambda x:x['item']!=i,self.data)
        return

    # agreement value on a given item between two coders
    def agr(self,cA,cB,i):
        kA = filter(lambda x:x['coder']==cA and x['item']==i,self.data)[0]
        kB = filter(lambda x:x['coder']==cB and x['item']==i,self.data)[0]
        return 1.0 - float(self.distance(kA['labels'],kB['labels']))

    # implements the "N"-notation
    def N(self,k="",i="",c=""):
        if(k!="" and i=="" and c==""):
            return len(filter(lambda x:k==x['labels'],self.data))
        elif(k!="" and i!="" and c==""):
            return len(filter(lambda x:k==x['labels'] and i==x['item'],self.data))
        elif(k!="" and c!="" and i==""):
            return len(filter(lambda x:k==x['labels'] and c==x['coder'],self.data))
        else:
            print "You must pass either i or c, not both!"

    # observed agreement between two coders
    def Ao(self,cA,cB):
        return float(sum(map(lambda x:self.agr(cA,cB,x),self.I)))/float(len(self.I))

    # average observed agreement between all coders
    def avg_Ao(self):
        s = self.C.copy()
        counter = 0.0
        total = 0.0
        for cA in self.C:
            s.remove(cA)
            for cB in s:
                total += self.Ao(cA,cB)
                counter += 1.0
        return total/counter

    # Agreement Coefficients

    def S(self):
        """Bennett, Albert and Goldstein 1954

        """
        Ae = 1.0/float(len(self.K))        
        return (self.avg_Ao() - Ae)/(1.0 - Ae)

    def pi(self):
        """Scott 1955

        """
        total = 0.0
        for k in self.K:
            total += self.N(k=k)**2
        Ae = (1.0/(4.0*float(len(self.I)**2)))*total
        return (self.avg_Ao()-Ae)/(1-Ae)

    def pi_avg(self):
        pass

    def kappa(self,cA,cB):
        """Cohen 1960

        """
        Ae = 0.0
        for k in self.K:
            Ae += (float(self.N(c=cA,k=k))/float(len(self.I))) * (float(self.N(c=cB,k=k))/float(len(self.I)))
        ret = (self.Ao(cA,cB)-Ae)/(1.0-Ae)
        return ret

    def kappa_avg(self):
        vals = {}
        for a in self.C:
            for b in self.C:
                if(a==b or "%s%s"%(b,a) in vals):
                    continue
                vals["%s%s"%(a,b)] = self.kappa(a,b)
        return sum(vals.values())/float(len(vals))

    def alpha(self):
        """Krippendorff 1980

        """
        De = 0.0
        for j in self.K:
            for l in self.K:
                De += float(self.N(k=j)*self.N(k=l))*self.distance(j,l)
        De = (1.0/(len(self.I)*len(self.C)*(len(self.I)*len(self.C)-1)))*De
        ret = 1.0 - (self.Do()/De)
        return ret

    # observed disagreement (for alpha)
    def Do(self):
        total = 0.0
        for i in self.I:
            for j in self.K:
                for l in self.K:
                    total += float(self.N(i=i,k=j)*self.N(i=i,k=l))*self.distance(l,j)
        ret = (1.0/float((len(self.I)*len(self.C)*(len(self.C)-1))))*total
        return ret

# Distance Metrics

def masi_distance(set1,set2):
    """Passonneau

    """
    return 1 - float(len(set1.intersection(set2)))/float(max(len(set1),len(set2)))

# straightforward equality test
def binary_distance(set1,set2):
    """

    """
    if(set1==set2):
        return 0
    else:
        return 1

def example_distance(set1,set2):
    """

    """
    try:
        all = set1.union(set2)
    except:
        all = sets.Set([set1,set2])
    if('stat' in all and 'ireq' in all):
        return 1.0
    elif('stat' in all and 'chck' in all):
        return 0.5
    elif('ireq' in all and 'chck' in all):
        return 0.5
    return 0.0

def interval_distance(set1,set2):
    """Krippendorff

    """
    try:
        return pow(len(set1.difference(set2)),2)
    except:
        return binary_distance(set1,set2)

def semrel_distance(set1, set2):
    '''Distance between semantic relatedness scores 
    '''
#    return math.fabs( float( list(set1)[0] ) - float( list(set2)[0] ) )
    return ( float( list(set1)[0] ) - float( list(set2)[0] ) )**2

# What's Available

# distance hash
distance = {
    'masi':masi_distance,
    'binary':binary_distance,
    'interval':interval_distance,
    'example':example_distance,
    'semrel':semrel_distance
    }

# agreement hash
agreement = [
    'kappa',
    'alpha',
    'pi',
    'S',
    ]

def readfile(fname):
    '''Default file format'''
    vals = []
    for l in open(fname):
        coder,item,labels = l.split()
        if(not coder in exclude and (include[0]=='' or coder in include)):
            vals.append("%s_%s_%s"%(coder,item,labels))
    return vals

def semrel_readfile(fname):
    '''
    File format for semantic relatedness scores:
        a list of score pairs (one for each coder)
    '''
    vals = []
    item = 1
    codera, coderb = 'a', 'b'
    for l in open(fname):
        labela, labelb = l.split()
        labela = str( math.floor(float(labela)) )
        labelb = str( math.floor(float(labelb)) )
        vals.append("%s_%s_%s"%(codera,item,labela))
        vals.append("%s_%s_%s"%(coderb,item,labelb))
        item += 1
    return vals


if(__name__=='__main__'):

    import optparse

    # process command-line arguments
    parser = optparse.OptionParser()
    parser.add_option("-d","--distance",dest="distance",default="binary",help="|".join(distance.keys()))
    parser.add_option("-a","--agreement",dest="agreement",default="",help="|".join(agreement))
    parser.add_option("-e","--exclude",dest="exclude",default="",help="coder names to exclude (comma-separated), e.g. jane,mike")
    parser.add_option("-i","--include",dest="include",default="",help="coder names to include, same format as exclude")
    parser.add_option("-f","--file",dest="file",help="file to read labelings from, each line with three columns: labeler, item, labels'")
    parser.add_option("-s","--sep",dest="sep",default=":",help="char/string that separates the three columns in the file")
    (options,remainder) = parser.parse_args()

    include = options.include.split(',')
    exclude = options.exclude.split(',')

    # print average distance between all specified sets
    if(options.distance in distance and options.agreement==""):
        vals = []
        for a in range(len(remainder)-1):
            for b in range(a+1,len(remainder)):
                setA = sets.Set(remainder[a].split(':'))
                setB = sets.Set(remainder[b].split(':'))
                vals.append(distance[options.distance](setA,setB))
        print "average %s distance: %f"%(options.distance,float(sum(vals))/float(len(vals)))

    # print average agreement over all specified labelings
    elif(options.agreement in agreement and options.distance in distance):
        data = DataSet(distance[options.distance])
        if(options.file):
#            vals = []
            vals = semrel_readfile(options.file)
#            for l in open(options.file):
#                coder,item,labels = l.split()
#                if(not coder in exclude and (include[0]=='' or coder in include)):
#                    vals.append("%s_%s_%s"%(coder,item,labels))
            data.load_array(vals,options.sep)
        else:
            data.load_array(remainder)
        data.sanitize()
        print "%s\t%s value: %f"%(options.include,options.agreement,data.get(options.agreement))

    # invalid option
    else:
        print "You specified an invalid distance or agreement metric, see 'agreement_metrics.py -h'"
