'''
Created on 13-Nov-2012

@author: gtholpadi
'''
import regex

def read_word_list_file(wlfile):
    '''
    Input: a file containing a list of words, one per line (ignore empty or '#'-prefixed lines)
    Output: a list containing the words
    '''
    wl = []
    for line in open(wlfile):
        line = line.decode('utf-8').strip()
        if line and not line.startswith('#'): wl.append(line)
    return wl

class TokenFilter:
    '''
    Input: a list of tokens
    Output: same list, with some tokens removed 
    '''
    def __init__(self):
        '''
        Constructor
        '''
        pass
    def filter(self, tokens):
        return tokens

class StopWord(TokenFilter):
    '''
    Input: a list of tokens
    Output: same list, with stopwords removed
    '''
    def __init__(self, stlist=None, stfiles=None):
        '''
        Input: stopword list, or list of stopword files
        '''
        self.stoplist = set([])
        if stlist:
            self.stoplist = self.stoplist.union(set(stlist))
        if stfiles:
            for stfile in stfiles:
                if stfile:
                    wl = read_word_list_file(stfile)
                    self.stoplist = self.stoplist.union(set(wl)) 
    def filter(self, tokens):
        return [t for t in tokens if not t in self.stoplist]

class SmallWord(TokenFilter):
    '''
    Input: a list of tokens
    Output: same list, with small words removed    
    '''
    def __init__(self, minsize=3):
        '''
        Input: stopword list, or list of stopword files
        '''
        self.minsize = minsize
    def filter(self, tokens):
        return [t for t in tokens if len(t)>=self.minsize]

class KeepWord(TokenFilter):
    '''
    Input: a list of tokens
    Output: same list, with all words removed, except those present in a list
    '''
    def __init__(self, kplist=None, kpfiles=None):
        '''
        Input: keepword list, or list of keepword files
        '''
        self.keeplist = set([])
        if kplist: 
            for w in kplist: self.keeplist[w] = True
        if kpfiles:
            for kpfile in kpfiles:
                wl = read_word_list_file(kpfile)
                self.keeplist = self.keeplist.union(set(wl)) 
    def filter(self, tokens):
        return [t for t in tokens if t in self.keeplist]

class KeepScript(TokenFilter):
    '''
    Input: a list of tokens
    Output: same list, with all words removed, except those 
    whose script matches given script.
    '''
    lascmap = {
               'en' : 'Latin',
               'de' : 'Latin',
               'fr' : 'Latin',
               'bn' : 'Bengali',
               'hi' : 'Devanagari',
               'mr' : 'Devanagari',
               'kn' : 'Kannada',
               'ml' : 'Malayalam',
               'ta' : 'Tamil',
               'te' : 'Telugu' 
               }
    @staticmethod
    def get_script(lang):
        if lang in KeepScript.lascmap:
            return KeepScript.lascmap[lang]
        else:
            return 'ALL'
    def __init__(self, lang=None, script=None):
        '''
        Input: language or Unicode script 
        '''
        if script: 
            self.script = script
        elif lang:
            self.script = KeepScript.get_script(lang)
        if self.script != 'ALL':
            pat = "(?u)\\P{Script=%s}" % self.script
            self.p = regex.compile(pat)
    def filter(self, tokens):
        if self.script == 'ALL':
            return tokens
        else:
            # new implementation using regex module; NOT TESTED YET!!
            return [t for t in tokens if not self.p.search(t)]
            # old implementation calling perl; works, but _very_ slow
#            from subprocess import Popen, PIPE, STDOUT
#            args = "perl -e 'use Encode; while(<STDIN>){chomp; $_ = decode(\"UTF-8\", $_); m/\\P{%s}/ or print encode(\"UTF-8\", $_).\"\\n\";}'" % (self.script) 
#            p = Popen(args, stdin=PIPE, stdout=PIPE, stderr=STDOUT, shell=True)
#            op, err = p.communicate(input=('\n'.join(tokens)).encode('utf-8'))
#            if not err:
#                return op.decode('utf-8').split()
#            else:
#                print >> sys.stderr, err
#                return tokens
                    

if __name__ == "__main__":
    import sys
    docfile, lang = sys.argv[1:]
    ks = KeepScript(lang=lang)
    tokens = ks.filter(open(docfile).read().decode('utf-8').split())
    print ('\n'.join(tokens)).encode('utf-8')