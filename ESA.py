"""Extended Suffix Array Class"""
from utils import *
import time
def lcp_functional(x, y):
    """Return the length of the longest common prefix of x and y.
    Because Python does not implement TCO, this will occasionally blow
    the stack (!) and therefore is deprecated.  It is only included as
    a reference implementation; use lcp instead."""
    if not x or not y or x[0] != y[0]:
        return 0
    else:
        return 1 + lcp(x[1:], y[1:])

def lcp(x, y):
    """Return the length of the longest common prefix of x and y"""
    #Because Python does not implement TCO, a functional
    #implementation will occasionally blow the stack (!).  That's why
    #this looks so ridiculous
    length = min(len(x), len(y))
    for i, (x, y) in enumerate(zip(x, y)):
        if x != y:
            return i
    return length

def lcp2(x, y, word):
    """Return the length of the longest common prefix of x and y"""
    #Because Python does not implement TCO, a functional
    #implementation will occasionally blow the stack (!).  That's why
    #this looks so ridiculous
    i = 0
    n = len(word) - max(x,y)
    while i < n:
        if word[x + i] == word[y + i]:
#            print word[x + i],word[y + i],i
            i += 1
        else:
            break
    return i

class ESA_deprecated(object):
    def __init__(self, word):
        self.word = word
        self._suffixes = suffixes_of(word)
        self.suf = [len(word) - len(self._suffixes[i]) + 1
                    for i in range(len(word)+1)]
        self.lcp = [0] + [lcp(x, y) for (x, y) in pairs(self._suffixes)]
        print len(self.lcp)
        self.skp = [min(gen_head((j for j in range(i+1, len(word) + 1)
                         if self.lcp[j] < self.lcp[i])),len(word) + 1)
                         for i in range(len(word) + 1)]

    def suffixes(self,i):
        return self._suffixes[i]
    
    def skipchain(self, i, d):
        n = len(self.word)
        print "skipchain"
        j = i + 1
        if i < n:
            while((j <= n) and (self.lcp[j] > d)):
                print "+1"
                j = self.skp[j]#should be +1? Tue Aug 21 13:27:34 EDT 2012
        else:
            j = n
        return j
    
class ESA(object):
    def __init__(self, word):
        self.word = word
        # self.suffixes = suffixes_of(word)
        # self.suf = [len(word) - len(self.suffixes[i]) + 1
        #             for i in range(len(word)+1)]
        n = len(word)
        print "sufs"
        self.suf = fast_sufs(word)
        print "lcp"
        self.lcp = [0] + [lcp2(self.suf[i], self.suf[j], word)
                          for (i, j) in gen_pairs_range(n + 1)]
        print len(self.lcp)
        print "skp"
        self.skp = [min(gen_head((j for j in xrange(i+1, n + 1)
                         if self.lcp[j] < self.lcp[i])),n + 1)
                         for i in xrange(n + 1)]

    def suffixes(self,i):
        return self.word[self.suf[i]:]

    def skipchain(self, i, d):
        n = len(self.word)
        print "skipchain"
        j = i + 1
        if i < n:
            while((j <= n) and (self.lcp[j] > d)):
                j = self.skp[j]#should be +1? Tue Aug 21 13:27:34 EDT 2012
        else:
            j = n
        return j
    
def alphabet(word):
    return ['$'] + sorted(list(set(word)))

S = 0
L = 1

def is_S(x):
    return x is None or x == S

def suffix_array(word,sigma = None):
    #1 First bucket the suffixes into array A, fill SL_array
    T = word + "$"
    N = len(T)
    sigma = alphabet(word)
    buckets = {c:[] for c in sigma}
    SL_array = [None for i in range(N)]
    S_yet = False
    dist = 0
    i = 0
    j = 0
    while i < N - 1:
        print "i:",i
        buckets[T[i]].append(i)
        if not T[i] == T[i + 1]:
            print "iffing"
            SL_array[i] = S if T[i] < T[i + 1] else L
            # print "before:",dist, S_yet
            # print "after:",dist, S_yet
        else:
            print "elsing"
            j = i + 1
            while T[j] == T[j + 1]:
                j += 1
                print "j:",j
            val = S if T[j] > T[i] else L
            print i,j
#            for k in range(i,j):
            while i < j:
                SL_array[i] = val
                i += 1
            i = j - 1        
        i += 1
    buckets['$'] = [N - 1]
    bucket = ([buckets[c] for c in sigma],[])#] if not c == '$'],[])
    dists = s_dist(SL_array)
    m = max(dists)
    #construct s-distance lists
    s_distance_lists = {j:{s:[b for b in buckets[s] if dists[b] == j]
                           for s in sigma} for j in range(1,m+1)}
    print "s distance list:",s_distance_lists
    #sort all type-S substrings
    print "buckets:",buckets
    type_s_lists = {c:[i for i in buckets[c] if is_S(SL_array[i])] for c in sigma}
    print "type_s_lists:",type_s_lists
    for j in range(1,m + 1):
        for c in sigma:
            for s_d in s_distance_lists[j][c]:
                k = s_d - j
                type_s_list = type_s_lists[c]
                if k in type_s_list:
                    type_s_list.remove(k)
                    type_s_list.insert(0,k)
    print "type_s_lists:",type_s_lists
                    
                    
    return SL_array,dists,s_distance_lists

def skew(word):
    n = len(word)
    word += "$" * ((3 - n) % 3)
    
    print word
def s_dist(xs):
    s_yet = False
    s_dists = []
    dist = 0
    for x in xs:
        dist += 1 * s_yet
        s_dists.append(dist)
        if is_S(x):
            s_yet = True
            dist = 0
    return s_dists
    
print ("loaded ESA")
