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

def gen_head(xs):
    """Return first element of xs if such an element, otherwise return
    empty list"""
    #We do this in order to avoid computing the entire list
    #comprehension in the assignment of ESA.skp.  This leads to a 3x
    #speedup at the expense of slightly uglier code.
    try:
        return xs.next()
    except StopIteration:
        return []
        
class ESA(object):
    def __init__(self, word):
        self.word = word
        self.suffixes = suffixes_of(word)
        self.suf = [len(word) - len(self.suffixes[i]) + 1
                    for i in range(len(word)+1)]
        self.lcp = [0] + [lcp(x, y) for (x, y) in pairs(self.suffixes)]
        print len(self.lcp)
        self.skp = [min(gen_head((j for j in range(i+1, len(word) + 1)
                         if self.lcp[j] < self.lcp[i])),len(word) + 1)
                         for i in range(len(word) + 1)]
        
        
class ESA2(object):
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
        
        

print ("loaded ESA")
