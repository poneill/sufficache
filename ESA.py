"""Extended Suffix Array Class"""
from utils import *

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
    i = 0
    length = min(len(x), len(y))
    for i, (x, y) in enumerate(zip(x, y)):
        if x != y:
            return i
    return length


class ESA(object):
    def __init__(self, word):
        self.word = word
        self.suffixes = suffixes_of(word)
        self.suf = [len(word) - len(self.suffixes[i]) + 1
                    for i in range(len(word)+1)]
        self.lcp = [0] + [lcp(x, y) for (x, y) in pairs(self.suffixes)]
        self.skp = [min([j for j in range(i+1, len(word) + 1)
                         if self.lcp[j] < self.lcp[i]] + [len(word) + 1])
                         for i in range(len(word) + 1)]
