"""Various utility functions for sufficache project"""

def normalize(xs):
    return map(lambda(x): x/float(sum(xs)),xs)

def lexicographic_cmp(xs,ys):
    i = 0
    for i, (x,y) in enumerate(itertools.izip(xs,ys)):
        if x != y:
            return cmp(xs[i:],ys[i:])
    return cmp(ys[i:],xs[i:])

def suffixes_of(word):
    suffixes = sorted([word[i:] for i in range(len(word))],
                      cmp=lexicographic_cmp)
    return [suf + "$" for suf in (suffixes + [''])]

def pprint(x):
    for row in x:
        print row
        
def nmers(n):
    if n == 1:
        return ["A","C","G","T"]
    else:
        return sum([map(lambda(b):b+c,nmers(n-1)) for c in base_pair_ordering],[])
    
def safe_log2(x):
    """Implements log2, but defines log2(0) = 0"""
    return math.log(x,2) if x > 0 else 0

def complement(base):
    return {"A":"T","T":"A","G":"C","C":"G"}[base]
    
def wc(word):
    return map(complement, word[::-1])

def split_on(xs, pred):
    """Split xs into a list of lists each beginning with the next x
    satisfying pred, except possibly the first"""
    indices = [i for (i,v) in enumerate(xs) if pred(v)]
    return [xs[i:j] for (i,j) in zip([0]+indices,indices+[len(xs)]) if i != j]

def tail(xs):
    return xs[1:] if xs else []

def partial_sums(xs,acc = 0):
    return [sum(xs[:i+1]) for i in range(len(xs))]

def pairs(xs):
    return zip([None] + xs, xs + [None])[1:-1]


