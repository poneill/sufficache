"""Various utility functions for sufficache project"""
import itertools,math

BASE_PAIR_ORDERING = "acgt"

def normalize(xs):
    return map(lambda(x): x/float(sum(xs)),xs)

def lexicographic_cmp(xs,ys):
    i = 0
    for i, (x,y) in enumerate(itertools.izip(xs,ys)):
        if x != y:
            return cmp(xs[i:],ys[i:])
    return cmp(ys[i:],xs[i:])

def suffixes_of(word,n=None):
    suffixes = sorted([word[i:] for i in range(len(word))],
                      cmp=lexicographic_cmp)
    return [suf + "$" for suf in (suffixes + [''])]

def enum_suffixes_of(word,n=None):
    suffixes = sorted([(i,word[i:]) for i in range(len(word))],
                      cmp=lexicographic_cmp,key = lambda (e,w):w)
    sufs = [(enum,suf + "$") for (enum,suf) in (suffixes + [(len(word),'')])]
    return sufs

def lexico(i,j,word):
    """Given indices i and j in word, determine whether the ith or jth
    suffix is lexically prior"""
    n = len(word)
    if not (i < n and j < n):
        return cmp(i,j)
    wi = word[i]
    wj = word[j]
    if not wi == wj:
        return cmp(wi, wj)
    else:
        return lexico(i + 1,j+1,word)

def lexico2(i,j,word):
    """Given indices i and j in word, determine whether the ith or jth
    suffix is lexically prior"""
    n = len(word)
    while (i < n and j < n):
        wi = word[i]
        wj = word[j]
        if not wi == wj:
            return cmp(wi, wj)
        else:
            i += 1
            j += 1
    return cmp(i,j)

def test_lexico(n):
    xs = "".join([random.choice("atgc") for i in range(10)])
    print xs
    for trial in range(n):
        i = random.randint(0,9)
        j = random.randint(0,9)
#        print xs[i:],xs[j:],i,j,"\t",lexico(i,j,xs)==lexicographic_cmp(xs[i:],xs[j:])
        print xs[i:],xs[j:],i,j,"\t",lexico(i,j,xs)==lexico2(i,j,xs)
    
def fast_sufs(word,n=None):
    sufs = sorted(range(len(word)),
                  cmp=lambda i,j:lexico2(i,j,word)) + [len(word)]
    return sufs

def faster_sufs(word):
    
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

def pairs(xs,gen=False):
    "given a list [x,y,z,...], return a list [(x,y),(y,z),...]"
    z = gen_zip if gen else zip
    return z([None] + xs, xs + [None])[1:-1]

def gen_pairs_range(n):
    i = 0
    while i < n - 1:
        yield (i,i + 1)
        i += 1
        
def gen_zip(xs,ys):
    n = min(len(xs),len(ys))
    for i in xrange(n):
        print i
        yield (xs[i],ys[i])

def separate(pred, lst):
    """separates lst into a list of elements satisfying pred and a list of 
    elements not satisfying it.
    """
    sheep = []
    goats = []
    for elem in lst:
        if pred(elem):
            sheep.append(elem)
        else:
            goats.append(elem)
    return (sheep, goats)

def transpose(xxs):
    return zip(*xxs)
