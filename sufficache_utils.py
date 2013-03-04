"""Various utility functions for sufficache project"""
import itertools,math
import random

#BASE_PAIR_ORDERING = "acgt"

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
    while (i < n and j < n):
        # wi = word[i]
        # wj = word[j]
        if not word[i] == word[j]:
            return cmp(word[i], word[j])
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
    def inline(i,j):
        return lexico(i,j,word)
    sufs = sorted(range(len(word)),
                  cmp=inline) + [len(word)]
    return sufs

def faster_sufs(word):
    pass

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
    return {"A":"T","T":"A","G":"C","C":"G",
            "a":"t","t":"a","g":"c","c":"g"}[base]
    
def wc_reference(word):
    return [complement(c) for c in word[::-1]]

def wc(word):
    """Reverse complement function"""
    # see wc_ref for non-terrible implementation
    new_word = ""
    #~3x speedup by inlining
    for c in word:
        new_word = {"A":"T","T":"A","G":"C","C":"G",
                     "a":"t","t":"a","g":"c","c":"g"}[c] + new_word
    return new_word

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
        
def random_site(n):
    return "".join([random.choice("acgt") for i in range(n)])

def choose2(xs,gen=False):
    """return list of choose(xs, 2) pairs, retaining ordering on xs"""
    if gen:
        return ((x1, x2) for i, x1 in enumerate(xs) for x2 in xs[i+1:])
    else:
        return [(x1, x2) for i, x1 in enumerate(xs) for x2 in xs[i+1:]]

def partition(pred, xs):
    part = []
    appended = False
    for x in xs:
        appended = False
        for p in part:
            if pred(x,p[0]):
                p.append(x)
                appended = True
                break
        if not appended:
            part.append([x])
    return part

def foldl(f,z,xs):
    if not xs:
        return z
    else: 
        return foldl(f,f(z,xs[0]),xs[1:])

def foldl1(f,xs):
    return foldl(f,xs[0],xs[1:])

def truncate(x,n=3):
    s = str(x)
    try:
        dec_point = s.index('.')
    except ValueError:
        return x
    return float(s[:dec_point + n + 1])

def concat(xxs):
    return sum(xxs,[])

print "loaded utils"

