import random
from utils import *
from math import log

kB  = 0.0019872041 #kcal/mol (!)
temp = 310.15 #37C
beta = -1/(kB*temp)

def verbose_gen(xs,n=1):
    for i,x in enumerate(xs):
        if i % n == 0:
            print i
        yield x
        
def contains_binding_sites(data):
    return all([[c in BASE_PAIR_ORDERING for c in site] for site in data])

def verbose_gen(xs,n=1):
    for i,x in enumerate(xs):
        if i % n == 0:
            print i
        yield x

        
class PSSM(list):
    bpo = {"a":0,"c":1,"g":2,"t":3}
    """Implements a position-specific scoring matrix.  The primary
    data representation is a list of lists of log-likelihoods for each
    character, for each column, of the pssm.  The characters are
    assumed to be ordered according to BASE_PAIR_ORDERING,
    i.e. alphabetically."""

    bpo = {"a":0,"c":1,"g":2,"t":3} #base-pair ordering
    
    # def bpo(self,c):
    #     if c in "ac":
    #         if c == 'a':
    #             return 0
    #         else:
    #             return 1
    #     elif c == "g":
    #         return 2
    #     else:
    #         return 3

    def __init__(self,data,background_probs = (0.25,)*4):
        """Given a representation of a binding motif and an optional
        tuple of background probabilities, initialize a PSSM object.
        PSSMs can be initialized from a raw binding motif, i.e. a list
        of strings, or a row-major matrix of counts.  The latter
        option assumes that the matrix is of the form:
        [[count_A0,count_C0,count_G0,count_T0],
         [count_A1,count_C1,count_G1,count_T1]
         ...
         [count_An,count_Cn,count_Gn,count_Tn]].

         The proper interpretation (raw motif or counts) is be inferred
         from the type of the list."""
        def convert_column(col):
            return [safe_log2(c/p)
                    for (c, p) in zip(normalize(col), background_probs)]
        def count(column):
            return [column.count(c) for c in BASE_PAIR_ORDERING]
        if not type(data[0][0]) is str:
            self.columns = [convert_column(col) for col in data]
        else:
            assert(contains_binding_sites(data))
            counts = [count(col) for col in transpose(data)]
            self.counts = counts
            self.columns = [convert_column(col) for col in counts]
            self.motif = data
        self.length = len(self.columns)
        self.consensus = self.get_consensus()
        
    def __len__(self):
        return len(self.columns)

    def __iter__(self):
        for col in self.columns:
            yield col

    def __getitem__(self,i):
        return self.columns[i]

    def __getslice__(self,i,j):
        return self.columns[i:j]

    def base_pair_ordering(self,base):
        if base in "ac":
            if base == 'a':
                return 0
            else:
                return 1
        else:
            if base == 'g':
                return 2
            else:
                return 3

    def base_pair_ordering2(self,base):
        if base == 'a':
            return 0
        elif base == 'c':
            return 1
        elif base == 'g':
            return 2
        else:
            return 3

    
    def base_pair_ordering3(self,base):
        return {"a":0,"c":1,"g":2,"t":3}

        
    def score(self,word):
        """Return log-odds score for word"""
        #return sum([col[BASE_PAIR_ORDERING.index(base)]
        #            for col, base in zip(self.columns, word)])
        return sum([col[PSSM.bpo[base]]
                    for col, base in zip(self.columns, word)])

    def partial_thresholds(self,theta):
        """Returns a list of partial thresholds to be interpreted as
        follows: After having read position i, you must have scored at
        least pt[i], hence pt[len(pssm)] >= theta if the word is a
        positive"""
        rev_maxes = map(max, self.columns[::-1])
        rev_thetas = reduce(lambda ths, x: ths + [ths[-1] - x], rev_maxes, [theta])
        return rev_thetas[::-1][1:]

    def list_all_scores(self,columns = None):
        """Enumerate all possible scores.  Do not call this method for
        non-trivial motifs: use sample_scores instead"""
        if not columns:
            return [0]
        else:
            return [w + x for w in columns[0]
                    for x in self.list_all_scores(columns[1:])]

    def sample_scores(self, n):
        return [sum(random.choice(w) for w in self.columns) for i in range(n)]

    def cutoff_bootstrap_ci(self,alpha,n):
        """Compute a bootstrap confidence interval for the cutoff
        values theta such that P(pssm.score(w) > theta) < alpha.  That
        is, find interval (a,b) such that P(theta in (a,b)) = .95"""
        cutoffs = sorted([self.cutoff(alpha,n) for i in range(200)])
        a = cutoffs[5] # 2.5%
        b = cutoffs[195] # 2.5%
        return (a,b)
        
    def cutoff(self, alpha, n):
        """Return a Monte Carlo estimate (over n trials) for the cutoff
        value theta such that P(pssm.score(w) > theta) < alpha"""
        samples = sorted(self.sample_scores(n), reverse=True)
        return samples[int(alpha*n)]

    def search_esa_reference(self,esa,theta):
        """Reference implementation of search esa"""
        scores = [(i,self.score(esa.word[i:i+self.length]))
                  for i in xrange(len(esa.word))]
        return filter(lambda (i,score): score > theta,scores)
    
    def search_esa(self,esa,theta):
        """Implements algorithm 1 from Beckstette et al."""
        #The pseudo-code given in Beckstette is incorrect.  I have
        #corrected it here.
        #As for style: sorry mom, sorry dad, sorry Cormen.
        matches = []
        C = {}
        thetas = self.partial_thresholds(theta)
        suf, lcp, skp, suffixes = esa.suf, esa.lcp, esa.skp, esa.suffixes
        depth = 0
        i = 0
        n = len(esa.word)
        if 'n' in esa.word:
            return [(-1, -1)] #if there exist ns in the string, deal with
                             #it later
        m = len(self)
        M = lambda d, char: self[d][BASE_PAIR_ORDERING.index(char)]
        

        while (i < n):
            if n - m < suf[i]: #if too far in to match
                while(n - m < suf[i] and (i < n)):
                    i += 1
                    depth = min(depth, lcp[i])
                if i >= n:
                    return matches
            if depth == 0:
                score = 0
            else:
                score = C[depth - 1]
            d = depth - 1
            sentinel = True #hacks around the do-while loop
            while(sentinel or (d < m -1 and score >= thetas[d])):
                sentinel = False
                d = d + 1
                score = score + M(d, suffixes(i)[d])
                C[d] = score
            if(d == m - 1 and score >= theta):
                matches.append((suf[i], score))
                while(i < n):
                    i += 1
                    if lcp[i] >= m:
                        matches.append((suf[i], score))
                    else:
                        break
            else:
                i = esa.skipchain(i, d)
            depth = lcp[i]
        return matches

    def get_consensus(self):
        return "".join([max(col,key=lambda c:col.count(c))
                        for col in transpose(self.motif)])

    def het_index(self,seq):
        n = len(self.motif)
        return sum(log((self.counts[i][PSSM.bpo[self.consensus[i]]]+0.5)/
                       (self.counts[i][PSSM.bpo[seq[i]]] + 0.5))
                   for i in range(len(seq)))

    def trap(self,seq,beta=beta):
        """Return the binding affinity as given by the TRAP model.
        See Manke 2008, Roider 2007."""
        n = len(self.motif)
        w = len(self.motif[0])
        lamb = 0.7
        ln_R_0 = 0.585 * w - 5.66
        E = 1/lamb * sum(log((self.counts[i][PSSM.bpo[self.consensus[i]]]+1)/
                             (self.counts[i][PSSM.bpo[seq[i]]] + 1))
                for i in range(len(seq)))
        #we define beta = -1/kBT whereas Manke defines b = 1/kBT,
        #hence change of sign
        return E + ln_R_0

    def slide_trap(self,genome):
        w = len(self.motif[0])
        return [self.trap(genome[i:i+w])
                           for i in verbose_gen(range(len(genome) - w + 1),100000)]
    def slide_score(self,genome):
        w = len(self.motif[0])
        return [self.score(genome[i:i+w])
                           for i in verbose_gen(range(len(genome) - w + 1),100000)]

    def slide_trap(self,genome):
        w = len(self.motif[0])
        return [self.trap(genome[i:i+w])
                           for i in verbose_gen(range(len(genome) - w + 1),100000)]
        
print("loaded PSSM")
