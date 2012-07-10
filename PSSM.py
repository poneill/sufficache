import random
from utils import *

def contains_binding_sites(data):
    return all([[c in BASE_PAIR_ORDERING for c in site] for site in data])

class PSSM(list):
    """Implements a position-specific scoring matrix.  The primary
    data representation is a list of lists of log-likelihoods for each
    character, for each column, of the pssm.  The characters are
    assumed to be ordered according to BASE_PAIR_ORDERING,
    i.e. alphabetically."""
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
        print data
        if not type(data[0][0]) is str:
            self.columns = [convert_column(col) for col in data]
        else:
            assert(contains_binding_sites(data))
            counts = [count(col) for col in transpose(data)]
            self.columns = [convert_column(col) for col in counts]
            self.motif = data
            
    def __len__(self):
        return len(self.columns)

    def __iter__(self):
        for col in self.columns:
            yield col

    def __getitem__(self,i):
        return self.columns[i]

    def __getslice__(self,i,j):
        return self.columns[i:j]

    def score(self,word):
        """Return log-odds score for word"""
        return sum([col[BASE_PAIR_ORDERING.index(base)]
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

    def cutoff(self, alpha, n):
        """Return a Monte Carlo estimate (over n trials) for the cutoff
        value theta such that P(pssm.score(w) > theta) < alpha"""
        samples = sorted(self.sample_scores(n), reverse=True)
        return samples[int(alpha*n)]

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
        def skipchain(lcp, skp, n, i, d):
            j = i + 1
            if i < n:
                while((j <= n) and (lcp[j] > d)):
                    j = skp[j]
                else:
                    j = n
            return j

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
                score = score + M(d, suffixes[i][d])
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
                i = skipchain(lcp, skp, n, i, d)
            depth = lcp[i]
        return matches


