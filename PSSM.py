import random
from sufficache_utils import *
from math import log,exp

kB  = 0.0019872041 #kcal/mol (!)
temp = 310.15 #37C
beta = -1/(kB*temp)

BASES = "acgt"
BASE_PAIR_ORDERING = {"a":0,"c":1,"g":2,"t":3}
INV_BASE_PAIR_ORDERING = {0:"a",1:"c",2:"g",3:"t"}

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
    """Implements a position-specific scoring matrix.  The primary
    data representation is a list of lists of log-likelihoods for each
    character, for each column, of the pssm.  The characters are
    assumed to be ordered according to BASE_PAIR_ORDERING,
    i.e. alphabetically."""

    bpo = {"a":0,"c":1,"g":2,"t":3} #base-pair ordering
    inv_bpo = {0:"a",1:"c",2:"g",3:"t"}
    
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
            return [column.count(c) for c in BASES]
        if not type(data[0][0]) is str: #"if data is a matrix of counts..."
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
    
    def upper_partial_thresholds(self,theta):
        """Returns a list of partial thresholds to be interpreted as
        follows: After having read position i, if you have scored at
        least pt[i] then the word will be a positive no matter what
        the rest of the word contains"""
        mins = map(min, self.columns)
        n = len(self.columns)
        thetas = [theta - sum(mins[i+1:]) for i in range(n)]
        return thetas


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
        """Implements algorithm 1 from Beckstette et al. (2006)"""
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
        M = lambda d, char: self[d][BASE_PAIR_ORDERING[char]]
        while (i < n):
            if i % 1000 == 0:
                print i
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
                #print "skipping"
                #print "d:",d
                #i_before = i
                i = esa.skipchain(i, d)
                #print "increased i by:",i - i_before
            depth = lcp[i]
        return matches

    def search_esa_for_top_k(self,esa,theta,k):
        """Implements algorithm 1 from Beckstette et al. (2006)"""
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
        M = lambda d, char: self[d][BASE_PAIR_ORDERING[char]]
        while (i < n):
            if i % 1000 == 0:
                print i
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
                print "adding"
                matches.append((suf[i], score))
                if len(matches) > k:
                    print "revising theta"
                    matches = sorted(matches,
                                     key = lambda(suf,score):score,
                                     reverse=True)[:k]
                    index,least_score = matches[-1]
                    theta = least_score
                    print "theta:",theta
                    thetas = self.partial_thresholds(theta)
                    print "thetas:",thetas
                while(i < n):
                    i += 1
                    if lcp[i] >= m:
                        matches.append((suf[i], score))
                    else:
                        break
            else:
                #print "skipping"
                #print "d:",d
                #i_before = i
                i = esa.skipchain(i, d)
                #print "increased i by:",i - i_before
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

    def info_score(self,seq):
        return (sum())
        
    def trap(self,seq,beta=beta,strand_correction=True):
        """Return the binding affinity as given by the TRAP model.
        See Manke 2008, Roider 2007."""
        if strand_correction:
            e_f = self.trap(seq,    strand_correction=False)
            e_b = self.trap(wc(seq),strand_correction=False)
            return log(exp(beta * e_f) + exp(beta * e_b))/beta
        else:
            n = len(self.motif)
            w = len(self.motif[0])
            lamb = 0.7
            ln_R_0 = 0.585 * w - 5.66
            E = 1/lamb * sum(log((self.counts[i][PSSM.bpo[self.consensus[i]]]+1)/
                                 float((self.counts[i][PSSM.bpo[seq[i]]] + 1)))
                             for i in range(len(seq)))
        #we define beta = -1/kBT whereas Manke defines b = 1/kBT,
        #hence change of sign
        return E + ln_R_0

    def slide_trap(self,genome,strand_correction=True):
        w = len(self.motif[0])
        return [self.trap(genome[i:i+w],strand_correction=strand_correction)
                           for i in verbose_gen(range(len(genome) - w + 1),100000)]
    def slide_score(self,genome):
        w = len(self.motif[0])
        return [self.score(genome[i:i+w])
                           for i in verbose_gen(range(len(genome) - w + 1),100000)]

    def enumerate_high_scoring_sites(self,theta):
        """Return sites over DNA alphabet scoring better or equal to theta"""
        distance = 1
        sites = enumerate_mutant_sites(self.get_consensus(),distance)
        while any(map(lambda site:self.score(site) > theta,sites)):
            distance += 1
            print distance
            sites = enumerate_mutant_sites(self.get_consensus(),distance)
        print "distance:",distance
        return [site for site in sites if self.score(site) >= theta]

    def min_score_given(self,partial_word):
        l = len(partial_word)
        return self.score(partial_word) + sum(map(min,self.columns[l:]))

    def max_score_given(self,partial_word):
        l = len(partial_word)
        return self.score(partial_word) + sum(map(max,self.columns[l:]))

    def next_best_word_naive(self,word):
        """Given a word, return the next best scoring 1-hamming dist. word"""
        bpo = PSSM.bpo
        inv_bpo = PSSM.inv_bpo
        best_column = None
        best_char = None
        best_diff = -1e6
        for j,char in enumerate(word):
            col = self.columns[j]
            val = col[bpo[char]]
            possible_vals = [(i,v) for (i,v) in enumerate(col)
                             if v < val or (v == val and i > bpo[char])] #later
                                                                       #in
                                                                       #lexico
                                                                       #order
            #print j,char,val,possible_vals
            if possible_vals:
                best_possible_index_val = max(possible_vals,key=lambda(x,y):y)
                best_possible_index,best_possible_val = best_possible_index_val
                best_possible_diff = best_possible_val - val
                if best_possible_diff > best_diff:
                    #print "accepting:",j,inv_bpo[best_possible_index]
                    best_diff = best_possible_diff
                    best_column = j
                    best_char = inv_bpo[best_possible_index]
        if best_column is None:
            raise Exception("Reached worst word")
        new_word = word[:best_column] + best_char + word[best_column+1:]
        #print word
        #print new_word
        return new_word
        
    def next_best_word(self,word):
        word_score = self.score(word)
        best_score = [self.score(self.next_best_word_naive(word))] # for mutability
        def walk_tree(partial_word):
            if random.random() < 0.0001:
                print partial_word,best_score[0]
            if (self.max_score_given(partial_word) < best_score[0] or
                self.min_score_given(partial_word) > word_score):
                return None
            elif len(partial_word) == len(self):
                best_score[0] = self.score(partial_word)
                return (best_score[0],partial_word)
            else:
                children = filter(lambda x:x,[walk_tree(partial_word + b)
                                              for b in BASES])
                if children:
                    return max(children, key = lambda(score,word):word)
                else:
                    return None
        result = walk_tree("")
        print best_score
        return result
                
        
    def fast_search(self,genome,theta):
        sites = self.enumerate_high_scoring_sites(theta)
        print "num sites:",len(sites)
        matches = []
        for i,site in enumerate(sites):
            print i,site
            score = self.score(site)
            matches.extend([(score,site) for site in re.findall(site,genome)])
            print len(matches)
        return matches
    
    
print("loaded PSSM")
