import random
from sufficache_utils import *
from math import log,exp

kB  = 0.0019872041 #kcal/mol (!)
temp = 310.15 #37C
beta = 1/(kB*temp)

BASES = "ACGT"
BASE_PAIR_ORDERING = {"A":0,"C":1,"G":2,"T":3}
INV_BASE_PAIR_ORDERING = {0:"A",1:"C",2:"G",3:"T"}

def verbose_gen(xs,n=1):
    for i,x in enumerate(xs):
        if i % n == 0:
            print i
        yield x
        
def contains_binding_sites(data):
    return all([[c in BASE_PAIR_ORDERING for c in site] for site in data])

class PSSM(list):
    """Implements a position-specific scoring matrix.  The primary
    data representation is a list of lists of log-likelihoods for each
    character, for each column, of the pssm.  The characters are
    assumed to be ordered according to BASE_PAIR_ORDERING,
    i.e. alphabetically."""


    bpo = {"A":0,"C":1,"G":2,"T":3} #base-pair ordering
    
    def __init__(self,data,background_probs = (0.25,)*4,pseudo_counts=False):
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
        def convert_column(column):
            """Return a dictionary of log odds scores for the column"""
            epsilon = 1e-10
            col_probs = normalize([column.count(b) + pseudo_counts
                                   for b in BASES])
            return {b:log(c/p + epsilon,2)
                    for (c, p, b) in zip(col_probs, background_probs,BASES)}
        def count(column):
            """Return a dictionary of counts for the column"""
            return {c:(column.count(c) + pseudo_counts) for c in BASES}
        if not type(data[0][0]) is str: #"if data is a matrix of counts..."
            self.columns = [convert_column(col) for col in data]
        else:
            assert(contains_binding_sites(data))
            counts = [count(col) for col in transpose(data)]
            self.counts = counts
            self.columns = [convert_column(col) for col in transpose(data)]
        self.motif = data
        self.length = len(self.columns)
        self.consensus = self.get_consensus()
        pseudocount = 1/float(len(data))
        print pseudocount
        self.fd_trap_columns = [{b:log((self.counts[i][self.consensus[i]]+pseudocount)/
                              float((self.counts[i][b] + pseudocount)))
                                 for b in BASES}
                                for i in range(self.length)]
        wc_motif = map(wc,self.motif)
        wc_counts = [count(col) for col in transpose(wc_motif)]
        wc_consensus = "".join([max(col,key=lambda c:col.count(c))
                        for col in transpose(wc_motif)])
        self.bk_trap_columns = [{b:log((wc_counts[i][wc_consensus[i]]+1)/
                              float((wc_counts[i][b] + 1)))
                                 for b in BASES}
                                for i in range(self.length)]
        
        
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
        if base in "AC":
            if base == 'A':
                return 0
            else:
                return 1
        else:
            if base == 'G':
                return 2
            else:
                return 3

    def base_pair_ordering2(self,base):
        if base == 'A':
            return 0
        elif base == 'C':
            return 1
        elif base == 'G':
            return 2
        else:
            return 3
    
    def base_pair_ordering3(self,base):
        return {"A":0,"C":1,"G":2,"T":3}

        
    def score(self,word,both_strands=True):
        """Return log-odds score for word"""
        fd_score = sum([col[base]
                    for col, base in zip(self.columns, word)])
        bk_score = (sum([col[base]
                    for col, base in zip(self.columns, wc(word))])
                    if both_strands else None)
        return max(fd_score,bk_score)

    def partial_thresholds(self,theta):
        """Returns a list of partial thresholds to be interpreted as
        follows: After having read position i, you must have scored at
        least pt[i], hence pt[len(pssm)] >= theta if the word is a
        positive"""
        rev_maxes = map(lambda col:max(col.values()), self.columns[::-1])
        rev_thetas = reduce(lambda ths, x: ths + [ths[-1] - x], rev_maxes, [theta])
        return rev_thetas[::-1][1:]
    
    def upper_partial_thresholds(self,theta):
        """Returns a list of partial thresholds to be interpreted as
        follows: After having read position i, if you have scored at
        least pt[i] then the word will be a positive no matter what
        the rest of the word contains"""
        mins = map(lambda col:min(col.values()), self.columns)
        n = len(self.columns)
        thetas = [theta - sum(mins[i+1:]) for i in range(n)]
        return thetas

    def sample_scores(self, n):
        return [sum(random.choice(col_dict.values())
                    for col_dict in self.columns) 
                for i in range(n)]

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
        return sum(log((self.counts[i][self.consensus[i]]+0.5)/
                       (self.counts[i][seq[i]] + 0.5))
                   for i in range(len(seq)))

    def trap(self,seq,beta=beta,both_strands=True,ns_binding=False):
        """Return the binding affinity as given by the TRAP model.
        See Manke 2008, Roider 2007.  Note that in Roider 2007, the
        trap score defined as E(s)*beta = TRAP(s), hence the TRAP
        score is dimensionless.  For consistency with other scoring
        methods, we wish to rescale the trap score by 1/beta = kbT so
        that the score will have units of kBT, and can be interpreted
        as a free energy of binding.  Therefore, we divide the final
        result by beta.  Also, we optionally implement non-specific
        binding by ensuring that -8kbt is the maximum binding score."""
        n = len(self.motif)
        w = len(self.motif[0])

        # TRAP values are reported up to an arbitrary constant, such
        # that the binding energy of the consensus sequence is
        # approximately zero for a 10 bp site. Assuming that the
        # \delta G of the consensus site is 2kbt/base, zero_delta_g =
        # 2*w, and the ns_threshold is 8kbt below that.
        ln_R_0 = 0.585 * w - 5.66
        consensus_energy_magnitude = 2 * w
        zero_delta_g = ln_R_0 + consensus_energy_magnitude
        ns_threshold = -8 #kbt
        ep_ns = consensus_energy_magnitude + ns_threshold
        ns_contrib = exp(-beta*ep_ns) if ns_binding else 0
        lamb = 0.7
        e_f = (1/lamb * sum([self.fd_trap_columns[i][seq[i]]
                          for i in range(len(seq))]) + ln_R_0)/beta 
        if both_strands:
            wc_seq = wc(seq)
            e_b = (1/lamb * sum([self.bk_trap_columns[i][seq[i]]
                          for i in range(len(seq))]) + ln_R_0)/beta
            return log(exp(-beta * e_f) + exp(-beta * e_b) + ns_contrib)/(-beta)
        else:
            return log(exp(-beta * e_f) + ns_contrib)/(-beta)

    def slide_trap(self,genome,both_strands=True,ns_binding=False):
        w = len(self.motif[0])
        return [self.trap(genome[i:i+w],both_strands=both_strands)
                           for i in verbose_gen(range(len(genome) - w + 1),100000)]
    
    def slide_score(self,genome,both_strands=True):
        w = len(self.motif[0])
        return [self.score(genome[i:i+w],both_strands=both_strands)
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
        return self.score(partial_word) + sum(map(lambda col:min(col.values()),
                                                  self.columns[l:]))

    def max_score_given(self,partial_word):
        l = len(partial_word)
        return self.score(partial_word) + sum(map(lambda col:max(col.values()),
                                                  self.columns[l:]))

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
