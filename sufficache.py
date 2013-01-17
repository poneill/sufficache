"""
Main module for sufficache library.

Usage

(1) Construct a PSSM object from a list of binding sites or pssm
matrix using sufficache.PSSM.

(2) Construct an ESA object from a genome to be searched

(3) Search the ESA with the PSSM, returning a list of
(site_index,pssm_score) tuples.

Example:

> import sufficache

Generate a fake motif of ten random octamers.

> motif = [sufficache.utils.random_site(8) for i in range(10)]

> pssm = sufficache.PSSM(motif)

Generate a fake genome:

> genome = sufficache.utils.random_site(1000)

> esa = sufficache.ESA(genome)

Search the genome:

> cutoff = 3.8 #let's say

> putative_binding_sites = pssm.search_esa(esa,cutoff)

> putative_binding_sites
# [(100, 3.916428337230777), (165, 3.916428337230777), (262, 3.916428337230777), (768, 4.086353338673089), (85, 3.916428337230777)]
"""

import sufficache_utils
from ESA import ESA
from PSSM import PSSM


