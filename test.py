at_lab = True
import time
import ESA,PSSM

def get_ecoli_genome():
    lab_file = "/home/poneill/ecoli/NC_000913.fna"
    with open(lab_file if at_lab else home_file) as f:
        genome = "".join([line.strip() for line in f.readlines()[1:]])
    return "".join(g for g in genome if g in "ATGC") # contains other iupac symbols

def get_crp_motif():
    lab_file = "/home/poneill/euksites/crp.csv"
    home_file = "/home/pat/Dropbox/entropy/crp.txt"
    with open(lab_file if at_lab else home_file) as f:
        lines = [line.strip() for line in f.readlines()[1:]]
    return [line.replace(",","") for line in lines]

def main():
    esa = ESA.ESA(get_ecoli_genome())
    motif = get_crp_motif()
    pssm = PSSM.PSSM(motif)
    results = pssm.search_esa(esa,20)
    print results
    
