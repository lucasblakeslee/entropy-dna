#!/usr/bin/env/python3

import math
from collections import Counter
import numpy as np
from itertools import combinations

seq = open("populus_deltoides_seq.txt", "r")
seq = seq.read()

#it's worth noting that this is an exon; it's actively expressed, so
#we would expect the entropy to be lower than an intron which aren't
#subject to evolutionary pressures.

randomseq = "".join(np.random.choice(["A","T","C","G"], size=len(seq)))

seq_list = list(seq)

seq_len = len(seq)

def main():
    n = find_n(seq_len)
    koslicki_definition(seq, n)

def find_n(seq_len):
    # let seq be a finite sequence of length |seq| (here len(seq)),
    # let n be the unique integer such that:
    # 4**(n) + n - 1 <= len(seq) < 4**(n+1) + (n + 1) - 1
    
    for n in range(seq_len): 
        if 4**n + n - 1 <= seq_len and 4**(n+1) + (n+1) -1 > seq_len:
            return n
    
def koslicki_definition(seq, n):
    P = complexity_function(seq, n)
    topological_entropy = ((math.log(P**(4**n +n -1),4))/(n))
    return topological_entropy

def complexity_function(seq, blocksize):
    #for a given sequence w, the complexity function pw:N -> N is defined as
    #p(n) = |{u:|u|=n and u appears as a subword of w}|
    #p(n) should be a representation of the number of n-length subwrds that appear in w

    combos = [seq[i:(i+blocksize)] for i in range(len(seq)+1-blocksize)]
    unique_combos = set(combos)
    return(len(unique_combos))

if __name__ == "__main__":
    # execute only if run as a script
    main()
