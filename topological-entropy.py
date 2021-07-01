#!/usr/bin/env/python3

import math
from collections import Counter
import numpy as np

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
    # highly non-optimized at this point, range should be narrower,
    # probably can be done with some kind of calculation to
    # approximate n first and then a local range can be searched.

def koslicki_definition(seq, n):
    Htop = (math.log(complexity_function(seq**(4**n +n -1)(n)), 4)/n)
    return Htop

def complexity_function(seq, n):
    #for a given sequence w, the complexity function pw:N -> N is defined as
    #pw(n) = |{u:|u|=n and u appears as a subword of w}|

main()
