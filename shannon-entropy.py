#!/usr/bin/env/python3

import math
import numpy as np
import collections
from collections import Counter

seq = "tatgagtcatcacataaaagctaatgtagtaaagcatccataccaattgatttaaaattaaactaatccgttgataaaacaaaaaatcaagagccttgagtcactccagactgagaagattgactcaagaacaaattttatttttattttattagaggctccattggaaaaataagacctaaacattaatcgagaagtgatgggaacgatggaacctgtaaatacataagatgttactgaaaacaaatcttaatgatttattgggaggatagcgaaatgaaccagaagacaatcgatttattttattatgagagatcatgggttacctctacaaataaaataactattgaaagactaaatatgcaagaggaaatattcgcccgcgaagattttttttatattctttttgattgctaaatttgatcaatctgatatcctaattcgaatatataaaaaaatttgttacaaccttataacaaataacaaaaacactattattgaaaatcaaaaaataaaatagaaaaaattctctataattataattaatgtgtataactagaataattttttctatttctatctataactattagatctagaactagtagaactctaaaatagaatatagattctaatttatatatattagatacaaatttatattctactaatattctattctacctaatatcctattctaataatctaagattctaatactaataaatagatcgaataagtaagaaataaattaaaataaatagatttaactaaattaagtgaaatctcaaagaatacgatgatttaatatattattttattcgtaaaagacatggatatttttttttaatcatttcattcgcgaggagctggatgagaagaaattctcatgtccggttctgtagtagagatggaattgagaaataaccatcaactataaccccaaaagaacccgattccgtaaacaacatagaggaagaatgaaaggaatagcttttcgaggaaatcgtatttgttttggaagatatgctcttcaagcacttgaatccgcttggattacatctaga"
#from Euphorbia iharanae genome, arbitrary choice

seq_list = list(seq)
#splits the string into a list

def Shannon_entropy_func(sequence, blocksize):
    N=len(sequence) -blocksize +1
    Block=[sequence[i:i +blocksize] for i in range(N)]
    entropy = sum([math.log(float(N)/Block.count(b)) for b in Block])/N
    print(entropy)

# maybe this should be return(entropy)

for blocksize in range(500):
    Shannon_entropy_func(seq_list, blocksize)
    blocksize += 1

#goes through every block length 1 - 500 (somewhat arbitrary max
#limit, need to look into where is the proper place to stop) and
#calculates entropy for that block.

#Still need to sum these all up








#----------------some old portions--------------------------

#n = len(seq)
# n = math.log((len(seq)), 4)

# n = int(n)
# results in approximately log base-4 of len(s) Here is
# where things feel illegal...the below arguments require an integer
# for n, whereas this n is a float. Converting i to the integer n gets
# rid of all the spicy information contained within, so we'll need to
# find a different way to do this.

# counts = Counter([s[i:i+n] for i in range(len(seq) - n)])
# iterates through every possible 5 letter "word" that appears in
# the original sequence, with overlaps
# in this case, 'aaata' appears the greatest number of times (15).
# shouldn't this iterate through words of every length rather than just 5 letter ones?

# shannon_entropy = Sfunction(sum([count * math.log(count) for count in counts.values()]))
