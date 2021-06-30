#!/usr/bin/env/python3

import math
import numpy as np
import collections
from collections import Counter

seq = "tatgagtcatcacataaaagctaatgtagtaaagcatccataccaattgatttaaaattaaactaatccgttgataaaacaaaaaatcaagagccttgagtcactccagactgagaagattgactcaagaacaaattttatttttattttattagaggctccattggaaaaataagacctaaacattaatcgagaagtgatgggaacgatggaacctgtaaatacataagatgttactgaaaacaaatcttaatgatttattgggaggatagcgaaatgaaccagaagacaatcgatttattttattatgagagatcatgggttacctctacaaataaaataactattgaaagactaaatatgcaagaggaaatattcgcccgcgaagattttttttatattctttttgattgctaaatttgatcaatctgatatcctaattcgaatatataaaaaaatttgttacaaccttataacaaataacaaaaacactattattgaaaatcaaaaaataaaatagaaaaaattctctataattataattaatgtgtataactagaataattttttctatttctatctataactattagatctagaactagtagaactctaaaatagaatatagattctaatttatatatattagatacaaatttatattctactaatattctattctacctaatatcctattctaataatctaagattctaatactaataaatagatcgaataagtaagaaataaattaaaataaatagatttaactaaattaagtgaaatctcaaagaatacgatgatttaatatattattttattcgtaaaagacatggatatttttttttaatcatttcattcgcgaggagctggatgagaagaaattctcatgtccggttctgtagtagagatggaattgagaaataaccatcaactataaccccaaaagaacccgattccgtaaacaacatagaggaagaatgaaaggaatagcttttcgaggaaatcgtatttgttttggaagatatgctcttcaagcacttgaatccgcttggattacatctaga"
randomseq = "".join(np.random.choice(["a","t","c","g"], size=len(seq)))
#from Euphorbia iharanae genome, arbitrary choice

seq_list = list(seq)
#splits the string into a list

def Shannon_entropy_func(sequence, blocksize):
    N=len(sequence) -blocksize +1
    Block=[sequence[i:i +blocksize] for i in range(N)]
    entropy = sum([math.log(float(N)/Block.count(b)) for b in Block])/N
    return (entropy)

def Shannon_entropy_dmp(sequence, blocksize):
    N = len(sequence) - blocksize + 1
    counter = Counter([sequence[i:i+blocksize] for i in range(N)])
    counts = np.array(list(counter.values()))

    # The following normalizations are to match Shannon_entropy_func
    # and are not necessarily correct
    normalization = -1/N
    normalization_add = np.log(N)

    entropy = np.sum(counts * np.log(counts)) * normalization + normalization_add
    return entropy


for blocksize in range(1, 20):
    # Using seq instead of seq_list lets you treat the subsequences as
    # strings instead of lists, which is faster
    slowvalue = Shannon_entropy_func(seq, blocksize)
    dmpvalue = Shannon_entropy_dmp(seq, blocksize)
    randvalue = Shannon_entropy_dmp(randomseq, blocksize)
    print(f"{blocksize:10d} {slowvalue:8.3f}  {dmpvalue:8.3f}     {randvalue:8.3f}")

    #    blocksize += 1

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
