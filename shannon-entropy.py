#!/usr/bin/env/python3

import math
import numpy as np
import collections
from collections import Counter
import matplotlib.pyplot as plt


seq = open("populus_deltoides_seq.txt", "r")
seq = seq.read()

#it's worth noting that this is an exon; it's actively expressed, so
#we would expect the entropy to be lower than an intron which aren't
#subject to evolutionary pressures.

randomseq = "".join(np.random.choice(["A","T","C","G"], size=len(seq)))

dmpvalue_list = []
randvalue_list = []

def main():
    for blocksize in range(1, 20):
        # Using seq instead of seq_list lets you treat the subsequences as
        # strings instead of lists, which is faster

        """
        The block size should probably be around the base-4 log of the
        total length. This for a random sequence would give about 1 copy on
        average of each possible substring. Therefore if you have a few
        strings that have many more than one copy, the sequence is not
        random. You can also see how many copies of strings you get that are
        several times that.
        """
        
        dmpvalue = Shannon_entropy_dmp(seq, blocksize)
        randvalue = Shannon_entropy_dmp(randomseq, blocksize)
        
        dmpvalue_list.append(dmpvalue)
        randvalue_list.append(randvalue)

        print(f"{blocksize:10d}     {dmpvalue:8.3f}     {randvalue:8.3f}")
    plt.plot(dmpvalue_list)
    plt.plot(randvalue_list)
    plt.show()


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


def Shannon_entropy_func(sequence, blocksize):
    """
    inefficient, not currently in use
    """
    N=len(sequence) -blocksize +1
    Block=[sequence[i:i +blocksize] for i in range(N)]
    entropy = sum([math.log(float(N)/Block.count(b)) for b in Block])/N
    return (entropy)



main()
