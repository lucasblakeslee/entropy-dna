#!/usr/bin/env/python3

import math
import numpy as np
import collections
from collections import Counter
import matplotlib.pyplot as plt
from pathlib import Path
from Bio import SeqIO

#it's worth noting that this is an exon; it's actively expressed, so
#we would expect the entropy to be lower than an intron which aren't
#subject to evolutionary pressures.


dmpvalue_list = []
randvalue_list = []
filename_list = []


tests_directory = Path("test_genomes")

def main():
    results = {}
    for subdirectory in tests_directory.iterdir():
        for filename in subdirectory.iterdir():
            unknown_org_seq = []
            for record in SeqIO.parse(filename, 'fasta'):
                unknown_org_seq.append(str(record.seq))
                break  # we only want the first one
            unknown_org_seq = "".join(unknown_org_seq).upper()
            randomseq = "".join(np.random.choice(["A","T","C","G"], size=len(unknown_org_seq)))
        for blocksize in range(1, 15):
            # Using seq instead of seq_list lets you treat the subsequences as
            # strings instead of lists, which is faster
            dmpvalue = Shannon_entropy_dmp(unknown_org_seq, blocksize)
            randvalue = Shannon_entropy_dmp(randomseq, blocksize)
            dmpvalue_list.append(dmpvalue)
            randvalue_list.append(randvalue)
            filename_list.append(filename.name)

#    print(f"{filename}    {blocksize:10d}     {dmpvalue:8.3f}     {randvalue:8.3f}")

#    print(filename_list)
    print(dmpvalue_list)
#    plt.scatter(dmpvalue_list, randvalue_list)
    plt.plot(dmpvalue_list, randvalue_list)
    plt.xlabel("Block size")
    plt.ylabel("Measured entropy")
    plt.show()


def Shannon_entropy_dmp(sequence, blocksize):
    N = len(sequence) - blocksize + 1
    counter = Counter([sequence[i:i+blocksize] for i in range(N)])
    counts = np.array(list(counter.values()))

    normalization = 1/N
    normalization_add = np.log(N)

    # it seems the way to normalize entropy is by dividing it by
    # information length, the ratio of which is metric entropy

    entropy = - np.sum(counts * np.log(counts)) * normalization + normalization_add
    return entropy


def Shannon_entropy_func(sequence, blocksize):
    """
    inefficient, not currently in use
    """
    N=len(sequence) -blocksize +1
    Block=[sequence[i:i +blocksize] for i in range(N)]
    entropy = sum([math.log(float(N)/Block.count(b)) for b in Block])/N
    return (entropy)


if __name__ == "__main__":
    # execute only if run as a script
    main()
