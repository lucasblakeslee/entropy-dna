#!bin/bash/usr/env python3
import math
import numpy as np
from pathlib import Path
from Bio import SeqIO
from collections import Counter

database_files = list(Path().glob("*databases_averages.db_txt"))
test_sequence = "xanthomonas-citri.fasta"


def main():
    unknown_org_seq = []
    for record in SeqIO.parse("test_genomes/gamma_test_genomes/"+test_sequence, 'fasta'):
        unknown_org_seq.append(str(record.seq))
        break  # we only want the first one
    unknown_org_seq = "".join(unknown_org_seq).upper()

    blocksize = 15
    unknown_counts = count_all_subsequences(unknown_org_seq, blocksize)
    for db_file in database_files:
        expected_class_counts = read_count_database(db_file, blocksize)
        result = func(unknown_counts, expected_class_counts)
        print(result, db_file.name)


def func(unknown_counts, expected_class_counts):
    """Given the (C_i, E_i) values of (counts in unknown sequence, expected counts
averaged from a class) we can define a function:"""
    unknown_counts_values = unknown_counts.values()
    expected_counts_values = expected_class_counts.values()
    logs_of_unknown_values = []
    logs_of_counts_values = []
    for entry in unknown_counts_values:
        logs_of_unknown_values.append(math.log(entry))
    for entry in expected_counts_values:
        logs_of_counts_values.append(math.log(entry))
    logs_unknown_arr = np.array(logs_of_unknown_values)
    logs_counts_arr = np.array(logs_of_counts_values)
    normalization = math.sqrt(np.sum(logs_unknown_arr)**2* np.sum(logs_counts_arr)**2)
    result = np.sum(logs_unknown_arr * logs_counts_arr)/normalization
    return result

def count_all_subsequences(seq, blocksize):
    counter = Counter(seq[i:(i+blocksize)] for i in range(len(seq)+1-blocksize))
    return counter

def read_count_database(filename, blocksize):
    result = {}
    for line in open(filename, 'r').readlines():
        subsequence, counts = line.split()
        subsequence = subsequence.strip()
        if len(subsequence) == blocksize:
            result[subsequence] = float(counts)
    return result


if __name__ == "__main__":
    main()

