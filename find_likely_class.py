#!/usr/bin/env python3

from pathlib import Path
from Bio import SeqIO
from collections import Counter
import numpy as np
import random
import re

# gammaproteo_database = Path("gammaproteobacteria_databases_averages.db_txt")
# epsilonproteo_database = Path("epsilonproteobacteria_databases_averages.db_txt")
#
# database_files = [gammaproteo_database, epsilonproteo_database]

database_files = list(Path().glob("*databases_averages.db_txt"))


database_files = list(Path().glob("*databases_averages.db_txt"))
tests_directory = Path("test_genomes")

def main():
    for subdirectory in tests_directory.iterdir():
        for filename in subdirectory.iterdir():
            unknown_org_seq = []
            for record in SeqIO.parse(filename, 'fasta'):
                unknown_org_seq.append(str(record.seq))
                break  # we only want the first one
            unknown_org_seq = "".join(unknown_org_seq).upper()
            # first_half, second_half = unknown_org_seq[:len(unknown_org_seq)//2],
                                      # unknown_org_seq[len(unknown_org_seq)//2:]
            # scales with sequence length
            blocksize = 15
            unknown_counts = count_all_subsequences(unknown_org_seq, blocksize)
            # can input random_seq(2000000)
            # also first_half and second_half
            print(subdirectory, filename.name)
            for db_file in database_files:
                db_counts = read_count_database(db_file, blocksize)
                prob = find_probability_vs_candidate(unknown_counts, db_counts)
                print(prob, db_file.name)

def find_probability_vs_candidate(test_counts, candidate_counts):
    test_sum = sum(test_counts.values())
    candidate_sum = sum(candidate_counts.values())
    expectation_factor = test_sum/candidate_sum
    missing_expectation = expectation_factor*1/2
    counts_and_expectations = [(count,
                                (expectation_factor*candidate_counts[subsequence])
                                if subsequence in candidate_counts
                                else missing_expectation)
                                for subsequence, count in test_counts.items()]
    log_probabilities = [counts*np.log(expectation) for counts, expectation in counts_and_expectations]
    log_probability = -np.sum(log_probabilities)
    return log_probability

def read_count_database(filename, blocksize):
    result = {}
    for line in open(filename, 'r').readlines():
        subsequence, counts = line.split()
        subsequence = subsequence.strip()
        if len(subsequence) == blocksize:
            result[subsequence] = float(counts)
    return result


def count_all_subsequences(seq, blocksize):
    counter = Counter(seq[i:(i+blocksize)] for i in range(len(seq)+1-blocksize))
    return counter


def find_all_subsequences(seq, blocksize):
    """Finds a list of all subsequences contained within the sequence
    'seq' of length 'blocksize.'

    """
    combos = [seq[i:(i+blocksize)] for i in range(len(seq)+1-blocksize)]
    return combos

def random_seq(seq_len):
    return ''.join(random.choice('CGTA') for _ in range(seq_len))

if __name__ == "__main__":
    main()
