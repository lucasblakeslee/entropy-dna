#!bin/bash/usr/env python3
import math
import numpy as np
from pathlib import Path
from Bio import SeqIO
from collections import Counter
import matplotlib.pyplot as plt

database_files = list(Path().glob("*databases_averages.db_txt"))
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

            blocksize = 15
            unknown_counts = count_all_subsequences(unknown_org_seq, blocksize)
            for db_file in database_files:
                expected_class_counts = read_count_database(db_file, blocksize)
                result = func(unknown_counts, expected_class_counts)
                results[f"result_{db_file}_{filename}"] = result
                print(result, db_file.name)
    for key,value in results.items():
        if key.startswith('result_epsilonproteobacteria'):
            x = value
            plt.scatter(x, 1, c='b')
        else:
            y = value
            plt.scatter(y, 1, c='g')
    plt.show()


def func(unknown_counts, expected_class_counts):
    """Given the (C_i, E_i) values of (counts in unknown sequence, expected counts
        averaged from a class) we can define a function:"""
    unknown_values, expected_values = find_counts_and_expectations(unknown_counts, expected_class_counts)
    unknown_log = np.log(unknown_values)
    expected_log = np.log(expected_values)
    normalization = math.sqrt(np.sum(unknown_log ** 2) * np.sum(expected_log**2))
    result = np.sum(unknown_log * expected_log)/normalization
    return result


def find_counts_and_expectations(test_counts, candidate_counts):
    test_sum = sum(test_counts.values())
    candidate_sum = sum(candidate_counts.values())
    expectation_factor = test_sum / candidate_sum
    missing_expectation = expectation_factor * 1 / 2
    counts_and_expectations = [(count,
                                (expectation_factor * candidate_counts[subsequence])
                                if subsequence in candidate_counts
                                else missing_expectation)
                               for subsequence, count in test_counts.items()]
    counts_and_expectations_np = np.array(counts_and_expectations)
    return counts_and_expectations_np[:,0], counts_and_expectations_np[:,1]


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

def plot_values(epsilon, gamma):
    plot

if __name__ == "__main__":
    main()

