#!/usr/bin/env python3

from pathlib import Path
from Bio import SeqIO
from collections import Counter
import numpy as np

gammaproteo_database = Path("gammaproteobacteria_averages.db_txt")
epsilonproteo_database = Path("epsilonproteobacteria_averages.db_txt")

database_files = [gammaproteo_database, epsilonproteo_database]

pseudoalteromonas_atlantica_seq = Path("pseudoalteromonas-atlantica.fasta")
# Pseudoalteromonas atlantica is in the Gammaproteobacteria class


def main():
    # unknown_org_seq = SeqIO.parse("pseudoalteromonas-atlantica.fasta", "fasta")[0]
    # 'FastaIterator' object is not scriptable
    unknown_org_seq = []
    for record in SeqIO.parse("pseudoalteromonas-atlantica.fasta", 'fasta'):
        unknown_org_seq.append(str(record.seq))
        break  # we only want the first one
    unknown_org_seq = str(unknown_org_seq)
    unknown_org_seq = unknown_org_seq.upper()
    blocksize = 10
    unknown_counts = count_all_subsequences(unknown_org_seq, blocksize)
    for db_file in database_files:
        db_counts = read_count_database(db_file, blocksize)
        prob = find_probability_vs_candidate(unknown_counts, db_counts)
        print(prob, db_file.name)
    


def find_probability_vs_candidate(test_counts, candidate_counts):
    test_sum = np.sum(test_counts.values())
    candidate_sum = np.sum(candidate_counts.values())
    print(test_sum)
    expectation_factor = test_sum/candidate_sum
    missing_expectation = expectation_factor*1/2
    counts_and_expectations = [(count,
                                (expectation_factor*candidate_counts[subsequence])
                                if subsequence in candidate_counts
                                else missing_expectation)
                                for subsequence, count in test_counts.items()]
    log_probabilities = [counts*np.log(expectation) for counts, expectation in counts_and_expectations]
    log_probability = np.sum(log_probabilities)
    return log_probability

def read_count_database(filename, blocksize):
    result = {}
    for line in open(filename, 'r').readlines():
        subsequence, counts = line.split()
        subsequence = subsequence.strip()
        if len(subsequence) == blocksize:
            result[subsequence] = counts
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


    """Apply Poisson distribution
    log(P(ci = ci_expected_t)) = sum(ci *ln(ci_expected_t) +k)
    ci is the count of ith subsequence ci
    t is the type of organism
    """

    

if __name__ == "__main__":
    main()
