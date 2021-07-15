#!/usr/bin/env python3

from pathlib import Path
from Bio import SeqIO

γproteo_database = Path("gammaproteobacteria_averages.db_txt")
εproteo_database = Path("epsilonproteobacteria_averages.db_txt")

pseudoalteromonas_atlantica_seq = Path("pseudoalteromonas-atlantica.fasta")
# Pseudoalteromonas atlantica is in the Gammaproteobacteria class


def main():
    unknown_org_seq = SeqIO.parse("pseudoalteromonas-atlantica.fasta", "fasta")
    for record in unknown_org_seq:
        unknown_org_seq = record
    unknown_org_seq = str(unknown_org_seq.seq)
    find_probability_T(unknown_org_seq)
    org_type = likely_class(probability_gamma, probability_epsilon)
    print("the sequence in {pseudoalteromonas_atlantica_seq} most likely belongs to an organism of type {org_type}".format)

def count_all_subsequences(seq, blocksize):
    counter = Counter(seq[i:(i+blocksize)] for i in range(len(seq)+1-blocksize))
    return counter


def find_all_subsequences(seq, blocksize):
    """Finds a list of all subsequences contained within the sequence
    'seq' of length 'blocksize.'

    """
    combos = [seq[i:(i+blocksize)] for i in range(len(seq)+1-blocksize)]
    return combos


def find_probability_T(unknown_org_seq):
    """Apply Poisson distribution

    log(P(ci = ci_expected_t)) = sum(ci *ln(ci_expected_t) +k)
    ci is the count of ith subsequence ci
    t is the type of organism
    """

    unknown_seq_list = []
    for blocksize in range(20):
        unknown_seq_list.append(find_all_subsequences(unknown_org_seq, blocksize))

    return probability_gamma, probability_epsilon

def likely_class(probability_gamma, probability_epsilon):
    gamma = "Gammaproteobacteria"
    epsilon = "Epsilonproteobacteria"
    unknown = "Unknown"
    if probability_gamma > probability_epsilon:
        return gamma
    if probability_epsilon > probability_gamma:
        return epsilon
    if probability_epsilon == probability_gamma:
        return unknown
    

if __name__ == "__main__":
    main()
