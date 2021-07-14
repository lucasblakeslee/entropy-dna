#!/usr/bin/env python3

from pathlib import Path
from bio import SeqIO
γproteo_database = Path("gammaproteobacteria_averages.db_txt")
εproteo_database = Path("epsilonproteobacteria_averages.db_txt")

pseudoalteromonas_atlantica_seq = Path("pseudoalteromonas-atlantica.fasta")
#Pseudoalteromonas atlantica is in the Gammaproteobacteria class


def main():
    unknown_org_seq = seqIO.parse("pseudoalteromonas-atlantica.fasta", "fasta")
    for record in sequences:
        unknown_org_seq = record
    unknown_org_seq = str(unknown_org_seq.seq)
    find_probability_T(unknown_org_seq):
    org_type = likely_class(probability_gamma, probability_epsilon):
    print(f"the sequence in {pseudomonas_atlantica_seq} most likely belongs to an organism of type {org_type}")

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
    """Apply Poisson distribution-esque equation

    log(P(ci = ci_expected_t)) = sum(ci *ln(ci_expected_t) +k)
    ci is the count of ith subsequence ci
    t is the type of organism
    """

    unknown_seq_list = []
    for blocksize in range(20):
        unknown_seq_list.append(find_all_sequences(unknown_org_seq, blocksize))

    return probability_gamma, probability_epsilon

def likely_class(probability_gamma, probability_epsilon):
    gamma = "Gammaproteobacteria"
    epsilon = "Epsilonproteobacteria"
    unknown = "Unknown"
    if probability_gamma > probability_epsilon:
        return gamma
    if probability_epsilon > probability_gamma:
        return epsilon
    if probability_epsilon = probability_gamma:
        return unknown
    

if __name__ == "__main__":
    main()









def main():
    """Gets FASTA file from directory, reads it using SeqIO.parse, and
    then obtains sequences (outside of comments). The
    Bio.SeqRecord.SeqRecord class is very useful for storing other
    information about the sequence, but for now we're just interested
    in the sequence itself.

    """
    
    #fixme: eliminate non AGTC values
    
    out_filename = filename.with_suffix(".db_txt").name
    entropy = make_database(out_filename, seq, 10)
    print(entropy)                

        
def make_database(out_filename, sequence, max_blocksize=20, min_count=5):
    result = []
    with open(out_filename, 'w') as out_file:
        for blocksize in range(1, max_blocksize+1):
            my_counts = count_all_subsequences(sequence, blocksize)
            for item, count in my_counts.most_common():
                if count < min_count:
                    break
                out_file.write(f"{item:s} {count}\n")
            entropy = np.sum(count*np.log(count) for count in my_counts.values()) *1/len(sequence)
            result.append((blocksize, entropy))
    return np.array(result)



