#!/bin/usr/env/Python3

from collections import Counter
import numpy as np
from pathlib import Path
from Bio import SeqIO

# Part 1 -- obtaining organisms' sequences

#1a. gammaproteobacteria
#"gammaproteobacteria/legionella-pneumophila.fasta"
#"gammaproteobacteria/amphritea-japonica.fasta"
#"gammaproteobacteria/haemophilus-influenzae.fasta"
#"gammaproteobacteria/xanthomonas-citri.fasta"

#1b. epsilonproteobacteria
#"epsilonproteobacteria/caminibacter-mediatlanticus.fasta"
#"epsilonproteobacteria/helicobacter-pylori.fasta"
#"epsilonproteobacteria/campylobacter-jejuni.fasta"
#"epsilonproteobacteria/sulfurovum-lithotrophicum.fasta"

# Part 2 -- making databases

def main():
    """Gets FASTA file from directory, reads it using SeqIO.parse, and
    then obtains sequences (outside of comments). The
    Bio.SeqRecord.SeqRecord class is very useful for storing other
    information about the sequence, but for now we're just interested
    in the sequence itself.

    """
    in_filename = "gammaproteobacteria/legionella-pneumophila.fasta"
    filename = Path('gammaproteobacteria/legionella-pneumophila.fasta')
    sequences = SeqIO.parse(in_filename, 'fasta')
    for record in sequences:
        example = record
    seq = str(example.seq)
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


def count_all_subsequences(seq, blocksize):
    counter = Counter(seq[i:(i+blocksize)] for i in range(len(seq)+1-blocksize))
    return counter


def find_all_subsequences(seq, blocksize):
    """Finds a list of all subsequences contained within the sequence
    'seq' of length 'blocksize.'

    """
    combos = [seq[i:(i+blocksize)] for i in range(len(seq)+1-blocksize)]
    return combos

if __name__ == "__main__":
    main()
