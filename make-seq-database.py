#!/bin/usr/env python3

from collections import Counter
import numpy as np
from pathlib import Path
from Bio import SeqIO
import os

# Part 1 -- obtaining organisms' sequences

#1a. gammaproteobacteria
#"gammaproteobacteria/thiohalobacter-thiocyanaticus.fasta"
#"gammaproteobacteria/allofrancisella-guangzhouensis.fasta"
#"gammaproteobacteria/haemophilus-influenzae.fasta"
#"gammaproteobacteria/francisella-opportunistica.fasta"
#"gammaproteobacteria/pasteurella-multocida.fasta"

#1b. epsilonproteobacteria
#"epsilonproteobacteria/caminibacter-mediatlanticus.fasta"
#"epsilonproteobacteria/helicobacter-pylori.fasta"
#"epsilonproteobacteria/campylobacter-jejuni.fasta"
#"epsilonproteobacteria/sulfurovum-lithotrophicum.fasta"
#"epsilonproteobacteria/arcobacter-venerupis.fasta"

# Part 2 -- making databases

def main():
    """Gets FASTA file from directory, reads it using SeqIO.parse, and
    then obtains sequences (outside of comments). The
    Bio.SeqRecord.SeqRecord class is very useful for storing other
    information about the sequence, but for now we're just interested
    in the sequence itself.
    """

    directory_names = (Path("epsilonproteobacteria"), Path("gammaproteobacteria"))
    for directory_name in directory_names:
        dbdirname = directory_name.parent.joinpath(directory_name.name + "_databases")
        for filename in directory_name.glob('*.fasta'):
            with open(filename) as in_filename:
                sequences = SeqIO.parse(in_filename, 'fasta')
                for record in sequences:
                    example = record
                seq = str(example.seq).upper()
                #fixme: eliminate non AGTC values
                #currently done in terminal with sed -e '/^[^>]/s/[^ATGCatcg]/N/g' filename.fasta
                out_filename = dbdirname.joinpath(filename.with_suffix(".db_txt").name)
                entropy = make_database(out_filename, seq, max_blocksize=15, min_blocksize=15, min_count=0)
                print(f"{filename.name}:")
                print(entropy)


def make_database(out_filename, sequence, max_blocksize=20, min_blocksize=1, min_count=5):
    result = []
    with open(out_filename, 'w') as out_file:
        for blocksize in range(min_blocksize, max_blocksize+1):
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
