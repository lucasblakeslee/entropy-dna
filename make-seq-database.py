#!/bin/usr/env/Python3

from collections import Counter
import numpy as np
from pathlib import Path


"""goal of this program is to:

1. get a representative sample of genomes for a certain organism class
from the NCBI Genome website (I plan to do this in the future,
currently just using downloaded FASTA files)

2. find the most commonly occuring sequences occuring in those genomes
and rank them (including those that appear 0 times)

3. write that information to a file (still need to decide on format
but leaning towards csv with column 1 being the sequence and column 2
being the number of times it occurs

"""

# Part 1 -- obtaining organisms' sequences

#1a. gammaproteobacteria

#i. Legionella pneumophila sequence
legionella_pneumophila_seq = open("gammaproteobacteria/legionella-pneumophila.fasta", "r")
next(legionella_pneumophila_seq)                  # advanced file pointer to next line
legionella_pneumophila_seq = legionella_pneumophila_seq.read()
legionella_pneumophila_seq = legionella_pneumophila_seq.replace("\n", "")

#ii. Amphritea japonica sequence
amphritea_japonica_seq = open("gammaproteobacteria/amphritea-japonica.fasta", "r")
next(amphritea_japonica_seq)
amphritea_japonica_seq = amphritea_japonica_seq.read()
amphritea_japonica_seq = amphritea_japonica_seq.replace("\n", "")

#iii. Haemophilus influenzae sequence
haemophilus_influenzae_seq = open("gammaproteobacteria/haemophilus-influenzae.fasta", "r")
next(haemophilus_influenzae_seq)
haemophilus_influenzae_seq = haemophilus_influenzae_seq.read()
haemophilus_influenzae_seq = haemophilus_influenzae_seq.replace("\n", "")

#iv. Xanthomonas citri sequence
xanthomonas_citri_seq = open("gammaproteobacteria/xanthomonas-citri.fasta", "r")
next(xanthomonas_citri_seq)
xanthomonas_citri_seq = xanthomonas_citri_seq.read()
xanthomonas_citri_seq = xanthomonas_citri_seq.replace("\n", "")

#-----------------------------------------------------------------
#1b. epsilonproteobacteria

#i. Caminibacter mediatlanticus sequence
caminibacter_mediatlanticus_seq = open("epsilonproteobacteria/caminibacter-mediatlanticus.fasta", "r")
next(caminibacter_mediatlanticus_seq)
caminibacter_mediatlanticus_seq = caminibacter_mediatlanticus_seq.read()
caminibacter_mediatlanticus_seq = caminibacter_mediatlanticus_seq.replace("\n", "")

#ii. Helicobacter pylori sequence
helicobacter_pylori_seq = open("epsilonproteobacteria/helicobacter-pylori.fasta", "r")
next(helicobacter_pylori_seq)
helicobacter_pylori_seq = helicobacter_pylori_seq.read()
helicobacter_pylori_seq = helicobacter_pylori_seq.replace("\n", "")

#iii. Campylobacter jejuni sequence
campylobacter_jejuni_seq = open("epsilonproteobacteria/campylobacter-jejuni.fasta", "r")
next(campylobacter_jejuni_seq)
campylobacter_jejuni_seq = campylobacter_jejuni_seq.read()
campylobacter_jejuni_seq = campylobacter_jejuni_seq.replace("\n", "")

#iv. Sulfurovum lithotrophicum sequence
sulfurovum_lithotrophicum_seq = open("epsilonproteobacteria/sulfurovum-lithotrophicum.fasta", "r")
next(sulfurovum_lithotrophicum_seq)
sulfurovum_lithotrophicum_seq = sulfurovum_lithotrophicum_seq.read()
sulfurovum_lithotrophicum_seq = sulfurovum_lithotrophicum_seq.replace("\n", "")

###################################################################################################
# Part 2: Find all subsequences #


def main():
    """
    """
    filename = Path("epsilonproteobacteria/sulfurovum-lithotrophicum.fasta")
    seq = read_sequence(filename)
    out_filename = filename.with_suffix(".db_txt").name
    entropy = make_database(out_filename, seq, 10)
    print(entropy)
    
def read_sequence(filename):
    seqfile = open(filename, "r")
    _=seqfile.readline()
    seq = seqfile.read()
    seq = seq.replace("\n", "")
    seq = seq.upper()
    #fixme: eliminate non ATCG values
    return seq
                

        
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
