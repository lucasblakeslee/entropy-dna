#!/bin/usr/env/Python3

from collections import Counter

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
    Makes a huge file, do not run
    """
    with open("epsilonproteobacteria_database.txt", 'a') as file1:
        for blocksize in range(3):
            a = find_all_subsequences(sulfurovum_lithotrophicum_seq, blocksize)
            my_dict = {i:a.count(i) for i in a}
            for item in my_dict:
                file1.write("%s\n" % item)

def find_all_subsequences(seq, blocksize):
    """Finds a list of all subsequences contained within the sequence
    'seq' of length 'blocksize.'

    """
    combos = [seq[i:(i+blocksize)] for i in range(len(seq)+1-blocksize)]
    return combos

if __name__ == "__main__":
    main()
