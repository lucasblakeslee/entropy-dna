# Entropy in DNA

## Motivation

In information theory, *entropy* can be thought of as a measure of the
"information" contained within a system. Oftentimes this can be interpreted as a
measure of disorder, or of the amount of surprise regarding the next piece of
information. For example, if one has the letter A repeated 100 times, the
sequence does not contain much information...there also isn't much surprise if
you were to have to guess the next letter from a random point.

The entropy of a sequence of DNA can give insight into a variety of information
about that sequence, namely it can serve as an estimator of coding and
non-coding regions of DNA, what are called *exons* and *introns* respectively.
Intron regions, not being expressed, are not subject to the same evolutionary
pressures under which exon regions are held (though recent literature suggests
that intron regions do play an integral role in gene expression regulation and
have been linked to some genetic diseases), so one might expect intron regions
to have a higher entropy (implying greater disorder) than exon regions.

Approximations for the *Shannon* and *topological* entropy of DNA sequences have
been described by Schmitt & Herzel (1992) and Koslicki (2011) respectively.
*Shannon* and *topological* entropy are calculated slightly differently, the
major difference being that topological entropy ignores the frequency with which
different sequences occur, attempting to get around the fact that with Shannon
entropy there is a positive correlation between the length of a sequence and its
entropy.

Previous methods for calculating the entropy of DNA sequences assume that
nucleotides and subsequences 1. are identically distributed, and
2. occur independently of one another. However, DNA is decidedly not random, and
even a cursory look at the structures of DNA will tell you that certain regions
will not follow the above assumptions, such as promoter regions which occur
upstream of a coding dna region, which lets proteins know where to bind to in
order to begin translation. This project seeks to calculate the entropy of a
genetic sequence using a novel method that takes into account common sequences
in other organisms of the same class, creating an entropy reduction to be
subtracted from the raw entropy (which can be calculated using either Shannon or
topological entropy). This should give a more accurate depiction of the entropy
of a DNA sequence, as it takes into account the expected common subsequences for
its class, reducing the amount of surprise.


## Code

WIP, but there currently exist programs to calculate the raw Shannon and
topological entropies of a sequence, as well as a program to make a database of
recurring subsequences for individual sequences, and then a separate program
that parses through a directory to find the average number of times a
subsequence occurs in a class (including those that occur 0 times).

Sequences have been obtained from the NCBI website's genome database
