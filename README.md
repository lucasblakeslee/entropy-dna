# Entropy in DNA

See paper at https://github.com/lucasblakeslee/entropy-dna/blob/master/entropy-paper.pdf

## Code

Running the code requires that some programs be executed before others.
First:

```
python3 make-seq-database.py
```

in order to make the sequence count databases

Then, you can run

```
python3 find_average_sequences.py
```

in order to create averaged class databases

Lastly,

```
python3 func_test.py
```

in order to test the probabilities of the test genomes falling into
either class

Feel free to try this with different classes and species!

## Acknowledgements
I would like to acknowledge my mentor, David Palmer, for the
invaluable help and guidance he provided me over the course of this
project.

Additionally, I thank the Institute for Computing in Research for
providing this opportunity

Sequences have been obtained from the NCBI website's genome database

Finally, thanks to all the maintainers of the open-source scientific
libraries I used, namely NumPy, BioPython, and Matplotlib!
