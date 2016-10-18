#!/usr/bin/env python

"""
Python course EBC 2016
Day 5 - Exercise instructions
Lucas Sinclair <lucas.sinclair@me.com>

An exercise designed to practice writing classes in python and parsing complex files.

This current file, `exercise_day5.py`, contains commands that are not to be modified.
Instead, the student should create a file called `fasta_parser.py` and place it in the same directory as this file. He should then edit and create the appropriate class in the `fasta_parser.py` file.
If the `fasta_parser.py` is made just right, then this current file should execute until the end without problems. Until the `fasta_parser.py` is not correct, this script will produce an error somewhere when you try to run it.

So to solve this exercise, it basically goes like this:

* Write some code in `fasta_parser.py`
* Run the current file `exercise_day5.py`, for instance by using the "execfile()" function.
* Observe and analyze the Exception produced.
* Rinse and repeat.

In summary, the current script will import the other file called `fasta_parser.py` and will test that the object you made provides a certain specific functionality and behavior.

The class to write will be called `FastaParser` and will parse "FASTA" type files. Two small example fasta files are provided. FASta files contain DNA sequences, one after the other. But not necessarily all on the same line.

There is one constraint: do not import any modules that are not built-in.
Modules such as `os` and `sys are okay but not the "Biopython" module for instance.

If you want to benefit from this exercise, restrain from searching for a fasta parser online and copying the code, obviously.

Once the exercise is completed, this file should be uploaded to your github repository called "python_homework" in a directory called "day5".
"""

# Step 1: Importing code from other files.
# Once the class exists in a file called "fasta_parser.py", it
# should be possible to import it and the following line should
# work. It just needs to be in the same directory as your current directory.
from fasta_parser import FastaParser

# Step 2: One should be able to make new instances from this class
# In fact, one should be able to make as many new objects as one wants
# Here we will just make two in this example.
# The class initialization should take one argument: the path of the fasta
# file to parse.
# Use the included fasta example files, or it won't work.
contigs = FastaParser("all_contigs.fasta")
genes = FastaParser("predicted_genes.fasta")

# Step 3: What if we give a path, but there is no file there ?
# Then your class should complain ! It must throw an exception
# of type IOError. To check this, we will use a function from
# the pytest module. It's like assert but for Exceptions.
# If you don't have pytest just install it with "pip install --user pytest"
import pytest
with pytest.raises(IOError): not_found = FastaParser('/file_does_not_exist.fasta')

# Step 4: What if we don't give a file path at all when making
# a new instance ? Then your class should complain !
# It must throw an exception of type TypeError.
with pytest.raises(TypeError): missing = FastaParser()

# Step 5: OK, we are ready to start reading data from the file.
# The first thing the class should do is report how many sequences
# there are in the fasta file
assert contigs.count == 5

# Step 6: We also want this to be reported as the true length property
# itself of the object.
assert len(genes) == 8

# Step 7: Now we need to access the sequences in the file. Let's say
# we want to use a simple indexing syntax that is natural for python.
# For instance asking for the second sequence would be as easy as genes[1]
assert genes[1] == "CGAGACTTATTCCTGAGATACTGTCCTTTCTCA"

# Step 8: What if we ask for a sequence number that is too high ?
# Then your class should complain !
# It must throw an exception of type IndexError.
with pytest.raises(IndexError): print contigs[93558]

# Step 9: What if the user doesn't ask for a certain sequence
# by using a number, but by using its ID ? This should work too.
# Check that we complain when we should complain by raising the
# right exception.
assert contigs['contig_C'] == "GACTACCAGGGTATCTAATCCCGTTTGCTCCCTTGGCTTTCGTGC"
with pytest.raises(KeyError): print contigs['lorem_ipsum']

# Step 10: We can do all sorts of things with a fasta file.
# Maybe the user wants to get all sequences that are shorter
# than a given length. Write a method called 'extract_length'
# that takes only one argument: the length above which we will
# filter sequences. Return the new sequences in a list.
assert len(genes.extract_length(30)) == 3
assert contigs.extract_length(2) == []

# Step 11: Write a method "length_dist" that takes one argument:
# the path at which a PDF file containing a graph of the length
# distribution of the sequences should be created. Use any plotting
# library such as matplotlib or ggplot.
genes.length_dist("~/test/genes_lengths.pdf")

# Step 12: Document your code. Each method should have a docstring
# that would enable someone to quickly understand your code.

# Step 13: Add style to your code. Google the term "PEP8" and make
# your code conform to it.
