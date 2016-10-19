#!/usr/bin/python2

# by Dmytro Kryvokhyzha dmytro.kryvokhyzha@ebc.uu.se

# Home work Day 5


class FastaParser(object):
    def __init__(self, filename):
        self.filename = filename
        self.name = []
        self.sequence = []
        self.count = 0

        # read fasta file
        fasFile = open(self.filename, 'r')
        seq = ''
        for line in fasFile:
            words = line.split()
            if words[0].startswith('>'):  # find sequence names
                self.count += 1
                self.name.append(words[0][1:])  # append sequence names
                if seq:
                    self.sequence.append(seq)
                seq = ''
            else:
                seq += words[0]  # append sequences
        self.sequence.append(seq)  # append the last sequence.

    def __len__(self):
        """Enables length property for number of sequences in a file """
        return self.count

    def __getitem__(self, i):
        """ Enables iteration through the sequences by index and names"""
        if isinstance(i, int):  # if index
            if i > self.count or i < 0:
                raise IndexError("Index is out of range")
            return self.sequence[i]
        else:   # if name
            if i not in self.name:
                raise KeyError("No gene with such name")
            seqIndex = self.name.index(i)
            return self.sequence[seqIndex]

    def extract_length(self, a):
        """ Function gets all sequences that are shorter than a given length"""
        if not(isinstance(a, int)):
            raise TypeError("Length should be integer")
        shortSeq = []
        for s in self.sequence:
            if len(s) < a:
                shortSeq.append(s)
        return shortSeq

    def length_dist(self, pdffile):
        """Plot the distribution of sequence length in fasta.
        Takes one argument: the path at which to output a PDF"""
        import matplotlib.pyplot as plt
        n = [len(i) for i in self.sequence]
        # Plot a histogram
        plt.hist(n, color="grey", bins=len(n))
        plt.title("Sequence length distribution", size=15)
        plt.xlabel("Sequence length")
        plt.ylabel("Frequency")
        plt.savefig(pdffile, dpi=90)
