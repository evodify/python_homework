#!/usr/bin/python2

# by Dmytro Kryvokhyzha dmytro.kryvokhyzha@ebc.uu.se

# Home work Day 5

class FastaParser(object):
	def __init__(self, filename):
		self.filename = filename
		fasFile = open(self.filename, 'r')
		
		# read file in
		seqDic = {}
		for line in fasFile:
			words = line.split()
			# print type(words[0].startswith('>'))
			if ">" in words[0]:
				seqkey = words
				seq = []
			else:
				seq.append(words[0])
			# seqDic[seqkey] = seq
		print seq


contigs = FastaParser("all_contigs.fasta")
# print contigs
# genes = FastaParser("predicted_genes.fasta")
# print genes