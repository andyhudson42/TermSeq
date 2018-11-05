#!/usr/bin/env python

import sys
import os
from os.path import basename
import re

#Terminator input keys and sequences to search against input file
term_keys = open('C:\\Users\\Andy\\Desktop\\TermSeq\\Term_keys.txt').read().split()
term_seqs = open('C:\\Users\\Andy\\Desktop\\TermSeq\\Term_seqs.txt').read().split()

if len(term_keys) != len(term_seqs):
	print("Unequal number of keys and seqs\n", file=sys.stderr)
	exit(1)

#Generate a dictionary of terminator keys and sequences
term_dict = dict(zip(term_keys, term_seqs))



#For each terminator in the library, make a new fasta file, then search all .txt files for a match and write each match to the file in fasta format
for key, term in term_dict.items():

	files_list = []

	for root, dirs, files in os.walk('./Term_Hits/'):
		for filename in files:
			if filename.endswith("T" + str(key) + ".txt"):
				files_list.append(os.path.join(root, filename))

		new_file_name = str("T" + key + "_matches.fasta")
		
		for file in files_list:
			seq_counter = 1
			with open('.\\MSA_Cons\\Matches\\' + new_file_name, "a") as fh:
				term_hits = []

				# for each terminator, search through each line of each fastq file for the first BamHI site, then an intervening sequence (poly-A sequence), terminator hairpin, then the next 8 nt (if available)
				with open(file, "r") as txt_fh:
					for line in txt_fh:
						m = re.search("GGATCC" + ".{8}" + term + ".{0,10}", line)
						if m:
							term_hits.append(m.group())

				for hit in term_hits:
					fasta_header = ">T" + str(key) + "_Seq" + str(seq_counter)
					fh.write(str(fasta_header) + "\n" + str(hit) +"\n")
					seq_counter += 1