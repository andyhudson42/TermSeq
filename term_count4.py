#!/usr/bin/env python

import sys
import os
from os.path import basename
from collections import defaultdict
import re

#Terminator input keys and sequences to search against input file
term_keys = open('C:\\Users\\Andy\\Desktop\\TermSeq\\Term_keys.txt').read().split()
term_seqs = open('C:\\Users\\Andy\\Desktop\\TermSeq\\Term_seqs.txt').read().split()

if len(term_keys) != len(term_seqs):
	print("Unequal number of keys and seqs\n", file=sys.stderr)
	exit(1)

#Add a "T" to each terminator key value
term_keys = ["T" + k for k in term_keys]

#Generate a dictionary of terminator keys and sequences
term_dict = dict(zip(term_keys, term_seqs))

files_list = []

for root, dirs, files in os.walk('./'):
	for filename in files:
		if filename.endswith(".fastq"):
			files_list.append(os.path.join(root, filename))


for file_name in files_list:
	file_title = file_name

	with open('.\\Term_Hits\\' + str(file_title) + '.out', "w") as summ_fh:
		summ_fh.write("Term ID\t" + "Terminator Sequence (5' to 3')\t" + "Total Hits\t" + "Terminated\t" + "Readthrough\t" + "%T termination" + "\n")

	#For each terminator in the library, make a new file, then search the input file and write each matching line in the file
	for key, term in term_dict.items():
		new_file_name = str(key + ".txt")
		#Check to see that the appropriate directory exists for dumping hits
		if not os.path.exists('.\\Term_Hits\\' + str(file_title) + "\\"):
			os.makedirs('.\\Term_Hits\\' + str(file_title) + "\\")

		with open('.\\Term_Hits\\' + str(file_title) + "\\" + new_file_name, "w") as fh:
			term_count = 0
			term_hits = []
			term_t = 0
			term_rt = 0
			with open(file_name, "r") as fastq_fh:
				for line in fastq_fh:
					if re.search(term, line):
						term_count += 1
						term_hits.append(line)
			for hit in term_hits:
				if len(hit) <= 75:
					term_t += 1
				elif len(hit) > 75:
					term_rt += 1

			if term_count == 0:
				term_eff = str("I")
			else:
				term_eff = str((term_t/term_count) * 100) + "%"
			fh.write(key + " terminator sequence: 5'-" + term + "-3'" + "\n" + "File name: " + file_name + "\n" + "Total Count: " + str(term_count) + "\n" + "Terminated: " + str(term_t) + "\n" + "Readthrough: " + str(term_rt) + "\n" + "Term. Eff. = " + str(term_eff) + "\n\n" + "Hits:\n\n")
			fh.write("".join(term_hits))
			with open('.\\Term_Hits\\' + str(file_title) + '.out', "a") as summ_fh:
				summ_fh.write(key + "\t" + "AAAAAAAA" + term + "TTTTTTTT" + "\t" + str(term_count) + "\t" + str(term_t) + "\t" + str(term_rt) + "\t" + str(term_eff) + "\n")