#!/usr/bin/env python

import sys
import os
# from Bio import SeqIO
from collections import defaultdict
import re


# get all input files from a directory or a file
# specify the input file or directory from the command line
def get_files_to_analyze(file_or_directory):
    """
    :param file_or_directory: name of either file or directory of files
    :return: a list of all files (even if there is only one)

    I am pasting the Python docs for some of the built-in functions so you can see what they do.
    We will automatically create a list of all files, so we don't need to specify them.

        sys.argv:

        The list of command line arguments passed to a Python script.
        argv[0] is the script name (it is operating system dependent whether this
        is a full pathname or not). If the command was executed using the -c command line
        option to the interpreter, argv[0] is set to the string '-c'. If no script name was
        passed to the Python interpreter, argv[0] is the empty string.
        ---

        os.walk(top, topdown=True, onerror=None, followlinks=False):

        Generate the file names in a directory tree by walking the tree either top-down or bottom-up.
        For each directory in the tree rooted at directory top (including top itself),
        it yields a 3-tuple (dirpath, dirnames, filenames).
        ---

    """

    # Python has some built in functions that we can use to see if a variable
    # refers to either a file or a directory. If the user specified a file
    # we will use that as the only entry in our list. If not, we will gather all the files from the
    # directory into the list.

    # This will be the list we send back
    files_list = []

    if os.path.isdir(file_or_directory):
        # Using a loop, we will go through every file in every directory
        # within the specified directory and create a list of the file names
        for root, dirs, files in os.walk(file_or_directory):
            for filename in files:
                files_list.append(os.path.join(root, filename))
    else:
        # We will just use the file given by the user
        files_list.append(os.path.abspath(file_or_directory))

    # We will sort the list when we send it back, so it always returns the
    # same order
    return sorted(files_list)


# def get_all_terminators_biopython(list_of_files):
#     """
#     Since you are looking at fastq files (which can be large) this function
#     iterates through them line by line, rather than concatenating them all at once
#     and storing them in memory. This should allow processing of an arbitrary number
#     of fastq files of any size.
#
#     This uses BioPython for reading through the sequencing files, but you could
#     do the same thing without BioPython (see the get_all_terminators() function)
#
#     :param list_of_files: List of all fastq files
#     :return: terminator_dictionary: A dictionary of all terminator sequences and counts
#     """
#
#     # This initializes our dictionary, and will set a default value
#     # of 0 if a terminator has never been encountered before.
#     # https://docs.python.org/3/library/collections.html?highlight=defaultdict#collections.defaultdict
#     terminator_dictionary = defaultdict(int)
#     # go through every file
#     for file in list_of_files:
#         # go through every sequence record in the file
#         # this creates a filehandle that can be iterated over
#         with open(file, "r") as fh:
#             # you need to specify what type of sequence file
#             # they are listed here: http://biopython.org/wiki/SeqIO
#             # now we can iterate over the fastq file record by record
#             for record in SeqIO.parse(fh, "fastq"):
#                 # now that we have the record http://biopython.org/wiki/SeqRecord
#                 # we can use your regular expression to see if it contains
#                 # a terminator. The sequence is found in record.seq
#                 m = re.search('A{8}.{32}T{4}', str(record.seq))
#
#                 # check if there was a match. If there was, add it to our
#                 # dictionary of terminators. If not, do nothing.
#                 if m:
#                     terminator_dictionary[m.group(0)] += 1
#
#     return terminator_dictionary


def get_all_terminators(list_of_files):
    """
    This is the straight-up Python version of the same thing, that does
    not require BioPython. It will look through all lines of the file for
    a match, not just the sequence containing ones. Otherwise it is identical.

    :param list_of_files: List of all fastq files
    :return: A dictionary of all terminator sequences and counts
    """

    # This initializes our dictionary, and will set a default value
    # of 0 if a terminator has never been encountered before.
    # https://docs.python.org/3/library/collections.html?highlight=defaultdict#collections.defaultdict
    terminator_dictionary = defaultdict(int)
    # go through every file
    for file in list_of_files:
        # go through every sequence record in the file
        # this creates a filehandle that can be iterated over
        with open(file, "r") as fh:
            for line in fh:
                m = re.search('A{8}.{32}T{4}', line)

                # check if there was a match. If there was, add it to our
                # dictionary of terminators. If not, do nothing.
                if m:
                    terminator_dictionary[m.group(0)] += 1

    return terminator_dictionary


def print_terminators(terminator_dict):
    """
    Prints out the terminators and their counts, sorted low to high
    :param terminator_dict: Dictionary of terminators and their counts
    :return: Success
    """

    # build a sorted dict using a list comprehension
    sorted_terminator_dict = [(k, terminator_dict[k]) for k in sorted(terminator_dict, key=terminator_dict.get)]
    for t,c in sorted_terminator_dict:
        # this just prints out the terminators and their counts separated by tab
        # https://docs.python.org/3/library/string.html#format-string-syntax
        print("{}\t{}".format(t, c))

        # if you want a new file per terminator, this will create them
        new_file_name = str(t + "_" + str(c) + ".txt")
        with open(new_file_name, "w") as fh:
            fh.write(t)