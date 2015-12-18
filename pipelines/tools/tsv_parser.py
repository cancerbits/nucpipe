#!/usr/bin/env python
"""
This is just a simple parser to grab, from a tsv file, a certain column.
I wrote it to pull output from biseq methcalling results, but it could
be used for anything... just give a tsv file and a key column to pull.
"""

__author__ = "Nathan C. Sheffield"

import sys, os, csv
import subprocess
from argparse import ArgumentParser


parser = ArgumentParser(description='TSV Parser')

parser.add_argument('-i', '--input-file', dest='input_file',
	help="Supply tsv file.",
	required=True)

parser.add_argument('-c', '--col', dest='column',
	default="uniqueSeqMotifCount", nargs="+",
	help="Choose column(s) to extract",
	required=True)

parser.add_argument('-k', '--keys', dest='include_keys',
 	default=False, action="store_true",
	help="Print keys?",
	required=False)



args = parser.parse_args()

input_open = open(args.input_file, 'rb')
f = csv.DictReader(input_open, delimiter="\t")

for row in f:  # iterates the rows of the file in order
	#print(row)
	for key in args.column:
		if args.include_keys:
			print("\t".join([key, row[key]]))
		else:
			print(row[key])
