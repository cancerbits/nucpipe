#!/usr/bin/python
"""
This is just a simple parser to grab, from a tsv file, a certain column.
I wrote it to pull output from biseq methcalling results, but it could
be used for anything... just give a tsv file and a key column to pull.
AUTHOR: Nathan Sheffield
"""

import sys, os, csv
import subprocess
from argparse import ArgumentParser


parser = ArgumentParser(description='TSV Parser')

parser.add_argument('-i', '--input-file', dest='input_file', help="Supply tsv file.")
parser.add_argument('-c', '--col', dest='column', default="uniqueSeqMotifCount",
					help="Choose column to extract")
args = parser.parse_args()

input_open = open(args.input_file, 'rb')
f = csv.DictReader(input_open, delimiter="\t")

try:
		for row in f:  # iterates the rows of the file in orders
			#print(row)
			print(row[args.column])
finally:
	pass