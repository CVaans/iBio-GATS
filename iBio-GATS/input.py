#!/usr/bin/python3
from Bio.SeqIO
import os
import sys

infile = input("Insert full file name including the fasta extension: ")
# remove any white space
infile = infile.strip()

# check that the file exists
if not os.path.is_file(infile):
    print("file not found")
    sys.exit()

count = 0
total_len = 0
for seq_record in SeqIO.parse(open(infile), "fasta"):
    count += 1
    total_len += len(seq_record.seq)

print("%i records with total sequence length %i" % (count, total_len))