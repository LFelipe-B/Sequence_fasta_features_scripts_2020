#!/usr/bin/python

# This script was writend by L. Felipe Benites in September 2020

# usage: python get_sequence_length.py yoursequence.fasta

import sys
from Bio import SeqIO

fasta_file = sys.argv[1]
csv_file = sys.argv[1]+ '_length.csv'


length = {}

with open(fasta_file, 'r') as fh:
    for seq_record in SeqIO.parse(fh, "fasta"):
        length[seq_record.id] = len(seq_record.seq)

# Create output file
with open(csv_file, 'w') as out_file:
    out_file.write("Sequence_id,Length\n")
    for seq in length:
        out_file.write(','.join([seq, str(length[seq])]) + '\n')
