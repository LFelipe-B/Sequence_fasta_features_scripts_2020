#!/usr/bin/env python

"""
This script was writen in September 2020 to calculate sequence GC and lenght stats using multifasta file. Usage python Get_GC_length.py sequence.fasta
"""

# Usage: python get_GC_length.py <sequences.fasta> <sequences_stats.csv>

import sys
from Bio import SeqIO


fasta_file = sys.argv[1]
csv_file = sys.argv[1]+ '_GC_length.csv'

distribution = {}

# Iterate over sequences: get GC content and length
with open(fasta_file, 'r') as fh:
    for seq_record in SeqIO.parse(fh, "fasta"):
        a=0; t=0; c=0; g=0; n=0; r=0; y=0; x=0
        for i in str(seq_record.seq):
            if "A" in i:
                a += 1
            elif "T" in i:
                t += 1
            elif "C" in i:
                c += 1
            elif "G" in i:
                g += 1
            elif "N" in i:
                n += 1
            elif "R" in i:
                r += 1
            elif "Y" in i:
                y += 1
            elif "X" in i:
                x += 1
            gc_content = (g + c + 0.) * 100 / (a + t + g + c + n + r + y+ x + 0.)
            distribution[seq_record.id] = gc_content,len(seq_record.seq)

# Create output file
with open(csv_file, 'w') as out_file:
    out_file.write("Sequence_id,GC%,Length\n")
    for seq in distribution:
        out_file.write(','.join([seq, str(distribution[seq][0]), str(distribution[seq][1])]) + '\n')
