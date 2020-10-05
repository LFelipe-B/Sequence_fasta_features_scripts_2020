#! /usr/local/bin/python3
# author: max-emil.schon@icm.uu.se
# co-author: l.felipebenites@gmail.com

# date: 16.02.2018
# SINGEK workshop Banyuls, France
# Notes: modified by L. Felipe Benites 21/2/2018 to account for older python versions (import division) and removed seaborn and matplot lib function; plus other stuffs (data frame names)
# Max wrote down in dt.frame the seq names "(...) I think it is best to store the information of distribution in a dictionary then and use the sequence names as the index".
#If need count the 'X' and 'Ns' use this dict (gc = {0:{'A':0,'G':0,'C':0,'T':0,'N':0,'X':0}, 1:{'A':0,'G':0,'C':0,'T':0,'N':0,'X':0}, 2:{'A':0,'G':0,'C':0,'T':0,'N':0,'X':0}}"
#I add the counts for "R", "N", and "Y".

# This script calculates the GC content of each codon from a nucl (CDS) file and output a .csv table with each sequence GC123 counts
# Cleaned V4 in September 2020 and changed name to get_sequence_GC3.py

#Usage python get_sequence_GC3.py <file.fasta>

from __future__ import print_function
from __future__ import division
from Bio import SeqIO
import pandas as pd
import sys

# seq by seq:
def gc_percent(counts):
    num_nucl = sum(counts.values())
    num_gc = sum([val for key,val in counts.items() if key == 'G' or key == 'C'])
    return num_gc/num_nucl*100

def gc_content(sequence):
    gc = {0:{'A':0,'G':0,'C':0,'T':0,'N':0,'X':0,'Y':0,'R':0}, 1:{'A':0,'G':0,'C':0,'T':0,'N':0,'X':0,'Y':0,'R':0}, 2:{'A':0,'G':0,'C':0,'T':0,'N':0,'X':0,'Y':0,'R':0}}
    for i, nucl in enumerate(sequence):
        try:
            gc[i % 3][nucl] += 1
        except KeyError:
            print("unknown nucleotide {}".format(nucl), end="\r")
    return {key:gc_percent(val) for key,val in gc.items()}

distribution = {}
for rec in SeqIO.parse(sys.argv[1], 'fasta'):
    distribution[rec.id] = gc_content(rec.seq)

df = pd.DataFrame.from_dict(distribution, 'index')
df.index.name='Sequence_id'
df.columns = [
              'GC1%',
              'GC2%',
              'GC3%'
              ]
df

# If print as a single file: df.to_csv("GC123_distribution.csv")
# To print with the input name "file.fasta"
df.to_csv(sys.argv[1] + '_GC123_distrib.csv')
