import os
from Bio import SeqIO
from skbio import DNA
from collections import Counter
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-f',required=True,help="Degenerate motif fasta file.")

myargs=parser.parse_args()

file = myargs.f
motifs = list(SeqIO.parse(file,'fasta'))

degenerates = [Counter(DNA(str(s.seq)).degenerates())[True] for s in motifs]

average_degen_positions = np.mean(degenerates)

mycounts = Counter()

for m in motifs:
    mycounts.update(Counter(str(m.seq)))
    
degen_total = 0
for c in mycounts.keys():
    if c == 'R':
        degen_total+=2*mycounts[c]
    elif c == 'Y':
        degen_total+=2*mycounts[c]
    elif c == 'W':
        degen_total+=2*mycounts[c]
    elif c == 'S':
        degen_total+=2*mycounts[c]
    elif c == 'M':
        degen_total+=2*mycounts[c]
    elif c == 'K':
        degen_total+=2*mycounts[c]
    elif c == 'B':
        degen_total+=3*mycounts[c]
    elif c == 'D':
        degen_total+=3*mycounts[c]
    elif c == 'H':
        degen_total+=3*mycounts[c]
    elif c == 'V':
        degen_total+=3*mycounts[c]
    elif c == 'N':
        degen_total+=4*mycounts[c]

average_degens = degen_total/len(motifs)

with open(file+"_degenerate_summary.txt",'w') as outfile:
    outfile.write("Average number of degenerate positions:"+'\t'+str(average_degen_positions))
    outfile.write('\n')
    outfile.write("Average number of bases represented in degenerate motifs:"+'\t'+str(average_degens))
    outfile.write('\n')
