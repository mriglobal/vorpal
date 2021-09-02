import os
import pandas as pd
import numpy as np
import re
from Bio import SeqIO
from skbio import Protein
from collections import Counter
import json
import argparse

parser = argparse.ArgumentParser(description="Takes CDS file of coding regions, a json for an amino alphabet, and a model coefficients table and maps the motifs back to the records. Outputs a table.")

#command line arguments
parser.add_argument('--seqs',required=True,help="File containing CDS data.")
parser.add_argument('-o',default='',help="Prefix for output files.")
parser.add_argument('-j',default=None,help="Amino acid translation dictionary in json format. Default: No re-encoding")
parser.add_argument('-m',required=True,help="Predictor motif coefficients table.")

myargs=parser.parse_args()

if not myargs.o:
    out_prefix = os.path.basename(myargs.m.split('.')[0])
else:
    out_prefix = myargs.o


if myargs.j:
    with open(myargs.j,'r') as infile:
         amino_alpha = json.load(infile)
         amino_alpha['X'] = 'X'

reduced = set(amino_alpha.values())

reverse_dict = {r:[] for r in reduced}

for r in reduced:
    for a in amino_alpha:
        if amino_alpha[a] == r:
            reverse_dict[r].append(a)

alpha_lookup = {}

for r in reverse_dict.keys():
    if len(reverse_dict[r]) > 1:
        alpha_lookup[r] = '[' + ''.join([a for a in reverse_dict[r]])+']'
    else:
        alpha_lookup[r] = reverse_dict[r][0]

def make_pattern(reduced_string):
    return re.compile(''.join([alpha_lookup[c] for c in reduced_string]))
seqs= list(SeqIO.parse(myargs.seqs,'fasta'))

coefs = pd.read_table(myargs.m,names=['motif','coef'])

#coefs.columns = ['motif','coef']

print("Finding motifs.")
data = {'amino_match':[],'start':[],'end':[],'motif':[],'name':[],'accession':[]}
for m in coefs['motif']:
    for s in seqs:
        results = re.finditer(make_pattern(m),str(s.seq.translate()).replace('*',''))
        if results:
            for r in results:
                data['amino_match'].append(r[0])
                data['start'].append(r.span()[0])
                data['end'].append(r.span()[1])
                data['motif'].append(m)
                data['name'].append(s.description)
                data['accession'].append(s.id.split(':')[0])

output = pd.DataFrame(data)
output['accession'] = output['accession'].str.replace('join\(','')

output = pd.merge(output,coefs)

output.to_csv(out_prefix+'_feature_table.tsv',sep='\t')
