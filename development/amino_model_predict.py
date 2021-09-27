import os
import pandas as pd
import numpy as np
from Bio import SeqIO
from skbio import Protein
import joblib
from sklearn.preprocessing import Binarizer
from sklearn.metrics import accuracy_score
from collections import Counter
import json
import argparse

parser = argparse.ArgumentParser(description="Takes CDS file of coding regions, a json for an amino alphabet, and a persistent model object and predicts on this data.")

#command line arguments
parser.add_argument('--seqs',required=True,help="File containing CDS data.")
parser.add_argument('-m', help="Meta data and groups table for genomic records.")
parser.add_argument('--model',required=True,help='Joblib file for sklearn classifier object.')
parser.add_argument('-f',required=True, help='Feature Array object for accompanying model.')
parser.add_argument('-b',action='store_true',default=False,help='Flag for feature vector binarization.')
parser.add_argument('-o',default='',help="Prefix for output files.")
parser.add_argument('-k',default=6,type=int,help="Amino word K size. Default:6")
parser.add_argument('-j',default=None,help="Amino acid translation dictionary in json format. Default: No re-encoding")
parser.add_argument('-p',type=int,default=os.cpu_count(),help="Number of processors to use. Default: Max available")
myargs=parser.parse_args()

cwd = os.getcwd()
if myargs.m:
  metafile = myargs.m
  meta = pd.read_table(metafile)
else:
  metafile = myargs.j
if not myargs.o:
    out_prefix = os.path.basename(metafile.split('.')[0])
else:
    out_prefix = myargs.o

clf = joblib.load(myargs.model)
feature_vector = pd.read_pickle(myargs.f)
cpus = myargs.p

seqs = list(SeqIO.parse(myargs.seqs,'fasta'))
accession_list = [s.id.split(":")[0].replace("join(",'') for s in seqs]

accession_dict = {a:[] for a in accession_list}

for s in seqs:
    if len(s.seq)%3 == 0 and len(s.seq)//3 > myargs.k:
        accession_dict[s.id.split(":")[0].replace("join(",'')].append(str(s.seq.translate()).replace('*',''))

if myargs.j:
    with open(myargs.j,'r') as infile:
         amino_alpha = json.load(infile)
         amino_alpha['X'] = 'X'
else:
    amino_alpha = {'F':'F',
              'Y':'Y',
              'W':'W',
              'M':'M',
              'L':'L',
              'I':'I',
              'V':'V',
              'A':'A',
              'T':'T',
              'S':'S',
              'N':'N',
              'H':'H',
              'Q':'Q',
              'E':'E',
              'D':'D',
              'R':'R',
              'K':'K',
              'C':'C',
              'P':'P',
              'G':'G',
              'X':'X'}

def amino_encode(x):
    return ''.join([amino_alpha[i] for i in x])

feature_counter = {a:Counter({f:0 for f in feature_vector}) for a in accession_list}

print("Counting {}mers.".format(myargs.k))
for a in accession_dict.keys():
    for f in accession_dict[a]:
        if 'J' not in f and 'B' not in f and 'Z' not in f:
            feature_counter[a].update(Protein(amino_encode(f)).kmer_frequencies(myargs.k))

data = pd.DataFrame(feature_counter)

data = data.fillna(0.0).T
print(data)
data = data[feature_vector]

data.index.name = "accession"

data.reset_index(inplace=True)

if myargs.m:
    complete_table = pd.merge(data,meta[['accession','label','groups']])
    complete_table = complete_table[complete_table['label'] > -1]
    complete_table = complete_table[complete_table['accession'].isin(meta['accession'])]
    X = complete_table.drop(['accession','label','groups'],axis=1).copy()
    if myargs.b:
        print("Binarizing features.")
        transformer = Binarizer()
        X = transformer.fit_transform(X)
    y = complete_table['label']
    print("Model accuracy:")
    print(accuracy_score(y,clf.predict(X)))

else:
    X = data.drop(['accession'],axis=1).values
    if myargs.b:
        print("Binarizing features.")
        transformer = Binarizer()
        X = transformer.fit_transform(X)

data['predict'] = clf.predict(X)
data['predict_proba'] = clf.predict_proba(X)[:,1]

data[['accession','predict','predict_proba']].to_csv(out_prefix+"_predictions.tsv",sep='\t',index=False)
