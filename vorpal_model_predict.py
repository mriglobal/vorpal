import os
import pandas as pd
import numpy as np
from Bio import SeqIO
from skbio import Protein
import joblib
from sklearn.preprocessing import Binarizer
from sklearn.metrics import accuracy_score
from collections import Counter
import argparse

parser = argparse.ArgumentParser(description="Takes beds files containing vorpal motifs and a persistent model object and predicts on this data.")

#command line arguments
parser.add_argument('--beds',required=True,help="Directory containing .bed files to predict on.")
parser.add_argument('-m', help="Meta data and groups table for genomic records. (Required for accuracy output. Optional.)")
parser.add_argument('--model',required=True,help='Joblib file for sklearn classifier object.')
parser.add_argument('-f',required=True, help='Feature Array object for accompanying model.')
parser.add_argument('-o',default='',help="Prefix for output files.")
parser.add_argument('--RVDB',action='store_true',default=False,help="Flag for RVDB fasta headers.")
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
os.chdir(myargs.beds)

bed_files = [file for file in os.listdir() if '.bed' in file]

beds = [pd.read_table(f,names=['chr','start','end','name','score']) for f in bed_files]

if myargs.RVDB:
    for i in range(len(beds)):
        beds[i]['chr'] = [a[2] for a in beds[i]['chr'].str.split('|').values]


feature_dict = {}

print("Building feature dict.")
for b in beds:
    accession_num = b['chr'].unique()[0]
    feature_dict[accession_num]={f:0 for f in feature_vector}
    feature_dict[accession_num].update(b['name'].value_counts().to_dict())

feature_table = pd.DataFrame(feature_dict).T.fillna(0.0)
feature_table.index.name = "accession"
feature_table.reset_index(inplace=True)
test_table = feature_table
data.index.name = "accession"

data.reset_index(inplace=True)

if myargs.m:
    complete_table = pd.merge(data,meta[['accession','label','groups']])
    complete_table = complete_table[complete_table['label'] > -1]
    complete_table = complete_table[complete_table['accession'].isin(meta['accession'])]
    X = complete_table[feature_vector].copy()
    y = complete_table['label']
    print("Model accuracy:")
    print(accuracy_score(y,clf.predict(X)))

else:
    X = data[feature_vector].values
    if myargs.b:
        print("Binarizing features.")
        transformer = Binarizer()
        X = transformer.fit_transform(X)

data['predict'] = clf.predict(X)
data['predict_proba'] = clf.predict_proba(X)[:,1]

data[['accession','predict','predict_proba']].to_csv(out_prefix+"_predictions.tsv",sep='\t',index=False)
