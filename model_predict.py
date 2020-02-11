import os
import pandas as pd
import joblib
import argparse

parser = argparse.ArgumentParser(description="Takes features mapped to new instances and predicts on phenotype using provided classifier object and feature dictionary.")

#command line arguments
parser.add_argument('--beds',required=True,help="Directory containing .bed files.")
parser.add_argument('-m',default='',help="Optional meta data table for data labels (used when scoring a test set.)")
parser.add_argument('-o',default=os.getcwd(),help="Output directory")
parser.add_argument('-k',required=True,help="Pickled K-mer motif feature array.")
parser.add_argument('-c',required=True,help="Joblib serialized classifier object.")
parser.add_argument('--RVDB',action='store_true',default=False,help="Flag for RVDB fasta headers.")
myargs=parser.parse_args()

cwd = os.path.abspath(myargs.o)
metafile = myargs.m
clf = joblib.load(myargs.c)
features = pd.read_pickle(myargs.k)
# with open('Coronavirus_complete_features.txt','r') as infile:
#     features = [r.strip() for r in infile.readlines()]
if metafile:
    meta = pd.read_table(metafile)
os.chdir(cwd+os.path.join('/',myargs.beds))

bed_files = [file for file in os.listdir() if '.bed' in file]

beds = [pd.read_table(f,names=['chr','start','end','name','score']) for f in bed_files]

if myargs.RVDB:
    for i in range(len(beds)):
        beds[i]['chr'] = [a[2] for a in beds[i]['chr'].str.split('|').values]

kmer_index = {f:0 for f in features}

feature_dict = {}

print("Building feature dict.")
for b in beds:
    accession_num = b['chr'].unique()[0]
    feature_dict[accession_num]=b['name'].value_counts().to_dict()

for a in feature_dict.keys():
    for k in kmer_index.keys():
        if k not in feature_dict[a].keys():
            feature_dict[a].update({k:0})

feature_table = pd.DataFrame(feature_dict).T.fillna(0.0)
feature_table.index.name = "accession"

X = feature_table.values

clf.predict(X)