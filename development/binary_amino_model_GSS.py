import os
import pandas as pd
import numpy as np
from Bio import SeqIO
from skbio import Protein
import joblib
from sklearn.model_selection import GridSearchCV
from sklearn.preprocessing import Binarizer
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupShuffleSplit
from sklearn.metrics import accuracy_score
from sklearn.utils import resample
from collections import Counter
import json
import argparse

parser = argparse.ArgumentParser(description="Takes CDS file of coding regions and fits a LASSO bag-of-words model with a reduced alphabet amino dictionary to a binary phenotype.")

#command line arguments
parser.add_argument('--seqs',required=True,help="File containing CDS data.")
parser.add_argument('-m', required=True,help="Meta data and groups table for genomic records.")
parser.add_argument('-o',default='',help="Prefix for output files.")
parser.add_argument('-k',default=6,type=int,help="Amino word K size. Default:6")
parser.add_argument('-j',default=None,help="Amino acid translation dictionary in json format. Default: No re-encoding")
parser.add_argument('-q',default=None,help="Quantile cutoff for reducing K-mer space. Default: None")
parser.add_argument('-b',action='store_true',default=False,help='Flag for feature vector binarization')
parser.add_argument('-s',type=float,default=0.10,help="Fraction size for group splits. Default: 0.10.")
parser.add_argument('-n',type=int,default=100,help="Number of splits for groups splits. Default: 100")
parser.add_argument('-i',type=int,default=500,help="Number of iterations for coordinate descent.")
parser.add_argument('-p',type=int,default=os.cpu_count(),help="Number of processors to use. Default: Max available")
parser.add_argument('-t',type=float,default=.00000001,help="Min loss tolerance for stopping. Default: .00000001")
parser.add_argument('-r',type=int,default=0,help="Number of resampled rows using stratified groups. (Default is no resample)")
myargs=parser.parse_args()

cwd = os.getcwd()
metafile = myargs.m
if not myargs.o:
    out_prefix = os.path.basename(metafile.split('.')[0])
else:
    out_prefix = myargs.o
split_size = myargs.s
iterations = myargs.i
tolerance = myargs.t
cpus = myargs.p

meta = pd.read_table(metafile)


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

feature_counter = {a:Counter() for a in accession_list}

print("Counting {}mers.".format(myargs.k))
for a in accession_dict.keys():
    for f in accession_dict[a]:
        if 'J' not in f and 'B' not in f and 'Z' not in f:
            feature_counter[a].update(Protein(amino_encode(f)).kmer_frequencies(myargs.k))

data = pd.DataFrame(feature_counter)

data = data.fillna(0.0).T
if myargs.q:
    print("Filtering out K-mers below {} quantile cutoff.".format(myargs.q))
    counts = data.sum()
    data = data[counts[counts > counts.quantile(float(myargs.q))].index]

data.index.name = "accession"
data.reset_index(inplace=True)

complete_table = pd.merge(data,meta[['accession','label','groups']],left_on='accession',right_on='accession')

complete_table = complete_table[complete_table['label'] > -1]
complete_table = complete_table[complete_table['accession'].isin(meta['accession'])]

if myargs.r:
    complete_table = resample(complete_table,n_samples=myargs.r,stratify=complete_table['groups'])

features = complete_table.drop(['accession','label','groups'],axis=1).copy()
if not myargs.b:
    X = features.values
else:
    print("Binarizing features.")
    transformer = Binarizer().fit(features.values)
    X = transformer.transform(features.values)

y = complete_table['label']
groups = complete_table['groups']

gss = GroupShuffleSplit(n_splits=100,test_size=myargs.s)
parameters = {'C':[.01,.1,1,10,100,1000,10000]}
logit = LogisticRegression(penalty='l1',verbose=1,solver='liblinear',max_iter=500,tol=.00000001,fit_intercept=False)
clf = GridSearchCV(logit,parameters,scoring='brier_score_loss',n_jobs=cpus,cv=gss,return_train_score=False)
print("Fitting model.")
clf.fit(X,y,groups=groups)

print("Training complete.")
print("Trained model score (accuracy):",accuracy_score(y,clf.predict(X)))

cv_results = pd.DataFrame(clf.cv_results_).sort_values('mean_test_score')
print("Max mean CV score: {}".format(cv_results['mean_test_score'].max()))
print("CV results:")
print(cv_results)

model_coef = pd.Series(dict(zip(features.columns[(clf.best_estimator_.coef_ !=0)[0]],clf.best_estimator_.coef_[(clf.best_estimator_.coef_ != 0)])))

print("Model predictors with coefficients greater than 0:")
print(model_coef)

os.chdir(cwd)
#output total feature numpy array as serialized pickle file
features.columns.values.dump(out_prefix+"_feature_array.pickle")
#output sparse model coefficients
model_coef.to_csv(out_prefix+"_model_coefficients.tsv",sep='\t')
cv_results.to_csv(out_prefix+"_cv_results.tsv",sep='\t')
#output serialized classifier object for persistence
joblib.dump(clf,out_prefix+"_CLF.joblib")
