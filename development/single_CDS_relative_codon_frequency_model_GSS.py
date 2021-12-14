import os
import pandas as pd
import numpy as np
from Bio import SeqIO
from collections import Counter
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GroupShuffleSplit
from sklearn.metrics import accuracy_score
from sklearn.utils import resample
import argparse

parser = argparse.ArgumentParser(description="Takes a file containing a collections of CDS sequences for single gene and fits a LASSO relative codon frequency model with a reduced alphabet amino dictionary to a binary phenotype.")

#command line arguments
parser.add_argument('--seqs',required=True,help="File containing CDS data.")
parser.add_argument('-m', required=True,help="Meta data and groups table for genomic records.")
parser.add_argument('-c',default=[.01,.1,1,10,100,1000,10000],nargs='+',type=float,help="List of Cs to search over. Default: [.01,.1,1,10,100,1000,10000]")
parser.add_argument('-o',default='',help="Prefix for output files.")
parser.add_argument('-s',type=float,default=0.10,help="Fraction size for group splits. Default: 0.10.")
parser.add_argument('-n',type=int,default=100,help="Number of splits for groups splits. Default: 100")
parser.add_argument('-i',type=int,default=500,help="Number of iterations for coordinate descent.")
parser.add_argument('-p',type=int,default=os.cpu_count(),help="Number of processors to use. Default: Max available")
parser.add_argument('-t',type=float,default=.00000001,help="Min loss tolerance for stopping. Default: .00000001")
parser.add_argument('-r',type=int,default=0,help="Number of resampled rows using stratified groups. (Default is no resample)")
myargs=parser.parse_args()

codons = pd.read_table(os.path.join(os.path.dirname(os.__file__),"codon_table_human.txt"),index_col='codon')

cwd = os.getcwd()
metafile = myargs.m
if not myargs.o:
    out_prefix = os.path.basename(metafile.split('.')[0])
else:
    out_prefix = myargs.o
split_size = myargs.s
num_splits = myargs.n
iterations = myargs.i
tolerance = myargs.t
cpus = myargs.p

meta = pd.read_table(metafile)

seqs = list(SeqIO.parse(myargs.seqs,'fasta'))
accession_list = [s.id.split(":")[0].replace("join(",'') for s in seqs]

accession_dict = {a:[] for a in accession_list}

feature_dict = {}

for a in accession_dict.keys():
    counter = Counter()
    for f in accession_dict[a]:
        counter.update([f[r:r+3] for r in range(0,len(f),3) if f[r:r+3] in codons.index])
    feature_dict[a] = pd.Series(counter)

data = pd.DataFrame(feature_dict)

data = data.T[~data.isna().all().values].fillna(0.0)

def get_relative(x):
    output_dict = {}
    for c in x[1].keys():
        output_dict[c] = x[1][c] / x[1][codons[codons['amino'] == codons.loc[c]['amino']].index].sum()
    return output_dict

print("Getting relative codon frequencies.")
relative_counts = [get_relative(row) for row in data.iterrows()]

data = pd.DataFrame(relative_counts,index=data.index)

data.index.name= 'accession'
model_data = pd.merge(data.reset_index(),meta).copy()

if myargs.r:
    print("Resampling training data {} times.".format(myargs.r))
    model_data = resample(model_data,n_samples=myargs.r,stratify=model_data['groups'])

gss = GroupShuffleSplit(n_splits=num_splits,test_size=split_size)

features = model_data[data.columns.values].copy()
X = features.values
y = model_data['label']
groups = model_data['groups']

parameters = {'C':myargs.c}

logit = LogisticRegression(penalty='l1',verbose=1,solver='liblinear',max_iter=iterations,tol=tolerance,fit_intercept=False)
clf = GridSearchCV(logit,parameters,scoring='brier_score_loss',cv=gss,n_jobs=-1,return_train_score=False)
print("Fitting model.")
clf.fit(X,y,groups=groups)

print("Training complete.")
print("Trained model score (accuracy):",accuracy_score(y,clf.predict(X)))

cv_results = pd.DataFrame(clf.cv_results_).sort_values('mean_test_score')
print("Max mean CV score: {}".format(cv_results['mean_test_score'].max()))
print("CV results:")
print(cv_results)

model_coef = pd.Series(dict(zip(features.columns[(clf.best_estimator_.coef_ !=0)[0]],clf.best_estimator_.coef_[(clf.best_estimator_.coef_ != 0)])))

print("Model predictors with non-zero coefficients:")
print(model_coef)

os.chdir(cwd)
#output total feature numpy array as serialized pickle file
features.columns.values.dump(out_prefix+"_feature_array.pickle")
#output sparse model coefficients
model_coef.to_csv(out_prefix+"_model_coefficients.tsv",sep='\t')
cv_results.to_csv(out_prefix+"_cv_results.tsv",sep='\t')
#output serialized classifier object for persistence
joblib.dump(clf,out_prefix+"_CLF.joblib")

