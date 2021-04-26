import os
import pandas as pd
import joblib
from sklearn.linear_model import LogisticRegressionCV
from sklearn.model_selection import GridSearchCV
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupShuffleSplit
from sklearn.metrics import accuracy_score
from sklearn.utils import resample
import argparse

parser = argparse.ArgumentParser(description="Takes feature-labeled .bed files and corresponding meta data as inputs to generate a sparse model of phenotype predictors.")

#command line arguments
parser.add_argument('--beds',required=True,help="Directory containing .bed files.")
parser.add_argument('-m', required=True,help="Meta data and groups table for genomic records.")
parser.add_argument('-o',default=os.getcwd(),help="Output directory")
parser.add_argument('-s',type=float,default=0.10,help="Fraction size for group splits. Default: 0.10.")
parser.add_argument('-n',type=int,default=100,help="Number of splits for groups splits. Default: 200")
parser.add_argument('--RVDB',action='store_true',default=False,help="Flag for RVDB fasta headers.")
parser.add_argument('-i',type=int,default=500,help="Number of iterations for coordinate descent.")
parser.add_argument('-p',type=int,default=os.cpu_count(),help="Number of processors to use. Default: Max available")
parser.add_argument('-t',type=float,default=.00000001,help="Min loss tolerance for stopping. Default: .00000001")
parser.add_argument('-r',type=int,default=0,help="Number of resampled rows using stratified groups. (Default is no resample)")
myargs=parser.parse_args()

cwd = os.path.abspath(myargs.o)
metafile = myargs.m
split_size = myargs.s
iterations = myargs.i
tolerance = myargs.t
cpus = myargs.p

meta = pd.read_table(metafile)
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
    feature_dict[accession_num]=b['name'].value_counts().to_dict()

feature_table = pd.DataFrame(feature_dict).T.fillna(0.0)
feature_table.index.name = "accession"
feature_table.reset_index(inplace=True)

accession_set = set(feature_table['accession'])

complete_table = pd.merge(feature_table,meta,left_on='accession',right_on='accession')
print("Dropping ambiguous labels.")
complete_table = complete_table[complete_table['label'] > -1]
complete_table = complete_table[complete_table['accession'].isin(accession_set)]
if myargs.r:
    print("Resampling training data {} times.".format(myargs.r))
    complete_table = resample(complete_table,n_samples=myargs.r,stratify=complete_table['groups'])

labels = complete_table['label']
groups = complete_table['groups']
features = complete_table.drop(['accession','label','groups','species'],axis=1).copy()
print("Assigning variables.")
X = features.values
y = labels
print(features)
gss = GroupShuffleSplit(n_splits=myargs.n,test_size=.1)

#parameters = {'C':[.00001,.0001,.001,.01,.1,1,10,100,1000,10000]} #Cs:10
parameters = {'C':[.01,.1,1,10,100,1000,10000]}
logit = LogisticRegression(penalty='l1',verbose=1,solver='liblinear',max_iter=iterations,tol=tolerance,fit_intercept=False)
clf = GridSearchCV(logit,parameters,scoring='brier_score_loss',cv=gss,n_jobs=cpus,return_train_score=False)
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
features.columns.values.dump(metafile.split('.')[0]+"_feature_array.pickle")
#output sparse model coefficients
model_coef.to_csv(metafile.split('.')[0]+"_model_coefficients.tsv",sep='\t')
cv_results.to_csv(metafile.split('.')[0]+"_cv_results.tsv",sep='\t')
#output serialized classifier object for persistence
joblib.dump(clf,metafile.split('.')[0]+"_CLF.joblib")
