import os
import pandas as pd
import joblib
from sklearn.linear_model import LogisticRegressionCV
from sklearn.model_selection import GridSearchCV
from xgboost import XGBClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupShuffleSplit
from sklearn.metrics import accuracy_score
import argparse
import json

parser = argparse.ArgumentParser(description="Takes feature-labeled .bed files and corresponding meta data as inputs to generate a sparse model of phenotype predictors.")

#command line arguments
parser.add_argument('--beds',required=True,help="Directory containing .bed files.")
parser.add_argument('-m', required=True,help="Meta data and groups table for genomic records.")
parser.add_argument('-o',default=os.path.abspath(os.getcwd()),help="Output directory")
parser.add_argument('-s',type=float,default=0.10,help="Fraction size for group splits. Default: 0.10.")
parser.add_argument('-n',type=int,default=100,help="Number of splits for groups splits. Default: 100")
parser.add_argument('--RVDB',action='store_true',default=False,help="Flag for RVDB fasta headers.")
parser.add_argument('-j',required=True,help="JSON file containing Grid Search parameters.")
parser.add_argument('-p',type=int,default=-1,help="Number of processors to use. Default: Max available")


cwd = myargs.o
out_prefix = os.path.join(myargs.o,metafile.split('.')[0]+os.path.basename(myargs.j).split('.')[0])
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
labels = complete_table['label']
groups = complete_table['groups']
features = complete_table.drop(['accession','label','groups','species'],axis=1).copy()
print("Assigning variables.")
X = features.values
y = labels
print(features)
gss = GroupShuffleSplit(n_splits=myargs.n,test_size=.1)

#parameters = {'C':[.00001,.0001,.001,.01,.1,1,10,100,1000,10000]} #Cs:10
with open(myargs.j,'r') as infile:
	parameters = json.load(infile)

tree = XGBClassifier(booster='gbtree',tree_method='gpu_hist',verbosity=2,objective="binary:logistic")
clf = GridSearchCV(tree,parameters,scoring='brier_score_loss',cv=gss,n_jobs=cpus,return_train_score=False)
print("Fitting model.")
clf.fit(X,y,groups=groups)
print("Training complete.")
print("Trained model score (accuracy):",accuracy_score(y,clf.predict(X)))

cv_results = pd.DataFrame(clf.cv_results_).sort_values('mean_test_score')
print("Max mean CV score: {}".format(cv_results['mean_test_score'].max()))
print("CV results:")
print(cv_results)

feat_import = pd.Series(dict(zip(features.columns[(clf.best_estimator_.feature_importances_ !=0)],clf.best_estimator_.feature_importances_[(clf.best_estimator_.feature_importances_ != 0)])))
print("Model feature importances:")
print(feat_import)

os.chdir(cwd)
#output total feature numpy array as serialized pickle file
features.columns.values.dump(out_prefix+"_feature_array.pickle")
#output sparse model coefficients
feat_import.to_csv(out_prefix+"_feat_importances.tsv",sep='\t',header=False)
cv_results.to_csv(out_prefix+"_cv_results.tsv",sep='\t')
#output serialized classifier object for persistence
joblib.dump(clf,out_prefix+"_CLF.joblib")
