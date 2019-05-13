import os
import pandas as pd
import joblib
from sklearn.model_selection import GridSearchCV
from sklearn.linear_model import SGDClassifier
import argparse

parser = argparse.ArgumentParser(description="Takes feature-labeled .bed files and corresponding meta data as inputs to generate a sparse model of phenotype predictors.")

#command line arguments
parser.add_argument('--beds',required=True,help="Directory containing .bed files.")
parser.add_argument('-m', required=True,help="Meta data table for genomic records.")
parser.add_argument('-o',default=os.getcwd(),help="Output directory")
parser.add_argument('-f',type=int,default=5,help="Number of folds for cross validation.")
myargs=parser.parse_args()

cwd = myargs.o
metafile = myargs.m
folds = myargs.f
# with open('Coronavirus_complete_features.txt','r') as infile:
#     features = [r.strip() for r in infile.readlines()]

meta = pd.read_table(metafile)
os.chdir(os.getcwd()+os.path.join(myargs.beds))

bed_files = [file for file in os.listdir() if '.bed' in file]

beds = [pd.read_table(f,names=['chr','start','end','name','score']) for f in bed_files]

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

complete_table = pd.merge(feature_table,meta)
print("Dropping ambiguous labels.")
complete_table = complete_table[complete_table['label'] > -1]
complete_table = complete_table[complete_table['accession'].isin(accession_set)]
labels = complete_table['label']
features = complete_table.drop(['accession','label'],axis=1).copy()
print("Assigning variables.")
X = features.values
y = labels

params = {'alpha':(.1,.01,.001,.0001,.00001,.000001)}
log_reg = SGDClassifier(n_jobs=-1,loss='log',verbose=3,penalty='elasticnet',max_iter=1000,tol=.00000001)
cv_clf = GridSearchCV(log_reg,params,cv=folds)
print("Fitting model.")
cv_clf.fit(X,y)
print("Training complete.")
print("Trained model score (accuracy):",cv_clf.score(X,y))
clf=cv_clf.best_estimator_
model_coef = pd.Series(dict(zip(features.columns[(clf.coef_ !=0)[0]],clf.coef_[(clf.coef_ != 0)])))

print("Model predictors with coefficients greater than 0:")
print(model_coef)

os.chdir(cwd)
#output total feature numpy array as serialized pickle file
features.columns.values.dump(metafile.split('.')[0]+"_feature_array.pickle")
#output sparse model coefficients
model_coef.to_csv(metafile.split('.')[0]+"_model_coefficients.tsv",sep='\t')
#output serialized classifier object for persistence
joblib.dump(clf,metafile.split('.')[0]+"_CLF.joblib")
