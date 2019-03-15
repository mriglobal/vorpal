import os
import pandas as pd
import joblib
from sklearn.linear_model import LogisticRegressionCV
import argparse

parser = argparse.ArgumentParser(description="Takes feature-labeled .bed files and corresponding meta data as inputs to generate a sparse model of phenotype predictors.")

#command line arguments
parser.add_argument('--beds',required=True,help="Directory containing .bed files.")
parser.add_argument('-m', required=True,help="Meta data table for genomic records.")
parser.add_argument('-o',default=os.getcwd(),help="Output directory")
myargs=parser.parse_args()

cwd = myargs.o
metafile = myargs.m
# with open('Coronavirus_complete_features.txt','r') as infile:
#     features = [r.strip() for r in infile.readlines()]

meta = pd.read_table(metafile)
os.chdir(os.getcwd()+os.path.join(myargs.beds))

bed_files = [file for file in os.listdir() if '.bed' in file]

beds = [pd.read_table(f,names=['chr','start','end','name','score']) for f in bed_files]

for i in range(len(beds)):
    beds[i]['chr'] = [a[2] for a in beds[i]['chr'].str.split('|').values]


feature_dict = {}

for b in beds:
    accession_num = b['chr'].unique()[0]
    feature_dict[accession_num]=b['name'].value_counts().to_dict()
        
feature_table = pd.DataFrame(feature_dict).T.fillna(0.0)
feature_table.index.name = "accession"
feature_table.reset_index(inplace=True)

complete_table = pd.merge(feature_table,meta)

features = complete_table.drop(['accession','taxid','pathogenic'],axis=1).copy()

X = features.values
y = meta['pathogenic']

clf = LogisticRegressionCV(penalty='l1',solver='liblinear',tol=.000000001)
clf.fit(X,y)

print("Trained model score (accuracy):",clf.score(X,y))

model_coef = pd.Series(dict(zip(features.columns[(clf.coef_ !=0)[0]],clf.coef_[(clf.coef_ != 0)])))

print("Model predictors with coefficients greater than 0:")
print(model_coef)

os.chdir(cwd)
#output sparse model coefficients
model_coef.to_csv(metafile.split('.')[0]+"_model_coefficients.tsv",sep='\t')
#output serialized classifier object for persistence
joblib.dump(metafile.split('.')[0]+"_CLF.joblib")