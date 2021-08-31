# -*- coding: utf-8 -*-
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import DBSCAN
import pandas as pd
#import hdbscan
import umap
import matplotlib.pyplot as plt
import numpy as np
import json
import argparse
import os
import logging

def is_file(string):
	if os.path.isfile(string):
		return string
	else:
		raise ValueError("Not a file")

parser = argparse.ArgumentParser(description="read 'amino_acid_data.txt' and attempt to reduce the amino alphabet using HDBscan")

# Params for aminos file/UMAP
parser.add_argument('-c', '--cols_to_inc', type=str, nargs='+', default=None,
                    help='list of columns in amino_acid_data.txt to use')
parser.add_argument('-d', '--drop_aminos', type=str, nargs='+', default=None,
                    help='one-lettered amino to hold out during clustering')
parser.add_argument('-a', '--aminos_file', type=is_file, default='amino_acid_data.txt',
                    help='amino_acid_data.txt file')
parser.add_argument('-o', '--out_prefix', type=str, default='amino_DBscan_clustering',
                    help='prefix of output .json file')
parser.add_argument('-u', '--umap_neighbors', type=int, default=2,
                    help='number of neighbors for UMAP to consider (default:2)')

# Params for DBScan
parser.add_argument('-m', '--min_samples', type=int, default=2,
                    help='min_samples for DBscan (points to be considered core point)')
parser.add_argument('-e', '--eps', type=float, default=.7,
                    help='epsilon parameter for DBscan (default .7)')
parser.add_argument('-g', '--algorithm', type=str, default='auto',
                    help='algorithm parameter for DBscan')
parser.add_argument('-p', '--metric_params', default=None,
                    help='Additional kwargs for DBscan')
parser.add_argument('-l', '--leaf_size', type=int, default=30,
                    help='Leaf size for DBscan (default:30)')
parser.add_argument('-x', '--metric', type=str, default='euclidean',
                    help='distance metric for DBscan')


args = parser.parse_args()
cols_to_inc = args.cols_to_inc
drop_aminos = args.drop_aminos
aminos_file = args.aminos_file
out_prefix = args.out_prefix
umap_neighbors = args.umap_neighbors
min_samples = args.min_samples
eps = args.eps
algorithm = args.algorithm
metric_params = args.metric_params
leaf_size = args.leaf_size
metric = args.metric


def run_analysis(cols_to_inc=cols_to_inc,
                 drop_aminos=drop_aminos,
                 aminos_file=aminos_file,
                 out_prefix=out_prefix,
                 umap_neighbors=umap_neighbors,
                 min_samples=min_samples,
                 eps=eps,
                 algorithm=algorithm,
                 metric_params=metric_params,
                 leaf_size=leaf_size,
                 metric=metric):
    
    
    AA_df = pd.read_table(aminos_file)
    AA_df.set_index('amino', inplace=True)
    
    if cols_to_inc:
        AA_df=AA_df.filter(cols_to_inc, axis=1)
    
    if drop_aminos:
        AA_df=AA_df.drop(drop_aminos, axis=0)
    
    # Agglomerative Heuristics of Amino Hierarchies AHAH
    for objcol in ['prop1','prop2','side_chain']:
        if objcol in AA_df.columns:
            AA_df = pd.get_dummies(AA_df, columns=[objcol])
    
    if 'full_name' in AA_df.columns:
        AA_df.drop('full_name', axis=1, inplace=True)
    AA_df.replace('', 0, inplace=True)
    AA_df.astype(dtype=float)
    AA_df.fillna(0, inplace=True)
    
    scaler = StandardScaler()
    
    AA_df_scaled = scaler.fit_transform(AA_df)
    AA_df_scaled = pd.DataFrame(AA_df_scaled, columns=AA_df.columns, index=AA_df.index)
    
    reducer = umap.UMAP(n_neighbors=umap_neighbors)
    AA_umap_embedding = reducer.fit_transform(AA_df_scaled)
    clusterer = DBSCAN(eps=eps, 
                       min_samples=min_samples, 
                       metric=metric, 
                       metric_params=metric_params,
                       algorithm=algorithm,
                       leaf_size=leaf_size)
    
    clusterer.fit(AA_umap_embedding)
    
    
    fig, ax = plt.subplots()
    #colorstr = "bgrcmykw"
    colorstr = ['tab:blue','tab:orange','tab:green', 'tab:red', 'tab:purple',
                'tab:brown', 'tab:pink', 'tab:grey', 'tab:olive', 'tab:cyan',
                'black','salmon','gold','lime','darkslategray',
                'navy','mediumslateblue','indigo','fuchsia','maroon','chocolate']

    cdict = {}
    for i in range(len(colorstr)):
        cdict[i-1] = colorstr[i]
    
    assignment_dict = {}
    for g in np.unique(clusterer.labels_):
        
        ix = np.where(clusterer.labels_ == g)
        
        if g != -1:
            lab = AA_df_scaled.index[ix[0][0]]
            
        else:
            lab = g
        
        for j in ix[0]:
            assignment_dict[AA_df_scaled.index[j]] = str(lab)
        
        ax.scatter(AA_umap_embedding[ix[0],0], 
                   AA_umap_embedding[ix[0],1], 
                   c=cdict[g],
                   label = lab)
    
    ax.legend()
    plt.show()
    
    
    if drop_aminos:
        for exclusion in drop_aminos:
            assignment_dict[exclusion] = exclusion
    
    for key, value in assignment_dict.items():
        if value == '-1':
            assignment_dict[key] = key
    
    write_dict = {"amino_alpha":assignment_dict}
    
    with open(out_prefix+'.json', "w") as outfile:
        json.dump(write_dict, outfile)
        
    with open(out_prefix+'_params.json', 'w') as outfile:
        json.dump(vars(args), outfile)
        
        

if __name__ == "__main__":
    
    for i in range(10):
        run_analysis(cols_to_inc=cols_to_inc,
                 drop_aminos=drop_aminos,
                 aminos_file=aminos_file,
                 out_prefix=out_prefix,
                 umap_neighbors=umap_neighbors,
                 min_samples=min_samples,
                 eps=eps,
                 algorithm=algorithm,
                 metric_params=metric_params,
                 leaf_size=leaf_size,
                 metric=metric)
    
