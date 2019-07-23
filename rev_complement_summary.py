import os
import pandas as pd
from skbio import DNA
import argparse

parser = argparse.ArgumentParser(description="Analysis and extraction of motifs in the feature array that are reverse complements of each other.")
parser.add_argument('-f',required=True,help="Pickled feature array producted by linear model.")
myargs = parser.parse_args()
kmers = pd.read_pickle(myargs.f)

def rev_comp(x):
    return str(DNA(x).reverse_complement())
    
feature_set = set(kmers)

reverse_set = set(map(rev_comp,kmers))

rev_comp_set = feature_set.intersection(reverse_set)

rev_comp_ratio = len(rev_comp_set)/len(feature_set)

with open(myargs.f.replace(".pickle","_complement_stats.txt"),'w') as outfile:
    outfile.write("Ratio of reverse complement motifs to total motifs in feature array: {}".format(rev_comp_ratio))

pd.Series(list(rev_comp_set)).to_csv(myargs.f.replace(".pickle","_rev_comp.txt",sep='\t',index=False))