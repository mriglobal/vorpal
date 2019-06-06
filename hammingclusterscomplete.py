import os
import pandas as pd
import time
from scipy.spatial.distance import pdist,squareform
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster, cophenet
from collections import Counter
import skbio.alignment
import skbio.sequence
import numpy as np
from skbio import DNA
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import pickle
import argparse

parser = argparse.ArgumentParser(description="Create Degenerate Primers from Kmer Clusters")
parser.add_argument('-p', required=True, help="Pickled Kmer Sparse DataFrame")
parser.add_argument('-n', default=4.0, type=float, help="Average Number of Allowed Degenerate Bases")
parser.add_argument('-q', default=0.95, type=float, help="Quantile of clustered Kmers")
parser.add_argument('-c', default=0, type=int, help="Number of chunks to split matrix into for count processing. Default 0")
myargs=parser.parse_args()

os.chdir(os.getcwd())
def dna_to_numeric(s):
    new = s.replace('G','0')
    new = new.replace('T','1')
    new = new.replace('C','2')
    new = new.replace('A','3')
    return new

def numeric_to_dna(s):
    new = str(s).replace('0','G')
    new = str(new).replace('1','T')
    new = str(new).replace('2','C')
    new = str(new).replace('3','A')
    return(new)
    
def get_chunks(x, num_chunks):
     chunks=[]
     position=0
     for r in range(num_chunks):
             if r != num_chunks-1:
                     chunks.append(x[position:position+x.shape[0]//num_chunks])
             else:
                     chunks.append(x[position:])
             position+=x.shape[0]//num_chunks
     return chunks

degen_base_num = float(myargs.n)
pickle_file = myargs.p
quantile=myargs.q
chunks=myargs.c

print("Reading in kmer table matrix. {}".format(time.asctime()))
kmers = pd.read_pickle(pickle_file)
basekey = {tuple(sorted(values)):keys for keys,values in DNA.degenerate_map.items()}
basekey.update({('A',):'A',('G',):'G',('C',):'C',('T',):'T'})
if chunks:
    print("Making {} data chunks. {}".format(chunks,time.asctime()))
    kmer_splits = get_chunks(kmers.to_coo().tocsr(),chunks)
    total_counts = []
    index = 0
    print("Getting kmer counts. {}".format(time.asctime()))
    for split in kmer_splits:
        split_counts = np.count_nonzero(split.todense(),axis=1).A1
        total_counts.append(pd.Series(split_counts,index=kmers.index[index:index+len(split_counts)]))
        index+=len(split_counts)
    counts = pd.concat(total_counts)
else:
    print("Getting kmer counts. {}".format(time.asctime()))
    counts = pd.Series(np.count_nonzero(kmers.to_dense(),axis=1),index=kmers.index)
#removing very large objects from memory
del(kmers)
if quantile != 0.0:
    print("Finding high frequency kmers at quantile: {}. {}".format(quantile,time.asctime()))
    hifreqkmers = counts[counts > counts.quantile(quantile)]
else:
    hifreqkmers = counts
print("Converting kmers to integer arrays.{}".format(time.asctime()))
df = pd.DataFrame(hifreqkmers.index.map(dna_to_numeric))

vectordf = df[0].apply(lambda x:list(map(int,x)))

final = pd.DataFrame(vectordf)[0].apply(pd.Series)
#removing more garbage from memory
del(df)
del(vectordf)
print("Performing hamming distance analysis on high frequency kmers. {}".format(time.asctime()))

kmerdist = pdist(final,'hamming')


print("Building kmer tree using average linkage with an average number of allowed based of: {} {}".format(degen_base_num,time.asctime()))
Z = linkage(kmerdist, 'average')
kmer_length=final.shape[1]
maxdist=round((degen_base_num/kmer_length) + .01 , 2)
clusters = fcluster(Z,maxdist,criterion='distance')
myclusters = {key:[] for key in set(clusters)}
for index, clust in enumerate(clusters):
    myclusters[clust].append(index)

clustergroups = []
for amp in Counter(clusters).keys():
    clustergroups.append(final.iloc[myclusters[amp]])

#no idea what this was for
# output = []
# for c in clustergroups:
#     output.append((list(c.index),list(c.sum()[c.sum() > c.sum().quantile(.995)].index)))

print("Building alignments for kmer motifs. {}".format(time.asctime()))
alignments = []
for c in clustergroups:
    group = [DNA(''.join(c.loc[i].map(numeric_to_dna))) for i in c.index]
    alignments.append(skbio.alignment.TabularMSA(group))

oligos = []
for n,a in enumerate(alignments):
    position_vars = [tuple(set(str(x))) for x in a.iter_positions()]
    degenseq = ''.join([basekey[tuple(sorted(p))] for p in position_vars])
    oligos.append(SeqRecord(Seq(degenseq,IUPAC.ambiguous_dna),id=pickle_file.split(sep='_')[0]+str(n),description=''))

# kmercoverage = {}
# for alignment in alignments:
#     for sequence in alignment:
#         kmercoverage[str(sequence)] = kmers.T[kmers.T[str(sequence)] > 0][str(sequence)].index
print("Writing {} degenerate kmers as fasta file. {}".format(len(oligos),time.asctime()))
# pickle.dump(kmercoverage,open(pickle_file.split(sep='.')[0]+'_'+str(degen_base_num)+'degenerate_coverage.pickle',mode='wb'))
SeqIO.write(oligos,pickle_file.split(sep='.')[0]+str(degen_base_num)+'_degenerate_primers.fasta','fasta')
# with open(pickle_file.split(sep='.')[0]+'_summary.txt','w') as outfile:
#     outfile.writelines(['Mean cosine distance of references :\t', str(np.mean(pdist(kmers.T,'cosine'))),'\n'])
#     outfile.writelines(['Mean cophenetic hamming distance of hifreq ',str(kmer_length), 'mers:\t', str(np.mean(cophenet(Z))),'\n'])
# ax = sns.heatmap(distmat)
#pd.Series(kmerdist).plot.hist(bins=18)

