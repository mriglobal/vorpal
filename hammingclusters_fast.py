import os
import pandas as pd
import time
from scipy.spatial.distance import pdist,squareform
from scipy.special import comb
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster, cophenet
from sklearn.metrics import pairwise_distances_chunked
import fastcluster
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
parser.add_argument('--temp',default=None,help="Location of temp directory for memory mapped matrix. (Default: None)")
parser.add_argument('--mem',default=0,type=int,help="Amount of memory allocated to process distance matrix chunks. (In MiB)")
myargs=parser.parse_args()
if not myargs.temp and myargs.mem:
    parser.error("Specifying memory allocation for chunked processing of distance matrix is only implemented with memory mapping. Give a path to a temp directory.")


os.chdir(os.getcwd())
def dna_to_numeric(s):
    '''integer encode dna sequence'''
    new = s.replace('G','0')
    new = new.replace('T','1')
    new = new.replace('C','2')
    new = new.replace('A','3')
    return new

def numeric_to_dna(s):
    '''decode integer encoding to dna sequence'''
    new = str(s).replace('0','G')
    new = str(new).replace('1','T')
    new = str(new).replace('2','C')
    new = str(new).replace('3','A')
    return(new)
    
def get_chunks(x, num_chunks):
    '''generic chunking recipe'''
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
#set up IUPAC degenerate character look up
basekey = {tuple(sorted(values)):keys for keys,values in DNA.degenerate_map.items()}
basekey.update({('A',):'A',('G',):'G',('C',):'C',('T',):'T'})
if chunks:
    print("Making {} data chunks. {}".format(chunks,time.asctime()))
    #coo converted to csr for row-wise slicing
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
    #fast numpy function for getting number of kmer appearances
    counts = pd.Series(np.count_nonzero(kmers.to_dense(),axis=1),index=kmers.index)
#removing very large objects from memory
del(kmers)
print("Finding high frequency kmers at quantile: {}. {}".format(quantile,time.asctime()))
hifreqkmers = counts[counts > counts.quantile(quantile)]
print("Converting kmers to integer arrays.{}".format(time.asctime()))
df = pd.DataFrame(hifreqkmers.index.map(dna_to_numeric))

vectordf = df[0].apply(lambda x:list(map(int,x)))

final = pd.DataFrame(vectordf)[0].apply(pd.Series)
#removing more garbage from memory
del(df)
del(vectordf)
print("Performing hamming distance analysis on high frequency kmers. {}".format(time.asctime()))

if myargs.temp:
    #writes distance matrix to specified location
    filename = os.path.join(myargs.temp,"distance.dat")
    kmerdist=np.memmap(filename, dtype='float32',mode='w+',shape=(comb(final.shape[0],2,exact=True),))
    if myargs.mem:
        dist_gen = pairwise_distances_chunked(final,metric='hamming',n_jobs=-1,working_memory=myargs.mem)
        position = final.shape[0]-1
        remaining = 1
        total = 0
        del(kmerdist)
        for temp in dist_gen:
            print("Processing distance matrix chunk. {}".format(time.asctime()))
            kmerdist=np.memmap(filename,dtype='float32',mode='r+')
            for r in range(len(temp)):
                kmerdist[total:total+position] = temp[r][remaining:]
                total+=position
                position-=1
                remaining+=1
            del(kmerdist)
        kmerdist=np.memmap(filename,dtype='float32',mode='r')
    else:
        print("Writing distance matrix to disk. {}".format(time.asctime()))
        kmerdist[:] = pdist(final,'hamming')
else:
    #if neither memory saving strategies are selected
    kmerdist = pdist(final,'hamming')

print("Building kmer tree using average linkage with an average number of allowed based of: {} {}".format(degen_base_num,time.asctime()))
Z = fastcluster.average(kmerdist)
kmer_length=final.shape[1]
maxdist=round((degen_base_num/kmer_length), 2)
clusters = fcluster(Z,maxdist,criterion='distance')
myclusters = {key:[] for key in set(clusters)}
for index, clust in enumerate(clusters):
    myclusters[clust].append(index)

clustergroups = []
for amp in Counter(clusters).keys():
    clustergroups.append(final.iloc[myclusters[amp]])

print("Building alignments for kmer motifs. {}".format(time.asctime()))
#group resulting clusters into de facto alignment objects
alignments = []
for c in clustergroups:
    group = [DNA(''.join(c.loc[i].map(numeric_to_dna))) for i in c.index]
    alignments.append(skbio.alignment.TabularMSA(group))

oligos = []
#find representative IUPAC base of observed positional variance
for n,a in enumerate(alignments):
    position_vars = [tuple(set(str(x))) for x in a.iter_positions()]
    degenseq = ''.join([basekey[tuple(sorted(p))] for p in position_vars])
    oligos.append(SeqRecord(Seq(degenseq,IUPAC.ambiguous_dna),id=pickle_file.split(sep='_')[0]+str(n),description=''))

print("Writing {} degenerate kmers as fasta file. {}".format(len(oligos),time.asctime()))
SeqIO.write(oligos,pickle_file.split(sep='.')[0]+str(degen_base_num)+'_degenerate_primers.fasta','fasta')

