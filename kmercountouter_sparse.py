from Bio import SeqIO
from Bio.SeqUtils import GC
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
import scipy.sparse
import pandas as pd
import pickle
import argparse
import os

parser = argparse.ArgumentParser(description="Kmer Counter for RNA viruses")
parser.add_argument('-r',required=True,help='Reference genome')
parser.add_argument('--seqs', required=True,help='Folder for Sequences')
parser.add_argument('-k', nargs='?', default=18,help='K size for counting: default 18')
parser.add_argument('-p',default=.1, type=float, help='Percent Variance in reference length for replicon binning')

myargs = parser.parse_args()
ksize = int(myargs.k)
percent = myargs.p
def replicontest(ref, example):
    refindex = 0
    for i in range(len(ref)):
        if len(example) < ref[refindex][1]+(ref[refindex][1]*percent) and len(example) > ref[refindex][1]-(ref[refindex][1]*percent):
            return refindex
        refindex += 1
    return -1

def kmercount(replicon, k):
    kd = {}
    for j in range(len(replicon)-k):
        kmer = Seq(replicon[j:j+k],IUPAC.unambiguous_dna)
        if str(kmer) in kd.keys():
            kd[str(kmer)]+=1
            kd[str(kmer.reverse_complement())]+=1
        else:
            kd[str(kmer)]=1
            kd[str(kmer.reverse_complement())]=1
    return kd
            
os.chdir(os.getcwd())
reference = myargs.r
refrec = list(SeqIO.parse(reference,"fasta", IUPAC.unambiguous_dna))
refstats = [(GC(f.seq),len(f.seq)) for f in refrec]
refreplicon=[]
for s in refrec:
    refreplicon.append(str(s.seq))

refkmerdicts = {key:{} for key in range(len(refreplicon))}
for x,replicon in enumerate(refreplicon):
    refkmerdicts[x] = kmercount(replicon,ksize)

recdata = {key:[] for key in range(len(refstats))}
refdirectory = myargs.seqs
os.chdir(refdirectory)
filelist = os.listdir()
for file in filelist:
    if ".fna" in file or ".fasta" in file:
        recordname = file.rsplit(sep="_",maxsplit=1)[0]
        records = list(SeqIO.parse(file,format="fasta", alphabet=IUPAC.unambiguous_dna))
        #if len(records) <= len(refstats):
        for r in records:
            i = replicontest(refstats, r)
            print(i)
            if i >= 0:
                recdata[i].append(r)

kmerseries = {key:[] for key in refkmerdicts}
labels = {key:[] for key in refkmerdicts}
kmerset = {key:set(refkmerdicts[key].keys()) for key in refkmerdicts.keys()}
kmerindex = {key:[] for key in refkmerdicts}

for n in kmerseries.keys():
    kmerseries[n].append(refkmerdicts[n])
    labels[n].append(refrec[n].id)

for repliconexample in recdata:
    for i,sequence in enumerate(recdata[repliconexample]):
        kmers = kmercount(str(sequence.seq),ksize)
        kmerseries[repliconexample].append(kmers)
        labels[repliconexample].append(sequence.id)
        kmerset[repliconexample].update(kmers.keys())
        print(str(i),len(sequence))
    kmerindex[repliconexample] = {k: i for i, k in enumerate(sorted(list(kmerset[repliconexample])))}



data = {key:[] for key in refkmerdicts}
row = {key:[] for key in refkmerdicts}
column = {key:[] for key in refkmerdicts}

for n in kmerseries.keys():
    for c,k in enumerate(kmerseries[n]):
        for key in k.keys():
            data[n].append(k[key])
            column[n].append(kmerindex[n][key])
            row[n].append(c)

kmer_coo = {key:scipy.sparse.coo_matrix((data[key],(row[key],column[key])), shape=(len(labels[key]), len(kmerset[key]))) for key in kmerseries.keys()}

for n in kmerseries.keys():
    sparse_df = pd.SparseDataFrame(kmer_coo[n].tocsr()).T.fillna(0.0)
    sparse_df.index = kmerindex[n].keys()
    sparse_df.columns = labels[n]
    sparse_df = sparse_df[~sparse_df.index.str.contains('H|V|Y|W|D|K|B|N|M|S|R')]
    sparse_df.to_pickle(reference.split(sep='.')[0]+str(n)+"_"+str(ksize)+"mers_sparse.pickle")
# for e in refkmerdicts:
#     for r in recdata[e]:
#         refkmerdicts[e] = {qkmer:refkmerdicts[e][qkmer] for qkmer in refkmerdicts[e] if qkmer in r}
#         print(len(refkmerdicts[e]))