from Bio import SeqIO
from Bio.SeqUtils import GC
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from sklearn.utils import resample
import scipy.sparse
import pandas as pd
import pickle
import argparse
import os

parser = argparse.ArgumentParser(description="Kmer Counter for RNA viruses")
parser.add_argument('-r',required=True,help='Reference genome')
parser.add_argument('--seqs', required=True,help='Folder for Sequences')
parser.add_argument('-n',default=0, type=int, help='Number of samples for stratified resampling.')
parser.add_argument('-m',default=None,help='Metadata file containing group labels for resampling.')
parser.add_argument('-k', nargs='?', default=17,help='K size for counting: default 17')
parser.add_argument('-p',default=.1, type=float, help='Percent Variance in reference length for replicon binning')
myargs = parser.parse_args()

if not myargs.n and myargs.m:
    parser.error("Metadata table provided but no resampling was specified.")
if myargs.n and not myargs.m:
    parser.error("Stratified resampling not possible without labels for 'groups'.")

if myargs.m:
    meta = pd.read_table(myargs.m)

ksize = int(myargs.k)
percent = myargs.p

def get_accession(seq_id):
    '''accession parser necessary to handle RVDB input'''
    if "|" in seq_id:
        return seq_id.split("|")[2]
    else:
        return seq_id

def replicontest(ref, example):
    '''test to see if given sequence should be binned
with the reference sequence based on length'''
    refindex = 0
    for i in range(len(ref)):
        if len(example) < ref[refindex][1]+(ref[refindex][1]*percent) and len(example) > ref[refindex][1]-(ref[refindex][1]*percent):
            return refindex
        refindex += 1
    return -1

def kmercount(replicon, k):
    '''canonical k-mer counting function'''
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
output_prefix = reference.split(sep='.')[0]+"_"+str(ksize)+"mers_"

#reference genome for binning
refrec = list(SeqIO.parse(reference,"fasta", IUPAC.unambiguous_dna))
#reference genome stasistics. GC content as binning criteria never implemented
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
        #bins sequences into groups based on provided reference. -1 indicates null group
        for r in records:
            i = replicontest(refstats, r)
            print(i)
            if i >= 0:
                recdata[i].append(r)

if myargs.n:
    rec_dict = {}
    groups_dict = dict(meta[['accession','groups']].values)
    for i in recdata:
        rec_dict[i] = {get_accession(s.id):s for s in recdata[i]}
        group_labels = [groups_dict[seq_id] for seq_id in rec_dict[i]]
        stratified = resample(pd.Series([seq_id for seq_id in rec_dict[i]]),n_samples=myargs.n,stratify=group_labels)
        recdata[i] = [rec_dict[seq_id] for seq_id in stratified]


kmerseries = {key:[] for key in refkmerdicts}
labels = {key:[] for key in refkmerdicts}
kmerset = {key:set(refkmerdicts[key].keys()) for key in refkmerdicts.keys()}
kmerindex = {key:[] for key in refkmerdicts}

for n in kmerseries.keys():
    kmerseries[n].append(refkmerdicts[n])
    #labels[n].append(refrec[n].id)

for repliconexample in recdata:
    for i,sequence in enumerate(recdata[repliconexample]):
        kmers = kmercount(str(sequence.seq),ksize)
        kmerseries[repliconexample].append(kmers)
        #labels[repliconexample].append(sequence.id)
        kmerset[repliconexample].update(kmers.keys())
        print(str(i),len(sequence))
    kmerindex[repliconexample] = {k: i for i, k in enumerate(sorted(list(kmerset[repliconexample])))}


#recipe for building COO sparse data object
data = {key:[] for key in refkmerdicts}
row = {key:[] for key in refkmerdicts}
column = {key:[] for key in refkmerdicts}

for n in kmerseries.keys():
    for c,k in enumerate(kmerseries[n]):
        for key in k.keys():
            data[n].append(k[key])
            column[n].append(kmerindex[n][key])
            row[n].append(c)

kmer_coo = {key:scipy.sparse.coo_matrix((data[key],(row[key],column[key])), shape=(len(kmerseries[key]), len(kmerset[key]))) for key in kmerseries.keys()}

for n in kmerseries.keys():
    out_name = output_prefix+"_"+str(n)
    if myargs.n:
        out_name = out_name+"_"+str(myargs.n)+"resample"
    sparse_df = pd.SparseDataFrame(kmer_coo[n].tocsr()).T.fillna(0.0)
    sparse_df.index = kmerindex[n].keys()
    #sparse_df.columns = labels[n]
#remove k-mers containing ambiguous bases
    sparse_df = sparse_df[~sparse_df.index.str.contains('H|V|Y|W|D|K|B|N|M|S|R')]
    sparse_df.to_pickle(out_name+str(n)+"_sparse.pickle")
