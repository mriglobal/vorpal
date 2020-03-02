import os
import pandas as pd
import numpy as np
import time
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import nt_search
from multiprocessing import Pool
import pickle
from skbio import DNA
import argparse

parser = argparse.ArgumentParser(description="Mapping Degenerate Kmers to Reference Sequences")
parser.add_argument('-r',required=True, help="Concatenated References Fasta")
parser.add_argument('-k',required=True, help='Degenerate Kmers Fasta')
parser.add_argument('-t',default=os.cpu_count(),type=int,help="Number of threads. Default all.")

myargs=parser.parse_args()
def pscore(primer,primerset):
'''function that populates bed column number 5
currently this is set to number of records the motif appears in
divided by the number of total records'''
    setcompletion = len(primerset)/len(totalset)
    return {'setcompletion':setcompletion,'total_score':int(round((setcompletion* 1000),0))}
    

sfile = myargs.r
kfile = myargs.k
os.chdir(os.getcwd())
cores = myargs.t

def make_splits(x, num_splits):
'''generic chunking recipe'''
    splits = []
    position = 0
    for r in range(num_splits):
        if r!= num_splits-1:
            splits.append(x[position:position+len(x)//num_splits])
        else:
            splits.append(x[position:])
        position+=len(x)//num_splits
    return splits

print("Loading references and motifs. {}".format(time.asctime()))
myseqs = list(SeqIO.parse(sfile,'fasta',IUPAC.unambiguous_dna))

totalset = set([r.id for r in myseqs])

kmers = list(SeqIO.parse(kfile,'fasta',IUPAC.ambiguous_dna))

alignments = {}

def find_alignments(kmer_list):
'''alignment function using nt_search'''
    my_align = {}
    for k in kmer_list:
        my_align[str(k.seq)] = {r.id:nt_search(str(r.seq),str(k.seq))[1:] for r in myseqs}
    return my_align
    
def multi_map(func, data):
'''generic multiprocess func mapping'''
    with Pool(cores) as pool:
        kmer_splits = make_splits(kmers,cores)
        results = pool.map(func, kmer_splits)
        return results

print("Mapping motifs with {} cores. {}".format(cores,time.asctime()))
mapped_kmers = multi_map(find_alignments,kmers)

for m in mapped_kmers:
    alignments.update(m)

#there may be a more effecient way to do this ¯\_(ツ)_/¯
del(mapped_kmers)

print("Building score table and bed files. {}".format(time.asctime()))
score_table = {}
bed = []
for primer in alignments.keys():
    score = pscore(primer,set([k for k in alignments[primer].keys() if alignments[primer][k]]))
    score_table[primer] = score
    for records in alignments[primer].items():
        if records[1]:
            for pos in records[1]:
                bed.append({'chr':records[0],'start':pos,'end':pos+len(primer),'name':primer,'score':score['total_score']})

print(len(bed))
score_df = pd.DataFrame(score_table).T
bed_df = pd.DataFrame(bed)
bed_df = bed_df[['chr','start','end','name','score']]
beds = bed_df.groupby('chr')
print("Writing bed files. {}".format(time.asctime()))
for accession in beds.groups.keys():
    beds.get_group(accession).to_csv(accession.replace('|','_')+'_primers.bed',sep='\t',header=False,index=False)
score_df.to_csv(kfile+"_score_table.csv")
