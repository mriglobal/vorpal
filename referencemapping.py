import os
import pandas as pd
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import nt_search
from Bio.SeqUtils import MeltingTemp as mt
import pickle
from skbio import DNA
import argparse

parser = argparse.ArgumentParser(description="Mapping Degenerate Kmers to Reference Sequences")
parser.add_argument('-r',required=True, help="Concatenated References Fasta")
parser.add_argument('-k',required=True, help='Degenerate Kmers Fasta')
parser.add_argument('-t',default=54.0, help='Target Melting Temperature')

myargs=parser.parse_args()


binaries = [key for key,values in DNA.degenerate_map.items() if len(values) ==2]
trinaries = [key for key,values in DNA.degenerate_map.items() if len(values) > 2]

def tm_estimator(oligo):
    return mt.Tm_GC(oligo,strict=False) - (2*oligo.count('N'))

def pscore(primer,primerset,target_tm):
    setcompletion = len(primerset)/len(totalset)
    binarycomposition = (1-(len([base for base in primer if base in binaries])/len(primer)))
    trinarycomposition = (1-(len([base for base in primer if base in trinaries])/len(primer)))
    melting = (1 - (tm_optimal - tm_estimator(primer))/tm_optimal)
    return {'setcompletion':setcompletion,'binary_comp':binarycomposition,'trinary_comp':trinarycomposition,'tm':melting,'total_score':int(round((setcompletion * binarycomposition * trinarycomposition * melting * 1000),0))}
    

sfile = myargs.r
kfile = myargs.k
tm_optimal = myargs.t
os.chdir(os.getcwd())

myseqs = list(SeqIO.parse(sfile,'fasta',IUPAC.unambiguous_dna))

totalset = set([r.id for r in myseqs])

kmers = list(SeqIO.parse(kfile,'fasta',IUPAC.ambiguous_dna))

alignments = {}

for k in kmers:
    alignments[str(k.seq)] = {r.id:nt_search(str(r.seq),str(k.seq))[1:] for r in myseqs}

score_table = {}
bed = []
for primer in alignments.keys():
    score = pscore(primer,set([k for k in alignments[primer].keys() if alignments[primer][k]]),tm_optimal)
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
for accession in beds.groups.keys():
    beds.get_group(accession).to_csv(accession.replace('|','_')+'_primers.bed',sep='\t',header=False,index=False)
score_df.to_csv(kfile+"_score_table.csv")
# for seq in myseqs:
#     alignments[seq.id] = [nt_search(str(seq.seq),str(k.seq)) for k in kmers]
