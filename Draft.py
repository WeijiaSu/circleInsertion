import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
import pandas as pd
pd.set_option("display.max_columns",40)

f=pd.read_table("/data/zhanglab/Weijia_Su/2021_fly_ecc/Fig1/GFP_171107/211203/171107_LW1_aubago_eggs.fastq.chop.fastq-HMS-Beagle.fa.TE+GFP+circles.txt")
f=f.drop_duplicates(["Readname"],keep="first")
l=list(f["Readname"])





#f1=pd.read_table("171107_LW1_aubago_eggs.fastq.chop.fastq-HMS-Beagle.fa.TE+GFP+.fa-HMS-Beagle_circle.fa.bam_AligTable.tsv")
#f1=f1.sort_values(["Readname","ReadStart","ReadEnd"])
#f1=f1.loc[~f1["Readname"].isin(l)]
#sub1=f1.loc[(f1["RefStart"]<3533-100) & (f1["RefEnd"]>3533+100)] 
#f1=f1.loc[f1["Readname"].isin(sub1["Readname"])]
#print(f1[0:10])
#print(f1.shape)

#f2=pd.read_table("171107_LW1_aubago_eggs.fastq.chop.fastq-HMS-Beagle.fa.TE+GFP+.fa-HMS-Beagle_circle.fa.paf",header=None,sep="\t")
f2=pd.read_table("171107_LW1_aubago_eggs.fastq.chop.fastq-HMS-Beagle.fa.sort.fa-HMS-Beagle_circle.fa.paf",header=None,sep="\t")
f2=f2.sort_values([0,2,3])
f2=f2.loc[~f2[0].isin(l)]
sub2=f2.loc[(f2[7]<3277-100) & (f2[8]>3533+100)]
f2=f2.loc[f2[0].isin(sub2[0])]
f2=f2[range(0,9)]
reads=set(f2[0])
new_f=pd.DataFrame(columns=range(0,11))
for r in reads:
	sub_f=f2.loc[f2[0]==r]
	min_s=min(list(sub_f[2]))
	max_e=max(list(sub_f[3]))
	sub_f[9]=min_s
	sub_f[10]=max_e
	new_f=new_f.append(sub_f)

new_f=new_f.loc[(new_f[9]>=100) & (new_f[10]<=new_f[1]-100)]
new_f.to_csv("171107_LW1_aubago_eggs.fastq.chop.fastq-HMS-Beagle.fa.sort.fa-HMS-Beagle_circle.fa.paf.candi.tsv",index=None,sep="\t")
print(new_f[0:100])
#print(new_f)
print(new_f.shape)

