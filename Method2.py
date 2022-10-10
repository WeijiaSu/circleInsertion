import pandas as pd
import os
from Bio import SeqIO

TE="/data/zhanglab/Weijia_Su/CommonDataSet/TE_full/HMS-Beagle.fasta"
OutName="171107"


def Selecting(reads,circle,insertion):
	cir=pd.read_table(circle)
	ins=pd.read_table(insertion)

	cir_l=set(cir["Readname"])
	ins_l=set(ins["read"])
	
	print(len(cir_l))
	print(len(ins_l))
	records=SeqIO.parse(reads,"fasta")
	SeqIO.write((rec for rec in records if rec.id not in cir_l and rec.id not in ins_l),OutName+"_selectedReads.fa","fasta")

reads="/data/zhanglab/Weijia_Su/2021_fly_ecc/Fig2/NonGFP_171107/171107_LW1_aubago_eggs.fastq.chop.fastq-HMS-Beagle.fa.TE+GFP_.fa"
circle="/data/zhanglab/Weijia_Su/2021_fly_ecc/Fig2/NonGFP_171107/1206/171107_LW1_aubago_eggs.fastq.chop.fastq-TE_full.fa.TE+GFP_circles.txt"
insertion="/data/zhanglab/Weijia_Su/2021_fly_ecc/Fig2/NonGFP_171107/171107_LW1_aubago_eggs.fastq.chop.fastq-HMS-Beagle.fa.TE+GFP+_AllTEInsertion_Final_insertion.tsv"


#Selecting(reads,circle,insertion)


def Mapping(reads):
	mapping="minimap2 -x map-ont -t 4 %s %s -Y > %s"%(TE,reads,OutName+"_HMS.paf")
	os.system(mapping)

Mapping(OutName+"_selectedReads.fa")
