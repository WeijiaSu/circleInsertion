import pandas as pd
import os
from Bio import SeqIO
import numpy as np


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

#Mapping(OutName+"_selectedReads.fa")

def FilterPaf(paf):
	f=pd.read_table(paf,header=None)
	f=f.sort_values([0])
	f=f[range(11)]
	r=f.drop_duplicates([0],keep="first")
	linAli=f.groupby([0]).filter(lambda x: len(x)==1)
	f=f.loc[~f[0].isin(linAli[0])]
	f.to_csv(OutName+"_MultiAlig.tsv",header=None,index=None,sep="\t")
FilterPaf(OutName+"_HMS.paf")

def map_ratio(sub_f):
	r_len=int(list(sub_f[1])[0])
	l=zip(sub_f[2],sub_f[3])
	a=np.array([0]*r_len)
	for i in l:
		a[i[0]:i[1]+1]=1
	return list(a).count(1)

def GetMapping(MulAlig):
	f=pd.read_table(MulAlig,header=None)
	f["mLen"]=0
	for r in set(f[0]):
		sub_f=f.loc[f[0]==r]
		p=map_ratio(sub_f)
		f.loc[f[0]==r,"mLen"]=p
	f.to_csv(OutName+"_chemReads.tsv",header=None,index=None,sep="\t")
	print(f[0:10])
	print(f.shape)

#GetMapping(OutName+"_MultiAlig.tsv")

def Junction(list1,list2):
	n1,n2,n3,n4=list1[0],list1[1],list2[0],list2[1]
	if n3>n2 and n3-n2>100:
		return False
	elif abs(n3-n1)<100 or abs(n2-n4)<100:
		return False
	else:
		return True


def isCircle(list1,list2):

	n1,n2,n3,n4=list1[0],list1[1],list2[0],list2[1]
	t1,t2,t3,t4=list1[2],list1[3],list2[2],list2[3]
	n_max=max(n1,n2,n3,n4)
	n_min=min(n1,n2,n3,n4)
	d1,d2=list1[4],list2[4]
	t_min=min(t1,t2,t3,t4)
	t_max=max(t1,t2,t3,t4)
	j=Junction(list1,list2)
	readlength=list1[-1]
	if d1!=d2 or j==False:
		return False
	else:
		if d1=="+" and t_max==t2:
			return True
		if d1=="-" and t_max==t4:
			return True
		else:
			return False
	
def CircleType(list1,list2):

	n1,n2,n3,n4=list1[0],list1[1],list2[0],list2[1]
	t1,t2,t3,t4=list1[2],list1[3],list2[2],list2[3]
	t_max=max(t1,t2,t3,t4)
	t_min=min(t1,t2,t3,t4)
	n_max=max(n1,n2,n3,n4)
	n_min=min(n1,n2,n3,n4)
	Overlap=n3-n2
	reflength=list1[-2]
	readlength=list1[-1]
	if isCircle(list1,list2) ==False:
		return "NC"
	else:
		return str(t_min)+"_"+str(t_max)+"_"+str(Overlap)

def JunctionReads(chemReads):
	f=pd.read_table(chemReads,header=None)
	f[12]=f[1]-f[11]
	f=f.loc[f[12]>=100]
	f=f.sort_values([0,2,3])
	d={}
	for r in set(f[0]):
		d[r]="NC"
		sub=f.loc[f[0]==r]
		l=list(zip(sub[2],sub[3],sub[7],sub[8],sub[4],sub[6],sub[1]))
		i=0
		while i<len(l)-1:
			list1=l[i]
			j=i+1
			while j < len(l):
				list2=l[j]
				cirType=CircleType(list1,list2)
				if cirType!="NC":
					d[r]=cirType
					break
				else:
					j+=1
			else:
				i+=1
			break
	f[13]=f[0].apply(lambda x: d[x])
	f.to_csv(OutName+"_junction.tsv",header=None,index=None,sep="\t")	
	return f
	
#JunctionReads(OutName+"_chemReads.tsv")


