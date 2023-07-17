import pandas as pd
import os
from Bio import SeqIO
import numpy as np
import argparse


pd.set_option("display.max_columns",40)
parser=argparse.ArgumentParser()
parser.add_argument("-TEmap","--TEpaf")
parser.add_argument("-OutName","--OutName")
args=parser.parse_args()

TEmap=args.TEpaf
OutName=args.OutName
Jun_type="1LTR_FL"
reads="/data/zhanglab/Weijia_Su/eccDNA/230619_fly_F2egg_Lig4Aub_gDNA_Tn5/%s.fastq.pre.fastq"%(OutName)


Chromosome=["chr2L","chr2R","chr3L","chr3R","chr4","chrX","chrY"]
def getChimeric_reads(infile):
	f=pd.read_table(infile)
	f["0"]=0
	file_A=f.drop_duplicates(["QName"],keep="first")[["QName","0","QLen"]]
	file_B=f[["QName","QStart","QEnd"]]
	file_A.to_csv(OutName+"_A.tsv",header=None,index=None,sep="\t")
	file_B.to_csv(OutName+"_B.tsv",header=None,index=None,sep="\t")
	bedtools="bedtools coverage -a %s -b %s> %s"%(OutName+"_A.tsv",OutName+"_B.tsv",OutName+"_cov.tsv")
	os.system(bedtools)
	cov=pd.read_table(OutName+"_cov.tsv",header=None)
	cov=cov.loc[cov[6]<=0.9]
	f=f.loc[f["QName"].isin(list(cov[0]))]
	f.to_csv(OutName+"_MultiAlig.tsv",index=None,sep="\t")
	rm="rm %s %s %s"%(OutName+"_A.tsv",OutName+"_B.tsv",OutName+"_cov.tsv")
	os.system(rm)


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
	
def JunCoor(list1,list2):

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

def JunctionReads(mutiAlig):
	f=pd.read_table(mutiAlig)
	f["info"]=f["QName"]+"_"+f["RName"]
	f=f.sort_values(["QName","QStart","QEnd"])
	d={}
	for r in set(f["info"]):
		d[r]="NC"
		sub=f.loc[f["info"]==r]
		l=list(zip(sub["QStart"],sub["QEnd"],sub["RStart"],sub["REnd"],sub["Strand"],sub["RLen"],sub["QLen"]))
		i=0
		while i<len(l)-1:
			list1=l[i]
			j=i+1
			while j < len(l):
				list2=l[j]
				Jcoor=JunCoor(list1,list2)
				if Jcoor!="NC":
					d[r]=Jcoor
					break
				else:
					j+=1
			else:
				i+=1
			break
	f[13]=f["info"].apply(lambda x: d[x])
	f.to_csv(OutName+"_junction.tsv",header=None,index=None,sep="\t")	
	return f
	

def circleType(x):
    cirCle_S=int(x.split("_")[0])
    cirCle_E=int(x.split("_")[1])
    cirCle_D=int(x.split("_")[2])
    TElength=int(x.split("_")[3])
    if cirCle_S <100 and cirCle_E>TElength-100 and cirCle_D<-100:
        return "1LTR_FL"
    elif cirCle_S <100 and cirCle_E>TElength-100 and cirCle_D>-100:
        return "2LTR_FL"
    elif cirCle_S <100:
        return "1LTR5_Frg"
    elif cirCle_E > TElength-100:
        return "1LTR3_Frg"
    else:
        return "nonLTR_Frg"

def GetCirType(Junction_reads):
	f=pd.read_table(Junction_reads,header=None)
	f=f.loc[f[11]!="NC"]
	f[12]=f[11]+"_"+f[6].apply(str)
	f["Type"]=f[12].apply(lambda x:circleType(x))
	f.to_csv(OutName+"_Type.tsv",index=None,header=None,sep="\t")

def GenomeMaapping(Jun_reads,Jun_type,reads):
	f=pd.read_table(Jun_reads,header=None)
	f=f.loc[f[12]==Jun_type]
	s=set(f[0])
	gmap=pd.read_table("%s.fastq_genome.paf"%(OutName),header=None)
	gmap=gmap.loc[gmap[0].isin(list(f[0]))]
	print(gmap[0:10])
	gmap=gmap.loc[gmap[5].isin(Chromosome)]
	print(gmap[0:10])
	gmap.to_csv(OutName+"_"+Jun_type+".Gmap.tsv",header=None,index=None,sep="\t")
	tmap=pd.read_table("%s.fastq_TE.paf"%(OutName),header=None,sep="\t")
	tmap=tmap.loc[tmap[0].isin(list(f[0]))]
	tmap.to_csv(OutName+"_"+Jun_type+".Tmap.tsv",header=None,index=None,sep="\t")

def CombineMapping(Gmap,Tmap):
	g=pd.read_table(Gmap,header=None)
	g=g.sort_values([0,2,3])
	g=g[range(0,9)]
	t=pd.read_table(Tmap,header=None)
	t=t.loc[t[0].isin(g[0])]
	t=t.sort_values([0,2,3])
	t=t[range(0,9)]
	g_columns=["rName","rLen","rGenome_s","rGenome_e","strand","gName","gLen","genome_s","genome_e"]
	t_columns=["rName","rLen","rTE_s","rTE_e","strand","tName","tLen","t_s","t_e"]
	g.columns=g_columns
	t.columns=t_columns
	combined=g.merge(t,on=["rName","rLen"],how="inner")
	combined["rTE_min"]=0
	combined["rTE_max"]=0
	for r in set(combined["rName"]):
		sub=combined.loc[combined["rName"]==r]
		min_=sub["rTE_s"].min()
		max_=sub["rTE_e"].max()
		combined.loc[combined["rName"]==r,"rTE_min"]=min_
		combined.loc[combined["rName"]==r,"rTE_max"]=max_
	combined1=combined.loc[(combined["rGenome_e"]<=combined["rTE_min"]+100) | (combined["rGenome_s"]>=combined["rTE_max"]-100)]
	r2=set(combined1["rName"])
	combined=combined.loc[combined["rName"].isin(r2)]
	print(len(set(combined["rName"])))
	combined.to_csv(OutName+"_"+Jun_type+"_cirIns_filter1.tsv",index=None,sep="\t")


def GetInsertion(CombineFile):
	f=pd.read_table(CombineFile)
	a=f.loc[f["rName"]=="38144642-1e3e-4a72-b5ff-9705bdb76ec9"]
	r=list(set(list(f["rName"])))
	#print(len(r))
	print(r)
	f=f.loc[f["tName"]!="HMS-Beagle"]
	fm=f.loc[f["gName"]=="chrM"]
	r=list(set(list(fm["rName"])))
	print(len(r))
	for read in r[20:]:
		print(r.index(read))
		sub=f.loc[f["rName"]==read]
		print("##################################")
		print(sub)
		print("")
	f["rGen_min"]=0
	f["rGen_max"]=0
	for r in set(f["rName"]):
		sub=f.loc[f["rName"]==r]
		min_=sub["rGenome_s"].min()
		max_=sub["rGenome_e"].max()
		f.loc[f["rName"]==r,"rGen_min"]=min_
		f.loc[f["rName"]==r,"rGen_max"]=max_
	f=f.loc[(f["rGen_min"]<f["rTE_min"]) & (f["rGen_max"]>f["rTE_max"])]

	f["d1"]=f["rGenome_e"]-f["rTE_min"]
	f["d2"]=f["rGenome_s"]-f["rTE_max"]
	f["d1"]=f["d1"].apply(lambda x: abs(x))
	f["d2"]=f["d2"].apply(lambda x: abs(x))
	f=f.loc[(f["d1"]<=500) | (f["d2"]<=500)]
	print(f.shape)
	print(f)
	print(len(set(f["rName"])))

getChimeric_reads(TEmap)
JunctionReads(OutName+"_MultiAlig.tsv")
GetCirType(OutName+"_junction.tsv")
GenomeMaapping(OutName+"_Type.tsv",Jun_type,reads)
CombineMapping(OutName+"_"+Jun_type+".Gmap.tsv",OutName+"_"+Jun_type+".Tmap.tsv")
#GetInsertion(OutName+"_"+Jun_type+"_cirIns_filter1.tsv")

