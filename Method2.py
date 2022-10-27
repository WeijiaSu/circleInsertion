import pandas as pd
import os
from Bio import SeqIO
import numpy as np

pd.set_option("display.max_columns",40)

#TE="/data/zhanglab/Weijia_Su/CommonDataSet/TE_full/HMS-Beagle.fasta"
TE="/data/zhanglab/Weijia_Su/CommonDataSet/TE_full.fa"
Genome="/data/zhanglab/Weijia_Su/Genomes/Dro/DM6/dm6_RM_1004/dm6.fa.masked"
OutName="171107_allTE"
reads="/data/zhanglab/Weijia_Su/Nanopore_Raw_Data/171107_LW1/Fly_HMS-Beagle-GFP_TEactive_gDNA.fastq"
Jun_type="1LTR_FL"



def FilterPaf(paf):
	f=pd.read_table(paf,header=None,sep=" ")
	print(f[0:10])
	print(f.shape)
	linAli=f.groupby([0]).filter(lambda x: len(x)==1)
	f=f.loc[~f[0].isin(linAli[0])]
	print(f.shape)
	f.to_csv(OutName+"_MultiAlig.tsv",header=None,index=None,sep="\t")

#FilterPaf(OutName+"_TEmap.paf")

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
	f=pd.read_table(mutiAlig,header=None)
	f["info"]=f[0]+"_"+f[5]
	print(f[0:10])
	print(f.shape)
	print(len(set(f[0])))
	f["info"]=f[0]+"_"+f[5]
#	f[12]=f[1]-f[11]
#	f=f.loc[f[12]>=100]
	f=f.sort_values([0,2,3])
	d={}
	for r in set(f["info"]):
		d[r]="NC"
		sub=f.loc[f["info"]==r]
		l=list(zip(sub[2],sub[3],sub[7],sub[8],sub[4],sub[6],sub[1]))
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
	
#JunctionReads(OutName+"_MultiAlig.tsv")

def circleType(x):
    cirCle_S=int(x.split("_")[0])
    cirCle_E=int(x.split("_")[1])
    cirCle_D=int(x.split("_")[-1])
    if cirCle_S <100 and cirCle_E>7062-100 and cirCle_D<-100:
        return "1LTR_FL"
    elif cirCle_S <100 and cirCle_E>7062-100 and cirCle_D>-100:
        return "2LTR_FL"
    elif cirCle_S <100:
        return "1LTR5_Frg"
    elif cirCle_E > 7062-100:
        return "1LTR3_Frg"
    else:
        return "nonLTR_Frg"

def GetCirType(Junction_reads):
	f=pd.read_table(Junction_reads,header=None)
	f=f.loc[f[14]!="NC"]
	f["Type"]=f[14].apply(lambda x:circleType(x))
	f.to_csv(OutName+"_Type.tsv",index=None,header=None,sep="\t")

#GetCirType(OutName+"_junction.tsv")

def GenomeMaapping(Jun_reads,Jun_type,reads):
	f=pd.read_table(Jun_reads,header=None)
	f=f.loc[f[15]==Jun_type]
	print(f[0:10])
	print(f.shape)
	s=set(f[0])
	print(len(s))
	records=SeqIO.parse(reads,"fastq")
	SeqIO.write((rec for rec in records if rec.id in s),OutName+"_"+Jun_type+".fa","fasta")
	mapping="minimap2 -x map-ont -t 4 %s %s -Y | awk '{for(i=1;i<=12;i++) printf $i" "; print ""}' FS='\t' > %s"%(Genome,OutName+"_"+Jun_type+".fa",OutName+"_"+Jun_type+".paf")
	os.system(mapping)

#GenomeMaapping(OutName+"_Type.tsv",Jun_type,reads)


def CombineMapping(Gmap,Tmap):
	g=pd.read_table(Gmap,header=None,sep=" ")
	g=g.sort_values([0,2,3])
	g=g[range(0,12)]
	print(g[0:10])
	t=pd.read_table(Tmap,header=None,sep=" ")
	t=t.loc[t[0].isin(g[0])]
	t=t.sort_values([0,2,3])
	t=t[range(0,12)]
	print(t[0:10])
	g_columns=["rName","rLen","rGenome_s","rGenome_e","strand","gName","gLen","genome_s","genome_e","genome_match","genome_align","genome_score"]
	t_columns=["rName","rLen","rTE_s","rTE_e","strand","tName","tLen","t_s","t_e","te_match","te_align","te_score"]
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
	
	combined=combined.loc[(combined["rGenome_e"]<=combined["rTE_min"]+100) | (combined["rGenome_s"]>=combined["rTE_max"]-100)]
	print(combined[0:20])
	print(combined.shape)
	print(len(set(combined["rName"])))
	combined.to_csv(OutName+"_cirIns_filter1.tsv",index=None,sep="\t")

#CombineMapping(OutName+"_"+Jun_type+".paf",OutName+"_TEmap.paf")


def GetInsertion(CombineFile):
	f=pd.read_table(CombineFile)
	f=f.loc[f["tName"]!="HMS-Beagle"]
	f=f.loc[f["gName"]=="chrM"]
	r=list(set(f["rName"]))
	for read in r[0:2]:
		sub=f.loc[f["rName"]==read]
		print(sub)
#	f["rGen_min"]=0
#	f["rGen_max"]=0
#	for r in set(f["rName"]):
#		sub=f.loc[f["rName"]==r]
#		min_=sub["rGenome_s"].min()
#		max_=sub["rGenome_e"].max()
#		f.loc[f["rName"]==r,"rGen_min"]=min_
#		f.loc[f["rName"]==r,"rGen_max"]=max_
#	f=f.loc[(f["rGen_min"]<f["rTE_min"]) & (f["rGen_max"]>f["rTE_max"])]

#	f["d1"]=f["rGenome_e"]-f["rTE_min"]
#	f["d2"]=f["rGenome_s"]-f["rTE_max"]
#	f["d1"]=f["d1"].apply(lambda x: abs(x))
#	f["d2"]=f["d2"].apply(lambda x: abs(x))
#	f=f.loc[(f["d1"]<=500) | (f["d2"]<=500)]
#	print(f.shape)
#	print(f)
#	print(len(set(f["rName"])))
GetInsertion(OutName+"_cirIns_filter1.tsv")

