import pandas as pd

def MappingRes(Data,circle,insertion):
	cir=pd.read_table(circle)
	ins=pd.read_table(insertion)
	
	cir_l=set(cir["Readname"])
	ins_l=set(ins["read"])
	
	f=pd.read_table(Data,header=None)
	f=f[range(0,9)]
	sub=f.loc[(f[7]<3533-300) & (f[8]>3533+100)]
	f=f.loc[f[0].isin(sub[0])]
	f=f.loc[f[0].apply(lambda x: x not in cir_l and x not in ins_l)]
	print(f[0:10])
	print(f.shape)
	r=f.drop_duplicates([0],keep="first")
	print(r.shape)
	

MappingRes("171107_LW1_aubago_eggs.fastq.chop.fastq-HMS-Beagle.fa.sort.fa-HMS-Beagle_circle.fa.paf",
"/data/zhanglab/Weijia_Su/2021_fly_ecc/Fig2/NonGFP_171107/1206/171107_LW1_aubago_eggs.fastq.chop.fastq-TE_full.fa.TE+GFP__circleAnalyze.txt",
"/data/zhanglab/Weijia_Su/2021_fly_ecc/Fig2/NonGFP_171107/171107_LW1_aubago_eggs.fastq.chop.fastq-HMS-Beagle.fa.TE+GFP+_AllTEInsertion_Final_insertion.tsv")
