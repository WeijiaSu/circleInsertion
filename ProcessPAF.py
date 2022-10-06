import pandas as pd
import os



def SelectReads(outName,Data,circle,insertion):
	cir=pd.read_table(circle)
	ins=pd.read_table(insertion)
	
	cir_l=set(cir["Readname"])
	ins_l=set(ins["read"])
	print(len(cir_l))
	print(len(ins_l))
	f=pd.read_table(Data,header=None)
	f=f[range(0,9)]
	sub=f.loc[(f[7]<3533-300) & (f[8]>3533+100)]
	f=f.loc[f[0].isin(sub[0])]
	f=f.loc[f[0].apply(lambda x: x not in cir_l and x not in ins_l)]
	print(f[0:10])
	print(f.shape)
	r=f.drop_duplicates([0],keep="first")
	print(r.shape)
	r[[0]].to_csv(outName+".candi.tsv",header=None,sep="\t",index=None)
	

def Re_Map(outName,Reads):
	seqtk="seqtk subseq %s %s > %s"%(Reads,outName+".candi.tsv",outName+".candi.fa")
	os.system(seqtk)
	#genome="/data/zhanglab/Weijia_Su/Genomes/Dro/dm6.fa.masked"
	genome="/data/zhanglab/Weijia_Su/Genomes/Dro/DM6/dm6_RM_1004/dm6.fa.masked"
	TE="/data/zhanglab/Weijia_Su/CommonDataSet/TE_full/HMS-Beagle.fasta"
	mapping1="minimap2 -x map-ont -t 4 %s %s -Y > %s"%(genome,outName+".candi.fa",outName+".candi_Genome.paf")
	os.system(mapping1)
	mapping2="minimap2 -x map-ont -t 4 %s %s -Y > %s"%(TE,outName+".candi.fa",outName+".candi_TE.paf")
	os.system(mapping2)


def MappingRes(GenomeMapping,TEMapping):
	genomeM=pd.read_table(GenomeMapping,header=None)
	genomeM=genomeM[range(0,9)]
	g1=genomeM.drop_duplicates([0],keep="first")

	TEM=pd.read_table(TEMapping,header=None)
	TEM=TEM[range(0,9)]
	t1=TEM.drop_duplicates([0],keep="first")
	
	TEM=TEM.loc[TEM[0].isin(genomeM[0])]
	
	print(genomeM[0:20])
	print(g1.shape)

	print(TEM[0:20])
	print(t1.shape)

	TEM=TEM.loc[TEM[0].isin(genomeM[0])]
	TEM=TEM.groupby([0]).filter(lambda x: len(x)==2)
	print(TEM.shape)
	print(TEM[0:10])


	combine=genomeM.merge(TEM,on=[0,1],how="inner")
	print(combine[0:20])
	print(combine.shape)
	c=combine.drop_duplicates([0],keep="first")
	print(c.shape)

#SelectReads("171107_LW1_aubago_eggs_HMS","171107_LW1_aubago_eggs.fastq.chop.fastq-HMS-Beagle.fa.sort.fa-HMS-Beagle_circle.fa.paf",
#"/data/zhanglab/Weijia_Su/2021_fly_ecc/Fig2/NonGFP_171107/1206/171107_LW1_aubago_eggs.fastq.chop.fastq-TE_full.fa.TE+GFP_circles.txt",
#"/data/zhanglab/Weijia_Su/2021_fly_ecc/Fig2/NonGFP_171107/171107_LW1_aubago_eggs.fastq.chop.fastq-HMS-Beagle.fa.TE+GFP+_AllTEInsertion_Final_insertion.tsv")

#reads="/data/zhanglab/Weijia_Su/2021_fly_ecc/Fig2/NonGFP_171107/171107_LW1_aubago_eggs.fastq.chop.fastq-HMS-Beagle.fa.sort.fa"
#Re_Map("171107_LW1_aubago_eggs_HMS",reads)

MappingRes("171107_LW1_aubago_eggs_HMS.candi_Genome.paf","171107_LW1_aubago_eggs_HMS.candi_TE.paf")
