import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
import pandas as pd
pd.set_option("display.max_columns",40)


def Junction(list1,list2):
	d1,d2=list1[4],list2[4]
	n1,n2,n3,n4=list1[0],list1[1],list2[0],list2[1]
	t1,t2,t3,t4=list1[2],list1[3],list2[2],list2[3]
	l=list1[5]
	
	if d1!=d2:
		return False
	if max([t1,t2,t3,t4])-min([t1,t2,t3,t4])>=200000:
		return False
	
	if d1=="+" and d2=="+":
		if n2>=l-100 and n3<=100 and n1>100 and n4<=l-100 and t3<t2 and t3-t2>-1000 and t4-t2>=100:
			return True
		else:
			return False
	elif  d1=="-" and d2=="-":
		if n1<=100 and n4>=l-100 and n2<l-100 and n3>100 and t3<t2 and t3-t2>-1000 and t4-t2>=100:
			return True
		else:
			return False


def getJunction(blastn):
	f=pd.read_table(blastn,header=None,sep="\t",comment="#")
	f=f.groupby([0,1]).filter(lambda x: len(x)>1)
	f[11]=0
	f[12]=0
	f[13]=""
	f=f.sort_values([0,1,11,12])
	f.loc[f[6]<f[7],11]=f[6]
	f.loc[f[6]<f[7],12]=f[7]
	f.loc[f[6]>f[7],11]=f[7]
	f.loc[f[6]>f[7],12]=f[6]
	f.loc[f[6]<f[7],13]="+"
	f.loc[f[6]>f[7],13]="-"
	f[14]=f[0].apply(str)+"_"+f[1].apply(str)
	f[15]=False
	r=f.drop_duplicates([14],keep="first")
	for r in list(r[14]):
		sub=f.loc[f[14]==r]
		l=list(zip(sub[4],sub[5],sub[11],sub[12],sub[13],sub[9]))
		i=0
		while i<len(l)-1:
			list1=l[i]
			j=i+1
			while j < len(l):
				list2=l[j]
				Junc=Junction(list1,list2)
				if Junc==False:
					j=j+1
				else:
					t1,t2,t3,t4=list1[2],list1[3],list2[2],list2[3]		
					t_min=min([t1,t2,t3,t4])
					t_max=max([t1,t2,t3,t4])
					s_local=sub.loc[(sub[11]>=t_min-5000) & (sub[12]<=t_max+5000)]
					if s_local.shape[0]<3:
						f.loc[(f[14]==r)&(f[4]==list1[0])&(f[5]==list1[1])&(f[11]==list1[2])&(f[12]==list1[3]),15]=Junc
					
						f.loc[(f[14]==r)&(f[4]==list2[0])&(f[5]==list2[1])&(f[11]==list2[2])&(f[12]==list2[3]),15]=Junc
						break
					else:
						j+=1
			i+=1
	s=f.loc[f[15]!=False]
	f.to_csv(blastn+"sel.tsv",header=None,index=None,sep="\t")
	print(s.shape)
	print(s[0:20])

blast=sys.argv[1]
getJunction(blast)
