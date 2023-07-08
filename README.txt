Method 1. 
1.	Map HMS-Beagle(+) reads to HMS-Beagle circle (Reads set A, bam1)
2.	Screen for reads contain the junction (mapping 3533-300 < s < e < 3533+100) (Reads set B)
3.	Map Reads set B to transposon masked dm6 genome (bam2)
4.	Map Reads set B to consensus HMS-Beagle.fa to check the original structure.(bam3)
5. 	Merge bam2 and bam3

**** bam1 and bam3 are not consistant. Changing the reference may also change the mapping result. ****

