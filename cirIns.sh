#!/bin/bash

read=$1
readName=$(basename $read);
#NanoStat --fastq $read -n $readName".stats"
#porechop -i $read --extra_end_trim 0 --discard_middle -o $readName".pre.fastq"
ref=$2
refName=$(basename $ref);
#minimap2 -ax map-ont $ref $readName".pre.fastq" -Y -t 16 | samtools view -bS | samtools sort > $readName"-"$refName".bam";
#samtools index $readName"-"$refName".bam";
#samtools stats $readName"-"$refName".bam" | head -n 30 > $readName"-"$refName".bam.stats";
#minimap2 -x map-ont -c -t 4 $ref $read -Y >$readName"-"$refName".paf"
minimap2 -ax map-ont $ref $read -Y -t 16 | samtools view -bS | samtools sort > $readName"-"$refName".bam";
python3 /data/zhanglab/Weijia_Su/PythonScrip/refCircle_Package/S1_Convert.py -bam $readName"-"$refName".bam";
