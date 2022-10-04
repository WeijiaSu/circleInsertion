#!/bin/bash

## Modify this job script accordingly and submit with the command:
##    sbatch HPC.sbatch
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=1   # 16 processor core(s) per node
#SBATCH --job-name='cirIns'
#SBATCH --partition="all"
#SBATCH --mem=50000
#SBATCH --mail-user=weijia.su@duke.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --output="cirIns-%j.out"
#SBATCH --error="cirIns-%j.err"
## LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

#reads="/data/zhanglab/Weijia_Su/2021_fly_ecc/Fig1/GFP_171107/211203/171107_LW1_aubago_eggs.fastq.chop.fastq-HMS-Beagle.fa.TE+GFP+.fa"
reads="/data/zhanglab/Weijia_Su/2021_fly_ecc/Fig2/NonGFP_171107/171107_LW1_aubago_eggs.fastq.chop.fastq-HMS-Beagle.fa.sort.fa"
ref="/data/zhanglab/Weijia_Su/Git/circleInsertion/HMS-Beagle_circle.fa"

/data/zhanglab/Weijia_Su/anaconda3/bin/time -v bash cirIns.sh $reads $ref;
