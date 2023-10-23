#!/bin/bash -l
#SBATCH --job-name=S5-busco
#SBATCH --partition=cuPartition
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --nodes=1

fasta=$1
out=$2

conda activate busco

busco -i $fasta -l /home/wangj/result/public/eudicots_odb10/ \
	-o $out -m genome --force --offline --augustus --cpu 50

conda deactivate
