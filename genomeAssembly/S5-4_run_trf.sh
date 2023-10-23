#!/bin/bash -l
#SBATCH --job-name=S5-trf 
#SBATCH --partition=m256Partition
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --nodes=1

fasta=$1
out=$2

/home/wangj/softwave/TRF/build/bin/trf ${fasta} 2 7 7 80 10 90 2000 -d -m -l 2 -h
python /home/wangj/softwave/TRF2GFF-master/TRF2GFF.py -d ${out}.2.7.7.80.10.90.2000.dat -o ../${out}_trf.gff3
