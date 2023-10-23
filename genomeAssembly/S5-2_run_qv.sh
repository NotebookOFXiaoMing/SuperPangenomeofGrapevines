#!/bin/bash -l
#SBATCH --job-name=S5-QV
#SBATCH --partition=cuPartition
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --nodes=1

fasta=$1
reads=$2
out=$3

export MERQURY=/home/wangj/softwave/merqury-1.3/
export PATH=/home/wangj/softwave/meryl-1.3/bin/:$PATH

/home/wangj/softwave/meryl-1.3/bin/meryl k=21 count output ${out}.meryl ${reads}
/home/wangj/softwave/merqury-1.3/merqury.sh ${out}.meryl ${fasta} ${out}.merqury

