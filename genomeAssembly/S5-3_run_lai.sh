#!/bin/bash -l
#SBATCH --job-name=S5-lai
#SBATCH --partition=cuPartition
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --nodes=1

fasta=$1
copyfasta=$2

export PATH=/home/wangj/softwave/LTR_FINDER_parallel-1.1/bin/:$PATH
export PATH=/home/wangj/softwave/LTR_retriever-2.9.0/bin/:$PATH
conda activate annotation

#rename the seqnames for LAI to work corrrectly
/home/wangj/softwave/seqkit replace -p "scaffold" -r "" $fasta > $copyfasta

perl /home/wangj/softwave/LTR_FINDER_parallel-1.1/LTR_FINDER_parallel \
	-seq $copyfasta -threads 49 -harvest_out

mv ${copyfasta}.finder.combine.scn ${copyfasta}.finder.scn

/home/wangj/softwave/LTR_retriever-2.9.0/LTR_retriever -threads 50 -genome ${copyfasta} \
	-inharvest ${copyfasta}.finder.scn

/home/wangj/softwave/LTR_retriever-2.9.0/LAI -t 50 -genome ${copyfasta} \
	-intact ${copyfasta}.pass.list -all ${copyfasta}.out

conda deactivate


