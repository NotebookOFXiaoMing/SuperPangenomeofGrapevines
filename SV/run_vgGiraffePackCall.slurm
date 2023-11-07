#!/bin/bash


#SBATCH --job-name="VGcall"
#SBATCH -n 4 #threads
#SBATCH -N 1 #node number
#SBATCH --mem=4000
#SBATCH --partition=tcum256c128Partition
#SBATCH --mail-user=mingyan24@126.com
#SBATCH --mail-type=BEGIN,END,FAIL

source activate syri
snakemake --cluster "sbatch --output=/data/myan/raw_data/practice/svpan.pipeline/slurm.logs/%j.out \
--error=/data/myan/raw_data/practice/svpan.pipeline/slurm.logs/%j.out --cpus-per-task={threads} \
--mail-type=END,FAIL --mail-user=mingyan24@126.com --mem={resources.mem_mb} \
--partition=tcum256c128Partition" \
--jobs 20 -s vgGiraffePackCall.smk --configfiles config02.yaml -k -p