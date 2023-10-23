#!/bin/bash -l

#SBATCH --job-name=S3_run_3ddna
#SBATCH --partition=m256Partition
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --nodes=1

set -e

#change arguments here
maindir=/home/wangj/result/grape_assemblies/

metafile=/home/wangj/result/public/grape_file_paths/Grape_samples_and_datapaths.csv


#inputs the sample name etc from the metafile via SLURM_ARRAY_TASK_ID
input="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F "," ' $1 == var {print;} ' $metafile)"
name=`echo $input | awk -F "," '{ print $2}'`
longname=`echo $input | awk -F "," '{ print $3 , " cultivar:" , $4}'`
hifi=`echo $input | awk -F "," '{ print $6}'`
hicfiles=`echo $input | awk -F "," '{ print $7}'`
hic=${hicfiles::-7}


echo $(date)" -- Starting the 3D-DNA for sample "$name" ("$longname")."

asmprefix=$maindir$name"/cleanasm/"$name".hap"
asmsuffix=".p_ctg.clean.fna"

#load environment
export PATH=/home/wangj/softwave/parallel-20200122/bin/:$PATH
module load anaconda3/latest

#create working directories !DELETE FOLDERS FROM PREVIOUS ATTEMPTS!
mkdir $maindir$name"/3ddna"


for i in 1 2 #haplotype 1 and 2
do

mkdir $maindir$name"/3ddna/hap"$i
cd $maindir$name"/3ddna/hap"$i
bash /home/wangj/softwave/3d-dna-201008/run-asm-pipeline.sh -r 0 -q 0 --editor-repeat-coverage 10 $asmprefix$i$asmsuffix ${maindir}${name}/juicer/hap${i}/aligned/merged_nodups.txt

done


echo $(date)" -- First 3D-DNA run is completed for "$name" ("$longname"). Inspect the results in Juicebox on local and move forward with the second iteration of 3D-DNA."
