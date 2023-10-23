#!/bin/bash -l
#SBATCH --job-name=S5_postassembly_QC
#SBATCH --partition=cuPartition
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
## (this is a wrapper script, no need to run with large resources)

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

asmprefix=$maindir$name"/3ddna2/hap"
asmsuffix=".p_ctg.clean_HiC.fasta"

mkdir $maindir$name"/QC"

### BUSCO ###
echo $(date)" -- Starting BUSCO..."
mkdir $maindir$name"/QC/busco"
cd $maindir$name"/QC/busco/"
for i in 1 2 
do
sbatch ${maindir}S5-1_run_busco.sh ${asmprefix}${i}/${name}.hap${i}${asmsuffix} ${name}.hap${i} 
done
echo $(date)" -- BUSCO jobs have been submitted."

### QV ###
echo $(date)" -- Starting MERQURY..."
mkdir $maindir$name"/QC/qv"
cd $maindir$name"/QC/qv/"
for i in 1 2 
do
sbatch ${maindir}S5-2_run_qv.sh ${asmprefix}${i}/${name}.hap${i}${asmsuffix} ${hifi} ${name}.hap${i}
done
echo $(date)" -- MERQURY jobs have been submitted."

### LAI ###
echo $(date)" -- Starting LAI..."
mkdir $maindir$name"/QC/lai"
cd $maindir$name"/QC/lai/"
for i in 1 2 
do
sbatch ${maindir}S5-3_run_lai.sh ${asmprefix}${i}/${name}.hap${i}${asmsuffix} ${name}.hap${i}${asmsuffix}
done
echo $(date)" -- LTR Retriever jobs have been submitted."

### TRF ###
echo $(date)" -- Starting TRF..."
mkdir $maindir$name"/QC/trf"
cd $maindir$name"/QC/trf/"
for i in 1 2
do
mkdir $maindir$name"/QC/trf/hap"$i
cd $maindir$name"/QC/trf/hap"$i
sbatch ${maindir}S5-4_run_trf.sh ${asmprefix}${i}/${name}.hap${i}${asmsuffix} ${name}.hap${i}${asmsuffix}
done
echo $(date)" -- TRF jobs have been submitted."


echo $(date)" -- All post-assembly quality jobs have been submitted for the sample "$name" ("$longname")."
