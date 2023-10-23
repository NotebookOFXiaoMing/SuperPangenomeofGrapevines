#!/bin/bash -l

#SBATCH --job-name=S4_run_3ddnareview
#SBATCH --partition=m256Partition
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --nodes=1

set -e

#change arguments here
maindir=/home/wangj/result/grape_assemblies/

metafile=/home/wangj/result/public/grape_file_paths/Grape_samples_and_datapaths.csv
ref=/home/wangj/result/call-snp/grape/references/PN40024.v4_11_05_21/PN40024.v4.REF.fasta

#inputs the sample name etc from the metafile via SLURM_ARRAY_TASK_ID
input="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F "," ' $1 == var {print;} ' $metafile)"
name=`echo $input | awk -F "," '{ print $2}'`
longname=`echo $input | awk -F "," '{ print $3 , " cultivar:" , $4}'`
hifi=`echo $input | awk -F "," '{ print $6}'`
hicfiles=`echo $input | awk -F "," '{ print $7}'`
hic=${hicfiles::-7}


echo $(date)" -- Starting 3D-DNA of manually reviewed assembly for sample "$name" ("$longname")."

asmprefix=$maindir$name"/cleanasm/"$name".hap"
asmsuffix=".p_ctg.clean.fna"

#load environment
export PATH=/home/wangj/softwave/parallel-20200122/bin/:$PATH
module load anaconda3/latest
#export PATH=/usr/bin/:$PATH

for i in 1 2 #haplotype 1 and 2
do

mkdir $maindir$name"/3ddna2/hap"$i  ##### put the reviewed assembly files to ../3ddna2
cd $maindir$name"/3ddna2/hap"$i
/home/wangj/softwave/3d-dna-201008/run-asm-pipeline-post-review.sh --sort-output -i 10000 -q 0 -r ../${name}.hap${i}.p_ctg.clean.0.review.assembly $asmprefix$i$asmsuffix ${maindir}${name}/juicer/hap${i}/aligned/merged_nodups.txt

assembly-stats ${name}.hap${i}${asmsuffix::-4}_HiC.fasta >> ${maindir}${name}/${name}.asm.stats

done


echo $(date)" -- The second 3D-DNA run is completed for "$name" ("$longname"). Move forward with the quality check."

###

echo $(date)" -- Starting minimap for sample "$name" ("$longname")."

#load environment
module unload anaconda3/latest
export PATH=/home/wangj/miniconda3/envs/assembly/bin/:$PATH
module load anaconda3/latest

for i in 1 2 #haplotype 1 and 2
do

cd $maindir$name"/3ddna2/hap"$i
asm=$name".hap"$i".p_ctg.clean_HiC.fasta"

minimap2 -x asm5 -t 50 $ref $asm > ${asm}.PN40024.paf
/home/wangj/softwave/dgenies-tools.py -i $asm -n ${name}.hap${i}_HiC -o ${asm}.idx

done


echo $(date)" -- Minimap for cleaned assemblies is completed for "$name" ("$longname"). You can upload them to D-Genies (https://dgenies.toulouse.inra.fr/run) to plot alignements."

