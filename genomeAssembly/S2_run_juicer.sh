#!/bin/bash -l

#SBATCH --job-name=S2_run_juicer
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


echo $(date)" -- Starting the first scafolding with juicer for sample "$name" ("$longname")."

asmprefix=$maindir$name"/cleanasm/"$name".hap"
asmsuffix=".p_ctg.clean.fna"

#load environment
module load anaconda3/latest
export PATH=/home/wangj/miniconda3/envs/assembly/bin/:$PATH

#create working directories !DELETE FOLDERS FROM PREVIOUS ATTEMPTS!
mkdir $maindir$name"/juicer"

for i in 1 2 #haplotype 1 and 2
do

#create required folders for juicer
mkdir $maindir$name"/juicer/hap"$i
cd $maindir$name"/juicer/hap"$i
mkdir fastq
mkdir references
mkdir restriction_sites

#prepare input files
sname=${name}hap${i}
cd references
ln -s $asmprefix$i$asmsuffix ./${sname}.fasta
bwa index $sname.fasta
python /home/wangj/softwave/juicer-1.6/misc/generate_site_positions.py HindIII $sname ${sname}.fasta
cd $maindir$name"/juicer/hap"$i
mv references/${sname}_HindIII.txt restriction_sites/${sname}_HindIII.txt
awk 'BEGIN { OFS = "\t" } { print $1, $NF }' restriction_sites/${sname}_HindIII.txt > restriction_sites/${sname}.chrom.sizes
ln -s ${hic}1.fq.gz fastq/reads_R1.fastq.gz
ln -s ${hic}2.fq.gz fastq/reads_R2.fastq.gz
ln -s /home/wangj/softwave/juicer-1.6/SLURM/scripts/ scripts

echo $(date)" -- Running Juicer for sample "$name" ("$longname") Haplotype "$i"."

#run juicer with the exact parameters
/home/wangj/softwave/juicer-1.6/SLURM/scripts/juicer.sh -g ${sname} -z references/${sname}.fasta -p restriction_sites/${sname}.chrom.sizes -y restriction_sites/${sname}_HindIII.txt -s HindIII -t 50 -q tcum256c128Partition -l tcum256c128Partition -D . 

echo $(date)" -- Juicer has sumitted the jobs for sample "$name" ("$longname") Haplotype "$i"."

done


echo $(date)" -- Initial scaffolding jobs have started for sample "$name" ("$longname"). After they are done, inspect the results and move forward with 3D-DNA."
