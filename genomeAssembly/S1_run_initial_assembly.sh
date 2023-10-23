#!/bin/bash -l

#SBATCH --job-name=S1_initial_assembly
#SBATCH --partition=tcum256c128Partition
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --nodes=1

set -e

### change arguments here ###
maindir=/home/wangj/result/grape_assemblies/

metafile=/home/wangj/result/public/grape_file_paths/Grape_samples_and_datapaths.csv
mito=/home/wangj/result/public/Vv_chloroplast_mitochondrial/vv.mitochondrion.fa
chlo=/home/wangj/result/public/Vv_chloroplast_mitochondrial/vv.chloroplast.fa
contaminant=/home/wangj/result/public/contam_in_euks/contam_in_euks

#inputs the sample name etc from the metafile via SLURM_ARRAY_TASK_ID
input="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F "," ' $1 == var {print;} ' $metafile)"
name=`echo $input | awk -F "," '{ print $2}'`
longname=`echo $input | awk -F "," '{ print $3 , " cultivar:" , $4}'`
hifi=`echo $input | awk -F "," '{ print $6}'`
hicfiles=`echo $input | awk -F "," '{ print $7}'`
hic=${hicfiles::-7}

### HIFIASM ####

echo $(date)" -- Starting the hifiasm assembly for sample "$name" ("$longname")."

#load environment
export PATH=/home/wangj/miniconda3/envs/assembly/bin/:$PATH

#create working directories !DELETE FOLDERS FROM PREVIOUS ATTEMPTS!
mkdir $maindir$name
mkdir $maindir$name"/hifiasm"
cd $maindir$name"/hifiasm"

#run hifiasm with Hi-C reads [h1&h2] and 50 threads [t]
hifiasm -o ${name}.asm -t 50 --h1 ${hic}1.fq.gz --h2 ${hic}2.fq.gz $hifi

echo $(date)" -- Hifiasm assembly of "$name" has finished. Assembly statistics are being generated."

#convert GFA to FASTA (assembly-stats requires FASTA)
awk '/^S/{print ">"$2;print $3}' ${name}.asm.hic.hap1.p_ctg.gfa > ${name}.asm.hic.hap1.p_ctg.fna
awk '/^S/{print ">"$2;print $3}' ${name}.asm.hic.hap2.p_ctg.gfa > ${name}.asm.hic.hap2.p_ctg.fna

samtools faidx ${name}.asm.hic.hap1.p_ctg.fna
samtools faidx ${name}.asm.hic.hap2.p_ctg.fna

awk 'OFS="\t" {print $1,0,$2}' ${name}.asm.hic.hap1.p_ctg.fna.fai > ${name}.asm.hic.hap1.p_ctg.fna.bed
awk 'OFS="\t" {print $1,0,$2}' ${name}.asm.hic.hap2.p_ctg.fna.fai > ${name}.asm.hic.hap2.p_ctg.fna.bed

assembly-stats ${name}.asm.hic.hap1.p_ctg.fna >> ${maindir}${name}/${name}.asm.stats
assembly-stats ${name}.asm.hic.hap2.p_ctg.fna >> ${maindir}${name}/${name}.asm.stats

echo $(date)" -- Assembly statistics are written to "$maindir$name"/"$name".asm.stats for sample "$name"."

### RM NONNUCLEAR GENOME ###

echo $(date)" -- Starting the cleanup of hifiasm assemblies for sample "$name" ("$longname")."

asmprefix=$maindir$name"/hifiasm/"$name".asm.hic.hap"
asmsuffix=".p_ctg.fna"

#create working directories !DELETE FOLDERS FROM PREVIOUS ATTEMPTS!
mkdir $maindir$name"/cleanasm"
cd $maindir$name"/cleanasm"

for i in 1 2 #for haplotype 1 & 2
do

echo $(date)" -- Starting the nonnuclear genome sequences check for sample "$name" ("$longname") Haplotype "$i"."

#run minimap for long assembly to mitochondrial reference [x] with 50 threads [t]
minimap2 -x asm5 -t 50 $mito ${asmprefix}${i}${asmsuffix} > ${name}.hap${i}${asmsuffix::-3}mito.paf

#run minimap for long assembly to chloroplast reference [x] with 50 threads [t]
minimap2 -x asm5 -t 50 $chlo ${asmprefix}${i}${asmsuffix} > ${name}.hap${i}${asmsuffix::-3}chlo.paf

#identify the contigs
awk 'OFS="\t" {print $1,$3,$4}' ${name}.hap${i}${asmsuffix::-3}mito.paf > ${name}.hap${i}${asmsuffix::-3}mito.paf.bed
awk 'OFS="\t" {print $1,$3,$4}' ${name}.hap${i}${asmsuffix::-3}chlo.paf > ${name}.hap${i}${asmsuffix::-3}chlo.paf.bed

bedtools coverage -a ${asmprefix}${i}${asmsuffix}.bed -b ${name}.hap${i}${asmsuffix::-3}mito.paf.bed > ${name}.hap${i}${asmsuffix::-3}mito.paf.bed.cov
bedtools coverage -a ${asmprefix}${i}${asmsuffix}.bed -b ${name}.hap${i}${asmsuffix::-3}chlo.paf.bed > ${name}.hap${i}${asmsuffix::-3}chlo.paf.bed.cov

awk '{ if ($7>.5) print $1} ' ${name}.hap${i}${asmsuffix::-3}mito.paf.bed.cov > ${name}.hap${i}${asmsuffix::-3}mito.contigs.txt
awk '{ if ($7>.5) print $1} ' ${name}.hap${i}${asmsuffix::-3}chlo.paf.bed.cov > ${name}.hap${i}${asmsuffix::-3}chlo.contigs.txt

#remove them
cat ${name}.hap${i}${asmsuffix::-3}mito.contigs.txt ${name}.hap${i}${asmsuffix::-3}chlo.contigs.txt > ${name}.hap${i}${asmsuffix::-3}rm.contigs.txt

/home/wangj/softwave/seqkit grep -v -f ${name}.hap${i}${asmsuffix::-3}rm.contigs.txt ${asmprefix}${i}${asmsuffix} > ${name}.hap${i}${asmsuffix::-3}mchl.fna
rm ${name}.hap${i}${asmsuffix::-3}rm.contigs.txt

#create mito & chlo fasta
/home/wangj/softwave/seqkit grep -f ${name}.hap${i}${asmsuffix::-3}mito.contigs.txt ${asmprefix}${i}${asmsuffix} > ${name}.hap${i}${asmsuffix::-3}mito.fna
/home/wangj/softwave/seqkit grep -f ${name}.hap${i}${asmsuffix::-3}chlo.contigs.txt ${asmprefix}${i}${asmsuffix} > ${name}.hap${i}${asmsuffix::-3}chlo.fna

done

### RM CONTAMINATIONS ###

#load environment for BLAST
export PATH=/home/wangj/softwave/ncbi-blast-2.13.0+/bin/:$PATH

for i in 1 2 #for haplotype 1 & 2
do

echo $(date)" -- Starting the contamination check for sample "$name" ("$longname") Haplotype "$i"."

/home/wangj/softwave/ncbi-blast-2.13.0+/bin/blastn -task megablast -db $contaminant -query ${name}.hap${i}${asmsuffix::-3}mchl.fna \
    -outfmt 6 -evalue 1e-4 -out ${name}.hap${i}${asmsuffix::-3}mchl.blast.out

#identify the sequences with high similarity to non-plant sequences.    
awk '$3 >98 && $4 >=50  && $4 <100 {print $1"\t"$7"\t"$8}' ${name}.hap${i}${asmsuffix::-3}mchl.blast.out > ${name}.hap${i}${asmsuffix::-3}mchl.blast.outf
awk '$3 >94 && $4 >=100 && $4 <200 {print $1"\t"$7"\t"$8}' ${name}.hap${i}${asmsuffix::-3}mchl.blast.out >> ${name}.hap${i}${asmsuffix::-3}mchl.blast.outf
awk '$3 >90 && $4 >=200 {print $1"\t"$7"\t"$8}' ${name}.hap${i}${asmsuffix::-3}mchl.blast.out >> ${name}.hap${i}${asmsuffix::-3}mchl.blast.outf
bedtools sort -i ${name}.hap${i}${asmsuffix::-3}mchl.blast.outf | bedtools merge -i stdin > ${name}.hap${i}${asmsuffix::-3}mchl.blast.outf.contamination

#remove the contaminant sequences and store them seperately
awk '{print $1}' ${name}.hap${i}${asmsuffix::-3}mchl.blast.outf.contamination | sort | uniq > ${name}.hap${i}${asmsuffix::-3}mchl.blast.outf.contamination.txt
/home/wangj/softwave/seqkit grep -v -f ${name}.hap${i}${asmsuffix::-3}mchl.blast.outf.contamination.txt ${name}.hap${i}${asmsuffix::-3}mchl.fna > ${name}.hap${i}${asmsuffix::-3}clean.fna
/home/wangj/softwave/seqkit grep -f ${name}.hap${i}${asmsuffix::-3}mchl.blast.outf.contamination.txt ${name}.hap${i}${asmsuffix::-3}mchl.fna > ${name}.hap${i}${asmsuffix::-3}contamination.fna

done

assembly-stats ${name}.hap1${asmsuffix::-3}clean.fna >> ${maindir}${name}/${name}.asm.stats
assembly-stats ${name}.hap2${asmsuffix::-3}clean.fna >> ${maindir}${name}/${name}.asm.stats

echo $(date)" -- Cleanup is finished for sample "$name". Inspect and move forward with Juicer."
