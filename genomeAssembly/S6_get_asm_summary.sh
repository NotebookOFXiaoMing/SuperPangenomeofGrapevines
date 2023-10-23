#!/bin/bash -l
#SBATCH --job-name=S6_asmsummary
#SBATCH --partition=cuPartition
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1

#set -e

metafile=/home/wangj/result/public/grape_file_paths/Grape_samples_and_datapaths.csv

#inputs the sample name etc from the metafile via SLURM_ARRAY_TASK_ID
input="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F "," ' $1 == var {print;} ' $metafile)"
name=`echo $input | awk -F "," '{ print $2}'`

asmstat=${name}/${name}.asm.stats

## if you ran this before, comment out this line
#echo -ne "sample\tsize\tcontign\tlargest\tN50\tL50\tN90\tL90\tscafn\tgaps\tqv\tbusco\tlai\n" > assemblysummaries.tsv
##

for i in 1 2
do
qvfile=${name}/QC/qv/${name}.hap${i}.merqury.qv
laifile=${name}/QC/lai/${name}.hap${i}.p_ctg.clean_HiC.fasta.out.LAI
buscofile=${name}/QC/busco/${name}.hap${i}/short_summary.specific.eudicots_odb10.${name}.hap${i}.txt

cleanasmline=`grep -n "hap"$i".p_ctg.clean.fna" $asmstat | cut -d : -f 1`
awk -v a=`expr $cleanasmline + 1` 'FNR==a' $asmstat > temp.txt

gsize=`grep -P 'sum\s\=\s(\d+)' temp.txt -o | awk '{ print $3 }'`
contign=`grep -P 'n\s\=\s(\d+)' temp.txt -o | awk '{ print $3 }'`
contiglarge=`grep -P 'largest\s\=\s(\d+)' temp.txt -o | awk '{ print $3 }'`

awk -v a=`expr $cleanasmline + 2` 'FNR==a' $asmstat > temp.txt
N50=`grep -P 'N50\s\=\s(\d+)' temp.txt -o | awk '{ print $3 }'`
L50=`grep -P 'n\s\=\s(\d+)' temp.txt -o | awk '{ print $3 }'`

awk -v a=`expr $cleanasmline + 6` 'FNR==a' $asmstat > temp.txt
N90=`grep -P 'N90\s\=\s(\d+)' temp.txt -o | awk '{ print $3 }'`
L90=`grep -P 'n\s\=\s(\d+)' temp.txt -o | awk '{ print $3 }'`

scafline=`grep -n "hap"$i".p_ctg.clean_HiC.fasta" $asmstat | cut -d : -f 1`
awk -v a=`expr $scafline + 1` 'FNR==a' $asmstat > temp.txt
scafn=`grep -P 'n\s\=\s(\d+)' temp.txt -o | awk '{ print $3 }'`

awk -v a=`expr $scafline + 9` 'FNR==a' $asmstat > temp.txt
gaps=`grep -P 'Gaps\s\=\s(\d+)' temp.txt -o | awk '{ print $3 }'`

qv=`awk '{print $4}' $qvfile`

lai=`awk 'FNR==2 {print $7}' $laifile`

busco=`awk 'FNR==9' $buscofile | xargs`

echo -ne "${name}.hap${i}\t${gsize}\t${contign}\t${contiglarge}\t${N50}\t${L50}\t${N90}\t${L90}\t${scafn}\t${gaps}\t${qv}\t${busco}\t${lai}\n" >> assemblysummaries.tsv

done