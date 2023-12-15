### Pan gene family


run orthofinder
```
orthofinder -f 00.all.peps/ -M msa -S diamond -T fasttree -a 120 -t 120
```

two result files

https://figshare.com/s/f9b61becb80cd9fa5dfd


#### get single copy gene families 

```
python getWGDinput.py ../CoreFamilyIds.txt ../Orthogroups.GeneCount.tsv ../Orthogroups.tsv WGDinputCore.txt
python getWGDinput.py ../SoftCoreFamilyIds.txt ../Orthogroups.GeneCount.tsv ../Orthogroups.tsv WGDinputSoftCore.txt
python getWGDinput.py ../DispensableFamilyIds.txt ../Orthogroups.GeneCount.tsv ../Orthogroups.tsv WGDinputDispensable.txt

python pairwise_CDS_id.py WGDinputCore.txt Core_pairwise_geneid.txt
python pairwise_CDS_id.py WGDinputSoftCore.txt SoftCore_pairwise_geneid.txt
python pairwise_CDS_id.py WGDinputDispensable.txt Dispensable_pairwise_geneid.txt

### random select genes
python randomSelectGenePair.py Core_pairwise_geneid.txt coreGenePairs.txt 10000
python randomSelectGenePair.py SoftCore_pairwise_geneid.txt softcoreGenePairs.txt 10000
python randomSelectGenePair.py Dispensable_pairwise_geneid.txt dispensableGenePairs.txt 10000

### alignemnt
ParaAT.pl -h coreGenePairs.txt -a ../../all.peps.fa -n ../../allCDS01.fa -p proc.txt -o core.axt -m mafft -f axt -c 1 -g 1
ParaAT.pl -h softcoreGenePairs.txt -a ../../all.peps.fa -n ../../allCDS01.fa -p proc.txt -o softcore.axt -m mafft -f axt -c 1 -g 1
ParaAT.pl -h dispenGenePairs.txt -a ../../all.peps.fa -n ../../allCDS01.fa -p proc.txt -o dispen.axt -m mafft -f axt -c 1 -g 1

### calculate kaks
snakemake -s calculateKaKs.smk --configfiles config.yaml --config input_folder=core.axt output_folder=core.kaks --cores 52 -p -k
snakemake -s calculateKaKs.smk --configfiles config.yaml --config input_folder=softcore.axt output_folder=softcore.kaks --cores 52 -p -k
snakemake -s calculateKaKs.smk --configfiles config.yaml --config input_folder=dispen.axt output_folder=dispen.kaks --cores 52 -p -k

python mergeKaKs.py dispen.kaks dispenKaKs.value
python mergeKaKs.py softcore.kaks softcoreKaKs.value
python mergeKaKs.py core.kaks coreKaKs.value

### axt convert to fasta
snakemake -s convertAXT2Fasta.smk --configfiles config.yaml --config input_folder=core.axt output_folder=core.fa --cores 52 -p -k
snakemake -s convertAXT2Fasta.smk --configfiles config.yaml --config input_folder=softcore.axt output_folder=softcore.fa --cores 52 -p -k
snakemake -s convertAXT2Fasta.smk --configfiles config.yaml --config input_folder=dispen.axt output_folder=dispen.fa --cores 52 -p -k

### calculate nucleotide diversity

snakemake -s calculateNucDiv.smk --configfiles config.yaml --config input_folder=core.fa output_folder=core.nucdiv --cores 52 -pn
snakemake -s calculateNucDiv.smk --configfiles config.yaml --config input_folder=softcore.fa output_folder=softcore.nucdiv --cores 52 -p -k
snakemake -s calculateNucDiv.smk --configfiles config.yaml --config input_folder=dispen.fa output_folder=dispen.nucdiv --cores 52 -p -k

cat core.nucdiv/*.nucdiv > core_nucdiv.txt
cat softcore.nucdiv/*.nucdiv > softcore_nucdiv.txt
cat dispen.nucdiv/*.nucdiv > dispen_nucdiv.txt
```
