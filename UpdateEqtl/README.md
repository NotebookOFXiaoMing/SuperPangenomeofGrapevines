# call snps

```
snakemake -s fastp.smk --cores 128 -p -k
bwa index -p refindex ref.fa
snakemake -s bwa.smk --cores 128 -p -k
snakemake -s picard.smk --cores 128 -p -k
snakemake -s bcftools.smk --cores 128 -p -n
snakemake -s bgzip.smk --cores 128 -p -k
vcftools --vcf 08.merged.vcf/merged.all.vcf --minDP 20 --maf 0.05 --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out merged.all.filtered
# After filtering, kept 11004014 out of a possible 71823356 Sites
```
