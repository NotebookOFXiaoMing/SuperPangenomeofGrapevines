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

~/my_data/myan/biotools/VCF2PCACluster-1.40/bin/VCF2PCACluster -InVCF merged.all.filtered.recode.vcf -OutPut snp.PCA
python convertVcfTo012Matrix.py merged.all.filtered.recode.vcf merged.all.filtered.recode.matrix
```

# Express data processing

```
library(tidyverse)

load("exp.df.wider.Rdata")
exp.df.wider

read_lines("113samples.id") -> sample.ids

exp.df.wider %>% 
  select(c("gene_id",sample.ids)) %>% 
  column_to_rownames("gene_id") -> exp.df

exp.df[rowSums(exp.df > 0.5) > 113*0.2,] -> exp.df.filter

exp.df.filter %>% 
  as.matrix() %>% 
  preprocessCore::normalize.quantiles() %>% 
  as.data.frame() -> exp.df.filter.quant.norm

colnames(exp.df.filter.quant.norm)<-colnames(exp.df.filter)
rownames(exp.df.filter.quant.norm)<-rownames(exp.df.filter)

exp.df.filter.quant.norm[1:6,1:6]

exp.df.filter.quant.norm %>% dim()


exp.df.filter.quant.norm %>% 
  rownames_to_column() %>% 
  filter(rowname=="VHPh1chr03G003100")

exp.df %>% 
  rownames_to_column() %>% 
  filter(rowname=="VHPh1chr03G003100")


exp.df.filter.quant.norm %>% 
  rownames_to_column() %>% 
  write_tsv("exp.df.filter.quant.norm.tsv")
```

# The top 20 hidden and confounding factors

```
conda activate R35
Rscript run_peer.R exp.df.filter.quant.norm.tsv peerFactors.tsv
```
