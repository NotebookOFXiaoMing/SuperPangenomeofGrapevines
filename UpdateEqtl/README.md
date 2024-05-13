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

```
read_tsv("peerFactors.tsv") %>% 
  as.data.frame()-> peerFactors
rownames(peerFactors)<-paste0("PEER",1:20)
read_tsv("snp.PCA.eigenvec") %>% 
  select(-2,-3) %>% 
  column_to_rownames("SampleName") %>% 
  t() %>% 
  as.data.frame() %>% 
  bind_rows(peerFactors) %>% 
  rownames_to_column() %>% 
  write_tsv("covariates.tsv")
```

# SV

```
vcftools --gzvcf merged.vgCall.vcf.gz --keep ../../22.SNPeQTL/113samples.id --minDP 3 --maf 0.05 --max-missing 0.8 --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out merged.vgCall.filtered
## kept 51998 out of a possible 114672 Sites
~/my_data/myan/biotools/PopLDdecay-3.42/bin/PopLDdecay -InVCF merged.vgCall.filtered.recode.edited.vcf -OutStat LDdecay
perl ~/my_data/myan/biotools/PopLDdecay-3.42/bin/Plot_OnePop.pl -inFile LDdecay.stat.gz -output svFig


beagle gt=merged.vgCall.filtered.recode.vcf nthreads=24 out=merged.vgCall.filtered.impute
bgzip -d merged.vgCall.filtered.impute.vcf.gz
python convertVcfTo012Matrix.py merged.vgCall.filtered.impute.vcf merged.vgCall.filtered.impute.012.matrix SV

python editSVvcf.py merged.vgCall.filtered.impute.vcf merged.vgCall.filtered.impute.edited.vcf
~/my_data/myan/biotools/VCF2PCACluster-1.40/bin/VCF2PCACluster -InVCF merged.vgCall.filtered.impute.edited.vcf -OutPut sv.PCA

## merge PCA and PEER factors

read_tsv("peerFactors.tsv") %>% 
  as.data.frame()-> peerFactors
rownames(peerFactors)<-paste0("PEER",1:20)
read_tsv("sv.PCA.eigenvec") %>% 
  select(-2,-3) %>% 
  column_to_rownames("SampleName") %>% 
  t() %>% 
  as.data.frame() %>% 
  bind_rows(peerFactors) %>% 
  rownames_to_column() %>% 
  write_tsv("sv.covariates.tsv")

## adjust sample order
read_tsv("merged.vgCall.filtered.impute.012.matrix") %>% 
  select(c("snpid",read_lines("113samples.id"))) %>% 
  write_tsv("merged.vgCall.filtered.impute.sorted.012.matrix")

## SV eqtl

cat merged.vgCall.filtered.impute.vcf | grep -v "#" | awk '{print $1"_"$2"_SV\t"$1"\t"$2}' > svPos.txt
Rscript ../../22.SNPeQTL/eqtl/run_matrixEqtl.R merged.vgCall.filtered.impute.sorted.012.matrix ../../22.SNPeQTL/eqtl/exp.df.filter.quant.norm.tsv sv.covariates.tsv svPos.txt ../../22.SNPeQTL/eqtl/genePos.tsv
```

# Manhattan Plot

```
grep "contig" merged.vgCall.filtered.recode.vcf | awk 'gsub("##contig=<ID=","")' | awk 'gsub("length=","")' | awk 'gsub(">","")' > chr.len

library(tidyverse)
library(ggrastr)

path<-"D:/abc"

chr.len<-read_csv(paste0(path,"chr.len"),col_names =FALSE) %>% 
  mutate(X2=X2+5000000) %>% 
  arrange(X1) %>% 
  rename("CHR"="X1",
         "LEN"="X2")
chr.len
sv.gwas.results<-read_tsv(paste0(path,"SV.cisEQTL.output")) %>% 
  mutate(CHR=str_extract(SNP,pattern = "chr[0-9]+"),
         BP=str_extract(SNP,pattern = "_[0-9]+_") %>% 
           str_replace_all("_","") %>% 
           as.numeric(),
         P=FDR) %>% 
  select(CHR,BP,P) %>% 
  arrange(CHR)

snp.gwas.results<-read_tsv(paste0(path,"snp.cisEQTL.output")) %>% 
  mutate(CHR=str_extract(SNP,pattern = "chr[0-9]+"),
         BP=str_extract(SNP,pattern = "_[0-9]+_") %>% 
           str_replace_all("_","") %>% 
           as.numeric(),
         P=FDR) %>% 
  select(CHR,BP,P) %>% 
  arrange(CHR)

manhattanPlot<-function(){
  chr.len %>% pull(LEN) %>% cumsum() -> x1
  chr.len %>% pull(LEN) -> x2
  head(x1,-1)
  temp.df<-data.frame(chromo=chr.len %>% pull(CHR),
                      chr_len=c(0,head(x1,-1)))
  
  
  sv.gwas.results %>% 
    left_join(temp.df,by=c("CHR"="chromo")) %>% 
    mutate(new_pos=BP+chr_len) -> new.sv.gwas.results
  
  snp.gwas.results %>% 
    left_join(temp.df,by=c("CHR"="chromo")) %>% 
    mutate(new_pos=BP+chr_len) -> new.snp.gwas.results
  
  
  cols<-c("#bdbadb","#fdb363","#f98177","#80b3d4","gray")
  p<-ggplot()+
    geom_point_rast(data=new.snp.gwas.results,
                    aes(x=new_pos,y=-log10(P),color=CHR),show.legend = FALSE,
                    shape=16,color="gray",alpha=0.7)+
    geom_point_rast(data=new.sv.gwas.results,
                    aes(x=new_pos,y=-log10(P),color=CHR),show.legend = FALSE,
                    shape=17,color="red",alpha=0.7)+
    scale_x_continuous(breaks =c(0,head(x1,-1)) + x2/2 ,
                       labels =temp.df %>% pull(chromo))+
    theme_bw()+
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line())+
    #scale_color_manual(values = rep(cols,5))+
    scale_y_continuous(expand = expansion(mult = c(0,0)),
                       limits = c(0,40))+
    labs(x=NULL)
  return(p)
}

pdf("gwas_man01.pdf",width = 16,height = 8)
manhattanPlot()
dev.off()

```
