library(MatrixEQTL)
library(tidyverse)
library(readxl)
library(DESeq2)
library(patchwork)

SNP_file_name<-"D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/RNAseq/secondBatch/eQTL/SV_sb.txt"
snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

expression_file_name<-"D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/RNAseq/secondBatch/eQTL/TPM_exp_sb_trt.txt"

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);


### peer计算的协变量 ###
### 06.sv.gwas/20230803.reanalysis.add.VS01.sample/08.our.own.illumina/06.sv.eqtl.secondBatch/

# library(peer,lib.loc = "/yourpath/peer")
# expr<-read.table("TPM_exp_sb.txt",row.names=1,header=TRUE,sep="\t")
# model = PEER()
# PEER_setPhenoMean(model,as.matrix(t(expr)))
# PEER_setNk(model,20)
# PEER_getNk(model)
# PEER_update(model)
# factors = as.data.frame(t(PEER_getX(model)))
# write.table(factors,"./peer_covariates.tsv",quote=F,row.names=F,sep="\t",col.names=F)

# pca_cov<-read_delim("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/RNAseq/secondBatch/eQTL/pca_covariates.txt",
#                     delim = "\t")
# peer_cov<-read_delim("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/RNAseq/secondBatch/eQTL/peer_covariates.tsv",
#                      delim = "\t",col_names = FALSE)
# 
# dim(peer_cov)
# dim(pca_cov)
# 
# colnames(peer_cov)<-colnames(pca_cov)[-1]
# 
# bind_rows(pca_cov,peer_cov %>% 
#             mutate(id=paste0("peer1",1:20)) %>% 
#           select(id,colnames(pca_cov)[-1])) %>% 
#   write_delim(file = "D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/RNAseq/secondBatch/eQTL/covariates.tsv",
#               delim = "\t")



covariates_file_name<-"D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/RNAseq/secondBatch/eQTL/covariates.tsv"

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}

## Run the analysis

snps_location_file_name<-"D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/RNAseq/secondBatch/eQTL/SVsloc.txt"
gene_location_file_name<-"D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/RNAseq/secondBatch/eQTL/eQTLgeneloc.txt"

snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

genepos %>% head()

cisDist<-2000
cisDist<-100000
cisDist<-5000
pvOutputThreshold_cis<-0.1
pvOutputThreshold_tra<-0.00000000001

errorCovariance<-numeric()

output_file_name_cis = tempfile();
output_file_name_tra = tempfile();

useModel<-modelLINEAR

me = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name     = output_file_name_tra,
  pvOutputThreshold     = pvOutputThreshold_tra,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);



unlink(output_file_name_tra);
unlink(output_file_name_cis);

me$cis$eqtls %>% head()

me$cis$eqtls %>% 
  filter(-log10(pvalue)>=5) %>% 
  pull(gene) %>% unique() %>% 
  write_lines("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/RNAseq/secondBatch/eQTL_log10pvaluelargerthan5GenesId.txt")


me$cis$eqtls %>% 
  filter(-log10(FDR)>=5) %>% 
  pull(gene) %>% unique() %>% 
  write_lines("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/RNAseq/secondBatch/eQTL_log10FDRlargerthan5GenesId.txt")



### cis-QTL manhattan plot
###CBM1chr18G024990
me$cis$eqtls %>% 
  filter(gene=="CBM1chr18G024990")

chr.len<-read_delim("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/NLR/05.NLR.distribution.on.chrom/CBM.hap1.final.fa.fai",
                    delim = "\t",
                    col_names = FALSE)
chr.len %>% 
  mutate(X2=X2+10000000) -> chr.len
chr.len %>% pull(X2) %>% cumsum() -> x1
chr.len %>% pull(X2) -> x2
dat03<-data.frame(chromo=paste0("chr",str_pad(1:19,side = "left",width = 2,pad = 0)),
                  chr_len=c(0,x1[1:18]))

cols<-c("#bdbadb","#fdb363","#f98177","#80b3d4")

me$cis$eqtls %>% 
  mutate(chr=str_extract(snps,pattern = "chr[0-9]+"),
         pos=str_extract(snps,pattern = "[0-9]+$") %>% as.numeric()) %>% 
  left_join(dat03,by=c("chr"="chromo")) %>% 
  mutate(pos01=pos+chr_len) %>% 
  filter(-log10(FDR)>15)

me$cis$eqtls %>% 
  mutate(chr=str_extract(snps,pattern = "chr[0-9]+"),
         pos=str_extract(snps,pattern = "[0-9]+$") %>% as.numeric()) %>% 
  left_join(dat03,by=c("chr"="chromo")) %>% 
  mutate(pos01=pos+chr_len) %>% 
  ggplot(aes(x=pos01,y=-log10(FDR)))+
  geom_point(aes(color=chr),
             show.legend = FALSE,
             size=1.5)+
  scale_x_continuous(breaks =c(0,x1[1:18]) + x2/2 ,
                     labels =paste0("chr",str_pad(1:19,side = "left",width = 2,pad = 0)))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line())+
  scale_color_manual(values = rep(cols,5))+
  scale_y_continuous(expand = expansion(mult = c(0,0)),
                     limits = c(0,20),
                     breaks = seq(0,20,5))+
  labs(x=NULL)+
  annotate(geom = "point",x=71458983,y=-log10(1.54765e-17),
           size=10,shape=17,color="#f98177") -> p01

gene.structure<-read_delim("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/RNAseq/secondBatch/eQTL/CBM1chr3G003050.gff",
                           delim = "\t",
                           col_names = FALSE)
gene.structure %>% 
  filter(X3=="exon")

gene.structure %>% 
  filter(X3=="exon") %>% 
  ggplot()+
  annotate(geom="segment",x=2489000,xend=2503903,y=2,yend=2,
           #arrow=arrow(angle = 60,type="closed",length = unit(20,'mm')),
           linewidth=2,color="gray")+
  annotate(geom="segment",x=2489000,xend=2495000,y=2,yend=2,
           arrow=arrow(angle = 30,type="closed",length = unit(5,'mm')),
           linewidth=2,color="gray")+
  annotate(geom="segment",x=2489000,xend=2500000,y=2,yend=2,
           arrow=arrow(angle = 30,type="closed",length = unit(5,'mm')),
           linewidth=2,color="gray")+
  annotate(geom="segment",x=2489994,xend=2489994,y=2,yend=2.4,color="black",
           lty="dashed")+
  annotate(geom="rect",xmin=2489994-1000,xmax=2489994+1000,
           ymin=2.4,ymax=3.4,color="black",
           lty="dashed",fill=NA)+
  annotate(geom = "text",x=2489994,y=2.9,label="85bp Ins")+
  annotate(geom = "text",x=2498000,y=2.8,label="CBM1chr3G003050")+
  geom_rect(aes(xmin=X4,xmax=X5,ymin=1,ymax=3),
            fill="#f98177")+
  theme_void()+
  scale_y_continuous(limits = c(0,4)) ->p02

pdf(file = "Figure5.pdf",width = 8,height = 6)
p02+
  p01+
  plot_layout(ncol=1,heights = c(1,5))
dev.off()

###### 第二批处理和对照的表达量 CBM1chr3G003050
### gene expression TPM of ck
read_delim("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/RNAseq/secondBatch/eQTL/TPM_exp_sb_ck.txt",
           delim = "\t") -> gene.exp.ck

### V062样本转录组数据只有处理组，没有对照组
sv<-read_delim(SNP_file_name,
               delim = "\t")
gene.exp.trt<-read_delim(expression_file_name,
                     delim = "\t")
match(colnames(sv %>% select(-V062))[-1],colnames(gene.exp.ck)[-1])

identical(colnames(sv %>% select(-V062))[-1],colnames(gene.exp.ck)[-1])

gene.exp.ck %>% 
  filter(gene_id=="CBM1chr3G003050")  %>% 
  select(-1) %>% 
  t() %>% as.data.frame() %>% 
  rename("TPM"="V1") %>%
  rownames_to_column(var="sample_id") -> exp.ck.CBM1chr3G003050

exp.ck.CBM1chr3G003050

sv %>%
  select(-V062) %>% 
  filter(svid=="chr03_2489994") %>% 
  select(-1) %>% 
  t() %>% as.data.frame() %>% 
  rename("genotype"="V1") %>% 
  rownames_to_column(var="sample_id") -> sv.ck



gene.exp.trt %>% 
  filter(gene_id=="CBM1chr3G003050")  %>% 
  select(-1) %>% 
  t() %>% as.data.frame() %>% 
  rename("TPM"="V1") %>%
  rownames_to_column(var="sample_id") -> exp.trt.CBM1chr3G003050

exp.trt.CBM1chr3G003050


sv %>% 
  filter(svid=="chr03_2489994") %>% 
  select(-1) %>% 
  t() %>% as.data.frame() %>% 
  rename("genotype"="V1") %>% 
  rownames_to_column(var="sample_id") -> sv.trt


exp.ck.CBM1chr3G003050 %>% 
  left_join(sv.ck) %>% 
  pull(genotype) %>% 
  table()


exp.ck.CBM1chr3G003050 %>% 
  left_join(sv.ck) %>% 
  mutate(group="ck") %>% 
  bind_rows(exp.trt.CBM1chr3G003050 %>% 
              left_join(sv.trt) %>% 
              mutate(group="trt")) %>% 
  mutate(genotype=factor(genotype)) %>% 
  ggplot(aes(x=genotype,y=TPM))+
  geom_boxplot(aes(fill=group),
               show.legend = TRUE,
               width=0.4)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = c(0.2,0.8),
        legend.title = element_blank())+
  scale_fill_manual(values = c("#fb8073","#80b1d3")) -> p03

dwm<-read_excel("D:/Jupyter/IAAS/grape/metaInfo/Grapevine_GWAS_SR-20230723.xlsx",
                sheet = "List with data",
                na="NaN") %>% 
  select(SLNo,DWM2022,DWM2023) %>% 
  mutate(SLNo=paste0("V",str_pad(SLNo,side = "left",width = 3,pad="0")))

### 基因型和孢子数
dwm %>% head()

sv.trt %>% head()

sv.trt %>% 
  left_join(dwm,by=c("sample_id"="SLNo")) %>% 
  head()


sv.trt %>% 
  left_join(dwm,by=c("sample_id"="SLNo")) %>% 
  ggplot(aes(x=factor(genotype),y=log2(DWM2023+1)))+
  geom_boxplot(aes(fill=factor(genotype)),
               show.legend = FALSE,
               width=0.6)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_fill_manual(values = c("#ffb564","#f98075","#7eb1d2"))+
  labs(x="genotype") -> p04

p04
p03+p04
p03
p04

### 基因型和气孔开度

SC<-read_excel("D:/Jupyter/IAAS/grape/metaInfo/Grapevine_GWAS_SR-20230723.xlsx",
                sheet = "List with data",
                na="NaN") %>% 
  select(SLNo,SC10,SC4,SC22) %>% 
  mutate(SLNo=paste0("V",str_pad(SLNo,side = "left",width = 3,pad="0")))

SC

sv.trt %>% 
  left_join(SC,by=c("sample_id"="SLNo")) %>% 
  ggplot(aes(x=factor(genotype),y=SC10))+
  geom_boxplot(aes(fill=factor(genotype)),
               show.legend = FALSE,
               width=0.4)+
  labs(title="SC10")+
  theme_bw()+
  theme(panel.grid = element_blank()) -> sc10.plot

sv.trt %>% 
  left_join(SC,by=c("sample_id"="SLNo")) %>% 
  ggplot(aes(x=factor(genotype),y=SC4))+
  geom_boxplot(aes(fill=factor(genotype)),
               show.legend = FALSE,
               width=0.4)+
  labs(title="SC4")+
  theme_bw()+
  theme(panel.grid = element_blank()) -> sc4.plot

sv.trt %>% 
  left_join(SC,by=c("sample_id"="SLNo")) %>% 
  ggplot(aes(x=factor(genotype),y=SC22))+
  geom_boxplot(aes(fill=factor(genotype)),
               show.legend = FALSE,
               width=0.4)+
  labs(title="SC22")+
  theme_bw()+
  theme(panel.grid = element_blank()) -> sc22.plot

sc4.plot+sc10.plot+sc22.plot


### 处理组表达量和孢子数
head(dwm)
head(sv.trt)
exp.trt.CBM1chr3G003050 %>% 
  left_join(dwm,by=c("sample_id"="SLNo")) %>% 
  left_join(sv.trt,by=c("sample_id"="sample_id")) %>% 
  select(TPM,DWM2022,genotype) %>% 
  na.omit() %>% 
  ggplot(aes(x=TPM,y=log2(DWM2022+1)))+
  geom_point(aes(color=factor(genotype),
                 shape=factor(genotype)))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.2),
        legend.title = element_blank()) -> p05
p05
exp.trt.CBM1chr3G003050 %>% 
  left_join(dwm,by=c("sample_id"="SLNo")) %>% 
  left_join(sv.trt,by=c("sample_id"="sample_id")) %>% 
  select(TPM,DWM2023,genotype) %>% 
  na.omit() %>% 
  ggplot(aes(x=TPM,y=log2(DWM2023+1)))+
  geom_point(aes(color=factor(genotype),
                 shape=factor(genotype)))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.2),
        legend.title = element_blank()) -> p06

p03+
  p04+
  #p05+
  p06+
  plot_layout(ncol = 3,nrow = 1)



##### 以上是 CBM1chr3G003050 基因相关的内容

##### 以下是探索其他位点

### chr14 28158455 这个位点周围的4个基因互相有重叠，检查一下基因注释
### 这个位点基因型的比例

me$ci$eqtls %>% 
  filter(-log10(FDR)>10) %>% 
  filter(str_sub(snps,1,5)=="chr14")

read_delim("D:/Jupyter/IAAs/grape/20230726.reanalysis/scripts/RNAseq/secondBatch/chr14_28158455.vcf",
           delim = "\t") %>% 
  select(-c(1:9)) %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(V2=str_sub(V1,1,3)) %>% 
  pull(V2) %>% table()

#./. 0/0 0/1 1/0 1/1 
#11   1   8   2  98

### eQTL关联出来的位点大部分都是这种某个基因型的比例很高这种吗？

### 看下两个抗病相关基因

cbm.hap1.rga<-read_delim("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/RNAseq/secondBatch/eQTL/V127_CBM.candidate.NLR.id",
                         delim = "\t",
                         col_names = FALSE) %>% 
  pull(X1) %>% 
  str_replace(".mRNA1","")

cbm.hap1.rga %>% length()

cbm.hap1.rga
me$cis$eqtls %>% pull(gene)

me$cis$eqtls %>% filter(-log10(pvalue)>5) %>% 
  pull(gene) %>% 
  intersect(cbm.hap1.rga)

me$cis$eqtls %>% filter(gene=="CBM1chr19G017230")


genotypeProp<-function(x){
  return(
    read_delim(x,delim = "\t") %>% 
      select(-c(1:9)) %>% 
      t() %>% as.data.frame() %>% 
      mutate(V2=str_sub(V1,1,3)) %>% 
      pull(V2) %>% table()
  )
}


genotypeProp("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/RNAseq/secondBatch/chr19_28227366.vcf")
genotypeProp("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/RNAseq/secondBatch/chr19_28238081.vcf")


sv %>%
  filter(svid=="chr19_28238081")

source("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/RNAseq/secondBatch/myexpandpheno.R")

myexpandpheno(g_id = "CBM1chr19G017230",v_id = "chr19_28238081")  

myexpandpheno(g_id = "CBM1chr19G017230",v_id = "chr19_28227366")

me$cis$eqtls %>%
  filter(gene=="CBM1chr19G017230")

me$cis$eqtls %>%
  filter(gene=="CBM1chr19G017960")

me$cis$eqtls %>% 
  filter(-log10(FDR)>10) %>% 
  filter(str_sub(snps,1,5)=="chr14")

exp.trt.CBM1chr3G003050 %>% 
  left_join(dwm,by=c("sample_id"="SLNo")) %>%
  mutate(V1=case_when(
    TPM <= 0.1 ~ "0",
    TPM > 0.1 ~ ">0")) %>% 
  ggplot(aes(x=V1,y=DWM2023))+
  geom_boxplot()

exp.trt.CBM1chr3G003050 %>% 
  left_join(dwm,by=c("sample_id"="SLNo")) %>%
  pull(TPM) %>% summary()

exp.trt.CBM1chr3G003050 %>% 
  left_join(dwm,by=c("sample_id"="SLNo")) %>% 
  left_join(sv.trt,by=c("sample_id"="sample_id")) %>% 
  select(TPM,DWM2023,genotype,sample_id) %>% 
  na.omit() %>% 
  filter(TPM>20)
### 对照表达量和孢子数量
exp.ck.CBM1chr3G003050 %>% 
  left_join(sv.ck) %>% 
  left_join(dwm,by=c("sample_id"="SLNo")) %>% head()


exp.ck.CBM1chr3G003050 %>% 
  left_join(sv.ck) %>% 
  left_join(dwm,by=c("sample_id"="SLNo")) %>% 
  filter(DWM2023<5000) %>% 
  write_delim(file = "ck.txt",delim = "\t")

exp.ck.CBM1chr3G003050 %>% 
  left_join(sv.ck) %>% 
  left_join(dwm,by=c("sample_id"="SLNo")) %>%
  ggplot(aes(x=DWM2023,y=TPM))+
  geom_point()+
  xlim(0,4000)

### 处理表达量和孢子数

exp.trt.CBM1chr3G003050 %>% 
  left_join(sv.trt) %>%
  left_join(dwm,by=c("sample_id"="SLNo")) %>% 
  filter(DWM2023<5000) %>% 
  write_delim(file = "trt.txt",delim = "\t")

sv.trt

## 第一批的对照和处理 CBM1chr3G003050基因表达

fb.tpm.ck<-read_delim("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/RNAseq/firstBatch/TPM_exp_fb_ck.txt",
                      delim = "\t") %>% 
  filter(gene_id=="CBM1chr3G003050") %>% 
  select(-1) %>% 
  t() %>% as.data.frame() %>% 
  rename("TPM"="V1") %>%
  rownames_to_column(var="sample_id")

fb.tpm.trt<-read_delim("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/RNAseq/firstBatch/TPM_exp_fb_trt.txt",
                       delim = "\t") %>% 
  filter(gene_id=="CBM1chr3G003050") %>% 
  select(-1) %>% 
  t() %>% as.data.frame() %>% 
  rename("TPM"="V1") %>%
  rownames_to_column(var="sample_id")

sv %>%
  #select(-V062) %>% 
  filter(svid=="chr03_2489994") %>% 
  select(-1) %>% 
  t() %>% as.data.frame() %>% 
  rename("genotype"="V1") %>% 
  rownames_to_column(var="sample_id") -> tmp.b

fb.tpm.ck %>% 
  left_join(tmp.b,by=c("sample_id"="sample_id")) %>% 
  mutate(group="CK") %>% 
  bind_rows(fb.tpm.trt %>% 
              left_join(tmp.b,by=c("sample_id"="sample_id")) %>% 
              mutate(group="TRT")) %>% 
  mutate(genotype=factor(genotype)) %>%
  pull(genotype) %>% table()

fb.tpm.ck %>% 
  left_join(tmp.b,by=c("sample_id"="sample_id")) %>% 
  mutate(group="CK") %>% 
  bind_rows(fb.tpm.trt %>% 
              left_join(tmp.b,by=c("sample_id"="sample_id")) %>% 
              mutate(group="TRT")) %>% 
  mutate(genotype=factor(genotype)) %>% 
  ggplot(aes(x=genotype,y=TPM))+
  geom_boxplot(aes(fill=group))

dwm<-read_excel("D:/Jupyter/IAAS/grape/metaInfo/Grapevine_GWAS_SR-20230723.xlsx",
                sheet = "List with data",
                na="NaN") %>% 
  select(SLNo,DWM2022,DWM2023) %>% 
  mutate(SLNo=paste0("V",str_pad(SLNo,side = "left",width = 3,pad="0")))
dwm

tmp.b.trt %>% 
  left_join(dwm,by=c("sample_id"="SLNo")) %>% 
  mutate(genotype=factor(genotype)) %>% 
  ggplot(aes(x=genotype,y=log2(DWM2022)))+
  geom_boxplot()+
  ggtitle("DWM2022")

tmp.b.trt %>% 
  left_join(dwm,by=c("sample_id"="SLNo")) %>% 
  mutate(genotype=factor(genotype)) %>% 
  ggplot(aes(x=genotype,y=log2(DWM2023)))+
  geom_boxplot()+
  ggtitle("DWM2023")

tmp.b.trt %>% 
  left_join(dwm,by=c("sample_id"="SLNo")) %>% 
  mutate(genotype=factor(genotype)) %>% 
  pull(genotype) %>% table()

tmp.b.trt %>% 
  left_join(dwm,by=c("sample_id"="SLNo")) %>% 
  head()

fb.tpm.ck %>% 
  left_join(tmp.b,by=c("sample_id"="sample_id")) %>% 
  mutate(group="CK") %>% 
  bind_rows(fb.tpm.trt %>% 
              left_join(tmp.b,by=c("sample_id"="sample_id")) %>% 
              mutate(group="TRT")) %>% 
  head()




##CBM1chr19G017230

sv<-read_delim(SNP_file_name,
               delim = "\t")
sv
gene.exp<-read_delim(expression_file_name,
                     delim = "\t")
gene.exp %>% 
  filter(gene_id=="CBM1chr19G017960") %>% 
  select(-1) %>% 
  t() %>% as.data.frame() %>% 
  rename("TPM"="V1") %>%
  rownames_to_column(var="sample_id") -> tmp.a

sv %>% 
  filter(svid=="chr19_29465446") %>% 
  select(-1) %>% 
  t() %>% as.data.frame() %>% 
  rename("genotype"="V1") %>% 
  rownames_to_column(var="sample_id") -> tmp.b

tmp.a %>% 
  left_join(tmp.b) %>% 
  mutate(genotype=factor(genotype)) %>% 
  ggplot(aes(x=genotype,y=TPM))+
  geom_boxplot()

tmp.a %>% 
  left_join(tmp.b) %>% 
  pull(genotype) %>% 
  table()





##CBM1chr19G017230

sv<-read_delim(SNP_file_name,
               delim = "\t")
sv
gene.exp<-read_delim(expression_file_name,
                     delim = "\t")
gene.exp %>% 
  filter(gene_id=="CBM1chr19G017230") %>% 
  select(-1) %>% 
  t() %>% as.data.frame() %>% 
  rename("TPM"="V1") %>%
  rownames_to_column(var="sample_id") -> tmp.a

sv %>% 
  filter(svid=="chr19_28238081") %>% 
  select(-1) %>% 
  t() %>% as.data.frame() %>% 
  rename("genotype"="V1") %>% 
  rownames_to_column(var="sample_id") -> tmp.b

tmp.a %>% 
  left_join(tmp.b) %>% 
  mutate(genotype=factor(genotype)) %>% 
  ggplot(aes(x=genotype,y=TPM))+
  geom_boxplot()

tmp.a %>% 
  left_join(tmp.b) %>% 
  pull(genotype) %>% 
  table()

## most significant site

me$cis$eqtls %>% head()


me$cis$eqtls %>% arrange() %>% head()

sv<-read_delim(SNP_file_name,
               delim = "\t")
sv
gene.exp<-read_delim(expression_file_name,
                     delim = "\t")
gene.exp %>% 
  filter(gene_id=="CBM1chr3G003050") %>% 
  select(-1) %>% 
  t() %>% as.data.frame() %>% 
  rename("TPM"="V1") %>%
  rownames_to_column(var="sample_id") -> tmp.a

sv %>% 
  filter(svid=="chr03_2489994") %>% 
  select(-1) %>% 
  t() %>% as.data.frame() %>% 
  rename("genotype"="V1") %>% 
  rownames_to_column(var="sample_id") -> tmp.b

tmp.a %>% 
  left_join(tmp.b) %>% 
  mutate(genotype=factor(genotype)) %>% 
  ggplot(aes(x=genotype,y=TPM))+
  geom_boxplot()

tmp.a %>% 
  left_join(tmp.b) %>% 
  pull(genotype) %>% 
  table()


me$cis$eqtls %>% 
  select(snps,pvalue) %>% 
  mutate(chr=str_extract(snps,pattern = "chr[0-9]+"),
         pos=str_extract(snps,pattern = "[0-9]+$")) %>% 
  ggplot(aes(x=pos,y=-log10(pvalue)))+
  geom_point()+
  facet_wrap(~chr)+
  theme(axis.text.x = element_blank())

me$cis$eqtls %>% 
  filter(gene=="CBM1chr19G017230")


me$cis$eqtls %>% 
  filter(-log10(FDR)>20)












### 差异表达 V127

read_delim("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/RNAseq/secondBatch/sampleName_clientId.txt",
                    delim = "\t") %>% 
  filter(str_detect(clientId,"127")) %>% head(n=10)
  pull(sampleName) -> v127.sample.name


sb.counts<-read_delim("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/RNAseq/secondBatch/eQTL/gene_counts.csv",
                      delim = ",") %>% 
  select(1,v127.sample.name)
sb.counts %>% head()

colnames(sb.counts)<-c("gene_name","CK1","CK2","CK3","TRT1","TRT2","TRT3")

sb.counts %>% column_to_rownames("gene_name") -> sb.counts

sb.counts<-sb.counts[rowSums(sb.counts) != 0,]

sb.counts %>% dim()
sb.counts %>% 
  rownames_to_column() %>% 
  filter(rowname=="CBM1chr3G003050")

mymeta<-data.frame(id=c("CK1","CK2","CK3","TRT1","TRT2","TRT3"),
                   group=rep(c("CK","TRT"),each=3))
mymeta
colnames(sb.counts) == mymeta$id

library(DESeq2)
dds<-DESeqDataSetFromMatrix(countData = sb.counts,
                            colData = mymeta,
                            design = ~group)
dds<-DESeq(dds)
res<-results(dds)
head(res)
res %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  mutate(group=case_when(
    log2FoldChange >= 2 & padj <= 0.05 ~ "UP",
    log2FoldChange <= -2 & padj <= 0.05 ~ "DOWN",
    TRUE ~ "Not_Change"
  )) %>% 
  filter(rowname=="CBM1chr3G003050")
  
ggplot(aes(x=log2FoldChange,y=-log10(padj)))+
  geom_point(aes(color=group))+
  scale_color_manual(values = c("UP"="red",
                                "DOWN"="blue",
                                "Not_Change"="gray"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = "top")

tmp.b %>% 
  filter(genotype==2)

### 差异表达 V011 V104 V106

read_delim("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/RNAseq/secondBatch/sampleName_clientId.txt",
           delim = "\t") %>% 
  filter(str_detect(clientId,"104")) %>% 
pull(sampleName) -> v104.sample.name


sb.counts<-read_delim("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/RNAseq/secondBatch/eQTL/gene_counts.csv",
                      delim = ",") %>% 
  select(1,v104.sample.name)
sb.counts %>% head()

colnames(sb.counts)<-c("gene_name","CK1","CK2","CK3","TRT1","TRT2","TRT3")

sb.counts %>% column_to_rownames("gene_name") -> sb.counts

sb.counts<-sb.counts[rowSums(sb.counts) != 0,]

sb.counts %>% dim()
sb.counts %>% 
  rownames_to_column() %>% 
  filter(rowname=="CBM1chr3G003050")

mymeta<-data.frame(id=c("CK1","CK2","CK3","TRT1","TRT2","TRT3"),
                   group=rep(c("CK","TRT"),each=3))
mymeta
colnames(sb.counts) == mymeta$id

library(DESeq2)
dds<-DESeqDataSetFromMatrix(countData = sb.counts,
                            colData = mymeta,
                            design = ~group)
dds<-DESeq(dds)
res<-results(dds)
head(res)
res %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  mutate(group=case_when(
    log2FoldChange >= 2 & padj <= 0.05 ~ "UP",
    log2FoldChange <= -2 & padj <= 0.05 ~ "DOWN",
    TRUE ~ "Not_Change"
  )) %>% 
  filter(rowname=="CBM1chr3G003050")

###### V106
read_delim("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/RNAseq/secondBatch/sampleName_clientId.txt",
           delim = "\t") %>% 
  filter(str_detect(clientId,"106")) %>% head()
  pull(sampleName) -> v106.sample.name

v106.sample.name


sb.counts<-read_delim("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/RNAseq/secondBatch/eQTL/gene_counts.csv",
                      delim = ",") %>% 
  select(1,v106.sample.name)
sb.counts %>% head()

colnames(sb.counts)<-c("gene_name","CK1","CK2","CK3","TRT1","TRT2","TRT3")

sb.counts %>% column_to_rownames("gene_name") -> sb.counts

sb.counts<-sb.counts[rowSums(sb.counts) != 0,]

sb.counts %>% dim()
sb.counts %>% 
  rownames_to_column() %>% 
  filter(rowname=="CBM1chr3G003050")

mymeta<-data.frame(id=c("CK1","CK2","CK3","TRT1","TRT2","TRT3"),
                   group=rep(c("CK","TRT"),each=3))
mymeta
colnames(sb.counts) == mymeta$id


dds<-DESeqDataSetFromMatrix(countData = sb.counts,
                            colData = mymeta,
                            design = ~group)
dds<-DESeq(dds)
res<-results(dds)
head(res)
res %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  mutate(group=case_when(
    log2FoldChange >= 2 & padj <= 0.05 ~ "UP",
    log2FoldChange <= -2 & padj <= 0.05 ~ "DOWN",
    TRUE ~ "Not_Change"
  )) %>% 
  filter(rowname=="CBM1chr3G003050")


###### V011
read_delim("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/RNAseq/secondBatch/sampleName_clientId.txt",
           delim = "\t") %>% 
  filter(str_detect(clientId,"-11-")|str_detect(clientId,"-11ck-")) %>%
  pull(sampleName) -> v011.sample.name

v011.sample.name


sb.counts<-read_delim("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/RNAseq/secondBatch/eQTL/gene_counts.csv",
                      delim = ",") %>% 
  select(1,v011.sample.name)
sb.counts %>% head()

colnames(sb.counts)<-c("gene_name","CK1","CK2","CK3","TRT1","TRT2","TRT3")

sb.counts %>% column_to_rownames("gene_name") -> sb.counts

sb.counts<-sb.counts[rowSums(sb.counts) != 0,]

sb.counts %>% dim()
sb.counts %>% 
  rownames_to_column() %>% 
  filter(rowname=="CBM1chr3G003050")

mymeta<-data.frame(id=c("CK1","CK2","CK3","TRT1","TRT2","TRT3"),
                   group=rep(c("CK","TRT"),each=3))
mymeta
colnames(sb.counts) == mymeta$id


dds<-DESeqDataSetFromMatrix(countData = sb.counts,
                            colData = mymeta,
                            design = ~group)
dds<-DESeq(dds)
res<-results(dds)
head(res)
res %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  mutate(group=case_when(
    log2FoldChange >= 2 & padj <= 0.05 ~ "UP",
    log2FoldChange <= -2 & padj <= 0.05 ~ "DOWN",
    TRUE ~ "Not_Change"
  )) %>% 
  filter(rowname=="CBM1chr3G003050")



library(ggrastr)
pdf(file = "Rplot15.pdf",width = 10,height = 10)
me$trans$eqtls %>% 
  select(snps,pvalue) %>% 
  mutate(chr=str_extract(snps,pattern = "chr[0-9]+"),
         pos=str_extract(snps,pattern = "[0-9]+$")) %>% 
  filter(-log10(pvalue)>=10) %>% 
  ggplot(aes(x=pos,y=-log10(pvalue)))+
  geom_point()+
  facet_wrap(~chr)+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
dev.off()

me$trans$eqtls %>% head()
me$trans$eqtls %>% dim()

me$trans$eqtls %>% 
  select(snps,pvalue,gene) %>% 
  mutate(chr=str_extract(snps,pattern = "chr[0-9]+"),
         pos=str_extract(snps,pattern = "[0-9]+$")) %>% 
  filter(-log10(pvalue)>=50)
