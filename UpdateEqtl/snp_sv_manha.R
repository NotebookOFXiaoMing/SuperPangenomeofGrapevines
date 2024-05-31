setwd("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/RNAseq/secondBatchReanalysis/SNP_SV_NG第一次审稿意见的结果/")
library(tidyverse)
library(ggrastr)
chr.len<-read_csv("manhaPlot/chr.len",col_names = FALSE) %>% 
  arrange(X1)

colnames(chr.len)<-c("CHR","LEN")

chr.len %>% pull(LEN) %>% cumsum() -> x1
chr.len %>% pull(LEN) -> x2
## remove last element of vectors
head(x1,-1)
temp.df<-data.frame(CHR=chr.len %>% pull(CHR),
                    chr_len=c(0,head(x1,-1)))

snp.eqtl<-read_tsv("manhaPlot/snp.eqtl.output.with.pop.stra") %>% 
  mutate(CHR=str_extract(SNP,pattern = "chr[0-9]+"),
         BP=str_extract(SNP,pattern = "_[0-9]+_") %>% 
           str_replace_all("_","") %>% 
           as.numeric(),
         P=FDR) %>% 
  select(CHR,BP,P)

snp.eqtl.dat<-snp.eqtl %>% 
  left_join(temp.df,by=c("CHR"="CHR")) %>% 
  mutate(new_pos=BP+chr_len)

sv.eqtl<-read_tsv("manhaPlot/sv.eqtl.output.with.pop.stra") %>% 
  mutate(CHR=str_extract(SNP,pattern = "chr[0-9]+"),
         BP=str_extract(SNP,pattern = "_[0-9]+") %>% 
           str_replace_all("_","") %>% 
           as.numeric(),
         P=FDR) %>% 
  select(CHR,BP,P)
sv.eqtl

sv.eqtl.dat<-sv.eqtl %>% 
  left_join(temp.df,by=c("CHR"="CHR")) %>% 
  mutate(new_pos=BP+chr_len)

sv.eqtl.dat<-sv.eqtl.dat %>% 
  mutate(new_y=-log10(P)+1,
         group="SV")

snp.eqtl.dat<-snp.eqtl.dat %>% 
  mutate(new_y=log10(P)-1,
         group="SNP")

chr.pos<-data.frame(x=c(0,head(x1,-1)) + x2/2,
                    y=0,
                    label=paste0("chr",str_pad(1:19,"left",width = 2,pad = "0")))
chr.pos
cols<-c("#999999","#e59f01","#56b4e8","#009f73")

pdf("manhaPlot/sv_snp_eqtl.pdf",width=12,height = 8)

bind_rows(sv.eqtl.dat,snp.eqtl.dat) %>% 
  ggplot(aes(x=new_pos,y=new_y))+
  geom_point_rast(aes(color=CHR,shape=group))+
  geom_vline(xintercept = 51712434,lty="dashed",color="grey")+
  scale_color_manual(values = rep(cols,5))+
  guides(color="none")+
  theme_bw(base_size = 15)+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.position = c(0.9,0.2),
        legend.background = element_rect(fill="transparent"),
        legend.title = element_blank(),
        axis.line.y = element_line(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())+
  geom_text(data=chr.pos,aes(x=x,y=y,label=label))+
  scale_y_continuous(breaks = c(-61,-41,-21,0,21),
                     labels = c(60,40,20,0,20))+
  labs(y="-log10(FDR)",x=NULL)+
  geom_hline(yintercept = -log10(5.90e-5),lty="dashed",
             color="grey")+
  geom_hline(yintercept = log10(1.45e-7),lty="dashed",
             color="grey")
dev.off()
