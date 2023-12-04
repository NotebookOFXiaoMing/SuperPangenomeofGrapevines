library(tidyverse)
library(RIdeogram)
library(tidyquant)
help(package="RIdeogram")

palette_dark() %>% as.vector()




dat01<-read_delim("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/NLR/V127_CBM.candidate.NLR.group",
                  delim="\t",
                  col_names=FALSE)
dat01

dat02<-read_delim("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/NLR/05.NLR.distribution.on.chrom/geneloc.txt",
                  delim = "\t",
                  col_names = FALSE)
dat02

dat01 %>% 
  left_join(dat02,by=c("X1"="X4")) %>% 
  arrange(X1.y)

chr.len<-read_delim("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/NLR/05.NLR.distribution.on.chrom/CBM.hap1.final.fa.fai",
                    delim = "\t",
                    col_names = FALSE)

chr.len %>% 
  pull(X2) %>% 
  cumsum() -> x1


x2<-chr.len %>% 
  pull(X2)
x2


dat03<-data.frame(chromo=paste0("chr",str_pad(1:19,side = "left",width = 2,pad = 0)),
                  chr_len=c(0,x1[1:18]))

dat02
dat01 %>% 
  left_join(dat02,by=c("X1"="X4")) %>% 
  arrange(X1.y) %>%
  left_join(dat03,by=c("X1.y"="chromo")) -> dat04

dat04 %>% 
  pull(X2.x) %>% unique()


dat04 %>% 
  left_join(data.frame(x=dat04 %>% 
                         pull(X2.x) %>% unique(),
                       #color=c(palette_dark() %>% as.vector())[1:8]),
                       color=c("#d53e4f","#f46d43","#fdae61","#fee08b",
                               "#e6f598","#abdda4","#66c2a5","#3288bd")),
            by=c("X2.x"="x")) -> dat04

hap1.cen<-read_delim("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/NLR/05.NLR.distribution.on.chrom/GrapeT2T_hap1-centro.bed",
                     delim = "\t",col_names = FALSE)
hap1.cen

grape_karyotype<-data.frame(Chr=chr.len %>% pull(X1),
                            Start=0,
                            End=chr.len %>% pull(X2),
                            CE_start=hap1.cen %>% pull(X2),
                            CE_end=hap1.cen %>% pull(X3))
grape_karyotype

write_delim(grape_karyotype,
            file = "D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/NLR/05.NLR.distribution.on.chrom/grape_karyotype.txt",
            delim = "\t")

# read_delim("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/NLR/05.NLR.distribution.on.chrom/V127_CBM.hap1.gff",
#            delim = "\t",
#            col_names = FALSE) %>% 
#   mutate(X1=case_when(
#     str_length(X1) == 4 ~ str_replace(X1,"chr","chr0"),
#     TRUE ~ X1
#   )) %>% 
#   write_delim(file = "D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/NLR/05.NLR.distribution.on.chrom/V127_CBM.hap1.gff",
#               delim = "\t",
#               col_names = FALSE)

gene_density<-GFFex(input = "D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/NLR/05.NLR.distribution.on.chrom/V127_CBM.hap1.gff",
                    karyotype = "D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/NLR/05.NLR.distribution.on.chrom/grape_karyotype.txt",
                    feature = "gene",
                    window = 1000000)

gene_density

nlr.dis<-data.frame(Type=dat04 %>% pull(X2.x),
                    Shape="circle",
                    Chr=dat04 %>% pull(X1.y),
                    Start=dat04 %>% pull(X2.y),
                    End=dat04 %>% pull(X3),
                    
                    color=dat04 %>% pull(color) %>% 
                      str_replace("#",""))
nlr.dis %>% head()

ideogram(karyotype = grape_karyotype,
         overlaid = gene_density,
         colorset1 = c("gray","gray","gray"),
         label = nlr.dis,
         label_type = "marker",
         Lx=80,Ly=25)
convertSVG("chromosome.svg", device = "png",width = 10,height = 8,dpi = 300)
convertSVG("chromosome.svg", device = "pdf")

