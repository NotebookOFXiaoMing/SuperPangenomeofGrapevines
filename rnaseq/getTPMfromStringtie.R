library(tidyverse)
file.list<-list.files("03.expression",
                      pattern = "*abund.tsv",
                      full.names = TRUE,
                      recursive = TRUE)
names(file.list)<-str_extract(file.list,pattern = "Unknown_BL435-01T[0-9]+")

df.list<-file.list%>%map(read_delim,delim="\t",show_col_types=FALSE,progress=FALSE)

## 9 is TPM 8 is FPKM 7 is coverage
df<-df.list[[names(df.list)[1]]] %>% 
  select(1,9)

colnames(df)<-c("gene_id",names(df.list)[1])

for (i in 2:length(names(df.list))){
  sample_name=names(df.list)[i]
  df.list[[sample_name]]%>%
    select(1,9) -> new.df
    colnames(new.df)<-c("gene_id",sample_name)
  
  df<-df %>% 
    left_join(new.df,by=c("gene_id"="gene_id"))
  
}

write_delim(df,file="TPM.txt",delim = "\t")
