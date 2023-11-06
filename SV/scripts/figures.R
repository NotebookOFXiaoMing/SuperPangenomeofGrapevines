library(tidyverse)
library(ggridges)

dat<-read_delim(snakemake@input[[1]],
                col_names = FALSE,
                delim = "\t")



dat %>% 
  group_by(X1) %>% 
  summarise(counts=n()) %>% 
  arrange(counts) %>% 
  pull(X1) -> x.levels

dat %>% 
  group_by(X1,X2) %>% 
  summarise(counts=n()) %>% 
  ungroup() %>% 
  mutate(X2=factor(X2,levels = rev(c("INS","DEL","INV","TRANS")))) %>%
  mutate(X1=factor(X1,levels = x.levels)) %>%  
  ggplot(aes(x=X1,y=counts))+
  geom_bar(stat="identity",
           aes(fill=X2))+
  theme_bw(base_size = 20)+
  theme(panel.grid = element_blank(),
        legend.position = "top",
        legend.title = element_blank())+
  labs(x=NULL,y=NULL) -> p1

dat %>% 
  mutate(X2=factor(X2,levels = rev(c("INS","DEL","INV","TRANS")))) %>% 
  ggplot(aes(x=log10(X3),y=X2))+
  geom_density_ridges(aes(fill=X2),alpha=0.8)+
  theme_bw(base_size = 20)+
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle=30,hjust=1,vjust=1))+
  labs(x="The length of Structure variations (bp)",y=NULL)+
  scale_x_continuous(breaks = c(2,4,6),
                     labels = c(100,10000,"1000000")) -> p2

pdf(file = snakemake@output[[1]])
print(p1)
print(p2)
dev.off()
