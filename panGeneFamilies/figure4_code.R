library(tidyverse)
library(scatterpie)
library(ggpmisc)
library(patchwork)
library(ggtree)
library(ggrastr)
library(pegas)
library(cowplot)
library(ggpubr)
library(ggh4x)
library(readxl)
sum(c(3914,12442,33696,464))
sum(c(8223,9252,45988,1054))


dat<-read_delim("Orthogroups.GeneCount.tsv",
                delim = "\t",
                col_names = TRUE) %>% 
  select(-Total) %>% 
  column_to_rownames("Orthogroup")

dat %>% colnames()
dat %>% dim()

dat[1:6,1:6]

dat.backup<-dat
head(dat.backup)

for(i in colnames(dat.backup)){
  print(i)
  print(dat.backup %>% select(i) %>% filter(!!as.name(i) !=0) %>% nrow())
}


which(rowSums(dat)==0)

dat[dat>0] = 1
dat[1:6,1:6]


df<-tibble(samples=character(),
           sampleNum=numeric(),
           Pan=numeric(),
           Core=numeric())

total_genome <- 144
sim<-100
samples<-dat %>% colnames()
samples
for(i in 1:sim){
  if(i%%10==0){
    print(i)
  }
  sim_samples<-sample(samples,replace = F)
  sim_pav<-dat[,sim_samples]
  for(j in 1:ncol(dat)){
    subsim<-sim_pav[,1:j]
    
    if(j==1){
      Npan<-sum(subsim)
      Ncore<-Npan
    }else{
      sum<-rowSums(subsim)
      Ncore<-nrow(subsim[sum==j,])
      Npan<-nrow(subsim[sum>0,])
    }
    df<-add_row(df,samples=as.character(i),
                sampleNum=j,
                Pan=Npan,
                Core=Ncore)
  }
}

load("df.Rdata")

df %>% 
  select(-samples) %>% 
  pivot_longer(!sampleNum) %>% 
  mutate(sampleNum=factor(sampleNum,levels = 1:total_genome)) -> longer.df


p01<-ggplot(data=longer.df,aes(x=sampleNum,y=value))+
  geom_boxplot(aes(fill=name,color=name),
               outlier.alpha = 0,
               #width=1,
               linewidth=0.1)+
  theme_bw(base_size = 10)+
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.5),
        legend.title = element_blank())+
  scale_fill_manual(values = c("Core"="#d4101a",
                               "Pan"="#3082b3"))+
  scale_color_manual(values = c("Core"="#d4101a",
                               "Pan"="#3082b3"))+
  labs(x="Sample number",y="Family number")+
  scale_x_discrete(breaks = c(1,seq(9,135,9),144))
p01
pdf(file = "figure4/fig4_1.pdf",width = 6,height = 6)
print(p01)
dev.off()

dat %>% rowSums() %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  rename("Total"=".") %>%
  mutate(Class=case_when(
    Total == total_genome ~ "Core",
    Total >= round(total_genome*0.9) & Total < total_genome ~ "SoftCore",
    Total < round(total_genome*0.9) & Total >= 2 ~ "Dispensable",
    Total == 1 ~ "Private"
  )) %>%
  mutate(Class=factor(Class,levels = c("Core","SoftCore","Dispensable","Private")))-> freq.df

freq.df %>% head()
freq.df %>% pull(Class) %>% table() %>% as.data.frame()



##没有修改之前

# pie.df<-data.frame(x=1,y=1,region=1,
#                    Core=3897,
#                    SoftCore=12587,
#                    Dispensable=33619,
#                    Private=465)

p02<-ggplot(data=freq.df,aes(x=Total))+
  geom_bar(aes(fill=Class),width=0.5)+
  theme_bw(base_size = 10)+
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle=60,vjust=1,hjust=1))+
  scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
  #scale_x_continuous(breaks = c(1,seq(9,135,9),146))+
  labs(x="Frequeny",y="Family number")+
  scale_fill_manual(values = c("Core"="#DC0000FF",
                               "SoftCore"="#2e89be",
                               "Dispensable"="#f36c42",
                               "Private"="#179437"))+
  scale_x_continuous(breaks = c(1,seq(15,135,15),144))

p02

freq.df %>% 
  pull(Class) %>% table() %>% 
  as.data.frame() %>% 
  rename("Class"=".") %>% 
  pivot_wider(names_from = Class,
              values_from = Freq) %>%
  mutate(x=1,y=1,region=1) -> pie.df
pie.df

x<-c(8223,9252,45988,1054)
x/sum(x)
y<-paste0(round(x/sum(x),4)*100,"%")
core.label<-paste0("Core\n",x[1]," (",y[1],")")
core.label
softcore.label<-paste0("SoftCore\n",x[2]," (",y[2],")")
dispen.label<-paste0("Dispensable\n",x[3]," (",y[3],")")
private.label<-paste0("Private\n",x[4]," (",y[4],")")
p03<-ggplot()+
  geom_scatterpie(data=pie.df,
                  aes(x,y,group=region,r=1),
                  cols=c("Core","SoftCore","Dispensable","Private"),
                  color=NA)+
  theme_void()+
  scale_fill_manual(name=NULL,
                    values = c("Core"="#DC0000FF",
                               "SoftCore"="#2e89be",
                               "Dispensable"="#f36c42",
                               "Private"="#179437"))+
  theme(legend.position = "none",
        legend.key.size = unit(0.2,'cm'))+
  coord_equal()+
  guides(fill=guide_legend(ncol=1))+
  annotate(geom = "text",x=1.3,y=2.1,label=core.label,angle=-15,size=3)+
  annotate(geom = "text",x=1.5,y=1.2,label=softcore.label,angle=0,size=3)+
  annotate(geom = "text",x=0.5,y=1,label=dispen.label,angle=0,size=3)+
  annotate(geom = "text",x=0.6,y=2.15,label=private.label,angle=20,size=3)+
  annotate(geom = "segment",x=0.95,xend=0.7,y=2,yend=2.1)
p03


p04<-p02+
  geom_plot(data=tibble(x=10,y=9500,plot=list(p03)),
            aes(x=x,y=y,label=plot),
            vp.width=0.8,vp.height=0.8)
p04

pdf(file = "figure4/fig4_2.pdf",width = 6,height = 6)
print(p04)
dev.off()

tree<-read.tree("SpeciesTree_rooted.txt")
ggtree(tree,branch.length = "none")+
  geom_tiplab() -> tree.data

ggplot_build(tree.data)$data[[3]] %>% 
  select(y,label) %>% 
  arrange(y) %>% pull(label) -> y.level

## 按照不同组给热图排序
y.level

read_excel("Supptable_phylogeny_grouping.xlsx") %>% 
  select(1,4,5) %>% 
  mutate(hap1=factor(hap1),
         hap2=factor(hap2)) %>% 
  pivot_longer(!sample) %>% 
  mutate(sample_id=paste(sample,name,sep="_")) %>% 
  select(3,4) %>% 
  filter(str_sub(sample_id,1,4) !="V003")-> phylo_group
phylo_group

phylo_group %>% arrange(value) %>%  
  pull(sample_id) -> y.level

bind_cols(dat,freq.df) %>% 
  select(-Total) %>% 
  mutate(x=paste(Class,rowname,sep="_")) %>% 
  select(-Class,-rowname) %>% 
  #mutate(x1=1:nrow(.)) %>% 
  pivot_longer(!c(x)) %>% 
  mutate(X1=str_extract(x,pattern = "[A-z]+") %>% str_replace("_OG",""),
         X2=str_extract(x,pattern = "OG[0-9]+")) %>% 
  mutate(X1=factor(X1,levels = c("Core","SoftCore","Dispensable","Private")),
         X2=factor(X2,levels = freq.df$rowname)) %>% 
  arrange(X1,X2) %>% 
  mutate(x=factor(x,levels = x %>% unique())) %>% 
  mutate(name=factor(name,levels = y.level))-> df.heatmap

df.heatmap[1:5,1:5]

data.frame(x = df.heatmap %>% pull(x) %>% unique(),
           x1 = 1:length(df.heatmap %>% pull(x) %>% unique())) %>% 
  right_join(df.heatmap,by=c("x"="x")) -> df.heatmap

df.heatmap[1:5,1:5]
df.heatmap %>% pull(x1) %>% summary()

p05<-ggplot(data=df.heatmap,aes(x=x1,y=name))+
  geom_tile_rast(aes(fill=factor(value)),
                 color=NA,
                 show.legend=FALSE)+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "right",
        legend.key.size = unit(0.5,'cm'),
        legend.title = element_blank(),
        legend.text = element_text(angle=90),
        legend.spacing.y = unit(1,'cm'),
        legend.key.width = unit(0.2,'cm'),
        axis.text.y = element_blank(),
        panel.border = element_blank(),
        axis.title.y = element_text(size=15))+
  labs(x=NULL,y=NULL)+
  scale_fill_manual(values = c("0"="#80b1d3","1"="#fb8073"),
                    labels=c("0"="Absent","1"="Present"))+
  scale_x_continuous(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))+
  guides(fill=guide_legend(byrow = TRUE))

## p5_2垂直的颜色条
p5_2<-phylo_group %>% 
  as.data.frame() %>% 
  arrange(value) %>% 
  mutate(x=1,
         sample_id=factor(sample_id,levels = sample_id)) %>% 
  ggplot(aes(x=x,y=sample_id))+
  geom_tile(aes(fill=factor(value)),
            show.legend = FALSE)+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "right",
        legend.key.size = unit(0.5,'cm'),
        legend.title = element_blank(),
        legend.text = element_text(angle=90),
        legend.spacing.y = unit(1,'cm'),
        legend.key.width = unit(0.2,'cm'),
        panel.border = element_blank(),
        axis.title.y = element_text(size=15))+
  scale_x_continuous(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))+
  labs(x=NULL,y="144 Haplotype genomes")+
  scale_fill_manual(values = c("gray","#7cae00","#00bfc4","#c77cff"))

## p5_1水平的颜色条
p5_1<-df.heatmap %>% filter(name=="V127_CBM_hap2") %>% 
  ggplot(aes(x=x1,y=1))+
  geom_tile_rast(aes(fill=X1),show.legend=FALSE)+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "right",
        legend.key.size = unit(0.5,'cm'),
        legend.title = element_blank(),
        legend.text = element_text(angle=90),
        legend.spacing.y = unit(1,'cm'),
        legend.key.width = unit(0.2,'cm'),
        axis.text.y = element_blank(),
        panel.border = element_blank(),
        axis.title.y = element_text(size=15))+
  scale_x_continuous(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))+
  labs(x=NULL,y=NULL)+
  scale_fill_manual(values = c("gray","#7cae00","#00bfc4","#c77cff"))

pdf(file = "figure4/fig4_3.pdf",width = 8,height = 4)
#print(p05)
plot_spacer()+p5_1+p5_2+p05+
  plot_layout(ncol = 2,nrow = 2,widths = c(1,40),heights = c(1,20))
dev.off()


# ### p05_1是用来表示分组的颜色条
# data.frame(x=c(0.5,3914.5,3914.5+12442,3914.5+12442+33696),
#            xend=c(3914.5,3914.5+12442,3914.5+12442+33696,3914.5+12442+33696+464),
#            group=c("Core","SoftCore","Dispensable","Private")) %>%
#   ggplot()+
#   geom_rect(aes(xmin=x,xmax=xend,ymin=1,ymax=10,fill=group),
#             show.legend = FALSE)+
#   #theme_void()+
#   scale_fill_manual(values = c("Core"="#DC0000FF",
#                                "SoftCore"="#2e89be",
#                                "Dispensable"="#f36c42",
#                                "Private"="#179437"))+
#   theme_bw()+
#   theme(axis.text.x = element_blank(),
#         axis.text.y = element_blank(),
#         axis.ticks = element_blank(),
#         panel.grid = element_blank(),
#         legend.position = "bottom",
#         legend.title = element_blank(),
#         #axis.text.y = element_blank(),
#         panel.border = element_blank())+
#   scale_x_continuous(expand = c(0,0))+
#   scale_y_discrete(expand = c(0,0)) -> p05_1
# p05_1


dat.backup %>% 
  rownames_to_column() %>% 
  left_join(freq.df,by=c("rowname"="rowname")) %>% 
  select(-Total) %>% 
  pivot_longer(!c(rowname,Class)) %>% 
  group_by(Class,name) %>% 
  summarise(total_genes=sum(value)) %>% 
  ungroup() %>% 
  mutate(name=factor(name,levels = y.level),
         Class=factor(Class,levels = c("Private","Dispensable","SoftCore","Core"))) %>% 
  ggplot(aes(x=total_genes,y=name))+
  geom_bar(stat="identity",aes(fill=Class),
           width = 1)+
  scale_fill_manual(values = c("Core"="#DC0000FF",
                               "SoftCore"="#2e89be",
                               "Dispensable"="#f36c42",
                               "Private"="#179437"),
                    breaks = c("Core","SoftCore","Dispensable","Private"))+
  theme_bw()+
  scale_x_continuous(expand = c(0,0))+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        axis.line.x = element_line(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  labs(x=NULL,y=NULL) -> p06
pdf(file = "figure4/fig4_4.pdf",width = 6,height = 4)
print(p06)
dev.off()

## 计算不同基因家族在不同地域葡萄中的比例
dat %>% 
  rownames_to_column() %>% 
  left_join(freq.df,by=c("rowname"="rowname")) %>% 
  select(c("Class",phylo_group %>% filter(value==0) %>% pull(sample_id))) %>% 
  rowwise() %>% 
  mutate(total_family=sum(c_across(where(is.numeric)))) %>% 
  filter(total_family != 0) %>% 
  pull(Class) %>% 
  table() %>% 
  as.data.frame() %>% 
  rename("Class"=".") %>% 
  mutate(group="0") -> group.0

dat %>% 
  rownames_to_column() %>% 
  left_join(freq.df,by=c("rowname"="rowname")) %>% 
  select(c("Class",phylo_group %>% filter(value==1) %>% pull(sample_id))) %>% 
  rowwise() %>% 
  mutate(total_family=sum(c_across(where(is.numeric)))) %>% 
  filter(total_family != 0) %>% 
  pull(Class) %>% 
  table() %>% 
  as.data.frame() %>% 
  rename("Class"=".") %>% 
  mutate(group="1") -> group.1

dat %>% 
  rownames_to_column() %>% 
  left_join(freq.df,by=c("rowname"="rowname")) %>% 
  select(c("Class",phylo_group %>% filter(value==2) %>% pull(sample_id))) %>% 
  rowwise() %>% 
  mutate(total_family=sum(c_across(where(is.numeric)))) %>% 
  filter(total_family != 0) %>% 
  pull(Class) %>% 
  table() %>% 
  as.data.frame() %>% 
  rename("Class"=".") %>% 
  mutate(group="2") -> group.2

dat %>% 
  rownames_to_column() %>% 
  left_join(freq.df,by=c("rowname"="rowname")) %>% 
  select(c("Class",phylo_group %>% filter(value==3) %>% pull(sample_id))) %>% 
  rowwise() %>% 
  mutate(total_family=sum(c_across(where(is.numeric)))) %>% 
  filter(total_family != 0) %>% 
  pull(Class) %>% 
  table() %>% 
  as.data.frame() %>% 
  rename("Class"=".") %>% 
  mutate(group="3") -> group.3

#save(group.0,group.1,group.2,group.3,file = "group0123.Rdata")
bind_rows(group.0,group.1,group.2,group.3) %>% 
  ggplot(aes(x=group,y=Freq))+
  geom_bar(stat="identity",aes(fill=Class))+
  scale_x_discrete(labels=c("Muscadine","North American","East Asian","European"))-> p07.1
p07.1
bind_rows(group.0,group.1,group.2,group.3) %>% 
  ggplot(aes(x=group,y=Freq))+
  geom_bar(stat="identity",aes(fill=Class),
           position = "fill")+
  scale_x_discrete(labels=c("Muscadine","North American","East Asian","European"))+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values = c("Core"="#DC0000FF",
                               "SoftCore"="#2e89be",
                               "Dispensable"="#f36c42",
                               "Private"="#179437"))+
  theme_bw(base_size = 10)+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 60,hjust=1,vjust=1))+
  labs(x=NULL,y="Proportion")-> p07.2
p07.2
pdf(file = "figure4/fig4_5.pdf",width = 4,height = 4)
print(p07.2)
dev.off()


dat.list<-list()
for(i in 1:length(samples)){
  dat %>% 
    rownames_to_column() %>% 
    left_join(freq.df,by=c("rowname"="rowname")) %>% 
    select(-c("rowname","Total")) %>% 
    select(c(samples[i],"Class")) %>% 
    filter(!!as.name(samples[i])!=0) %>% 
    pull(Class) %>% table() %>% 
    as.data.frame() %>% 
    rename("Class"=".") %>% 
    mutate(sample_id=samples[i]) -> dat.list[[i]]
}

phylo_group %>% 
  pull(value) %>% 
  table()

phylo_group %>% 
  mutate(phylo_group=case_when(
    value == "0" ~ "Muscadine",
    value == "1" ~ "North American",
    value == "2" ~ "East Asian",
    value == "3" ~ "European"
  )) %>% 
  mutate(phylo_group=factor(phylo_group,levels = c("Muscadine","North American","East Asian","European"))) -> phylo_group


### softcore在不同类别中的数量
kruskal.test(Freq~phylo_group,data=dat.list  %>% 
               bind_rows() %>% 
               left_join(phylo_group,by=c("sample_id"="sample_id")) %>% 
               filter(Class=="SoftCore"))

PMCMRplus::kwAllPairsNemenyiTest(Freq~phylo_group,data=dat.list  %>% 
                                   bind_rows() %>% 
                                   left_join(phylo_group,by=c("sample_id"="sample_id")) %>% 
                                   filter(Class=="SoftCore"),dist="Chisquare")
dat.list  %>% 
  bind_rows() %>% 
  left_join(phylo_group,by=c("sample_id"="sample_id")) %>% 
  filter(Class=="SoftCore") %>% 
  pull(phylo_group) %>% table()
p07_1<-dat.list  %>% 
  bind_rows() %>% 
  left_join(phylo_group,by=c("sample_id"="sample_id")) %>% 
  filter(Class=="SoftCore") %>% 
  ggplot(aes(x=phylo_group,y=Freq))+
  stat_boxplot(geom = "errorbar",
               width=0.2,
               linewidth=0.1)+
  geom_boxplot(fill="#3182b3",linewidth=0.1,
               width=0.4,outlier.alpha = 0)+
  theme_bw(base_size = 10)+
  theme(panel.grid = element_blank())+
  labs(x=NULL,y="The number of SoftCore Family")+
  geom_text(data=data.frame(x=1:4,y=c(9000,9200,9200,9300),label=c("c","b","b","a")),
            aes(x=x,y=y,label=label),
            size=5)+
  ggtitle(expression(paste(italic(P)," = 1.11",''%*%'',"10",''^-'')^8))+
  scale_x_discrete(labels=c("Muscadine\n(n=8)","North American\n(n=43)",
                            "East Asian\n(n=50)","European\n(n=43)"))+
  theme(axis.text.x = element_text(angle = 60,hjust=1,vjust=1),
        plot.title = element_text(hjust=0.5))

## Dispensable在不同类别中的数量
kruskal.test(Freq~phylo_group,data=dat.list  %>% 
               bind_rows() %>% 
               left_join(phylo_group,by=c("sample_id"="sample_id")) %>% 
               filter(Class=="Dispensable"))

PMCMRplus::kwAllPairsNemenyiTest(Freq~phylo_group,data=dat.list  %>% 
                                   bind_rows() %>% 
                                   left_join(phylo_group,by=c("sample_id"="sample_id")) %>% 
                                   filter(Class=="Dispensable"),dist="Chisquare")

p07_2<-dat.list  %>% 
  bind_rows() %>% 
  left_join(phylo_group,by=c("sample_id"="sample_id")) %>% 
  filter(Class=="Dispensable") %>% 
  ggplot(aes(x=phylo_group,y=Freq))+
  stat_boxplot(geom = "errorbar",
               width=0.2,
               linewidth=0.1)+
  geom_boxplot(fill="#e66a44",linewidth=0.1,
               width=0.4,outlier.alpha = 0)+
  theme_bw(base_size = 10)+
  theme(panel.grid = element_blank())+
  labs(x=NULL,y="The number of Dispensable Family")+
  geom_text(data=data.frame(x=1:4,y=c(7500,9200,9100,8800),label=c("b","a","a","a")),
            aes(x=x,y=y,label=label),
            size=5)+
    ggtitle(expression(paste(italic(P)," = 1.97",''%*%'',"10",''^-'')^5))+
    scale_x_discrete(labels=c("Muscadine\n(n=8)","North American\n(n=43)",
                              "East Asian\n(n=50)","European\n(n=43)"))+
    theme(axis.text.x = element_text(angle = 60,hjust=1,vjust=1),
          plot.title = element_text(hjust=0.5))

### private在不同类别中的比例

kruskal.test(Freq~phylo_group,data=dat.list  %>% 
               bind_rows() %>% 
               left_join(phylo_group,by=c("sample_id"="sample_id")) %>% 
               filter(Class=="Private"))

PMCMRplus::kwAllPairsNemenyiTest(Freq~phylo_group,data=dat.list  %>% 
                                   bind_rows() %>% 
                                   left_join(phylo_group,by=c("sample_id"="sample_id")) %>% 
                                   filter(Class=="Private"),dist="Chisquare")

p07_3<-dat.list  %>% 
  bind_rows() %>% 
  left_join(phylo_group,by=c("sample_id"="sample_id")) %>% 
  filter(Class=="Private") %>% 
  ggplot(aes(x=phylo_group,y=Freq))+
  stat_boxplot(geom = "errorbar",
               width=0.2,
               linewidth=0.1)+
  geom_boxplot(fill="#238e3b",linewidth=0.1,
               width=0.4,outlier.alpha = 0)+
  theme_bw(base_size = 10)+
  theme(panel.grid = element_blank())+
  labs(x=NULL,y="The number of Private Family")+
  geom_text(data=data.frame(x=1:4,y=c(40,25,25,25),label=c("c","b","b","a")),
            aes(x=x,y=y,label=label),
            size=5)+
  ggtitle(expression(paste(italic(P)," = 5.35",''%*%'',"10",''^-'')^7))+
  scale_x_discrete(labels=c("Muscadine\n(n=8)","North American\n(n=43)",
                            "East Asian\n(n=50)","European\n(n=43)"))+
  theme(axis.text.x = element_text(angle = 60,hjust=1,vjust=1),
        plot.title = element_text(hjust=0.5))


pdf(file = "figure4/fig4_5_1.pdf",width = 12,height = 4)
p07_1+p07_2+p07_3+plot_layout(nrow = 1)
dev.off()
## 在cultivar 和wild中不同类别基因家族的差异

dat.list  %>% 
  bind_rows() %>% 
  left_join(read_tsv("grape_Cultivar_and_wild.txt",col_names = FALSE) %>% 
              mutate(X3="hap1",
                     X4="hap2") %>% 
              pivot_longer(!c(X1,X2)) %>% 
              mutate(X5=paste(X1,value,sep="_")) %>% 
              select(X5,X2),by=c("sample_id"="X5")) %>% 
  filter(Class=="Core") %>% 
  pull(X2) %>% table()

kruskal.test(Freq~X2,data=dat.list  %>% 
               bind_rows() %>% 
               left_join(read_tsv("grape_Cultivar_and_wild.txt",col_names = FALSE) %>% 
                           mutate(X3="hap1",
                                  X4="hap2") %>% 
                           pivot_longer(!c(X1,X2)) %>% 
                           mutate(X5=paste(X1,value,sep="_")) %>% 
                           select(X5,X2),by=c("sample_id"="X5")) %>% 
               filter(Class=="SoftCore"))
p07_4<-dat.list  %>% 
  bind_rows() %>% 
  left_join(read_tsv("grape_Cultivar_and_wild.txt",col_names = FALSE) %>% 
              mutate(X3="hap1",
                     X4="hap2") %>% 
              pivot_longer(!c(X1,X2)) %>% 
              mutate(X5=paste(X1,value,sep="_")) %>% 
              select(X5,X2),by=c("sample_id"="X5")) %>% 
  filter(Class=="SoftCore") %>% 
  ggplot(aes(x=X2,y=Freq))+
  stat_boxplot(geom = "errorbar",
               width=0.2,
               linewidth=0.1)+
  geom_boxplot(fill="#3182b3",linewidth=0.1,
               width=0.4,outlier.alpha = 0)+
  theme_bw(base_size = 10)+
  theme(panel.grid = element_blank())+
  ylim(8800,9200)+
  labs(x=NULL,y="The number of SoftCore Family")+
  ggtitle(expression(paste(italic(P)," = 1.78",''%*%'',"10",''^-'')^2))+
  scale_x_discrete(labels=c("Cultivar\n(n=94)","Wild\n(n=50)"))+
  theme(plot.title = element_text(hjust=0.5))

kruskal.test(Freq~X2,data=dat.list  %>% 
               bind_rows() %>% 
               left_join(read_tsv("grape_Cultivar_and_wild.txt",col_names = FALSE) %>% 
                           mutate(X3="hap1",
                                  X4="hap2") %>% 
                           pivot_longer(!c(X1,X2)) %>% 
                           mutate(X5=paste(X1,value,sep="_")) %>% 
                           select(X5,X2),by=c("sample_id"="X5")) %>% 
               filter(Class=="Dispensable"))

p07_5<-dat.list  %>% 
  bind_rows() %>% 
  left_join(read_tsv("grape_Cultivar_and_wild.txt",col_names = FALSE) %>% 
              mutate(X3="hap1",
                     X4="hap2") %>% 
              pivot_longer(!c(X1,X2)) %>% 
              mutate(X5=paste(X1,value,sep="_")) %>% 
              select(X5,X2),by=c("sample_id"="X5")) %>% 
  filter(Class=="Dispensable") %>% 
  ggplot(aes(x=X2,y=Freq))+
  stat_boxplot(geom = "errorbar",
               width=0.2,
               linewidth=0.1)+
  geom_boxplot(fill="#e66a44",linewidth=0.1,
               width=0.4,outlier.alpha = 0)+
  theme_bw(base_size = 10)+
  theme(panel.grid = element_blank())+
  labs(x=NULL,y="The number of Dispensable Family")+
  ggtitle(expression(paste(italic(P)," = 9.2",''%*%'',"10",''^-'')^1))+
  scale_x_discrete(labels=c("Cultivar\n(n=94)","Wild\n(n=50)"))+
  theme(plot.title = element_text(hjust=0.5))

kruskal.test(Freq~X2,data=dat.list  %>% 
               bind_rows() %>% 
               left_join(read_tsv("grape_Cultivar_and_wild.txt",col_names = FALSE) %>% 
                           mutate(X3="hap1",
                                  X4="hap2") %>% 
                           pivot_longer(!c(X1,X2)) %>% 
                           mutate(X5=paste(X1,value,sep="_")) %>% 
                           select(X5,X2),by=c("sample_id"="X5")) %>% 
               filter(Class=="Private"))
p07_6<-dat.list  %>% 
  bind_rows() %>% 
  left_join(read_tsv("grape_Cultivar_and_wild.txt",col_names = FALSE) %>% 
              mutate(X3="hap1",
                     X4="hap2") %>% 
              pivot_longer(!c(X1,X2)) %>% 
              mutate(X5=paste(X1,value,sep="_")) %>% 
              select(X5,X2),by=c("sample_id"="X5")) %>% 
  filter(Class=="Private") %>% 
  ggplot(aes(x=X2,y=Freq))+
  stat_boxplot(geom = "errorbar",
               width=0.2,
               linewidth=0.1)+
  geom_boxplot(fill="#e66a44",linewidth=0.1,
               width=0.4,outlier.alpha = 0)+
  theme_bw(base_size = 10)+
  theme(panel.grid = element_blank())+
  labs(x=NULL,y="The number of Private Family")+
  ggtitle(expression(paste(italic(P)," = 1.9",''%*%'',"10",''^-'')^1))+
  scale_x_discrete(labels=c("Cultivar\n(n=94)","Wild\n(n=50)"))+
  theme(plot.title = element_text(hjust=0.5))+
  ylim(0,25)

pdf(file = "figure4/fig4_5_2.pdf",width = 12,height = 4)
p07_4+p07_5+p07_6+plot_layout(nrow = 1)
dev.off()

freq.df %>% 
  filter(Class=="Core") %>% 
  pull(rowname) %>% 
  write_lines("CoreFamilyIds.txt")

freq.df %>% 
  filter(Class=="SoftCore") %>% 
  pull(rowname) %>% 
  write_lines("SoftCoreFamilyIds.txt")


freq.df %>% 
  filter(Class=="Dispensable") %>% 
  pull(rowname) %>% 
  write_lines("DispensableFamilyIds.txt")

freq.df %>% 
  filter(Class=="Private") %>% 
  pull(rowname) %>% 
  write_lines("PrivateFamilyIds.txt")

## 每个类别的基因家族基因的数量


read_tsv("CoreFamilyIds.txt",col_names = FALSE) %>% 
  left_join(read_tsv("Orthogroups.GeneCount.tsv") %>% 
              select(Orthogroup,Total),
            by=c("X1"="Orthogroup")) %>% 
  pull(Total) %>% sum()

1860296

read_tsv("SoftCoreFamilyIds.txt",col_names = FALSE) %>% 
  left_join(read_tsv("Orthogroups.GeneCount.tsv") %>% 
              select(Orthogroup,Total),
            by=c("X1"="Orthogroup")) %>% 
  pull(Total) %>% sum()

2206741

read_tsv("DispensableFamilyIds.txt",col_names = FALSE) %>% 
  left_join(read_tsv("Orthogroups.GeneCount.tsv") %>% 
              select(Orthogroup,Total),
            by=c("X1"="Orthogroup")) %>% 
  pull(Total) %>% sum()

1725948

read_tsv("PrivateFamilyIds.txt",col_names = FALSE) %>% 
  left_join(read_tsv("Orthogroups.GeneCount.tsv") %>% 
              select(Orthogroup,Total),
            by=c("X1"="Orthogroup")) %>% 
  pull(Total) %>% sum()

2861

(2861+1725948)/(2861+1725948+2206741+1860296)

## 注释到结构域的基因的比例

orthogroups<-read_delim("Orthogroups.tsv",
                        delim = "\t",
                        col_names = TRUE)
head(orthogroups)

identical(orthogroups$Orthogroup,freq.df$rowname)

head(freq.df)
bind_cols(orthogroups,freq.df) %>% 
  #mutate(X1=str_count(V127,"mRNA1")) %>% 
  #filter(X1>=2) %>% 
  select(-c(rowname,Total,Orthogroup)) %>% 
  pivot_longer(!Class) %>% 
  filter(Class=="Private") %>% 
  na.omit() %>% 
  pull(value) %>% 
  str_split(pattern = ", ") %>% 
  unlist() %>% 
  write_lines(file = "allPrivateID.txt")


bind_cols(orthogroups,freq.df) %>% 
  #mutate(X1=str_count(V127,"mRNA1")) %>% 
  #filter(X1>=2) %>% 
  select(-c(rowname,Total,Orthogroup)) %>% 
  pivot_longer(!Class) %>% 
  filter(Class=="Dispensable") %>% 
  na.omit() %>% 
  pull(value) %>% 
  str_split(pattern = ", ") %>% 
  unlist() %>% 
  write_lines(file = "allDispensableID.txt")


bind_cols(orthogroups,freq.df) %>% 
  #mutate(X1=str_count(V127,"mRNA1")) %>% 
  #filter(X1>=2) %>% 
  select(-c(rowname,Total,Orthogroup)) %>% 
  pivot_longer(!Class) %>% 
  filter(Class=="SoftCore") %>% 
  na.omit() %>% 
  pull(value) %>% 
  str_split(pattern = ", ") %>% 
  unlist() %>% 
  write_lines(file = "allSoftCoreID.txt")

bind_cols(orthogroups,freq.df) %>% 
  #mutate(X1=str_count(V127,"mRNA1")) %>% 
  #filter(X1>=2) %>% 
  select(-c(rowname,Total,Orthogroup)) %>% 
  pivot_longer(!Class) %>% 
  filter(Class=="Core") %>% 
  na.omit() %>% 
  pull(value) %>% 
  str_split(pattern = ", ") %>% 
  unlist() %>% 
  write_lines(file = "allCoreID.txt")


# myfun<-function(x){
#   return(read_delim(x,delim = "\t",comment = "##") %>% 
#            select(`#query`,PFAMs))
# }
# 
# list.files("06.emapper.collections/",pattern = "*.annotations$",full.names = TRUE)%>%
#   map(myfun)%>%
#   bind_rows() -> allemapperPfams
# 
# save(allemapperPfams,file = "allemapperPfams.Rdata")
load("allemapperPfams.Rdata")
dim(allemapperPfams)
head(allemapperPfams,n=4000) %>% DT::datatable()
tail(allemapperPfams)

read_delim("allCoreID.txt",
           delim = "\t",col_names = FALSE)%>%
  left_join(allemapperPfams,by=c("X1"="#query"))%>% 
  mutate(PFAMs=replace_na(PFAMs,"-")) %>% 
  mutate(group=case_when(
    PFAMs == "-" ~ "nodomain",
    TRUE ~ "domain"
  )) %>% pull(group) %>% table()
9182711/(9182711+1055866)

read_delim("allSoftCoreID.txt",
           delim = "\t",col_names = FALSE)%>%
  left_join(allemapperPfams,by=c("X1"="#query"))%>% 
  mutate(PFAMs=replace_na(PFAMs,"-")) %>% 
  mutate(group=case_when(
    PFAMs == "-" ~ "nodomain",
    TRUE ~ "domain"
  )) %>% pull(group) %>% table()

10459389/(10459389+1282457)

read_delim("allDispensableID.txt",
           delim = "\t",col_names = FALSE)%>%
  left_join(allemapperPfams,by=c("X1"="#query"))%>% 
  mutate(PFAMs=replace_na(PFAMs,"-")) %>% 
  mutate(group=case_when(
    PFAMs == "-" ~ "nodomain",
    TRUE ~ "domain"
  )) %>% pull(group) %>% table()

6828436/(6828436+1093913)

read_delim("allPrivateID.txt",
           delim = "\t",col_names = FALSE)%>%
  left_join(allemapperPfams,by=c("X1"="#query"))%>% 
  mutate(PFAMs=replace_na(PFAMs,"-")) %>% 
  mutate(group=case_when(
    PFAMs == "-" ~ "nodomain",
    TRUE ~ "domain"
  )) %>% pull(group) %>% table()
2810/(2810+1778)

p08<-data.frame(domain=c(9182711,10459389,6828436,2810),
           nodomain=c(1055866,1282457,1093913,1778),
           Class=c("Core","SoftCore","Dispensable","Private")) %>% 
  pivot_longer(!Class) %>% 
  mutate(Class=factor(Class,levels = c("Core","SoftCore","Dispensable","Private")),
         name=factor(name,levels = c("nodomain","domain"))) %>% 
  ggplot(aes(x=Class,y=value))+
  geom_bar(stat="identity",
           aes(fill=name),
           position="fill",show.legend = FALSE)+
  theme_bw(base_size = 10)+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        axis.text.x = element_text(angle = 60,hjust = 1,vjust=1))+
  scale_y_continuous(expand = expansion(mult = c(0,0)))+
  labs(x=NULL,y="Annotated domain")+
  scale_fill_manual(values = c("gray","#007AC1FF"))

pdf(file = "figure4/fig4_6.pdf",width=4,height = 4)
print(p08)
dev.off()


## KaKs value

read_tsv("coreKaKs.value",col_names = FALSE) %>% 
  select(3,4,5) %>% 
  mutate(Class="Core") -> core.kaks

read_tsv("softcoreKaKs.value",col_names = FALSE) %>% 
  select(3,4,5) %>% 
  mutate(Class="SoftCore") -> softcore.kaks

read_tsv("dispenKaKs.value",col_names = FALSE) %>% 
  select(3,4,5) %>% 
  mutate(Class="Dispensable") -> dispen.kaks

bind_rows(core.kaks,softcore.kaks,dispen.kaks) %>% 
  mutate(Class=factor(Class,levels = c("Core","SoftCore","Dispensable"))) %>% 
  filter(X5<1) %>% 
  pull(Class) %>% table()

kruskal.test(X5~Class,data=bind_rows(core.kaks,softcore.kaks,dispen.kaks) %>% 
               mutate(Class=factor(Class,levels = c("Core","SoftCore","Dispensable"))) %>% 
               filter(X5<1))

PMCMRplus::kwAllPairsNemenyiTest(X5~Class,data=bind_rows(core.kaks,softcore.kaks,dispen.kaks) %>% 
                                   mutate(Class=factor(Class,levels = c("Core","SoftCore","Dispensable"))) %>% 
                                   filter(X5<1))


p09<-bind_rows(core.kaks,softcore.kaks,dispen.kaks) %>% 
  mutate(Class=factor(Class,levels = c("Core","SoftCore","Dispensable"))) %>% 
  ggplot(aes(x=Class,y=X5))+
  geom_boxplot(linewidth=0.1,
               aes(fill=Class),
               width=0.4,
               show.legend = FALSE,
               outlier.alpha = 0)+
  scale_y_continuous(limits = c(0,1))+
  theme_bw(base_size = 10)+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle=60,hjust=1,vjust=1),
        plot.title = element_text(hjust=0.5))+
  labs(x=NULL,y="dn/ds")+
  scale_fill_manual(values = c("#cd131c","#2b83b1","#e66a44"))+
  geom_text(data=data.frame(x=c(1,2,3),y=c(0.57,0.61,0.65),
                            label=c("C","B","A")),
            aes(x=x,y=y,label=label))+
  ggtitle(expression(paste(italic(P)," < 2.2",''%*%'',"10",''^-'')^16))+
  scale_x_discrete(labels=c("Core\n(n=9,504)","SoftCore\n(n=9,388)","Dispensable\n(n=9,056)"))
p09
pdf(file = "figure4/fig4_7.pdf",width = 4,height = 4)
print(p09)
dev.off()

read_tsv("core_nucdiv.txt",col_names = FALSE) %>% 
  mutate(Class="Core") -> core.nucdiv


read_tsv("softcore_nucdiv.txt",col_names = FALSE) %>% 
  mutate(Class="SoftCore") -> softcore.nucdiv

read_tsv("dispen_nucdiv.txt",col_names = FALSE) %>% 
  mutate(Class="Dispensable") -> dispen.nucdiv


kruskal.test(X1~Class,data=bind_rows(core.nucdiv,softcore.nucdiv,dispen.nucdiv) %>% 
               mutate(Class=factor(Class,levels = c("Core","SoftCore","Dispensable"))))

PMCMRplus::kwAllPairsNemenyiTest(X1~Class,data=bind_rows(core.nucdiv,softcore.nucdiv,dispen.nucdiv) %>% 
                                   mutate(Class=factor(Class,levels = c("Core","SoftCore","Dispensable"))))

bind_rows(core.nucdiv,softcore.nucdiv,dispen.nucdiv) %>% 
  mutate(Class=factor(Class,levels = c("Core","SoftCore","Dispensable"))) %>% 
  pull(Class) %>% table()

bind_rows(core.nucdiv,softcore.nucdiv,dispen.nucdiv) %>% 
  mutate(Class=factor(Class,levels = c("Core","SoftCore","Dispensable"))) %>% 
  group_by(Class) %>% 
  summarise(mean_value=mean(X1))

p10<-bind_rows(core.nucdiv,softcore.nucdiv,dispen.nucdiv) %>% 
  mutate(Class=factor(Class,levels = c("Core","SoftCore","Dispensable"))) %>% 
  ggplot(aes(x=Class,y=X1))+
  geom_boxplot(linewidth=0.1,
               aes(fill=Class),
               width=0.4,
               show.legend = FALSE,
               outlier.alpha = 0)+
  scale_y_continuous(limits = c(0.4,0.8))+
  theme_bw(base_size = 10)+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle=60,hjust=1,vjust=1),
        plot.title = element_text(hjust=0.5))+
  labs(x=NULL,y="Nucleotide diversity")+
  scale_fill_manual(values = c("#cd131c","#2b83b1","#e66a44"))+
  geom_text(data=data.frame(x=c(1,2,3),y=c(0.7,0.71,0.72),
                            label=c("B","B","A")),
            aes(x=x,y=y,label=label))+
  ggtitle(expression(paste(italic(P)," = 2.6",''%*%'',"10",''^-'')^8))+
  scale_x_discrete(labels=c("Core\n(n=9,748)","SoftCore\n(n=9,709)","Dispensable\n(n=9,687)"))

p10
pdf(file = "figure4/fig4_8.pdf",width = 4,height = 4)
print(p10)
dev.off()
