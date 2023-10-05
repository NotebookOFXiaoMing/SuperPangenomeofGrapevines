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
#help(package="ggrastr")

dat<-read_delim("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/20230921.remove.V003/Orthogroups.GeneCount.tsv",
                delim = "\t",
                col_names = TRUE) %>% 
  select(-Total) %>% 
  column_to_rownames("Orthogroup")

dat %>% colnames()
dat %>% dim()

dat[1:6,1:6]

dat.backup<-dat
head(dat)

which(rowSums(dat)==0)

dat[dat>0] = 1

dat[1:6,1:6]

#dat[rowSums(dat)==0,]

df<-tibble(samples=character(),
           sampleNum=numeric(),
           Pan=numeric(),
           Core=numeric())

total_genome <- 144
sim<-100
samples<-dat %>% colnames()

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

#save(df,file = "D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/20230921.remove.V003/df.Rdata")
# write_delim(dat,file = "D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/20230815.use.all.hap/matrix01.txt",
#             delim = "",
#             col_names = FALSE)
load("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/20230921.remove.V003/df.Rdata")
df %>% 
  select(-samples) %>% 
  pivot_longer(!sampleNum) %>% 
  mutate(sampleNum=factor(sampleNum,levels = 1:total_genome)) -> longer.df

### p01是Pan Core箱线图
# p01<-ggplot(data=longer.df,aes(x=sampleNum,y=value))+
#   geom_boxplot(aes(fill=name),
#                outlier.alpha = 0,
#                #width=1,
#                linewidth=0.1,
#                color="gray")+
#   theme_bw(base_size = 20)+
#   theme(panel.grid = element_blank(),
#         legend.position = c(0.8,0.5),
#         legend.title = element_blank())+
#   scale_fill_manual(values = c("Core"="darkred",
#                                "Pan"="darkblue"))+
#   labs(x="Sample number",y="Family number")+
#   scale_x_discrete(breaks = c(1,seq(9,135,9),144))





p01<-ggplot(data=longer.df,aes(x=sampleNum,y=value))+
  # geom_boxplot(aes(fill=name),
  #              outlier.alpha = 0,
  #              width=1,
  #              linewidth=0.01,
  #              color="gray")+
  geom_jitter(aes(color=name),
              width = 0.2,
              alpha=0.8)+
  theme_bw(base_size = 20)+
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.5),
        legend.title = element_blank(),
        legend.text = element_text(size=20),
        axis.text.x = element_text(angle=60,vjust=1,hjust=1))+
  scale_color_manual(values = c("Core"="#DC0000FF",
                               "Pan"="#2e89be"))+
  guides(color = guide_legend(override.aes = list(size = 5)))+
  labs(x="Sample number",y="Family number")+
  scale_x_discrete(breaks = c(1,seq(15,135,15),144))

longer.df %>% 
  filter(name=="Pan") %>% 
  group_by(sampleNum) %>% 
  summarise(value=median(value)) %>% 
  pull(value) -> tmp.x

tmp.df<-data.frame(x=numeric(),
                   y=numeric())
for(i in 1:length(tmp.x)){
  if(i<length(tmp.x)){
    tmp.df<-add_row(tmp.df,x=i,y=tmp.x[i+1]-tmp.x[i])}
  
}

tmp.df %>% 
  ggplot(aes(x=x,y=y))+
  geom_point()+
  geom_line()

### n=125时差值开始出现个位数

dat[1:5,1:5]
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

## 每个类别里的基因id
freq.df %>% 
  select(-Total) %>% 
  left_join(read_delim("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/20230921.remove.V003/Orthogroups.tsv",
                       delim = "\t"),
            by=c("rowname"="Orthogroup")) %>% 
  rio::export("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/20230921.remove.V003/PanGeneFamily_id.xlsx")

freq.df %>% 
  select(-Total) %>% 
  left_join(read_delim("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/20230921.remove.V003/Orthogroups.GeneCount.tsv",
                       delim = "\t"),
            by=c("rowname"="Orthogroup")) %>% 
  rio::export("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/20230921.remove.V003/PanGeneFamily_genenumber.xlsx")

### NLR基因数量
freq.df %>% head()
freq.df %>% dim()
# dat %>% rowSums() %>% 
#   as.data.frame() %>% 
#   rownames_to_column() %>% 
#   rename("Total"=".") %>% 
#   mutate(Class=case_when(
#     Total == total_genome ~ "Core",
#     Total >= round(total_genome*0.9) & Total < total_genome ~ "SoftCore",
#     Total < round(total_genome*0.9) & Total >= 2 ~ "Dispensable",
#     Total == 1 ~ "Private"
#   )) %>% 
#   mutate(Class=factor(Class,levels = c("Core","SoftCore","Dispensable","Private")))-> freq.df
## p02是直方图
p02<-ggplot(data=freq.df,aes(x=Total))+
  geom_bar(aes(fill=Class),width=0.5)+
  theme_bw(base_size = 20)+
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
# pie.df<-data.frame(x=1,y=1,region=1,
#                    Core=3897,
#                    SoftCore=12587,
#                    Dispensable=33619,
#                    Private=465)
pie.df
x<-c(3914,12442,33696,464)
sum(x)

x/sum(x)

### p03是饼图
p03<-ggplot()+
  geom_scatterpie(data=pie.df,
                  aes(x,y,group=region,r=1),
                  cols=c("Core","SoftCore","Dispensable","Private"),
                  color=NA)+
  theme_void()+
  scale_fill_manual(labels=c("Core 3914 (7.75%)",
                             "SoftCore 12442 (24.63%)",
                             "Dispensable 33696 (66.70%)",
                             "Private 464 (0.92%)"),
                    name=NULL,
                    values = c("Core"="#DC0000FF",
                               "SoftCore"="#2e89be",
                               "Dispensable"="#f36c42",
                               "Private"="#179437"))+
  theme(legend.position = "none",
        legend.key.size = unit(0.2,'cm'))+
  coord_equal()+
  guides(fill=guide_legend(ncol=1))+
  annotate(geom = "text",x=1.3,y=2.1,label="Core\n3914 (7.75%)",angle=-15,size=3)+
  annotate(geom = "text",x=1.5,y=1.2,label="SoftCore\n12442 (24.63%)",angle=0,size=3)+
  annotate(geom = "text",x=0.5,y=1,label="Dispensable\n33696 (66.70%)",angle=0,size=3)+
  annotate(geom = "text",x=0.6,y=2.15,label="Private\n464 (0.92%)",angle=20,size=3)+
  annotate(geom = "segment",x=0.95,xend=0.7,y=2,yend=2.1)

p03

## p04是直方图组合饼图
p04<-p02+
  geom_plot(data=tibble(x=10,y=5500,plot=list(p03)),
            aes(x=x,y=y,label=plot),
            vp.width=0.8,vp.height=0.8)
p04
pdf(file = "D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/20230921.remove.V003/plot01_1.pdf",
    width=12,height = 5)
p01+p04+
  plot_layout(ncol=2)
dev.off()
identical(dat %>% rownames(),
          freq.df$rowname)

tree<-read.tree("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/20230921.remove.V003/SpeciesTree_rooted.txt")
ggtree(tree,branch.length = "none")+
  geom_tiplab() -> tree.data

ggplot_build(tree.data)$data[[3]] %>% 
  select(y,label) %>% 
  arrange(y) %>% pull(label) -> y.level

ggtree(tree,branch.length = "none")+
  scale_x_continuous(expand = c(0,0)) -> singlecopygenetree

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


### p05是基因家族存在缺失热图
p05<-ggplot(data=df.heatmap,aes(x=x1,y=name))+
  geom_tile_rast(aes(fill=factor(value)),
                 color=NA)+
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
        #axis.text.y = element_blank(),
        panel.border = element_blank(),
        axis.title.y = element_text(size=15))+
  labs(x=NULL,y="144 Haplotype genomes")+
  # geom_vline(xintercept = 3987.5,lty="dashed",color="gray")+
  # geom_vline(xintercept = 3987.5+12587,lty="dashed",color="gray")+
  # geom_vline(xintercept = 3987.5+12587+33619,lty="dashed",color="gray")+
  scale_fill_manual(values = c("0"="#80b1d3","1"="#fb8073"),
                    labels=c("0"="Absent","1"="Present"))+
  scale_x_continuous(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))+
  guides(fill=guide_legend(byrow = TRUE))
# annotate(geom = "rect",
#          xmin = 0,xmax=3987,ymin = 147,ymax=152,fill="#bdbadd")+
# annotate(geom = "text",label = "Core",x=3987/2,y=150,size=7)+
# annotate(geom = "rect",
#          xmin = 3987,xmax=3987+12587,ymin = 147,ymax=152,fill="#fdb363")+
# annotate(geom = "text",label = "SoftCore",x=3987+12587/2,y=150,size=7)+
# annotate(geom = "rect",
#          xmin = 3987+12587,xmax=3987+12587+33619,ymin = 147,ymax=152,fill="#f98278")+
# annotate(geom = "text",label = "Dispensable",x=3987+12587+33619/2,y=150,size=7)+
# annotate(geom = "rect",
#          xmin = 3987+12587+33619,xmax=3987+12587+33619+400,ymin = 147,ymax=152,fill="#80b3d5")

x<-c(3914,12442,33696,464)
### p05_1是用来表示分组的颜色条
data.frame(x=c(0.5,3914.5,3914.5+12442,3914.5+12442+33696),
           xend=c(3914.5,3914.5+12442,3914.5+12442+33696,3914.5+12442+33696+464),
           group=c("Core","SoftCore","Dispensable","Private")) %>%
  ggplot()+
  geom_rect(aes(xmin=x,xmax=xend,ymin=1,ymax=10,fill=group),
            show.legend = FALSE)+
  #theme_void()+
  scale_fill_manual(values = c("Core"="#DC0000FF",
                               "SoftCore"="#2e89be",
                               "Dispensable"="#f36c42",
                               "Private"="#179437"))+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        #axis.text.y = element_blank(),
        panel.border = element_blank())+
  scale_x_continuous(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0)) -> p05_1



# singlecopygenetree +
#   p05 +
#   plot_layout(ncol = 2,widths = c(1,20)) -> p05.1
# print(p05.1)


dat.backup[1:6,1:6]
freq.df %>% head()

### 不同类别基因家族额基因数量所占的比例

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
  group_by(name) %>% 
  summarise(total_value=sum(total_genes)) %>% 
  ungroup() %>% 
  left_join(dat.backup %>% 
              rownames_to_column() %>% 
              left_join(freq.df,by=c("rowname"="rowname")) %>% 
              select(-Total) %>% 
              pivot_longer(!c(rowname,Class)) %>% 
              group_by(Class,name) %>% 
              summarise(total_genes=sum(value)) %>% 
              ungroup() %>% 
              mutate(name=factor(name,levels = y.level),
                     Class=factor(Class,levels = c("Private","Dispensable","SoftCore","Core"))) %>% 
              filter(Class=="Dispensable"|Class=="Private") %>% 
              group_by(name) %>% 
              summarise(dp_value=sum(total_genes)) %>% 
              ungroup(),
            by=c("name"="name")) %>% 
  mutate(dp_prop=dp_value/total_value) %>% 
  pull(dp_prop) %>% mean()



### p06是基因数量的堆积柱形图
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


pdf(file = "D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/20230921.remove.V003/plot03_6.pdf",
    width = 12,height = 4)
p05_1+plot_spacer()+
  p05+p06+plot_layout(ncol = 2,heights = c(1,20),
                      widths = c(3,1))
dev.off()

### 看看不同类别基因家族在不同地域葡萄中的比例

read_excel("D:/Jupyter/IAAS/grape/metaInfo/Supptable_phylogeny_grouping.xlsx") %>% 
  select(1,4,5) %>% 
  mutate(hap1=factor(hap1),
         hap2=factor(hap2)) %>% 
  pivot_longer(!sample) %>% 
  mutate(sample_id=paste(sample,name,sep="_")) %>% 
  select(3,4) %>% 
  filter(str_sub(sample_id,1,4) !="V003")-> phylo_group

phylo_group
phylo_group %>% head()

phylo_group %>% filter(value==0) %>% pull(sample_id)

freq.df %>% head()

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


bind_rows(group.0,group.1,group.2,group.3) %>% 
  ggplot(aes(x=group,y=Freq))+
  geom_bar(stat="identity",aes(fill=Class))+
  scale_x_discrete(labels=c("Muscadine","North American","East Asian","European"))-> p07.1

bind_rows(group.0,group.1,group.2,group.3) %>% 
  ggplot(aes(x=group,y=Freq))+
  geom_bar(stat="identity",aes(fill=Class),
           position = "fill")+
  scale_x_discrete(labels=c("Muscadine","North American","East Asian","European"))+
  scale_fill_manual(values = c("Core"="#DC0000FF",
                               "SoftCore"="#2e89be",
                               "Dispensable"="#f36c42",
                               "Private"="#179437"))+
  theme_bw(base_size = 20)+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 60,hjust=1,vjust=1))+
  labs(x=NULL,y=NULL)-> p07.2

p07.2

# freq.df %>%
#   filter(Class=="Core") %>%
#   pull(rowname) %>%
#   write_lines(file = "D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/20230921.remove.V003/CoreFamilyId.txt")
# 
# freq.df %>%
#   filter(Class=="SoftCore") %>%
#   pull(rowname) %>%
#   write_lines(file = "D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/20230921.remove.V003/SoftCoreFamilyId.txt")
# 
# freq.df %>%
#   filter(Class=="Dispensable") %>%
#   pull(rowname) %>%
#   write_lines(file = "D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/20230921.remove.V003/DispensableFamilyId.txt")
# 
# freq.df %>%
#   filter(Class=="Private") %>%
#   pull(rowname) %>%
#   write_lines(file = "D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/20230921.remove.V003/PrivateFamilyId.txt")

read_delim("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/20230815.use.all.hap/Orthogroups.tsv",
           delim = "\t") %>% colnames()

### 挑选核心可变基因去做结构域注释，看不同的比例
orthogroups<-read_delim("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/20230921.remove.V003/Orthogroups.tsv",
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
  write_lines(file = "D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/20230921.remove.V003/allPrivatePepID.txt")


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
  write_lines(file = "D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/20230921.remove.V003/allDispensableID.txt")


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
  write_lines(file = "D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/20230921.remove.V003/allSoftCoreID.txt")

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
  write_lines(file = "D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/20230921.remove.V003/allCoreID.txt")


myfun<-function(x){
  return(read_delim(x,delim = "\t",comment = "##") %>% 
           select(`#query`,PFAMs))
}

list.files("06.emapper.collections/",pattern = "*.annotations",full.names = TRUE)%>%
  map(myfun)%>%
  bind_rows() -> allemapperPfams

save(allemapperPfams,file = "allemapperPfams.Rdata")
load("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/20230921.remove.V003/allemapperPfams.Rdata")
dim(allemapperPfams)
allemapperPfams %>%
  mutate(`#query`=str_replace(`#query`,".mRNA1","")) -> allemapperPfams
allemapperPfams %>% head()
## dim(allemapperPfams) 4700000

read_delim("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/20230921.remove.V003//allCoreID.txt",
           delim = "\t",col_names = FALSE)%>%
  left_join(allemapperPfams,by=c("X1"="#query"))%>% 
  mutate(PFAMs=replace_na(PFAMs,"-")) %>% 
  mutate(group=case_when(
    PFAMs == "-" ~ "nodomain",
    TRUE ~ "domain"
  )) %>% pull(group) %>% table()

## 965284    96545

965284/(965284+96545)

read_delim("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/20230921.remove.V003/allSoftCoreID.txt",
           delim = "\t",col_names = FALSE)%>%
  left_join(allemapperPfams,by=c("X1"="#query"))%>% 
  mutate(PFAMs=replace_na(PFAMs,"-")) %>% 
  mutate(group=case_when(
    PFAMs == "-" ~ "nodomain",
    TRUE ~ "domain"
  )) %>% pull(group) %>% table()

## 2331130   481700
2331130/(2331130+481700)

read_delim("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/20230921.remove.V003/allDispensableID.txt",
           delim = "\t",col_names = FALSE)%>%
  left_join(allemapperPfams,by=c("X1"="#query"))%>% 
  mutate(PFAMs=replace_na(PFAMs,"-")) %>% 
  mutate(group=case_when(
    PFAMs == "-" ~ "nodomain",
    TRUE ~ "domain"
  )) %>% pull(group) %>% table()

## 868463   763608

868463/(868463+763608)
read_delim("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/20230921.remove.V003/allPrivatePepID.txt",
           delim = "\t",col_names = FALSE)%>%
  left_join(allemapperPfams,by=c("X1"="#query"))%>% 
  mutate(PFAMs=replace_na(PFAMs,"-")) %>% 
  mutate(group=case_when(
    PFAMs == "-" ~ "nodomain",
    TRUE ~ "domain"
  )) %>% pull(group) %>% table()

#480      759
480/(480+759)


### p08是有结构域的基因的比例
data.frame(nodomain=c(96545,481700,763608,759),
           domain=c(965284,2331130,868463,480),
           Class=c("Core","SoftCore","Dispensable","Private")) %>% 
  pivot_longer(!Class) %>% 
  mutate(Class=factor(Class,levels = c("Core","SoftCore","Dispensable","Private")),
         name=factor(name,levels = c("nodomain","domain"))) %>% 
  ggplot(aes(x=Class,y=value))+
  geom_bar(stat="identity",
           aes(fill=name),
           position="fill",show.legend = FALSE)+
  theme_bw(base_size = 20)+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        axis.text.x = element_text(angle = 60,hjust = 1,vjust=1))+
  scale_y_continuous(expand = expansion(mult = c(0,0)))+
  labs(x=NULL,y="Annotated domain")+
  scale_fill_manual(values = c("gray","#007AC1FF")) -> p08

p08

###p09是核苷酸多样性的图
read_delim("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/20230921.remove.V003/nucleotide_diversity.txt",
           delim = "\t") %>% 
  mutate(x=factor(x,levels = c("Core","SoftCore","Dispensable","Private"))) %>% 
  ggplot(aes(x=x,y=log10(y)))+
  geom_boxplot(aes(fill=x),
               show.legend = FALSE,
               outlier.alpha = 0,
               width=0.5)+
  theme_bw(base_size = 20)+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle=60,hjust=1,vjust=1))+
  scale_fill_manual(values = c("Core"="#DC0000FF",
                               "SoftCore"="#2e89be",
                               "Dispensable"="#f36c42",
                               "Private"="#179437"))+
  labs(x=NULL,y="Nucleotide Diversity")+
  scale_y_continuous(labels = c(0.001,0.01,0.1),
                     breaks = c(-3,-2,-1)) -> p09
p09  
#ylim(0,1.5)+
# stat_compare_means(method = "wilcox.test",#label.y = 25,
#                    comparisons = list(c("Core","SoftCore"),
#                                       c("Dispensable","SoftCore"),
#                                       #c("Private","SoftCore"),
#                                       c("Private","Dispensable")),
#                    symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
#                                       symbols = c("****", "***", "**", "*", "ns")))
# scale_y_continuous(breaks = c(0,0.5,1),limits = c(0,1))-> p09

### p10是Ka/Ks的箱线图
#bind_rows(private.kaks,dispen.kaks,softcore.kaks,core.kaks) %>% 
load("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/20230921.remove.V003/kaks01.Rdata")
colnames(kaks01)<-c("Ka/Ks","Class")  
kaks01 %>% 
  filter(`Ka/Ks`<=2) %>% 
  mutate(Class=factor(Class,levels=c("Core","SoftCore","Dispensable","Private"))) %>% 
  ggplot(aes(x=Class,y=`Ka/Ks`))+
  geom_boxplot(aes(fill=Class),
               show.legend = FALSE,
               outlier.alpha = 0)+
  theme_bw(base_size = 20)+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle=60,hjust=1,vjust=1))+
  scale_fill_manual(values = c("Core"="#DC0000FF",
                               "SoftCore"="#2e89be",
                               "Dispensable"="#f36c42",
                               "Private"="#179437"))+
  labs(x=NULL,y="Ka/Ks")+
  scale_y_continuous(breaks = c(0,1,2),limits = c(0,2)) -> p10
  # stat_compare_means(method = "wilcox.test",#label.y = 25,
  #                    comparisons = list(c("Core","SoftCore"),
  #                                       c("SoftCore","Dispensable"),
  #                                       c("Private","Dispensable")),
  #                    symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
  #                                       symbols = c("****", "***", "**", "*", "ns"))) -> p10

p10

pdf(file = "D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/20230921.remove.V003/plot04.pdf",
    width = 8,height = 8)
p07.2+p08+p09+p10+
  plot_layout(ncol=2)
dev.off()
### 计算核苷酸多样性

save(orthogroups,freq.df,file = "D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/20230921.remove.V003/orthogroups_freq_df.Rdata")
bind_cols(orthogroups,freq.df) %>% 
  #mutate(X1=str_count(V127,"mRNA1")) %>% 
  #filter(X1>=2) %>% 
  select(-c(rowname,Total)) %>%
  #pivot_longer(!Class) %>% 
  filter(Class=="Core") %>% 
  select(-Class) %>% 
  pivot_longer(!Orthogroup) %>% 
  na.omit() -> df.core

df.core %>% head()
df.core %>% pull(Orthogroup) %>% unique() -> core.family.ids

for (i in 1:length(core.family.ids)){
  df.core %>% 
    filter(Orthogroup==core.family.ids[i]) %>% 
    pull(value) %>% 
    str_split(pattern = ", ") %>% 
    unlist() %>% 
    write_lines(paste0("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/20230921.remove.V003/geneIDS/Core/",
                       core.family.ids[i],".txt"))
    
}  



bind_cols(orthogroups,freq.df) %>% 
  #mutate(X1=str_count(V127,"mRNA1")) %>% 
  #filter(X1>=2) %>% 
  select(-c(rowname,Total)) %>%
  #pivot_longer(!Class) %>% 
  filter(Class=="SoftCore") %>% 
  select(-Class) %>% 
  pivot_longer(!Orthogroup) %>% 
  na.omit() -> df.softcore

df.softcore %>% pull(Orthogroup) %>% unique() -> softcore.family.ids

for (i in 1:length(softcore.family.ids)){
  df.softcore %>% 
    filter(Orthogroup==softcore.family.ids[i]) %>% 
    pull(value) %>% 
    str_split(pattern = ", ") %>% 
    unlist() %>% 
    write_lines(paste0("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/20230921.remove.V003/geneIDS/SoftCore/",
                       softcore.family.ids[i],".txt"))
  
}  


bind_cols(orthogroups,freq.df) %>% 
  #mutate(X1=str_count(V127,"mRNA1")) %>% 
  #filter(X1>=2) %>% 
  select(-c(rowname,Total)) %>%
  #pivot_longer(!Class) %>% 
  filter(Class=="Dispensable") %>% 
  select(-Class) %>% 
  pivot_longer(!Orthogroup) %>% 
  na.omit() -> df.dispensable

df.dispensable %>% pull(Orthogroup) %>% unique() -> dispensable.family.ids

for (i in 1:length(dispensable.family.ids)){
  df.dispensable %>% 
    filter(Orthogroup==dispensable.family.ids[i]) %>% 
    pull(value) %>% 
    str_split(pattern = ", ") %>% 
    unlist() %>% 
    write_lines(paste0("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/20230921.remove.V003/geneIDS/Dispensable/",
                       dispensable.family.ids[i],".txt"))
  
} 


bind_cols(orthogroups,freq.df) %>% 
  #mutate(X1=str_count(V127,"mRNA1")) %>% 
  #filter(X1>=2) %>% 
  select(-c(rowname,Total)) %>%
  #pivot_longer(!Class) %>% 
  filter(Class=="Private") %>% 
  select(-Class) %>% 
  pivot_longer(!Orthogroup) %>% 
  na.omit() -> df.private

df.private %>% pull(Orthogroup) %>% unique() -> private.family.ids

for (i in 1:length(private.family.ids)){
  df.private %>% 
    filter(Orthogroup==private.family.ids[i]) %>% 
    pull(value) %>% 
    str_split(pattern = ", ") %>% 
    unlist() %>% 
    write_lines(paste0("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/20230921.remove.V003/geneIDS/Private/",
                       private.family.ids[i],".txt"))
  
} 


myfun<-function(x){
  pai_list<-c()
  dna<-read.dna(x,format = "fasta")
  nseq<-length(labels(dna))
  
  if(nseq<2){
    return(pai_list)
  }else{
    for(i in 1:ncol(combn(1:nseq,2)))
    {
      j1<-combn(1:nseq,2)[1,i]
      j2<-combn(1:nseq,2)[2,i]
      pai<-nuc.div(dna[c(j1,j2),])
      pai_list<-append(pai_list,pai)
    }
    return(pai_list)
  }
  
}



### /home/wangj/my_data/myan/grape/06.sv.gwas/09.panGeneFamily/07.kaksPi
list.files("core_fasta/",pattern = "*.fasta",full.names = TRUE)%>%
  map(myfun)%>%unlist()->core.pai

list.files("softcore_fasta/",pattern = "*.fasta",full.names = TRUE)%>%
  map(myfun)%>%unlist()->softcore.pai

list.files("dispensable_fasta/",pattern = "*.fasta",full.names = TRUE)%>%
  map(myfun)%>%
  unlist()->dispensable.pai

list.files("private_fasta/",pattern = "*.fasta",full.names = TRUE)%>%
  map(myfun)%>%
  unlist() -> private.pai

bind_rows(data.frame(y=core.pai,x="Core"),
          data.frame(y=softcore.pai,x="SoftCore"),
          data.frame(y=dispensable.pai,x="Dispensable"),
          data.frame(y=private.pai,x="Private")) -> panGeneFamily.nucleotide.diversity

save(panGeneFamily.nucleotide.diversity,file="panGenefamilyNucleotideDiversity.Rdata")

#load("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/20230815.use.all.hap/panGenefamilyNucleotideDiversity.Rdata")






p01+p04+plot_layout(nrow = 1) -> p.row1
p05+p06+plot_layout(nrow = 1,widths = c(2,1)) -> p.row2
p07+p08+p09+plot_layout(nrow = 1)
pdf(file="Figure4.1.pdf",width = 12,height = 14)
plot_grid(p.row1,p.row2,p.row3,
          nrow = 3,ncol = 1,
          rel_heights = c(1,2,1))
dev.off()

read_delim("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/20230815.use.all.hap/kaks/private.kaks.txt",
           delim = "\t") %>% 
  select(`Ka/Ks`) %>% 
  mutate(Class="Private") %>% 
  na.omit()  -> private.kaks

private.kaks

read_delim("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/20230815.use.all.hap/kaks/dispen.kaks.txt",
           delim = "\t") %>% 
  select(`Ka/Ks`) %>% 
  mutate(Class="Dispensable") %>% 
  na.omit()  -> dispen.kaks

dispen.kaks

read_delim("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/20230815.use.all.hap/kaks/softcore.kaks.txt",
           delim = "\t") %>% 
  select(`Ka/Ks`) %>% 
  mutate(Class="SoftCore") %>% 
  na.omit()  -> softcore.kaks

softcore.kaks


read_delim("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/20230815.use.all.hap/kaks/core.kaks.txt",
           delim = "\t") %>% 
  select(`Ka/Ks`) %>% 
  mutate(Class="Core") %>% 
  na.omit()  -> core.kaks

private.kaks

### 不同类别的基因家族的表达，先用CBM的基因试试

ck.exp<-read_delim("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/RNAseq/secondBatchReanalysis/ckTPM.txt",
                   delim = "\t") %>% 
  pivot_longer(!gene_id)
ck.exp
read_delim("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/20230921.remove.V003/allCoreID.txt",
           col_names = FALSE,delim = "\t") %>% 
  filter(str_detect(X1,"CBM1")) %>% 
  mutate(X2="Core") -> df.core

read_delim("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/20230921.remove.V003/allSoftCoreID.txt",
           col_names = FALSE,delim = "\t") %>% 
  filter(str_detect(X1,"CBM1")) %>% 
  mutate(X2="SoftCore") -> df.softcore

read_delim("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/20230921.remove.V003/allDispensableID.txt",
           col_names = FALSE,delim = "\t") %>% 
  filter(str_detect(X1,"CBM1")) %>% 
  mutate(X2="Dispensable") -> df.dispensable

read_delim("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/20230921.remove.V003/allPrivatePepID.txt",
           col_names = FALSE,delim = "\t") %>% 
  filter(str_detect(X1,"CBM1")) %>% 
  mutate(X2="Private") -> df.private

bind_rows(df.core,df.softcore,df.dispensable,df.private) %>% 
  right_join(ck.exp,by=c("X1"="gene_id")) %>%
  na.omit() %>% 
  mutate(X2=factor(X2,levels = c("Core","SoftCore","Dispensable","Private"))) -> CoreSoftCoreDispenPrivate.exp
  
save(CoreSoftCoreDispenPrivate.exp,
     file = "grape/CoreSoftCoreDispenPrivate.Rdata")
  


##CBM hap1 抗病基因在染色体上的分布

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
                       color=c(palette_dark() %>% as.vector())[1:8]),
            by=c("X2.x"="x")) -> dat04


grape_karyotype<-data.frame(Chr=chr.len %>% pull(X1),
                            Start=0,
                            End=chr.len %>% pull(X2))
grape_karyotype



nlr.dis<-data.frame(Type=dat04 %>% pull(X2.x),
                    Shape="triangle",
                    Chr=dat04 %>% pull(X1.y),
                    Start=dat04 %>% pull(X2.y),
                    End=dat04 %>% pull(X3),
                    
                    color=dat04 %>% pull(color) %>% 
                      str_replace("#",""))
nlr.dis %>% head()

ideogram(karyotype = grape_karyotype,label = nlr.dis,label_type = "marker")
convertSVG("chromosome.svg", device = "png",width = 10,height = 8,dpi = 300)
convertSVG("chromosome.svg", device = "pdf",width = 8,height = 4)


### 还缺KaKs的图
### 还缺GO富集的图
### KEGG富集 结构域富集
## all samples private
orthogroups<-read_delim("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/20230815.use.all.hap/Orthogroups.tsv",
                        delim = "\t") %>% 
  select(Orthogroup,contains("V"))

bind_cols(orthogroups,freq.df) %>% head() %>% DT::datatable()

bind_cols(orthogroups,freq.df) %>% head()
#mutate(X1=str_count(V127,"mRNA1")) %>% 
#filter(X1>=2) %>% 
filter(Class=="Private") %>% 
  #select(V127_hap1,V127_hap2) %>% 
  mutate(new_col=paste(V127_hap1,V127_hap2,sep=", ") %>% str_replace_all(", ","\t")) %>% 
  select(new_col) %>% 
  filter(new_col != "NA\tNA") %>% 
  mutate(new_col=str_replace(new_col,"NA\t",""))%>% 
  mutate(new_col=str_replace(new_col,"\tNA","")) %>% 
  write_delim("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/20230815.use.all.hap/private.txt",
              delim = "\t",
              col_names = FALSE,
              quote = "none")

bind_cols(orthogroups,freq.df) %>% 
  #mutate(X1=str_count(V127,"mRNA1")) %>% 
  #filter(X1>=2) %>% 
  filter(Class=="Core") %>% 
  #select(V127_hap1,V127_hap2) %>% 
  mutate(new_col=paste(V127_hap1,V127_hap2,sep=", ") %>% str_replace_all(", ","\t")) %>% 
  select(new_col) %>% 
  filter(new_col != "NA\tNA") %>% 
  mutate(new_col=str_replace(new_col,"NA\t",""))%>% 
  mutate(new_col=str_replace(new_col,"\tNA","")) %>% 
  write_delim("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/20230815.use.all.hap/core.txt",
              delim = "\t",
              col_names = FALSE,
              quote = "none")

bind_cols(orthogroups,freq.df) %>% 
  #mutate(X1=str_count(V127,"mRNA1")) %>% 
  #filter(X1>=2) %>% 
  filter(Class=="SoftCore") %>% 
  #select(V127_hap1,V127_hap2) %>% 
  mutate(new_col=paste(V127_hap1,V127_hap2,sep=", ") %>% str_replace_all(", ","\t")) %>% 
  select(new_col) %>% 
  filter(new_col != "NA\tNA") %>% 
  mutate(new_col=str_replace(new_col,"NA\t",""))%>% 
  mutate(new_col=str_replace(new_col,"\tNA","")) %>% 
  write_delim("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/20230815.use.all.hap/softcore.txt",
              delim = "\t",
              col_names = FALSE,
              quote = "none")

bind_cols(orthogroups,freq.df) %>% 
  #mutate(X1=str_count(V127,"mRNA1")) %>% 
  #filter(X1>=2) %>% 
  filter(Class=="Dispensable") %>% 
  #select(V127_hap1,V127_hap2) %>% 
  mutate(new_col=paste(V127_hap1,V127_hap2,sep=", ") %>% str_replace_all(", ","\t")) %>% 
  select(new_col) %>% 
  filter(new_col != "NA\tNA") %>% 
  mutate(new_col=str_replace(new_col,"NA\t",""))%>% 
  mutate(new_col=str_replace(new_col,"\tNA","")) %>% 
  write_delim("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/20230815.use.all.hap/dispensable.txt",
              delim = "\t",
              col_names = FALSE,
              quote = "none")

# bind_cols(orthogroups,freq.df) %>% 
#   mutate(X1=str_count(V127,"mRNA1")) %>% 
#   filter(X1>=2) %>% 
#   filter(Class=="Core") %>% 
#   select(1,2) %>% 
#   write_delim("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/core.txt",
#               delim = "\t",
#               col_names = FALSE)
# 
# bind_cols(orthogroups,freq.df) %>% 
#   mutate(X1=str_count(V127,"mRNA1")) %>% 
#   filter(X1>=2) %>% 
#   filter(Class=="SoftCore") %>% 
#   select(1,2) %>% 
#   write_delim("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/softcore.txt",
#               delim = "\t",
#               col_names = FALSE)
# 
# bind_cols(orthogroups,freq.df) %>% 
#   mutate(X1=str_count(V127,"mRNA1")) %>% 
#   filter(X1>=2) %>% 
#   filter(Class=="Dispensable") %>% 
#   select(1,2) %>% 
#   write_delim("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/dispensable.txt",
#               delim = "\t",
#               col_names = FALSE)
# 
# 
# bind_cols(orthogroups,freq.df) %>%
#   filter(Class=="Core") %>% 
#   pull(V127) %>% 
#   str_split(", ") %>% unlist() %>% 
#   write_lines(file = "D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/CoreGeneID.txt")
# 
# bind_cols(orthogroups,freq.df) %>%
#   filter(Class=="SoftCore") %>% 
#   pull(V127) %>% 
#   str_split(", ") %>% unlist() %>% 
#   write_lines(file = "D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/SoftCoreGeneID.txt")
# 
# bind_cols(orthogroups,freq.df) %>%
#   filter(Class=="Dispensable") %>% 
#   pull(V127) %>% 
#   str_split(", ") %>% unlist() %>% 
#   write_lines(file = "D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/DispensableGeneID.txt")
# 
# bind_cols(orthogroups,freq.df) %>%
#   filter(Class=="Private") %>% 
#   pull(V127) %>% 
#   str_split(", ") %>% unlist() %>% 
#   write_lines(file = "D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/PrivateGeneID.txt")


## GO富集分析
## 以下代码在linux下使用，win系统下为啥报错暂时没有想明白

# library(tidyverse, quietly = TRUE)
# library(formattable, quietly = TRUE)
# #BiocManager::install("AnnotationForge")
# library(AnnotationForge, quietly = TRUE)
# library(seqinr, quietly = TRUE)
# library(clusterProfiler, quietly = TRUE)
# 
# emapper <- read_delim("D:/Jupyter/IAAS/grape/gwas/cbmEggnog.emapper.annotations", 
#                       "\t", escape_double = FALSE, col_names = FALSE, 
#                       comment = "#", trim_ws = TRUE) %>%
#   dplyr::select(GID = X1, 
#                 COG = X7,
#                 Gene_Name = X8,
#                 Gene_Symbol = X9,
#                 GO = X10,
#                 KO = X12,
#                 Pathway = X13
#   )
# 
# gene_info <- dplyr::select(emapper,  GID, Gene_Name) %>%
#   dplyr::filter(!is.na(Gene_Name)) %>%
#   dplyr::filter(Gene_Name != '-') 
# eggnog_anno <- length(gene_info$GID)
# 
# 
# gene2go <- dplyr::select(emapper, GID, GO) %>%
#   separate_rows(GO, sep = ',', convert = F) %>%
#   filter(!is.na(GO)) %>%
#   filter(str_detect(GO, '^GO')) %>%
#   mutate(EVIDENCE = 'IEA') 
# 
# 
# makeOrgPackage(gene_info=gene_info,
#                go=gene2go,
#                maintainer='MingYan <mingyan24@126.com>',
#                author='MingYan',
#                outputDir="./",
#                tax_id=29760,
#                genus='V',
#                species='v',
#                goTable="go",
#                version="1.0")
# 
# pkgbuild::build('./org.Vv.eg.db', dest_path = ".")

# geneid<-read_lines("CoreGeneID.txt")
# ego<-enrichGO(gene = geneid,OrgDb = org.Vv.eg.db::org.Vv.eg.db,keyType = 'GID',ont = "ALL",qvalueCutoff = 0.05,pvalueCutoff = 0.05)
# save(ego,file = "ego.Rdata")
#BiocManager::install("AnnotationDbi")

# library(clusterProfiler)
# 
# ?enricher
# 
# library(AnnotationHub)
# ah <- AnnotationHub()
# vv<-ah[["AH108013"]]
core_ego<-enrichGO(gene = geneid,
                   OrgDb = vv,
                   keyType = 'ENTREZID',
                   ont = "ALL",
                   qvalueCutoff = 0.05,
                   pvalueCutoff = 0.05)
load("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/core_ego.Rdata")

data.frame(core_ego) %>% 
  select(-geneID)
dotplot(core_ego)+
  labs(title="Core") -> p06

load("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/softcore_ego.Rdata")

data.frame(softcore_ego)
dotplot(softcore_ego) +
  labs(title="SoftCore")-> p07

load("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/dispen_ego.Rdata")

data.frame(dispen_ego)
dotplot(dispen_ego)+
  labs(title="Dispensable") -> p08

pdf(file = "D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/plot03.pdf",
    width = 12,height = 8)
p06+p07+p08+
  plot_layout(nrow = 2)
dev.off()
### KaKs


read_delim("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/private_kaks/CBM1chr13G019100.mRNA1-CBM1chr15G000090.mRNA1-CBM1chr15G000180.mRNA1-CBM1chr15G000220.mRNA1-CBM1chr15G000320.mRNA1-CBM1chr15G000380.mRNA1.cds_aln.kaks",
           delim = "\t") %>% 
  select(`Ka/Ks`)

list.files("private.kaks",
           pattern = "*.kaks",
           full.names = TRUE) %>% 
  map(read_delim,delim="\t",show_col_types=FALSE,progress=FALSE) %>% 
  bind_rows() %>% 
  select(`Ka/Ks`) %>% 
  mutate(Class="private") -> private.kaks


list.files("core.kaks/",
           pattern = "*.kaks",
           full.names = TRUE) %>% 
  map(read_delim,delim="\t",show_col_types=FALSE,progress=FALSE) %>% 
  bind_rows() %>% 
  select(`Ka/Ks`) %>% 
  mutate(Class="core") -> core.kaks

core.kaks

list.files("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/softcore_kaks/",
           pattern = "*.kaks",
           full.names = TRUE) %>% 
  map(read_delim,delim="\t",show_col_types=FALSE,progress=FALSE) %>% 
  bind_rows() %>% 
  select(`Ka/Ks`) %>% 
  mutate(Class="softcore") -> softcore.kaks

softcore.kaks

list.files("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/dispensable_kaks/",
           pattern = "*.kaks",
           full.names = TRUE) %>% 
  map(read_delim,delim="\t",show_col_types=FALSE,progress=FALSE) %>% 
  bind_rows() %>% 
  select(`Ka/Ks`) %>% 
  mutate(Class="dispensable") -> dispensable.kaks

dispensable.kaks


p09<-bind_rows(private.kaks,core.kaks,softcore.kaks,dispensable.kaks) %>% 
  mutate(Class=factor(Class,levels = c("core","softcore","dispensable","private"))) %>% 
  ggplot(aes(x=Class,y=`Ka/Ks`))+
  geom_boxplot(aes(fill=Class),
               outlier.alpha = 0)+
  ylim(NA,1.5)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = "none")

## diversity

library(tidyverse)
library(pegas)

myfun<-function(x){return(nuc.div(read.dna(x,format = "fasta")))}
list.files("private.fasta/",pattern = "*.fasta",full.names = TRUE)%>%map(myfun)%>%unlist() -> private.div
list.files("core.fasta/",pattern = "*.fasta",full.names = TRUE)%>%map(myfun)%>%unlist() -> core.div
list.files("softcore.fasta/",pattern = "*.fasta",full.names = TRUE)%>%map(myfun)%>%unlist() -> softcore.div
list.files("dispensable.fasta/",pattern = "*.fasta",full.names = TRUE)%>%map(myfun)%>%unlist() -> dispen.div

bind_rows(data.frame(div=private.div,group="private"),
          data.frame(div=core.div,group="core"),
          data.frame(div=softcore.div,group="softcore"),
          data.frame(div=dispen.div,group="dispensable")) %>% 
  write_delim(file = "nucle_diversity.txt",delim = "\t")

nuc_div<-read_delim("D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/20230815.use.all.hap/nucle_diversity.txt",
                    delim = "\t")
nuc_div %>% 
  mutate(y=factor(group,levels = c("core","softcore","dispensable","private")))
nuc_div %>% pull(group) %>% table()

p10<-ggplot(data=nuc_div%>% 
              mutate(group=factor(group,levels = c("core","softcore","dispensable","private"))) %>% 
              filter(div>0),
            aes(x=group,y=div))+
  geom_boxplot(aes(fill=group),outlier.alpha = 0)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  labs(x="Class",y="Nucleotide diversity")+
  theme(legend.title = element_blank(),
        legend.position = "none")

p10


pdf(file = "D:/Jupyter/IAAS/grape/20230726.reanalysis/scripts/PanGeneFamily/plot04.pdf",
    width = 8,height = 6)
p10+p09
dev.off()
