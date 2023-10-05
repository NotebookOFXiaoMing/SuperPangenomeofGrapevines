library(clusterProfiler)
library(tidyverse)
term2gene<-read_delim("term2gene.txt",delim = "\t",col_names = FALSE)
term2name<-read_delim("go.tb",delim = "\t")

core.gene.list<-read_lines("allCoreID.txt")
enricher(gene = core.gene.list,
         TERM2NAME = term2name,
         TERM2GENE = term2gene,
         pvalueCutoff = 0.05,
         qvalueCutoff = 0.05) -> core.enrich

softcore.gene.list<-read_lines("allDispensableID.txt")
enricher(gene = softcore.gene.list,
         TERM2NAME = term2name,
         TERM2GENE = term2gene,
         pvalueCutoff = 0.05,
         qvalueCutoff = 0.05) -> softcore.enrich

dispen.gene.list<-read_lines("allSoftCoreID.txt")
enricher(gene = dispen.gene.list,
         TERM2NAME = term2name,
         TERM2GENE = term2gene,
         pvalueCutoff = 0.05,
         qvalueCutoff = 0.05) -> dispen.enrich

private.gene.list<-read_lines("allPrivatePepID.txt")
enricher(gene = private.gene.list,
         TERM2NAME = term2name,
         TERM2GENE = term2gene,
         pvalueCutoff = 0.05,
         qvalueCutoff = 0.05) -> private.enrich


save(core.enrich,softcore.enrich,dispen.enrich,private.enrich,
     file="enrich.Rdata")
