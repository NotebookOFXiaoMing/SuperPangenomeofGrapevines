library(peer)

args <- commandArgs(trailingOnly = TRUE)

expr<-read.table(args[1],row.names=1,header=TRUE,sep="\t")
model = PEER()
PEER_setPhenoMean(model,as.matrix(t(expr)))
PEER_setNk(model,20)
PEER_getNk(model)
PEER_update(model)
factors = as.data.frame(t(PEER_getX(model)))
colnames(factors)<-colnames(expr)
write.table(factors,args[2],quote=F,row.names=F,sep="\t",col.names=TRUE)