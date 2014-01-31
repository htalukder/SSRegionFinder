#(on $ R-devel)
#Getting Methylation Data 

library(minfiLocal)
load("/cbcb/lab/hcorrada/methyl/tcga/colon_Mset.rda")

gr=granges(colon_Mset)
clust=bumphunter::clusterMaker(as.character(seqnames(gr)),pos=start(gr),maxGap=1000)

m=getM(colon_Mset)

# pheno data
pd =phenoData(colon_Mset)