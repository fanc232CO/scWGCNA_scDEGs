library(Seurat)
library(edgeR)
#fn_seu<-'./preproc/GSE157883_seu_obj.rds'
#seu_obj<-readRDS(fn_seu)
#pseudo-bulk samples
#y<-Seurat2PB(seu_obj,sample='cell_ontology',cluster='patient')
#save(y,file='GSE157883_y.RData')
load('GSE157883_y.RData')
#filtering and normlization
keep.samples<- y$samples$lib.size > 5e4
y<-y[,keep.samples]
keep.genes<-filterByExpr(y,group=y$samples$cluster,min.count=10,min.total.count=20)
y<-y[keep.genes,]
#TMM normalization
y<-normLibSizes(y)
#keep only Astrocytes
mat1<-y$counts
sele_cols<-grep('Astrocytes',colnames(mat1))
mat1<-mat1[,sele_cols]
df1<-y$samples
df1<-df1[df1$sample=='Astrocytes',]
df1$group<-gsub('[[:digit:]]','',df1$cluster)
df1$group<-factor(df1$group,levels=c('C','PD'))
y$samples<-df1
y$counts<-mat1
#design matrix
design<-model.matrix(~0+group,data=y$samples)
#disperse calculation
y<-estimateDisp(y,design,robust=T)
#differential expression analysis
et<-exactTest(y,pair=c('C','PD'))
save(y,et,file='res_exactTest.RData')
