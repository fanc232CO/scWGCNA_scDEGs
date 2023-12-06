library(Seurat)
options(header=T,stringsAsFactors=F)
#load and process files
mat_expr<-read.csv(gzfile('../yangsy_GSE157783/matrix.tsv.gz'),sep='\t')
df_genes<-read.csv(gzfile('../yangsy_GSE157783/features.tsv.gz'),sep='\t')
df_genes<-df_genes[,c('gene','row')]
df_genes$row<-as.numeric(df_genes$row)
df_genes<-df_genes[order(df_genes$row),]
#for running of scDRS, transfer to gene symbol
fn_map='/home/fancong/database_fan/ID_map/hs_hg38_102_gtf_ensembl2symbol.map'
df_map=read.table(fn_map,header=F,stringsAsFactors=F)
l_ensembl<-sapply(df_genes$gene,function(x){strsplit(x,split='[.]')[[1]][1]})
l_symbol<-df_map$V2[match(l_ensembl,df_map$V1)]
#ensemble to no symbol, delete
mat_expr<-mat_expr[!is.na(l_symbol),]
l_symbol<-na.omit(l_symbol)
#multi ensemble to one symbol, delete.
mat_expr<-mat_expr[-which(duplicated(l_symbol)),]
rownames(mat_expr)<-l_symbol[-which(duplicated(l_symbol))]
l1<-gsub('[.]','-',colnames(mat_expr))
colnames(mat_expr)<-l1
df_meta<-read.csv(gzfile('../yangsy_GSE157783/barcodes.tsv.gz'),sep='\t')
rownames(df_meta)<-df_meta$barcode
df_meta<-df_meta[,-1]
seu_obj<-CreateSeuratObject(counts=mat_expr,meta.data=df_meta)
#saveRDS(seu_obj,file='./GSE157783_seu_obj.rds')
saveRDS(seu_obj,file='./GSE157783_seu_obj_symbol.rds')
