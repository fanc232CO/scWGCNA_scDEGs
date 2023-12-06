library(Seurat)
library(rliger)
library(hdWGCNA)
library(Matrix)
CT<-'Astrocytes'
TS<-'Cerebrum'
#load files
matr_expr<-readRDS('/home/fancong/database_fan/pair_scATAC_scRNAseq/scRNAseq/Cerebrum_gene_count.RDS')
df_cell<-readRDS('/home/fancong/database_fan/pair_scATAC_scRNAseq/meta/df_cell.RDS')
#solve batch effect to generate seu_obj
get_liger_batch<-function(){
    ary_expr<-tapply(df_cell1$sample,df_cell1$Experiment_batch,function(x){matr_expr1[,intersect(x,colnames(matr_expr1))]}) #array format with name
    list_expr<-list('exp1'=ary_expr[1][[1]],'exp2'=ary_expr[2][[1]],'exp3'=ary_expr[3][[1]],'exp4'=ary_expr[4][[1]],'exp5'=ary_expr[5][[1]])
    l_obj<-createLiger(list_expr)
    l_obj<-normalize(l_obj)
    l_obj<-selectGenes(l_obj)
    l_obj<-scaleNotCenter(l_obj)
    l_obj<-online_iNMF(l_obj,k=20,max.epoch=5)
    l_obj<-quantile_norm(l_obj)
    l_obj<-runUMAP(l_obj)
    return(l_obj)
}
df_cell1<-df_cell[which(df_cell$Main_cluster_name==CT&df_cell$Organ==TS),c('sample','Fetus_id','Main_cluster_name','Experiment_batch')]
l_cell<-intersect(df_cell1$sample,colnames(matr_expr))
rownames(df_cell1)<-df_cell1$sample
df_cell1<-df_cell1[l_cell,]
matr_expr1<-matr_expr[,l_cell]
rm(matr_expr,df_cell)#save memory
seu_obj<-CreateSeuratObject(counts=matr_expr1,meta.data=df_cell1) #delete sample column.
liger_obj<-get_liger_batch()
rm(matr_expr1,df_cell1) #save memory
fn_liger<-paste0('liger_',CT,'_batch.rds')
saveRDS(liger_obj,file=fn_liger)
#liger_obj<-readRDS(fn_liger)
#transfer iMAF matrix to seurat objec
seu_obj@reductions$ctiNMF<-CreateDimReducObject(
    loadings=t(liger_obj@W),
    embeddings=liger_obj@H.norm[colnames(seu_obj),],
    key='ctiNMF_',
    assay='RNA'
)
VariableFeatures(seu_obj)<-liger_obj@var.genes
rm(liger_obj)#save memory
#Why scale again?
seu_obj<-ScaleData(seu_obj,features=VariableFeatures(seu_obj))
seu_obj<-RunUMAP(seu_obj,reduction='ctiNMF',dims=1:20)
saveRDS(seu_obj,file='seu_cerebrum_Astrocytes.rds')
seu_obj<-SetupForWGCNA(seu_obj,gene_select='fraction',fraction=0.05,wgcna_name='WGCNA')
seu_obj<-MetacellsByGroups(
       seurat_obj=seu_obj,
       group.by=c('Fetus_id'), #character vector required.
       ident.group='Fetus_id',
       reduction='umap'
)
seu_obj<-NormalizeMetacells(seu_obj)
meta_obj<-seu_obj@misc$WGCNA[[2]]
seu_obj <- SetDatExpr(
  seu_obj,
  group_name = unique(meta_obj@meta.data$Fetus_id), # the name of the group of interest in the group.by column
  group.by='Fetus_id', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  slot = 'data' # using normalized data
)
seu_obj <- TestSoftPowers(
  seu_obj,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)
saveRDS(seu_obj,file='seu_softPower.rds')
power_table <- GetPowerTable(seu_obj)
write.csv(power_table,file='power_table.csv')
