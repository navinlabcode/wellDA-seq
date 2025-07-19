# Inputs:
message(sample_name)
message(dir_res)
print(obja)
print(head(df_meta))

# Output:
# ident_atac_clusters

ident_atac_clusters <- NULL

DefaultAssay(obja) <- 'peaks'
library(clustree)
p <- clustree::clustree(obja, prefix='peaks_snn_res.')
ggsave(file.path(dir_res, 'clustree_peaks_snn_res.pdf'), p,
       width = 10, height = 10)
# for (x in colnames(df_meta)[tmp]) { df_meta[, x] <- as.factor(df_meta[, x])}
p <- clustree::clustree(obja, prefix='peaks_hdbscan_umap_minPts.')
ggsave(file.path(dir_res, 'clustree_peaks_hdbscan_umap_minPts.pdf'), p, 
       width = 10, height = 10)

tmp <- str_detect(colnames(df_meta), 'peaks_hdbscan_umap_minPts.')
print(colnames(df_meta)[tmp])
tmp <- str_detect(colnames(df_meta), 'peaks_snn_res.')
print(colnames(df_meta)[tmp])

#------------------- ~~~ Choices per sample ~~~ -------------------  

hdbscan_num_nb = 30; snn_res = 0.6
if (sample_name == 'DCIS22T') {hdbscan_num_nb = 30; snn_res = 0.6}
if (sample_name == 'DCIS28T') {hdbscan_num_nb = 30; snn_res = 0.6}
if (sample_name == 'DCIS41T') {hdbscan_num_nb = 20; snn_res = 0.6}
if (sample_name == 'DCIS51T') {hdbscan_num_nb = 10; snn_res = 0.6}
if (sample_name == 'DCIS67T') {hdbscan_num_nb = 20; snn_res = 0.6}
hdbscan_num_nb = 20; snn_res = 0.6
#------------------- ~~~ Viz & manually explore ~~~ -------------------  
# // additionally use the hdbscan to remove the possible doublets

hdbscan_str <- sprintf('peaks_hdbscan_umap_minPts.%s', hdbscan_num_nb)
snn_str <- sprintf('peaks_snn_res.%s', snn_res)
ident_use <- snn_str
ident_use <- hdbscan_str

Idents(obja) <- ident_use 
print(table(Idents(obja)))
UMAPPlot(obja, label=T, group.by=ident_use) + UMAPPlot(obja, label=T, group.by='celltypes')
ruok::qtable_scatter_hclust(
  pred  = obja@meta.data[, 'celltypes'],
  groundtruth = obja@meta.data[, ident_use])

#------------------- ~~~ Finalize a choice ~~~ -------------------  

ident_atac_clusters <- hdbscan_str
Idents(obja) <- ident_atac_clusters

p1 <- UMAPPlot(obja, label=T) + theme(aspect.ratio = 1) + 
  labs(title = ident_atac_clusters, caption = sprintf('%s cells', ncol(obja)))
p2 <- UMAPPlot(obja, group.by='celltypes', label=T) + theme(aspect.ratio = 1) + 
  labs(title = 'celltypes before refine', caption = sprintf('%s cells', ncol(obja)))
ggsave(file.path(dir_res, sprintf('dr.celltypes_before_refine.pdf')), 
       p1+p2, width = 12, height = 5, useDingbats = F)

message('Final choice: ident_atac_clusters = ', ident_atac_clusters)
