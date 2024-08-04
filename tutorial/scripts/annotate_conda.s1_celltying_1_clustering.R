# Inputs:
# sample_name
# dir_res
# f_a; obja
# f_b; objd
# df_meta; f_meta
# Output:
# ident_atac_clusters
update_coda(df_meta, obja, objd)
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

hdbscan_num_nb = 20; snn_res = 0.6

if (sample_name == 'DCIS66T_chip2') {hdbscan_num_nb = 20; snn_res = 0.6}
#------------------- ~~~ Viz & manually explore ~~~ -------------------  
# // additionally use the hdbscan to remove the possible doublets

hdbscan_str <- sprintf('peaks_hdbscan_umap_minPts.%s', hdbscan_num_nb)
snn_str <- sprintf('peaks_snn_res.%s', snn_res)
ident_use <- snn_str
ident_use <- hdbscan_str

Idents(obja) <- ident_use 
print(table(Idents(obja)))
UMAPPlot(obja, label=T, group.by=snn_str)  
UMAPPlot(obja, label=T, group.by=hdbscan_str)  
UMAPPlot(obja, label=T, group.by=ident_use)  
UMAPPlot(obja, label=T, group.by='celltypes')
ruok::qtable_scatter(
  pred = obja@meta.data[, snn_str], 
  groudtruth = obja@meta.data[, hdbscan_str])

#------------------- ~~~ Pick a choice ~~~ -------------------  

ident_atac_clusters <- hdbscan_str
Idents(obja) <- ident_atac_clusters

p <- UMAPPlot(obja, label=T) + theme(aspect.ratio = 1) + 
  labs(title = ident_atac_clusters, caption = sprintf('%s cells', ncol(obja)))
ggsave(file.path(dir_res, sprintf('dr.ident_atac_clusters.pdf')), 
       p, width = 6, height = 5, useDingbats = F)
