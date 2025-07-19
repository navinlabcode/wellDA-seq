
if (ncol(obja) < 100) {
  n.neighbors = 5
  k.param <- 5
} else {
  n.neighbors = 20
  k.param <- 20
}
message(n.neighbors); message(k.param)

DefaultAssay(obja) = 'peaks'
obja <- RunTFIDF(obja)
obja <- FindTopFeatures(obja, min.cutoff = 'q5')
obja <- RunSVD(obja, assay='peaks', 
               n = ifelse(ncol(obja) < 100, ncol(obja)-1, 100), 
               reduction.key = 'LSI_', reduction.name = 'lsi')
p <- DepthCor(obja)+geom_hline(yintercept = c(0.75, -0.75) , lty='dashed', color='red')
bad_lsi <- p$data$Component[abs(p$data$counts) > 0.75]
dim_lsi <- head(setdiff(c(1:100), bad_lsi), 50)
obja <- RunUMAP(
  obja, reduction = 'lsi', dims = dim_lsi, verbose = T,
  metric = 'correlation',
  n.neighbors=n.neighbors)
obja <- FindNeighbors(obja, reduction = "lsi", dims = dim_lsi, verbose = F,
                      k.param=k.param)
DefaultAssay(obja) = 'peaks'
tmp <- str_detect(colnames(obja@meta.data), pattern = 'peaks_snn|peaks_hdbscan_umap')
tmp <- colnames(obja@meta.data)[tmp]
for (x in tmp) {obja@meta.data[[x]] <- NULL}; rm(tmp); rm(x)
for (sr3_snn_res in rev(seq(0.1, 0.6, by=0.1))){
  obja <- FindClusters(object = obja, resolution = sr3_snn_res, 
                       verbose = F, n.start = 10, algorithm = 3)
  sr3_snn_str <- sprintf('peaks_snn_res.%s', sr3_snn_res)
}
for (hdbscan_num_nb in rev(seq(5, 50, by=5)) ){
  obja <- run_hdbscan_sr(obja, n_neighbors = hdbscan_num_nb, reduction = 'umap')
  hdbscan_str <- sprintf('peaks_hdbscan_umap_minPts.%s', hdbscan_num_nb)
}
