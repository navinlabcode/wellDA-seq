#--------------------------
# Initiate wellDA object
#--------------------------
if (T) {
  library(Seurat)
  library(Signac)
  library(tidyverse)
  library(fs)
  library(copykit)
  library(VennDiagram)
  library(SummarizedExperiment)
  library(ggpubr); library(ruok); library(scales)
  
  setwd("/volumes/USR1/yyan/project/coda")
  source('/volumes/USR1/yyan/project/coda/rsrc/utils.R')
}
options <- commandArgs(trailingOnly = TRUE)

#------ DCIS ------
if (F) {
  run_name <- 'DCIS66T_chip2'
  f_obja <- 'rds/DCIS66T_chip2_atac/ready.signac_default.rds'
  f_objd <- 'rds/DCIS66T_chip2_dna/ready_clean.copykit.rds'
}

if (length(options) > 0) {
  run_name = options[[1]]; message(run_name)
  f_obja   = options[[2]]; message(f_obja)
  f_objd   = options[[3]]; message(f_objd)
}


dir_res <- sprintf('./rds_coda/%s', run_name); message(dir_res)
dir_create(dir_res)
#------ coda cell dictionary ------
if (!file.exists(file.path(dir_res, 'objd.rds')) | !file.exists(file.path(dir_res, 'obja.rds'))) {
# if (T) {
  obja <- read_rds(f_obja)
  objd <- read_rds(f_objd)
  
  ## adhoc deal with older DNA objects
  if (!'clones' %in% colnames(colData(objd)))  {
    colData(objd)$clones <- colData(objd)$superclones
  }
  if (!'ploidy_class' %in% colnames(colData(objd))) {
    colData(objd)$ploidy_class <- ifelse(colData(objd)$is_aneuploid, 'aneuploid', 'diploid')
  }
  

  str(Cells(obja))
  str(colnames(objd))
  
  cnames_a <- Cells(obja)
  cnames_d <- colnames(objd)
  

  cnames_a <- Cells(obja); str(cnames_a)
  cnames_d <- colnames(objd); str(cnames_d)
  
  library(grDevices)
  p <- VennDiagram::venn.diagram(
    list(ATAC=cnames_a, DNA=cnames_d), filename = NULL,
    # col = c('blue', 'gold'), cat.col = c('blue', 'gold'),
    disable.logging=T, cat.default.pos = "text")
  pdf(file.path(dir_res, 'vennplot.cellnames.beforeinteresection.pdf'), width = 5, height = 5, useDingbats = F)
  grid.draw(p)
  dev.off()
  cnames_coda <- intersect(Cells(obja), colnames(objd)); str(cnames_coda)
  
  #------ Know what cells are removed ------
  pa <- DimPlot(obja, cells.highlight = setdiff(Cells(obja), cnames_coda)) + 
    rremove('legend') + 
    labs(title = 'ATAC', 
         caption = sprintf('%s / %s cells are to be removed', length(setdiff(Cells(obja), cnames_coda)), ncol(obja)))
  pa <- pa + theme(aspect.ratio = 1)
  colData(objd)$coda_exclude <- !colnames(objd) %in% cnames_coda
  pd <- copykit::plotUmap(objd, label = 'coda_exclude') + 
    scale_fill_manual(values = c(`FALSE` = 'grey', `TRUE`='red')) + 
    rremove('legend') + theme(aspect.ratio = 1) + 
    labs(title = 'CNA', 
         caption = sprintf('%s / %s cells are to be removed', sum(colData(objd)$coda_exclude), ncol(objd)))
  pdf(file.path(dir_res, 'dimplot.cellnames.beforeinteresection.pdf'), width = 5, height = 5, useDingbats = F, onefile = T)
  print(pa)
  print(pd)
  dev.off()
  colData(objd)$coda_exclude <- NULL # no more neede
  
  obja <- obja[, cnames_coda]
  objd <- objd[, cnames_coda]
  str(cnames_coda)
  
  table(colData(objd)$subclones)
  table(colData(objd)$superclones)
  
  colData(objd)$sum_bincounts <- colSums(bincounts(objd))

  #------ Create coda cell meta.data data frame ------
  if ('celltypes' %in% colnames(obja@meta.data)) {
    obja$seurat_clusters <- Idents(obja) <- obja$celltypes
  }
  if (!'clones' %in% colnames(colData(objd))) {
    colData(objd)$clones <- colData(objd)$subclones
  }
  
  
  df_d <- as.data.frame(colData(objd))
  if (any( str_detect(colnames(df_d), '^ATAC_') | str_detect(colnames(df_d), '_ATAC$') )) {
    # this object has been in coda once
    tmp <- str_detect(colnames(df_d), '^ATAC_') | str_detect(colnames(df_d), '^CNA_') | str_detect(colnames(df_d), '_ATAC$') | str_detect(colnames(df_d), '_CNA$')
    tmp <- which(tmp); length(tmp)
    df_d <- dplyr::select(df_d, -c(tmp))
  }
  # df_d <- df_d %>% dplyr::select(dplyr::where(function(x) any(!is.na(x))))
  # colnames(df_d) <- paste0('CNA_', colnames(df_d)) ## stop adding CNA_*prefix
  colData(objd) <- as(df_d, 'DFrame')
  # df_d$cellname <- rownames(df_d)
  
  
  df_a <- obja@meta.data
  if (any( str_detect(colnames(df_a), '^CNA_') | str_detect(colnames(df_a), '_CNA$') )) {
    # this object has been in coda once
    tmp <- str_detect(colnames(df_a), '^ATAC_') | str_detect(colnames(df_a), '^CNA_') | str_detect(colnames(df_a), '_ATAC$') | str_detect(colnames(df_a), '_CNA$')
    tmp <- which(tmp)
    df_a <- dplyr::select(df_a, -c(tmp))
  }
  # colnames(df_a) <- paste0('ATAC_', colnames(df_a)) ## stop adding ATAC_*prefix
  obja@meta.data <- df_a
  
  # df_a$cellname <- rownames(df_a)
  stopifnot( all.equal(as.character(Cells(obja)), rownames(df_a)) )
  stopifnot( all.equal(colnames(objd), rownames(colData(objd))) )
  stopifnot( all.equal(as.character(Cells(obja)), colnames(objd)) )
  stopifnot( all.equal(as.character(Cells(obja)), cnames_coda) )
  
  umap_a <- as.data.frame(obja@reductions$umap@cell.embeddings)
  colnames(umap_a) <- c('UMAP_1_ATAC', 'UMAP_2_ATAC')
  umap_d <- as.data.frame(reducedDim(objd, 'umap'))
  colnames(umap_d) <- c('UMAP_1_CNA', 'UMAP_2_CNA')
  stopifnot( all.equal( rownames(umap_a), rownames(umap_d) ) )
  
  coda_metadf <- cbind(df_d, df_a, umap_a, umap_d)
  coda_metadf$cellname <- cnames_coda
  coda_metadf <- left_update2(coda_metadf, df_a, df_d)
  
  obja@meta.data <- subset(coda_metadf, select = -c(UMAP_1_ATAC, UMAP_2_ATAC, UMAP_1_CNA, UMAP_2_CNA))
  colData(objd) <- as(subset(coda_metadf, select = -c(UMAP_1_ATAC, UMAP_2_ATAC, UMAP_1_CNA, UMAP_2_CNA)), 'DFrame')
  
  if (!identical(cnames_d, cnames_coda)) {
    objd <- runDistMat(objd, metric = 'manhattan', n_threads = 20)
    objd <- runPhylo(objd)
  }
  colData(objd)$sample <- colnames(objd)
  # export
  write_rds(objd, file.path(dir_res, 'objd.rds'))
  write_rds(obja, file.path(dir_res, 'obja.rds'))
  write.csv(coda_metadf, file.path(dir_res, 'metadata.csv'))
  write_rds(coda_metadf, file.path(dir_res, 'metadata.df.rds'))

} else {
  # load
  message('load the existing files')
  f_obja <- file.path(dir_res, 'obja.rds')
  f_objd <- file.path(dir_res, 'objd.rds')
  obja <- read_rds(f_obja)
  objd <- read_rds(f_objd)  
  coda_metadf <- read_rds(file.path(dir_res, 'metadata.df.rds'))
}

#------ Reclusering ATAC and DNA ------
objd <- copykit::runPca(objd)
# umap
n_neighbors = 30
objd <- runUmap(objd, n_neighbors=n_neighbors, n_threads=20, min_dist = 0)
if (T) {
  k_superclones <- 50; k_subclones <- 50
  findClusters_embedding <- 'umap'
  findClusters_method <- 'hdbscan'
  findClusters_ncomponents <- 2
}
objd  <- findClusters(
  objd, 
  k_superclones = k_superclones,
  k_subclones = k_subclones,
  embedding = findClusters_embedding, 
  method = findClusters_method, 
  ncomponents = findClusters_ncomponents)
# consensu phylo
if ('clones' %in% colnames(objd@colData) ) {
	objd <- calcConsensus(objd, consensus_by = 'clones')
} else {
	objd <- calcConsensus(objd, consensus_by = 'subclones')
}
objd <- runConsensusPhylo(objd)

DefaultAssay(obja) = 'peaks'
p <- DepthCor(obja)+geom_hline(yintercept = c(0.75, -0.75) , lty='dashed', color='red')
bad_lsi <- p$data$Component[abs(p$data$counts) > 0.75]
dim_lsi <- head(setdiff(c(1:100), bad_lsi), 50)
obja <- RunUMAP(
  obja, reduction = 'lsi', dims = dim_lsi, verbose = T,
  metric = 'correlation',
  n.neighbors=ifelse(ncol(obja) < 100, 5, 30))
obja <- FindNeighbors(obja, reduction = "lsi", dims = dim_lsi, verbose = F,
                          k.param=ifelse(ncol(obja) < 20, 5, 30))
DefaultAssay(obja) = 'peaks'
tmp <- str_detect(colnames(obja@meta.data), pattern = 'peaks_snn')
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

update_coda(coda_metadf, obja, objd)

write_rds(objd, file.path(dir_res, 'objd.rds'))
write_rds(obja, file.path(dir_res, 'obja.rds'))
write.csv(coda_metadf, file.path(dir_res, 'metadata.csv'))
write_rds(coda_metadf, file.path(dir_res, 'metadata.df.rds'))

## 
# stop('\nNext Run master.coda_refine.R\n')


#-------------------------- Viz --------------------------  
pal_ploidy_class = c('aneuploid'='#cd0bbc', 'diploid'='#61d04f', 'Unknown'='grey', `NA`='grey')
pal_ploidy_is_aneuploid = c(`TRUE`='#cd0bbc', `FALSE`='#61d04f', 'Unknown'='grey', `NA`='grey')
# pal_cna_clones = new_palette_D(levels(coda_metadf$clones), pal = 'ironMan')
pal_cna_clones = new_palette_D(sort(unique(coda_metadf$clones)), pal = 'ironMan')
pal_cna_superclones = ArchR::paletteDiscrete(coda_metadf$superclones)
pal_cna_subclones = ArchR::paletteDiscrete(coda_metadf$subclones, set = 'ironMan')
pal_atac_clusters = ArchR::paletteDiscrete(coda_metadf$seurat_clusters, set = 'circus')
pal_run = new_palette_D(unique(coda_metadf$run), pal = 'Dynamic')
pal_combo = list(
  `superclones` = pal_cna_superclones,
  `subclones` = pal_cna_subclones,
  `clones` = pal_cna_clones,
  `ploidy_class` = pal_ploidy_class, 
  `is_aneuploid` = pal_ploidy_is_aneuploid,
  `run` = pal_run
)


pdf(file.path(dir_res, 'atac.umap.ATAC_seurat_clusters.pdf'), 
    width = 6, height = 5, useDingbats = F)
p <- DimPlot(obja, reduction = 'umap', group.by='seurat_clusters', label=T) + theme(aspect.ratio = 1) +
    labs(caption = 'ATAC space') +
    scale_color_manual(values=pal_atac_clusters)
print(p)
dev.off()

z_opts <- c('superclones', 'subclones', 'clones', 
            'ploidy_class', 'run')
z_opts <- intersect(z_opts, colnames(colData(objd)))
for (z in z_opts) {
  message(z)
  pal_use <- pal_combo[[z]]
  pdf(file.path(dir_res, sprintf('cna.umap.%s.pdf', z)), width = 6, height = 5, useDingbats = F)
  p <- plotUmap(objd, label = z) + theme(aspect.ratio = 1) +
    labs(caption = 'CNA space', fill = z) +
    scale_fill_manual(values=pal_use) 
  print(p)
  dev.off()  
}

for (z in z_opts) {
  pal_use <- pal_combo[[z]]
  pdf(file.path(dir_res, sprintf('atac.umap.%s.pdf', z)), width = 6, height = 5, useDingbats = F)
  p <- DimPlot(obja, reduction = 'umap', group.by=z) + 
    theme(aspect.ratio = 1) +
    labs(caption = 'ATAC space') +
    scale_color_manual(values=pal_use) 
  print(p)
  dev.off()  
}

pdf(file.path(dir_res, 'cna.umap.ATAC_seurat_clusters.pdf'), width = 6, height = 5, useDingbats = F)
p <- plotUmap(objd, label = 'seurat_clusters') + theme(aspect.ratio = 1) +
  labs(caption = 'CNA space') +
  scale_fill_manual(values=pal_atac_clusters)
print(p)
dev.off()

if (any(str_detect(colnames(obja@meta.data), 'run'))) {
  try(colData(objd)$run <- as.factor(colData(objd)$run.ATAC))
  pd <- plotUmap(objd[, sample(colnames(objd), ncol(objd))], label = 'run') + theme(aspect.ratio = 1) + 
    scale_fill_manual(values = pal_run)
  pa <- UMAPPlot(obja, group.by = 'run', shuffle=T) + theme(aspect.ratio = 1) + 
    scale_color_manual(values = pal_run)
  ggsave(pd, filename = file.path(dir_res, 'cna.umap.run.pdf'), width = 6, height = 5, useDingbats = F)
  ggsave(pa, filename = file.path(dir_res, 'atac.umap.run.pdf'), width = 6, height = 5, useDingbats = F)
}

library(ComplexHeatmap)
for (y in z_opts) {
  
  confusionmat <- table(coda_metadf[[y]],
                        coda_metadf[['seurat_clusters']]) %>%
    as.data.frame.matrix() %>% as.matrix()
  rownames(confusionmat) <- paste0('DNA_', rownames(confusionmat))
  colnames(confusionmat) <- paste0('ATAC_', colnames(confusionmat))
  write.table(confusionmat, file=file.path(
    dir_res, sprintf('table.ATAC_clusters_x_%s.txt', y)))
  pdf(file.path(
    dir_res, sprintf('confusion_matrix.ATAC_clusters_x_%s.pdf', y)), 
    width = 6, height = 6, useDingbats = F)
  p <- Heatmap(
    confusionmat, name='cell number',
    cluster_rows = F, cluster_columns = T, 
    row_names_side  = 'left', 
    col = circlize::colorRamp2(quantile(confusionmat, c(0.01, 0.99)), hcl.colors(n=2, 'Purples 2', rev = T)),
    show_row_dend = T, show_column_dend = T, 
    width = unit(3, 'inch'),
    height = unit(3/ncol(confusionmat)*nrow(confusionmat), 'inch'),
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (confusionmat[i, j] > 0) {
        grid.text(sprintf("%d", confusionmat[i, j]), x, y, gp = gpar(fontsize = 10, col='black'))
      } 
    },
    use_raster = T
  )
  draw(p)
  dev.off()
}  

if (all(c('celltypes', 'clones') %in% colnames(coda_metadf))) {
  p1 <- coda_dimplot2(coda_metadf, 'celltypes', pal='circus', ncol=1, do_shuffle = T)
  p2 <- coda_dimplot2(coda_metadf, 'clones', pal='stallion', ncol=1, 
                      color_specific = c('D' = 'lightsteelblue'), do_shuffle = T)
  p <- wrap_plots(p1, p2, ncol=2, widths = c(1,1))
  ggsave(filename = file.path(dir_res, 'coda_dimplot4.pdf'), 
         plot = p,
         width = 10, height = 8, useDingbats = F)
}

if ('run' %in% colnames(df_meta)) {
  ggsave(filename = file.path(dir_res, 'coda_dimplot.run.pdf'), 
         plot = coda_dimplot2(coda_metadf, 'run', pal = 'bear', do_shuffle = 'yes'),
         width = 8, height = 5, useDingbats = F)  
}

try(colData(objd)$ploidy_confidence_gt_1 <- colData(objd)$ploidy_confidence>1)
try(
  colData(objd)$run <- as.numeric(as.factor(colData(objd)$run.ATAC)))
try(dev.off())
pdf(file.path(dir_res, 'cna.heatmap.combo_coda.%03d.pdf'), 
    width = 15, height = 10, useDingbats = F, onefile = F)
for (z in c('clones', 'subclones', 'superclones', 'seurat_clusters')){
  cat('\n', z, '...')
  tmp <- c('ploidy_class', 'is_aneuploid', 'clones', 'subclones', 'superclones', 
           'seurat_clusters', 'run')
  tmp <- intersect(tmp, colnames(colData(objd)))
  pal_list <- list(`superclones`=pal_cna_superclones,
                   `subclones`=pal_cna_subclones,
                   `clones` = pal_cna_clones,
                   `seurat_clusters`=pal_atac_clusters,
                   `ploidy_class`=pal_ploidy_class,
                   `run` = pal_run)
  p <- try( copykit::plotHeatmap(
    objd, 
    col = pal_copykit_ratio_clip1,
    order_cells = 'consensus_tree', 
    row_split = z,
    label=tmp, 
    label_colors = pal_list[tmp]))
  if ('try-error' %in% class(p)) {next()}
  draw(p)
}; cat('\n'); try(rm(tmp))
# colData(objd)$subclones <- NULL # tmp
dev.off()


cat('DONE')
#--------------------------  --------------------------  


#------ QC ------
library(ggpubr)
colnames(colData(objd))
# library(ggbeeswarm)
pdf(file.path(dir_res, 'boxplot.qc.pdf'), width = 8, height = 5, 
    onefile = T, useDingbats = F)
for (x in c('subclones', 'superclones', 'ploidy_class', 
            'clones', 'seurat_clusters')) {
  pal_use <- switch (x,
                     superclones = pal_cna_superclones,
                     subclones = pal_cna_subclones,
                     clones = pal_cna_clones,
                     ploidy_class = pal_ploidy_class, 
                     seurat_clusters = pal_atac_clusters
  )
  
  for (z in c('sum_bincounts', 'ReadsKept', 'TotalReads', 
              'MedianBinCount', 
              'reads_count', 'nFrags', 'TSS.enrichment')) {
    
    p1 <-  ggplot(coda_metadf, aes_string(x=x, y=z, color=x)) + 
      ggbeeswarm::geom_quasirandom(dodge.width=0.5) + 
      scale_color_manual(values = pal_use)+ 
      rremove('legend') + rotate_x_text() + 
      stat_compare_means()
    
    if ('try-error' %in% class(p1)) {next()}
    if (!z %in% c('TSS.enrichment', 'MedianBinCount')) {
      p1 <- p1 + scale_y_log10()
    }
    print(p1)
  }
}; rm(x); rm(z)
dev.off()
library(ggpubr)
for (z in c('clones', 'subclones')) {
  p <- ggscatter(as.data.frame(colData(objd)),
                 x='ploidy', y='ploidy_confidence',
                 facet.by= z) +
    geom_hline(yintercept = 1, lty='dashed', color='blue')  +
    geom_vline(xintercept = 2, lty='dashed', color='blue') 
  print(p)
  ggsave(file.path(dir_res, sprintf('qc.scquantum.ploidy_by_%s.pdf', z)), p, width = 10,height = 10, useDingbats = F)
}


# stop('Stop as expected')
cat('-- DONE -- ')
