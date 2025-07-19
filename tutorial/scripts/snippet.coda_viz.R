#---------------------------
# dir_snippet; df_meta; objd; 
# pal_list <- list(`superclones`=pal_cna_superclones,
# `subclones`=pal_cna_subclones,
# `clones` = pal_cna_clones,
# `seurat_clusters`=pal_atac_clusters,
# `ploidy_class`=pal_ploidy_class,
# `run` = pal_run)                 
#---------------------------
message(dir_snippet)
print(head(df_meta))
print(objd)
# z_tofig_opts <- c('clones',  'seurat_clusters', 'celltypes')
# z_tofig_opts <- c('subclones', 'superclones', 'seurat_clusters', 'celltypes')
print(z_tofig_opts)

pal_ploidy_class = c('aneuploid'='#cd0bbc', 'diploid'='#61d04f', 'Unknown'='grey', `NA`='grey')
pal_ploidy_is_aneuploid = c(`TRUE`='#cd0bbc', `FALSE`='#61d04f', 'Unknown'='grey', `NA`='grey')
pal_cna_clones = new_palette_D(levels(df_meta$clones), pal = 'stallion')
pal_cna_superclones = new_palette_D(levels(df_meta$superclones))
pal_cna_subclones = new_palette_D(levels(df_meta$subclones), pal = 'ironMan')
pal_atac_clusters = new_palette_D(levels(df_meta$seurat_clusters), pal = 'circus')
pal_celltypes = new_palette_D(levels(df_meta$celltypes), pal = 'circus')
pal_run = new_palette_D(unique(df_meta$run), pal = 'scales_hue') # consistent with copykit
pal_ploidy_type_final <- c('A' = '#cd0bbc', 'D'='#61d04f')
if ('D' %in% names(pal_cna_clones)) {pal_cna_clones['D'] <- '#61d04f'}
#------ confusion matrix ------
  
message('confusionmat')
for (y in c('ploidy_class', 'clones', 'subclones', 'superclones', 'celltypes')) {
  if (!y %in% colnames(df_meta)) {next()}
  
  confusionmat <- table(df_meta[[y]],
                        df_meta[['seurat_clusters']]) %>%
    as.data.frame.matrix() %>% as.matrix()
  rownames(confusionmat) <- paste0('DNA_', rownames(confusionmat))
  colnames(confusionmat) <- paste0('ATAC_', colnames(confusionmat))
  write.table(confusionmat, file=file.path(
    dir_snippet, sprintf('table.ATAC_clusters_x_%s.txt', y)))
  pdf(file.path(
    dir_snippet, sprintf('confusion_matrix.ATAC_clusters_x_%s.pdf', y)), 
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

#------ coda_dimplot2 ------
message('coda_dimplot2')
if ('celltypes' %in% colnames(df_meta)) {
  ggsave(filename = file.path(dir_snippet, 'coda_dimplot.celltypes.pdf'), 
         plot = coda_dimplot2(df_meta, 'celltypes', pal = 'circus'),
         width = 8, height = 5, useDingbats = F)  
}
ggsave(filename = file.path(dir_snippet, 'coda_dimplot.seurat_clusters.pdf'), 
       plot = coda_dimplot2(df_meta, 'seurat_clusters', pal = 'circus'),
       width = 8, height = 5, useDingbats = F)
ggsave(filename = file.path(dir_snippet, 'coda_dimplot.subclones.pdf'), 
       plot = coda_dimplot2(df_meta, 'subclones', pal = 'ironMan'),
       width = 8, height = 5, useDingbats = F)
ggsave(filename = file.path(dir_snippet, 'coda_dimplot.clones.pdf'), 
       plot = coda_dimplot2(df_meta, 'clones', pal = 'stallion', color_specific = c('D' = '#61d04f')),
       width = 8, height = 5, useDingbats = F)
ggsave(filename = file.path(dir_snippet, 'coda_dimplot.superclones.pdf'), 
       plot = coda_dimplot2(df_meta, 'superclones', pal = 'stallion'),
       width = 8, height = 5, useDingbats = F)

if (all(c('celltypes', 'clones') %in% colnames(df_meta))) {
  p1 <- coda_dimplot2(df_meta, 'celltypes', pal='circus', ncol=1)
  p2 <- coda_dimplot2(df_meta, 'clones', pal='stallion', ncol=1, 
                      color_specific = c('D' = '#61d04f'))
  p <- wrap_plots(p1, p2, ncol=2, widths = c(1,1))
  ggsave(filename = file.path(dir_snippet, 'coda_dimplot4.pdf'), 
         plot = p,
         width = 10, height = 8, useDingbats = F)
}
# stop()
z_viz_opts <- c('ploidy_class', 'clones', 'subclones', 'superclones', 
                'seurat_clusters', 'celltypes', 'run', 'ploidy_type_final')
z_viz_opts <- intersect(z_viz_opts, colnames(colData(objd)))
z_tofig_opts <- intersect(z_tofig_opts, colnames(colData(objd)))

pal_list <- list(`superclones`=pal_cna_superclones,
                 `subclones`=pal_cna_subclones,
                 `clones` = pal_cna_clones,
                 `seurat_clusters`=pal_atac_clusters,
                 `ploidy_class`=pal_ploidy_class,
                 `celltypes` = pal_celltypes,
                 `ploidy_type_final` = pal_ploidy_type_final,
                 `run` = pal_run)

#------ plotheatmap ratio ------
message('plotheatmap ratio')
for (z in z_tofig_opts){
  cat('\n', z, '...')
  if (!z %in% colnames(colData(objd))) {next('not present.')}
  
  objd <- calcConsensus(objd, consensus_by = z)
  n_z_opts <- ncol(consensus(objd))
  order_cells <- NULL
  if (n_z_opts > 1) {
    objd <- runConsensusPhylo(objd)
    order_cells <- 'consensus_tree'
  }
  pdf(file.path(dir_snippet, sprintf('cna.heatmap_consensus.combo_coda.%s.pdf', z)), 
      width = 15, height = 5, onefile = F, useDingbats = F)
  p <- try(plotHeatmap(objd, label = z, 
                       col = pal_copykit_ratio_clip1,
                       group = 'run', consensus=T, 
                       label_colors = pal_list[z]))
  if (!'try-error' %in% class(p)) {draw(p)}; rm(p)
  dev.off()
  
  pdf(file.path(dir_snippet, 
                paste0(sprintf('cna.heatmap.combo_coda.%s', z),  '.%03d.pdf')), 
      width = 15, height = 10, useDingbats = F, onefile = F)  
  p <- try( copykit::plotHeatmap(
    objd, 
    col = pal_copykit_ratio_clip1,
    order_cells = order_cells, 
    label=z_viz_opts, 
    label_colors = pal_list[z_viz_opts]))
  if (!'try-error' %in% class(p)) {draw(p)}; rm(p)
  
  p <- try( copykit::plotHeatmap(
    objd, 
    col = pal_copykit_ratio_clip1,
    order_cells = order_cells, 
    row_split = z,
    label=z_viz_opts, 
    label_colors = pal_list[z_viz_opts]))
  if (!'try-error' %in% class(p)) {draw(p)} ; rm(p)
  
  # 'ploidy_type_final' %in% colnames(colData(objd))
  # p <- try( copykit::plotHeatmap(
  #   objd, 
  #   col = pal_copykit_ratio_clip1,
  #   order_cells = order_cells,
  #   row_split = 'ploidy_type_final',
  #   label=z_viz_opts, 
  #   label_colors = pal_list[z_viz_opts]))
  # if (!'try-error' %in% class(p)) {draw(p)} ; rm(p)
  
  dev.off()
}; cat('\n')

if (ncol(objd) < 500) {
  # too slow and no need for big object
  pdf(file.path(dir_snippet, 'cna.heatmap_by_hclust.combo_coda.pdf'),
      width = 15, height = 10, useDingbats = F, onefile = F)
  p <- try( copykit::plotHeatmap(
    objd,
    col = pal_copykit_ratio_clip1,
    order_cells = 'hclust',
    label=z_viz_opts, n_threads = 40,
    label_colors = pal_list[z_viz_opts]) )
  if (!'try-error' %in% class(p)) {draw(p)} ; rm(p)
  dev.off()
}
#------ plotheatmap integer ------
message('plotheatmap integer')
if ('integer_scquantum' %in% assayNames(objd)) {
# if (T) {
  
  assay(objd, 'integer') <- assay(objd, 'integer_scquantum')
  for (z in z_tofig_opts){
    cat('\n', z, '...')
    
    objd <- calcConsensus(objd, assay = 'integer', 
                          consensus_by = z)
    
    n_z_opts <- ncol(consensus(objd))
    order_cells <- NULL
    if (n_z_opts > 1) {
      objd <- runConsensusPhylo(objd)
      order_cells <- 'consensus_tree'
    }
    
    pdf(file.path(dir_snippet, sprintf('cna.heatmap_consensus_int.combo_coda.%s.pdf', z)), 
        width = 15, height = 5, onefile = F, useDingbats = F)
    p <- try( draw( plotHeatmap(objd, label = z, 
                                assay = 'integer',
                                # col = pal_copykit_int_clip_0_2_6,
                                label_colors = pal_list[z],
                                group = 'run', consensus=T) ) )
    if (!'try-error' %in% class(p)) {draw(p)}; rm(p)
    dev.off()
    
    pdf(file.path(dir_snippet, 
                  paste0(sprintf('cna.heatmap_int.combo_coda.%s',z), '%03d.pdf')),
      width = 15, height = 10, useDingbats = F, onefile = F)
    p <- try( copykit::plotHeatmap(
      objd,
      assay = 'integer',
      # col = pal_copykit_int_clip_0_2_6,
      order_cells = order_cells,
      label=z_viz_opts,
      label_colors = pal_list[z_viz_opts]) )
    if (!'try-error' %in% class(p)) {draw(p)}; rm(p)
    
    p <- try( copykit::plotHeatmap(
      objd,
      assay = 'integer',
      # col = pal_copykit_int_clip_0_2_6,
      order_cells = order_cells,
      row_split = z,
      label=z_viz_opts,
      label_colors = pal_list[z_viz_opts]))
    if (!'try-error' %in% class(p)) {draw(p)}; rm(p)
    
    # p <- try( copykit::plotHeatmap(
    #   objd,
    #   assay = 'integer',
    #   # col = pal_copykit_int_clip_0_2_6,
    #   order_cells = order_cells,
    #   row_split = 'ploidy_type_final',
    #   label=z_viz_opts,
    #   label_colors = pal_list[z_viz_opts]))
    # if (!'try-error' %in% class(p)) {draw(p)}; rm(p)
  
    dev.off()
    
  }; cat('\n')
  
}

#------ plotheatmap fixed integer ------
message('plotheatmap fixed integer')
if ('integer_fixed' %in% assayNames(objd)) {
# if (T) {
  
  assay(objd, 'integer') <- assay(objd, 'integer_fixed')
  for (z in z_tofig_opts){
    cat('\n', z, '...')
    
    objd <- calcConsensus(objd, assay = 'integer', 
                          consensus_by = z)
    
    n_z_opts <- ncol(consensus(objd))
    order_cells <- NULL
    if (n_z_opts > 1) {
      objd <- runConsensusPhylo(objd)
      order_cells <- 'consensus_tree'
    }
    
    pdf(file.path(dir_snippet, sprintf('cna.heatmap_consensus_int_fixed.combo_coda.%s.pdf', z)), 
        width = 15, height = 5, onefile = F, useDingbats = F)
    p <- try( draw( plotHeatmap(objd, label = z, 
                                assay = 'integer',
                                # col = pal_copykit_int_clip_0_2_6,
                                label_colors = pal_list[z],
                                group = 'run', consensus=T) ) )
    if (!'try-error' %in% class(p)) {draw(p)}; rm(p)
    dev.off()
    
    pdf(file.path(dir_snippet, 
                  paste0(sprintf('cna.heatmap_int_fixed.combo_coda.%s',z), '%03d.pdf')),
        width = 15, height = 10, useDingbats = F, onefile = F)
    p <- try( copykit::plotHeatmap(
      objd,
      assay = 'integer',
      # col = pal_copykit_int_clip_0_2_6,
      order_cells = order_cells,
      label=z_viz_opts,
      label_colors = pal_list[z_viz_opts]) )
    if (!'try-error' %in% class(p)) {draw(p)}; rm(p)
    
    p <- try( copykit::plotHeatmap(
      objd,
      assay = 'integer',
      # col = pal_copykit_int_clip_0_2_6,
      order_cells = order_cells,
      row_split = z,
      label=z_viz_opts,
      label_colors = pal_list[z_viz_opts]))
    if (!'try-error' %in% class(p)) {draw(p)}; rm(p)
    
    # p <- try( copykit::plotHeatmap(
    #   objd,
    #   assay = 'integer',
    #   # col = pal_copykit_int_clip_0_2_6,
    #   order_cells = order_cells,
    #   row_split = 'ploidy_type_final',
    #   label=z_viz_opts,
    #   label_colors = pal_list[z_viz_opts]))
    # if (!'try-error' %in% class(p)) {draw(p)}; rm(p)
    
    dev.off()
    
  }; cat('\n')
  
}
message('done viz')
