## Visualization purpose only!

message(dir_snippet)
# print(head(df_meta))
print(objd)
# z_tofig_opts <- c('clones', 'celltypes')
print(z_tofig_opts)

pal_cna_clones = new_palette_D(levels(df_meta$clones), pal = 'stallion')
pal_celltypes = new_palette_D(levels(df_meta$celltypes), pal = 'circus')
# pal_celltypes <- PAL_DCIS_CELLTYPES ## adhoc !!!
pal_run = new_palette_D(unique(df_meta$run), pal = 'scales_hue') # consistent with copykit
if ('D' %in% names(pal_cna_clones)) {pal_cna_clones['D'] <- 'darkseagreen1'}

## helpful for seeing `run` on the heatmap
if (F) {
  ## Careful: here adds adhoc cell name prefix so
  ## never use the object after this to avoid overwriting original objects. 
  adhoc_cellprefix <- as.character(sample(x=1:ncol(objd), size=ncol(objd), replace = F))
  rownames(df_meta) <- paste0(adhoc_cellprefix, rownames(df_meta))
  obja <- RenameCells(obja, new.names = paste0(adhoc_cellprefix, Cells(obja)))
  colnames(objd) <- paste0(adhoc_cellprefix, colnames(objd))
}

#------ confusion matrix ------

message('confusionmat')
for (y in c('clones')) {
  if (!y %in% colnames(df_meta)) {next()}
  
  confusionmat <- table(df_meta[[y]],
                        df_meta[['celltypes']]) %>%
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

ggsave(filename = file.path(dir_snippet, 'coda_dimplot.clones.pdf'), 
       plot = coda_dimplot2(df_meta, 'clones', pal = 'stallion', color_specific = c('D' = 'lightsteelblue')),
       width = 8, height = 5, useDingbats = F)

if (all(c('celltypes', 'clones') %in% colnames(df_meta))) {
  p1 <- coda_dimplot2(df_meta, 'celltypes', pal='circus', ncol=1)
  p2 <- coda_dimplot2(df_meta, 'clones', pal='stallion', ncol=1, 
                      color_specific = c('D' = 'lightsteelblue'))
  p <- wrap_plots(p1, p2, ncol=2, widths = c(1,1))
  ggsave(filename = file.path(dir_snippet, 'coda_dimplot4.pdf'), 
         plot = p,
         width = 10, height = 8, useDingbats = F)
}


#------------------- ~~~ heatmap ~~~ -------------------  
para.plotHeatmap.group <- 'run'
if (length(unique(objd@colData[['run']])) == 1) {
  para.plotHeatmap.group = NULL
} 
z_viz_opts <- c('clones', 'celltypes')
z_viz_opts <- intersect(z_viz_opts, colnames(colData(objd)))
z_tofig_opts <- intersect(z_tofig_opts, colnames(colData(objd)))

# pal_list <- list(`superclones`=pal_cna_superclones,
#                  `subclones`=pal_cna_subclones,
#                  `clones` = pal_cna_clones,
#                  `seurat_clusters`=pal_atac_clusters,
#                  `ploidy_class`=pal_ploidy_class,
#                  `celltypes` = pal_celltypes,
#                  `ploidy_type_final` = pal_ploidy_type_final,
#                  `run` = pal_run)
pal_list <- list(`clones` = pal_cna_clones,
                 `celltypes` = pal_celltypes,
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
                       group = para.plotHeatmap.group, consensus=T, 
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
  
  # p <- try( copykit::plotHeatmap(
  #   objd, 
  #   col = pal_copykit_ratio_clip1,
  #   order_cells = order_cells, 
  #   row_split = z,
  #   label=z_viz_opts, 
  #   label_colors = pal_list[z_viz_opts]))
  # if (!'try-error' %in% class(p)) {draw(p)} ; rm(p)
  
  dev.off()
}; cat('\n')

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
                                group = para.plotHeatmap.group, consensus=T) ) )
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
    
    # p <- try( copykit::plotHeatmap(
    #   objd,
    #   assay = 'integer',
    #   # col = pal_copykit_int_clip_0_2_6,
    #   order_cells = order_cells,
    #   row_split = z,
    #   label=z_viz_opts,
    #   label_colors = pal_list[z_viz_opts]))
    # if (!'try-error' %in% class(p)) {draw(p)}; rm(p)
    
    dev.off()
    
  }; cat('\n')
  
}


#------ plotheatmap the fixed integer ------
message('plotheatmap the fixed integer')
# if ('integer_fixed' %in% assayNames(objd)) {
if (T) {
  
  try(assay(objd, 'integer') <- assay(objd, 'integer_fixed'))
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
                                group = para.plotHeatmap.group, consensus=T) ) )
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
    
    # p <- try( copykit::plotHeatmap(
    #   objd,
    #   assay = 'integer',
    #   # col = pal_copykit_int_clip_0_2_6,
    #   order_cells = order_cells,
    #   row_split = z,
    #   label=z_viz_opts,
    #   label_colors = pal_list[z_viz_opts]))
    # if (!'try-error' %in% class(p)) {draw(p)}; rm(p)
    
    dev.off()
    
  }; cat('\n')
  
}
message('done viz')