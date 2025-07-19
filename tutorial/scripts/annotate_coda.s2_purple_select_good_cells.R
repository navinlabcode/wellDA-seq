#---------------------------
# good cells pass_scquantum_QC=TRUE, clones!='Z'               ----    
#---------------------------
objd[['pass_scquantum_QC']] <- T
if (sample_name %in% c('DCIS28T', 'DCIS35T', 'DCIS67T')) {
  objd[['pass_scquantum_QC']] <- objd[['pass_scquantum_sc']]
}
if (sample_name %in% c('BCMDCIS69')) {
  objd[['pass_scquantum_QC']] <- T
}

is_good_aneu <- objd[['pass_scquantum_QC']]==TRUE & objd[['clones']] != 'Z'; print(table(good_aneuploid=is_good_aneu))
colData(objd)$is_good_aneu <- as.factor(is_good_aneu)

pal_list <- list(
  'subclones' = new_palette_D(levels(objd[['subclones']]), 'ironMan'), 
  'pass_scquantum' = c(`TRUE`='steelblue1', `FALSE`='plum1'), 
  'pass_scquantum_sc' = c(`TRUE`='steelblue1', `FALSE`='plum1'), 
  'pass_scquantum_QC' = c(`TRUE`='steelblue1', `FALSE`='plum1'), 
  'seurat_clusters' = new_palette_D(levels(objd[['seurat_clusters']]),'circus'), 
  'is_good_aneu' = new_palette_D(levels(objd[['is_good_aneu']]), 'scales_hue'))

pdf(file = file.path(dir_res, 'plotHeatmap.pass_QC.%03d.pdf'), onefile = F, useDingbats = F, 
    width = 15, height = 10)
z = 'subclones'; z_viz_opts <- c('subclones', 'pass_scquantum', 'pass_scquantum_sc', 'pass_scquantum_QC', 'seurat_clusters', 'is_good_aneu')
p <- try( copykit::plotHeatmap(
  objd, 
  col = pal_copykit_ratio_clip1,
  order_cells = 'consensus_tree', 
  row_split = z,
  label=z_viz_opts,
  label_colors = pal_list[z_viz_opts]))
if (!'try-error' %in% class(p)) {draw(p)}
# p <- try( copykit::plotHeatmap(
#   objd,
#   assay = 'integer',
#   order_cells = 'consensus_tree',
#   row_split = z,
#   label=z_viz_opts,
#   label_colors = pal_list[z_viz_opts]))
# if (!'try-error' %in% class(p)) {draw(p)}
dev.off()