#---------------------------
# Process the aneuploid cells compartments: 
# 1. [here] identify the C, DC:noise (diploid-like with random events), and Z:doublets (carry aneuploid pattern but lighter) populations
# 2. [next] merge similar subclones
# 
# The C group is ready for subclones
# The Z group is removed from analysis
# The DC group is further filtering
#---------------------------

message(sample_name)
print(objd)
print(dir_res)

library(magrittr)
#------ prepare assay ------
  
assay(objd, 'integer') <- assay(objd, 'integer_fixed')
objd <- calcConsensus(objd, consensus_by = 'subclones', assay = 'integer')
objd <- runConsensusPhylo(objd)

#------ prepare colors ------
pal_list <- list(
  'subclones' = new_palette_D(levels(objd[['subclones']]), 'ironMan'), 
  'pass_scquantum' = c(`TRUE`='steelblue1', `FALSE`='plum1'), 
  'pass_scquantum_sc' = c(`TRUE`='steelblue1', `FALSE`='plum1'), 
  'pass_scquantum_QC' = c(`TRUE`='steelblue1', `FALSE`='plum1'), 
  'seurat_clusters' = new_palette_D(levels(objd[['seurat_clusters']]),'circus'))
#------------------- ~~~ Identify C/DC/Z  by scquantum ~~~ -------------------  
library(ggpubr)
# ploidy_score (=abs(1-confidence_ratio)), confidence_ratio
# ploidy, ploidy_confidence (confidence_ratio)
#---------------------------
# Compare which scquantum-based strategy works better:               ----    
# simply filtering at single-cell level vs at subclone level
# Conclusion: It seemed subclone level works better. But I decided to remove 
# any cell if it cannot pass either one.  
#---------------------------
if (T) {
  subclones_qc_scquantum <- colData(objd)[, c('subclones', 'ploidy', 'ploidy_score')] %>%
    as.data.frame() %>%
    dplyr::group_by(subclones) %>%
    dplyr::summarise(
      avg_ploidy = median(ploidy), 
      avg_ploidy_score = median(ploidy_score), 
      pass_ploidy_pct = sum(ploidy >= 2)/length(ploidy), 
      pass_confidence_ratio_pct = sum(ploidy_score <= 0.05)/length(ploidy_score))
  subclones_qc_scquantum %<>% 
    dplyr::mutate(
      pass_avg_ploidy_flag = avg_ploidy >=2, 
      pass_avg_ploidy_score_flag = avg_ploidy_score <= 0.05, 
      pass_ploidy_pct_flag = pass_ploidy_pct >= 50/100, 
      pass_confidence_ratio_pct_flag = pass_confidence_ratio_pct >= 50/100
    )
  # subclones_qc_scquantum %<>% 
  #   dplyr::mutate(pass_scquantum = pass_ploidy_pct_flag & pass_confidence_ratio_pct_flag)
  subclones_qc_scquantum %<>% 
    dplyr::mutate(pass_scquantum = pass_avg_ploidy_flag & pass_avg_ploidy_score_flag)
  # subclones_qc_scquantum %<>% 
  #   dplyr::mutate(pass_scquantum = pass_avg_ploidy_flag)
  colData(objd)[['pass_scquantum']] <- deframe(subclones_qc_scquantum[, c('subclones', 'pass_scquantum')])[colData(objd)[['subclones']]]
  # view(subclones_qc_scquantum)
  table(colData(objd)[['pass_scquantum']])
}

if (T) {
  objd[['pass_scquantum_sc']] <- objd[['ploidy']] >= 2 & objd[['ploidy_score']] <= 0.05
  pass_scquantum_pct_per_subclone <- colData(objd) %>% as.data.frame() %>%
    dplyr::group_by(subclones) %>%
    dplyr::summarise(pass_scquantum_sc_pct = sum(pass_scquantum_sc)/length(pass_scquantum_sc)) %>%
    deframe()
  barplot(pass_scquantum_pct_per_subclone)
  print(table(objd[['pass_scquantum_sc']]))
}

library(ggpubr)


table(objd[['pass_scquantum']])
objd[['pass_scquantum_QC']] <- objd[['pass_scquantum']] & objd[['pass_scquantum_sc']]

table(objd[['pass_scquantum_QC']])

for (z in c('clones', 'subclones')) {
  p <- ggscatter(as.data.frame(colData(objd)),
                 x='ploidy', y='confidence_ratio', color='pass_scquantum_QC',
                 facet.by= z) +
    geom_hline(yintercept = c(0.95, 1.05), lty='dashed', color='blue')  +
    geom_vline(xintercept = 2, lty='dashed', color='blue') 
  print(p)
  ggsave(file.path(dir_res, sprintf('qc.scquantum.ploidy_by_%s.pdf', z)), p, width = 10,height = 10, useDingbats = F)
}
# stop()
pdf(file = file.path(dir_res, 'plotHeatmap.pass_scquantum.%03d.pdf'), onefile = F, useDingbats = F, 
    width = 15, height = 10)
z = 'subclones'; z_viz_opts <- c('subclones', 'pass_scquantum', 'pass_scquantum_sc', 'pass_scquantum_QC', 'seurat_clusters')
p <- try( copykit::plotHeatmap(
  objd, 
  col = pal_copykit_ratio_clip1,
  order_cells = 'consensus_tree', 
  row_split = z,
  label=z_viz_opts,
  label_colors = pal_list[z_viz_opts]))
if (!'try-error' %in% class(p)) {draw(p)}
p <- try( copykit::plotHeatmap(
  objd,
  assay = 'integer',
  # col = pal_copykit_int_clip_0_2_6,
  order_cells = 'consensus_tree',
  row_split = z,
  label=z_viz_opts,
  label_colors = pal_list[z_viz_opts]))
if (!'try-error' %in% class(p)) {draw(p)}

draw( plotHeatmap(objd, label = 'subclones', 
                  label_colors = pal_list['subclones'], 
                  col = pal_copykit_ratio_clip1,
                  group = 'pass_scquantum_QC', consensus=T) )
dev.off()
