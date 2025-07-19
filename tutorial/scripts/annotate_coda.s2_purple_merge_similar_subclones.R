#---------------------------
# Process the aneuploid cells compartments: 
# 1. [before] identify the C, DC:noise (diploid-like with random events), and Z:doublets (carry aneuploid pattern but lighter) populations
# 2. [here] merge similar subclones
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
  
#------ prepare colors ------
pal_list <- list(
  'subclones' = new_palette_D(levels(objd[['subclones']]), 'ironMan'), 
  'seurat_clusters' = new_palette_D(levels(objd[['seurat_clusters']]),'circus'), 
  'celltypes' = new_palette_D(levels(objd[['celltypes']]),'circus'))

#------------------- ~~~ Merge similar subclones ~~~ -------------------  
n_bin_diff_tolerate <- 3 # requires subclones showing >3 bins diff
#------ Specify parameters ------
# para for merging similar clones
merge_to_k = NULL

if (sample_name == 'DCIS22T') {merge_to_k = 4}
if (sample_name == 'DCIS28T') {merge_to_k = 6}
if (sample_name == 'DCIS35T') {merge_to_k = 2}
if (sample_name == 'DCIS36T') {merge_to_k = 8}
if (sample_name == 'DCIS41T') {merge_to_k = 15}
if (sample_name == 'DCIS51T') {merge_to_k = 3}
if (sample_name == 'DCIS68T') {merge_to_k = 7}
if (sample_name == 'DCIS67T') {merge_to_k = 11}
if (sample_name == 'BCMDCIS69') {merge_to_k = 7}
if (sample_name == 'DCIS66T_chip2') {merge_to_k = 4}
message(merge_to_k)
#------ Identify the C, DC, and z groups ------
merge_for_z = 'subclones'

try( assay(objd, 'integer') <- assay(objd, 'integer_fixed') )

pdf(file.path(dir_res, sprintf('plotHeatmap.%s_before_merge.pdf', merge_for_z)), 
    width = 12, height = 5, onefile = T, useDingbats = F)

objd <- calcConsensus(objd, consensus_by = merge_for_z, assay = 'integer')
objd <- runConsensusPhylo(objd)
draw( plotHeatmap(objd, label = merge_for_z, 
                  # group = 'run', 
                  consensus=T, label_colors = pal_list[merge_for_z], 
                  assay = 'integer') )

objd <- calcConsensus(objd, consensus_by = merge_for_z, assay = 'segment_ratios')
objd <- runConsensusPhylo(objd)
draw( plotHeatmap(objd, label = merge_for_z, 
                  # group = 'run', 
                  consensus=T, label_colors =  pal_list[merge_for_z], 
                  col = pal_copykit_ratio_clip1,
                  assay = 'segment_ratios') )
dev.off()

objd <- calcConsensus(objd, consensus_by = merge_for_z, assay = 'segment_ratios')
objd <- runConsensusPhylo(objd)
ck_css <- copykit::consensus(objd); print(range(ck_css))
if (!is.integer(ck_css)) {
  # using ratio
  gain_cutoff <- 3/2; loss_cutoff <- 1/2
  if (max(ck_css) <= gain_cutoff) {gain_cutoff <- 1.01}
  if (min(ck_css) >= loss_cutoff) {loss_cutoff <- 0.99}
  ck_css_gains <- 1 * (ck_css >= gain_cutoff)
  ck_css_loses <- -1 * (ck_css <= loss_cutoff)
} else {
  # using integer
  gain_cutoff <- 3; loss_cutoff <- 1
  if (max(ck_css) <= gain_cutoff) {gain_cutoff <- 2.01}
  if (min(ck_css) >= loss_cutoff) {loss_cutoff <- 1.99}
  ck_css_gains <- ck_css >= 3
  ck_css_loses <- ck_css <= 1
}

ck_css_events <- abs(ck_css_gains) + abs(ck_css_loses)

#------ which subclones are diploid ------ [deprecated]
# n_events_per_subclone <- colSums(abs(ck_css_events))
# n_max_events_diploid <- pmin(quantile(n_events_per_subclone, .75), 10) # <=#events are considered diploid subclone
# subclone_is_diploid <- n_events_per_subclone <= n_max_events_diploid
# subclone_is_diploid[!subclone_is_diploid]
#------ which subclones are similar based on manhattan distance ------
# sum(abs(ck_css_events[, 'c37'] - ck_css_events[, 'c35'])) # manhattan distance  
ident_css_dist <- dist(t(ck_css_events), method = 'manhattan')
ident_css_dist[ident_css_dist<=n_bin_diff_tolerate] <- 0
merge_dict <- dendextend::cutree(
  hclust(ident_css_dist, 'ward.D2'), k=merge_to_k, 
  order_clusters_as_data=F)
merge_dict <- merge_dict[colnames(ck_css_events)]
print(sort(merge_dict))
pdf(file.path(dir_res, sprintf('diagnose_merge_similar_%s.pdf', merge_for_z)), 
    width = 12, height = 12, onefile = F, useDingbats = F)
draw(Heatmap(as.matrix(ident_css_dist), name = 'distance',
             row_split = merge_dict, column_split = merge_dict, 
             cell_fun = function(j, i, x, y, width, height, fill) {
               if (as.matrix(ident_css_dist)[i, j] <= 20) {
                 grid.text(sprintf("%.0f", as.matrix(ident_css_dist)[i, j]), x, y, gp = gpar(fontsize = 10))
               }
             }
))
dev.off()

#------ which subclones are similar based on #bins showing diff events ------
# // this is to remove the noisy small events mostly located at chr-chr bourndry
get_bins_num_of_max_diff_event <- function(x, y) {
  # x: profile X
  # y: profile Y
  d <- x - y
  o <- rle(d)
  i <- o$values!=0
  if (sum(i) > 0) {
    return( max(o$lengths[i]) )
  } else {
    return(0)
  }
}

dim(ck_css_events)
ident_css_bin_dist <- proxy::dist(
  t(ck_css_events), method = get_bins_num_of_max_diff_event)
range(ident_css_bin_dist)

pdf(file.path(dir_res, sprintf('diagnose_merge_similar_%s.nbins.pdf', merge_for_z)), 
    width = 12, height = 12, onefile = F, useDingbats = F)
draw(Heatmap(as.matrix(ident_css_bin_dist), name = 'len of seg',
             # row_split = merge_dict, column_split = merge_dict, 
             cell_fun = function(j, i, x, y, width, height, fill) {
               if (as.matrix(ident_css_bin_dist)[i, j] <= 5) {
                 grid.text(sprintf("%.0f", as.matrix(ident_css_bin_dist)[i, j]), x, y, gp = gpar(fontsize = 10, col='white'))
               } else if (as.matrix(ident_css_bin_dist)[i, j] <= 20) {
                 grid.text(sprintf("%.0f", as.matrix(ident_css_bin_dist)[i, j]), x, y, gp = gpar(fontsize = 10))
               } else {}
             }
))
dev.off()

#------ propose subclones identification ------
f_newclones <- file.path(dir_res, sprintf('cna.manual_merge_%s.csv', merge_for_z))
write_csv(x = enframe(merge_dict, name = 'old', value = 'new') %>% arrange(new), file=f_newclones)
write_csv(x = data.frame(old = colnames(ck_css_events), new=colnames(ck_css_events)), file=f_newclones)
write_csv(x = data.frame(unique(objd@colData[, c('subclones', 'clones')])), file=f_newclones)
stop(
  'Manually merge the subclones.'
)

if (T) {
  new_clones <- read_csv(f_newclones, col_types = c('c', 'c'))
  new_clones <- deframe(new_clones[, c('new', 'old')])
} else {
  new_clones <- structure(
    names(merge_dict),
    names=paste0(clone_prefix, merge_dict))
}

colData(objd)[, paste0(merge_for_z, '_raw')] <- colData(objd)[, merge_for_z]
colData(objd)[, 'clones'] <- forcats::fct_recode(
  colData(objd)[[paste0(merge_for_z, '_raw')]], !!!new_clones)
colData(objd)[, 'clones'] <- forcats::fct_relevel(
  colData(objd)[, 'clones'], 
  mixedsort(unique(names(new_clones)))
)
table(colData(objd)[, 'clones'])
if (any(c('X', 'Z') %in% colData(objd)[['clones']])) {
  tmp <- colData(objd)[['clones']] %in% c('X', 'Z'); sum(tmp)
  objd <- objd[, !tmp]
  subset_coda(df_meta, obja, objd, colnames(objd))
  objd[['clones']] <- fct_drop(objd[['clones']])
  update_coda(df_meta, obja, objd)
}
table(colData(objd)[, 'clones'])
colData(objd)[, 'clones'] <- pretty_seamless_factor(
  colData(objd)[, 'clones'], prefix='C')
table(colData(objd)[, 'clones'])

