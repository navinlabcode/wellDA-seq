#---------------------------
# Process the diploid compartments: 
# 1. identify the D, CD (progenitor or non-cancer), and X (noise)
# 2. merge similar subclones
#---------------------------
message(sample_name)
print(objd)

n_bin_diff_tolerate <- 3 # requires subclones showing >3 bins diff
#------ Specify parameters ------
merge_to_k = NULL
if (sample_name == 'DCIS22T') {merge_to_k = 4}
if (sample_name == 'DCIS28T') {merge_to_k = 4}
if (sample_name == 'DCIS35T') {merge_to_k = 2}
if (sample_name == 'DCIS36T') {merge_to_k = 6}
if (sample_name == 'DCIS41T') {merge_to_k = 2}
if (sample_name == 'DCIS51T') {merge_to_k = 2}
if (sample_name == 'DCIS68T') {merge_to_k = 5}
if (sample_name == 'DCIS67T') {merge_to_k = 4}
if (sample_name == 'BCMDCIS69') {merge_to_k = 4}
if (sample_name == 'DCIS66T_chip2') {merge_to_k = 4}
#------ sub-cluster if needed ------
subclone_to_resubcluster <- NULL 
if (sample_name == c('DCIS41T') ) {subclone_to_resubcluster <- 'c1'; k_nb = 3}
if (!is.null(subclone_to_resubcluster)) {
  objd <- copykit_find_sub_clusters(
    objd, subcluster_by = 'subclones', cluster = subclone_to_resubcluster, 
    saveto_col_str = 'newsubclusters', k_nb = k_nb)
  table(objd[['newsubclusters']])
  plotUmap(objd, label = 'newsubclusters')
}
if ('newsubclusters' %in% colnames(colData(objd))) {
  objd[['subclones']] <- objd[['newsubclusters']]
  objd[['newsubclusters']] <- NULL
  objd[['subclones']] <- pretty_seamless_factor(
    objd[['subclones']], prefix='c')
  table(objd[['subclones']])
  update_coda(df_meta, obja, objd)
}


message(merge_to_k)
#------ Identify the D, A, and A/D groups ------
merge_for_z = 'subclones'
pal_cna_merge_z <- new_palette_D(levels(df_meta[, merge_for_z]), 'ironMan')
pal_cna_merge_z <- list(pal_cna_merge_z); names(pal_cna_merge_z) <- merge_for_z
assay(objd, 'integer') <- assay(objd, 'integer_fixed')
print(table(objd[[merge_for_z]]))


pdf(file.path(dir_res, sprintf('plotHeatmap.%s_before_merge.pdf', merge_for_z)), 
    width = 12, height = 5, onefile = T, useDingbats = F)

objd <- calcConsensus(objd, consensus_by = merge_for_z, assay = 'integer')
objd <- runConsensusPhylo(objd)
draw( plotHeatmap(objd, label = merge_for_z, 
                  # group = 'run', 
                  consensus=T, label_colors = pal_cna_merge_z, 
                  assay = 'integer') )

objd <- calcConsensus(objd, consensus_by = merge_for_z, assay = 'segment_ratios')
objd <- runConsensusPhylo(objd)
draw( plotHeatmap(objd, label = merge_for_z, 
                  # group = 'run', 
                  consensus=T, label_colors = pal_cna_merge_z, 
                  col = pal_copykit_ratio_clip1,
                  assay = 'segment_ratios') )
dev.off()


objd <- calcConsensus(objd, consensus_by = merge_for_z, assay = 'integer')
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
print(dim(ck_css_events))
#------ which subclones are diploid ------ [deprecated]
# n_events_per_subclone <- colSums(abs(ck_css_events))
# n_max_events_diploid <- pmin(quantile(n_events_per_subclone, .75), 10) # <=#events are considered diploid subclone
# subclone_is_diploid <- n_events_per_subclone <= n_max_events_diploid
# subclone_is_diploid[!subclone_is_diploid]
#------ which subclones are similar ------
# sum(abs(ck_css_events[, 'c37'] - ck_css_events[, 'c35'])) # manhattan distance  
ident_css_dist <- dist(t(ck_css_events), method = 'manhattan')
ident_css_dist[ident_css_dist<=n_bin_diff_tolerate] <- 0
merge_dict <- dendextend::cutree(
  hclust(ident_css_dist, 'average'), k=merge_to_k, 
  order_clusters_as_data=F)
merge_dict <- merge_dict[colnames(ck_css_events)]
print(sort(merge_dict))

pdf(file.path(dir_res, sprintf('diagnose_merge_similar_%s.pdf', merge_for_z)), 
    width = 12, height = 12, onefile = F, useDingbats = F)
draw(Heatmap(as.matrix(ident_css_dist), name = 'distance',
             row_split = merge_dict, column_split = merge_dict))
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
message(
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
print( table( colData(objd)[, paste0(merge_for_z, '_raw')] ))
colData(objd)[, 'clones'] <- forcats::fct_recode(
  colData(objd)[[paste0(merge_for_z, '_raw')]], !!!new_clones)
colData(objd)[, 'clones'] <- forcats::fct_relevel(
  colData(objd)[, 'clones'], 
  mixedsort(unique(names(new_clones)))
)
table(colData(objd)[, 'clones'])
# colData(objd)[, 'clones'] <- pretty_seamless_factor(colData(objd)[, 'clones'], prefix='c')
# table(colData(objd)[, 'clones'])

