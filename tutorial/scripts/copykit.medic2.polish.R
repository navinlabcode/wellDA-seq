#---------------------------
# Given a copykit object, yield the MEDIC2 tree.              ----    
#---------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(tidyverse); library(forcats)
  library(fs)
  library(copykit)
  library(VennDiagram)
  library(SummarizedExperiment)
  library(ggpubr); library(ruok); library(scales)
  library(ComplexHeatmap)
  library(ggtree)
  library(forcats)
})
setwd("/volumes/USR1/yyan/project/coda")
source('/volumes/USR1/yyan/project/coda/rsrc/utils.R')
source('/volumes/USR1/yyan/project/coda/rsrc/color_pal.R')

options <- commandArgs(trailingOnly = TRUE)
if (length(options)>0) {
  cat('Reading user parameters:\n')
  sample_name <- options[[1]]
} else {
  sample_name = 'DCIS66T_chip2'
}
message(sample_name)

f_objd <- file.path('/volumes/USR1/yyan/project/coda', 'rds_coda_ready', 
                    sample_name, 'aneuploid_epi', 'objd.rds')
if (F) {
  f_objd <- file.path('/volumes/USR1/yyan/project/coda', 'rds_coda_ready', 
                      'mda231wt_kaile', 'objd.rds')
}
objd <- read_rds(f_objd); print(objd)

dir_proj <- dirname(f_objd)
dir_res <- file.path(dir_proj, 'medic2_polish'); fs::dir_create(dir_res)


# *** MEDIC2 workflow is shared from Kaile Wang ***
# objd <- calcConsensus(objd, assay = 'integer_scquantum', 
#                       consensus_by = 'clones')
try(assay(objd, 'integer') <- assay(objd, 'integer_fixed'))
objd <- calcConsensus(objd, assay = 'integer', 
                      consensus_by = 'clones')
css_mat <- objd@consensus; dim(css_mat); print(range(css_mat))
head(rowRanges(objd))

#------ Identify the noisy segments ------
cna_neutral <- 2
bin_size_cutoff <- 10
# if (sample_name == 'DCIS35T') {bin_size_cutoff <- 20}
# clone_focus <- 'C7'
dict_cna_fix <- lapply(colnames(css_mat), function(clone_focus) {
  rle_o <- rle( css_mat[, clone_focus] )
  which(rle_o$lengths <= 10)
  rle_seg_end <- cumsum(rle_o$lengths)
  rle_seg_start <- c(1, head(rle_seg_end, length(rle_seg_end)-1)+1)
  cna_fix <- data.frame(s=rle_seg_start, e=rle_seg_end, cna=rle_o$values, len=rle_o$lengths)
  cna_fix$new_cna <- cna_fix$cna
  cna_fix$new_cna[cna_fix$len <= bin_size_cutoff] <- cna_neutral
  cna_fix$clones <- clone_focus
  return(cna_fix)
}); names(dict_cna_fix) <- colnames(css_mat)
# dict_cna_fix <- do.call('rbind', dict_cna_fix)

#------ Polish the consensus CNAs ------
for (clone_focus in names(dict_cna_fix)) {
  message(clone_focus)
  cna_fix <- dict_cna_fix[[clone_focus]]
  cna_fix <- cna_fix[cna_fix$cna!=cna_fix$new_cna, , drop=F]
  for (i in seq_len(nrow(cna_fix))) {
    css_mat[cna_fix$s[i]:cna_fix$e[i], clone_focus] <- cna_fix$new_cna[i]
  }
}
print(range(css_mat)); print(range(objd@consensus))
identical(objd@consensus, css_mat)
#------ make medicc input matrix ------

dummy_diploid <- rep.int(2, times = nrow(objd))
medicc_input <- cbind(
  objd@rowRanges %>% 
    dplyr::as_tibble() %>% 
    dplyr::select(seqnames, start, end), 
  css_mat
) %>%
  # mutate(diploid=dummy_diploid) %>%  ## medicc2 by default always has a 'diploid' name; assign diploid cell CN profile to be any other profiles as the root
  dplyr::rename(chrom=seqnames) %>% 
  gather('sample_id', 'CN', -chrom, -start, -end) %>% 
  dplyr::select(sample_id,chrom, everything())

#------ polish the edge of each chrom ------

medicc_input <- medicc_input %>%
  dplyr::group_by(sample_id, chrom) %>%
  dplyr::mutate(
    CN = case_when(
      row_number() %in% 1:bin_size_cutoff ~ CN[bin_size_cutoff+1],  # For the first 10 rows, use the CN value of the 11th row
      row_number() %in% (n() - 1):n() ~ CN[n() - bin_size_cutoff],  # For the last 10 rows, use the CN value of the 11th row from the bottom
      TRUE ~ CN  # For all other rows, keep the original CN value
    )
  )

# subclone chr start end copy_number_int
# sample_id	chrom	start	end	CN
# C1	chr1	977837	1200863	2
# C1	chr1	1200864	1455238	2
# C1	chr1	1455239	1758057	2
unique(medicc_input$sample_id)
write_tsv(medicc_input, file = file.path(dir_res, "medicc2_input.tsv"))
# medicc_input <- read_tsv(file = paste0("metrics/medicc_files/", pro_name_d, "/",pro_name_d, c("_medicc2_input.tsv")))

#---- run medicc2 ---
# conda activate medicc_env
# medicc2 -a CN --total-copy-numbers -j 40 -vv input_path/*_medicc2_input.tsv output_path
# if (!file.exists(file.path(dir_res, 'medicc2_input_final_tree.new'))) {
if (T) {
  medicc2_cmd <- paste0(
    '/volumes/USR1/yyan/anaconda3/envs/medicc2/bin/medicc2 ', 
    '-a CN ', 
    '--total-copy-numbers -j 40 -vv ', 
    file.path(dir_res, "medicc2_input.tsv"), ' ', 
    dir_res
  )
  cat(medicc2_cmd)
  system(medicc2_cmd)
} else {
  cat('medicc2 has been performed.\n')
}
#------ viz ------
pal_cna_clones = new_palette_D(levels(objd@colData$clones), pal = 'stallion')
size_clones <- c(table(fct_drop(objd@colData$clones)))
prop_clones <- c(prop.table(table(fct_drop(objd@colData$clones))))

# ggtree_misc_meta <- data.frame(
#   Newick_label = c(names(size_clones), 'diploid'),
#   n = c(size_clones, 0),
#   prop=c(prop_clones, 0))

ggtree_misc_meta <- data.frame(
  Newick_label = c(names(size_clones)),
  n = c(size_clones),
  prop=c(prop_clones))

#------ read medicc2 ------

medic <- ape::read.tree(
  file.path(dir_res, 'medicc2_input_final_tree.new'))

my_sub <- factor(levels(objd@colData$clones), 
                 levels = levels(objd@colData$clones))
list_samples <- split(my_sub, my_sub)

tree <- ggtree::groupOTU(medic, list_samples)
treeplt <- ggtree::ggtree(tree, ladderize=F, size = 0.2)
treeplt <- treeplt %<+% ggtree_misc_meta
treeplt <- treeplt + 
  ggtree::geom_tiplab(size=3, aes(color=group),hjust = -0.4, alpha=1)+
  ggtree::geom_tippoint(aes(color=group, size = prop), alpha=1)+
  scale_colour_manual(values = pal_cna_clones) +
  geom_text(aes(x=branch, 
                label=plyr::mapvalues(round(branch.length, 0),from = "0",to = ""), 
                vjust=-.5), size = 3) +
  theme(legend.position = "none") + 
  ggtree::geom_rootpoint() 
treeplt <- treeplt + ggtree::theme_tree() + guides(color = guide_legend(override.aes = list(size = 3)))


cowplot::ggsave2(file.path(dir_res, 'ggtree.medicc2.pdf'), treeplt, 
                 width = 5, height = 4)
cowplot::ggsave2(
  file.path(dir_res, 'paper.ggtree.medicc2.pdf'), 
  treeplt + scale_colour_manual(values = pal_cna_clones_paper), 
  width = 5, height = 4)

#------ viz medicc2 tree by ME ------
# MEDIC2 distance -> ME tree -> delete diploid
medic_dist <- read_tsv(
  file.path(dir_res, 'medicc2_input_pairwise_distances.tsv'))

medic_dist <- medic_dist %>% as.data.frame() %>% column_to_rownames(var="sample_id") %>% as.dist()
me_tree <- ape::fastme.bal(medic_dist)
me_tree <- ape::root.phylo(me_tree, outgroup = "diploid", resolve.root = TRUE)
# me_tree <- ape::drop.tip(me_tree, tip = 'diploid')

my_sub <- factor(levels(objd@colData$clones), 
                 levels = levels(objd@colData$clones))
# my_sub <- my_sub[!(my_sub %in% other_clst)] %>% droplevels()
list_samples <- split(my_sub, my_sub)

tree <- ggtree::groupOTU(me_tree, list_samples)
treeplt <- ggtree::ggtree(tree, ladderize=F, size = 0.2)
treeplt <- treeplt %<+% ggtree_misc_meta
 
treeplt <- treeplt + 
  ggtree::geom_tiplab(size=3, aes(color=group),hjust = -0.4, alpha=1)+
  ggtree::geom_tippoint(aes(color=group, size = prop), alpha=1)+
  scale_colour_manual(values = pal_cna_clones) +
  geom_text(aes(x=branch,
                label=plyr::mapvalues(round(branch.length, 0), from = "0",to = ""),
                vjust=-.5), size = 3) +
  theme(legend.position = "none") +
  ggtree::geom_rootpoint() +
  ggtree::theme_tree2() + 
  guides(color = guide_legend(override.aes = list(size = 3)))

treeplt
cowplot::ggsave2(file.path(dir_res, 'ggtree.medicc2_ME.pdf'), treeplt, 
                 width = 5, height = 4)



#------ plotheatmap integer ------
message('plotheatmap integer')
write_rds(css_mat, file.path(dir_res, 'css_mat.rds'))
objd@consensus <- css_mat

z <- 'clones'
para.plotHeatmap.group <- 'run'
if (length(unique(objd@colData[['run']])) == 1) {
  para.plotHeatmap.group = NULL
} 
n_z_opts <- ncol(objd@consensus)
order_cells <- NULL
if (n_z_opts > 1) {
  objd <- runConsensusPhylo(objd)
  order_cells <- 'consensus_tree'
}

pdf(file.path(dir_res, sprintf('cna.heatmap_consensus_int.combo_coda.%s.pdf', z)), 
    width = 15, height = 5, onefile = F, useDingbats = F)
p <- try( draw( plotHeatmap(
  objd, label = z, 
  assay = 'integer',
  order_cells = order_cells,
  # col = pal_copykit_int_clip_0_2_6,
  label_colors = list(clones=pal_cna_clones),
  group = para.plotHeatmap.group, consensus=T) ) )
if (!'try-error' %in% class(p)) {draw(p)}; rm(p)
dev.off()

pdf(file.path(dir_res, sprintf('paper.cna.heatmap_consensus_int.combo_coda.%s.pdf', z)), 
    width = 15, height = 5, onefile = F, useDingbats = F)
p <- try( draw( plotHeatmap(
  objd, label = z, 
  assay = 'integer',
  order_cells = order_cells,
  # col = pal_copykit_int_clip_0_2_6,
  label_colors = list(clones=pal_cna_clones_paper),
  group = para.plotHeatmap.group, consensus=T) ) )
if (!'try-error' %in% class(p)) {draw(p)}; rm(p)
dev.off()


message('done viz')
cat('Done: ', sample_name, '\n')