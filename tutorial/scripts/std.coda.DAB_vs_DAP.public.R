#---------------------------
# cis-/trans-effect: Assess diff bins vs diff peaks              ----    
#---------------------------
setwd('/volumes/USR1/yyan/project/coda')
if (T) {
  library(Signac); library(copykit); library(Seurat); library(tidyverse)
  library(MatrixGenerics); library(Matrix)
  library(ComplexHeatmap); library(ggpubr)
  source('./rsrc/utils.R')
  source('./rsrc/utils.copykit.R')
  
  find_mat_pos <- function(index, n_rows, n_cols, matrix.byrow=FALSE){
    c <- index %/% n_rows 
    r <- index %% n_rows
    
    if (r == 0) {
      r <- n_rows
    } else {
      c <- c + 1
    }
    return(c('i'=r, 'j'=c))
  }

  is_CNA_varbin.int <- function(csc, base=2) {
    csc[is.na(csc)] <- base
    o <- rowMin(csc) == base & rowMax(csc) == base
    return(!o)
  }
  
}
if (F) {
  library(future)
  plan("multisession", workers = 8)
  plan()
  options(future.globals.maxSize = 50 * 1024 ^ 3) # for 50 Gb RAM
}

# ~~~ 0. Input ~~~ --------------------------------------  
sample_name <- 'DCIS66T_chip2'
args  <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  sample_name = args[[1]]
}


message(sample_name)
dir_coda <- file.path("/volumes/USR1/yyan/project/coda/",
                      "rds_coda_ready", 
                      sample_name, "aneuploid_epi")
ident_str <- dna_grp_str <- 'clones'

ncells_min_thre <- 15 ## the final working
ncells_min_thre <- 3 ## take a step back to run all. 


#------ setup output ------
dir_dab <- file.path(dir_coda, 'dab_two_clones') ## dab results [reusable]
dir_dap <- file.path(dir_coda, 'dap_two_clones') ## dap results [reusable]
dir_res <- file.path(dir_coda, 'dosage_effect_two_clones') ## dab vs dap

library(fs)
dir_create(dir_dab); dir_create(dir_dap); dir_create(dir_res); 

load_coda(dir_coda)

print(head(df_meta))
print(obja)
print(objd)
print(coda_dimplot2(df_meta, 'clones'))

#------ prepare colors ------
colnames(df_meta)
z_str_list <- c('clones', 'celltypes', 'cellstates', 'run')
pal_list <- lapply(z_str_list, function(z) {
  if (!z %in% colnames(df_meta)) {return(NULL)}
  
  pal_ploidy_class = c('aneuploid'='#cd0bbc', 'diploid'='#61d04f', 'Unknown'='grey', `NA`='grey')
  pal_ploidy_is_aneuploid = c(`TRUE`='#cd0bbc', `FALSE`='#61d04f', 'Unknown'='grey', `NA`='grey')
  switch (z,
    clones = new_palette_D(levels(df_meta$clones), pal = 'stallion'), 
    superclones = new_palette_D(levels(df_meta$superclones)),
    subclones = new_palette_D(levels(df_meta$subclones), pal = 'ironMan'),
    run = new_palette_D(unique(df_meta$run), pal = 'scales_hue'),
    celltypes = new_palette_D(levels(df_meta$celltypes), pal = 'circus'),
    cellstates = new_palette_D(levels(df_meta$cellstates), pal = 'circus'),
    ploidy_is_aneuploid = pal_ploidy_is_aneuploid,
    ploidy_class = pal_ploidy_class
  )
}); names(pal_list) <- z_str_list; pal_list <- pal_list[!sapply(pal_list, is.null)]
pal_list[['varbin_is_dab']] <- c('no' = 'ghostwhite', 'yes'='goldenrod2')
pal_list[['varbin_has_dap']] <- c('no' = 'ghostwhite', 'yes'='dodgerblue2')
pal_list[['varbin_has_peak']] <- c('closed' = 'ghostwhite', 'open'='mistyrose4')
pal_list[['varbin_has_CNA']] <- c('no' = 'ghostwhite', 'yes'='mistyrose4')

# ~~~ 1. Propose which two clones to compare ~~~ --------------------------------------  
cat('# ~~~ 1. Propose which two clones to compare ~~~ --------------------------------------  \n')


#------ remove small clones ------
ncells_per_dna_grp <- table(df_meta[[dna_grp_str]]); ncells_per_dna_grp <- c(ncells_per_dna_grp); print(ncells_per_dna_grp)
idx <- ncells_per_dna_grp >= ncells_min_thre
dna_grp_opts_use <- names(ncells_per_dna_grp)[idx]
cat('Candidate clones to compare:\n')
print(dna_grp_opts_use)
if (is_empty(dna_grp_opts_use) | length(dna_grp_opts_use) < 2){
  stop('Less than 2 clones having >= ', ncells_min_thre, ' cells. \n')
}
#------ load medic2 pair-wise distance ------
mat_medic2_pairwise <- read.table(
  file.path(dir_coda, 'medic2/medicc2_input_pairwise_distances.tsv'), header = T, sep = ' ')
mat_medic2_pairwise <- read_tsv(
  file.path(dir_coda, 'medic2/medicc2_input_pairwise_distances.tsv'))
mat_medic2_pairwise <- column_to_rownames(mat_medic2_pairwise, 'sample_id')
mat_medic2_pairwise <- as.matrix(mat_medic2_pairwise)
# view(mat_medic2_pairwise)
#------ the most distant/similar/abundant clones ------
dist_to_diploid <- mat_medic2_pairwise[dna_grp_opts_use, 'diploid']
dist_to_diploid <- sort(dist_to_diploid, decreasing = T)
dna_grp_opts_srt_by_aneu <- names(dist_to_diploid)
mat_medic2_pairwise_use <- mat_medic2_pairwise[dna_grp_opts_use, dna_grp_opts_use, drop=F]

library(ComplexHeatmap)
pdf(file.path(dir_res, 'mat_distance.pdf'), width = 5.5, height = 5, useDingbats = F)
draw(Heatmap(mat_medic2_pairwise_use, name = 'distance', 
             col = pal_heatmap_auto_c(c(mat_medic2_pairwise_use), pal = 'Viridis', rev = F),
             heatmap_width = unit(4, 'inch'), heatmap_height = unit(4, 'inch')))
dev.off()

#------ greedy list all possible comparison ------
df_comparisons_all <- as.data.frame(t(combn(x=dna_grp_opts_srt_by_aneu, 2)))
colnames(df_comparisons_all) <-  c('clone_w', 'clone_v')
head(df_comparisons_all)
df_comparisons_all <- tidyr::unite(
  df_comparisons_all, 'cp_code', c('clone_w', 'clone_v'), 
  sep = "_VS_", remove = TRUE)

comparing_clones_strategy_opts <- c('distant', 'similar', 'abundant')
mat_medic2_pairwise_use_eval <- mat_medic2_pairwise_use
diag(mat_medic2_pairwise_use_eval) <- NA
mat_medic2_pairwise_use_eval[lower.tri(mat_medic2_pairwise_use_eval)] <- NA
# comparing_clones_strategy_opts <- c('similar')
list_comparisons <- lapply(comparing_clones_strategy_opts, function(co) {
  out <- NULL
  if (co == 'abundant') {
    v <- ncells_per_dna_grp[dna_grp_opts_use]
    out <- names(v)[order(v, decreasing = T)]
    out <- head(out, 2)
    out <- intersect(dna_grp_opts_srt_by_aneu, out)
    return(c(out, co))
  }
  
  if (co %in% c('distant', 'similar')) {
    co_fun <- switch (co,
      distant = max, 
      similar = min)
    out_list <- which(mat_medic2_pairwise_use_eval == co_fun(mat_medic2_pairwise_use_eval, na.rm = T))
    print(out_list)
    out <- sapply(out_list, function(mat_loc) {
      mat_loc_i_j <- find_mat_pos(
        mat_loc, 
        n_rows = nrow(mat_medic2_pairwise_use_eval), 
        n_cols = ncol(mat_medic2_pairwise_use_eval))
      mat_loc_out <- c(rownames(mat_medic2_pairwise_use_eval)[mat_loc_i_j['i']], 
                       colnames(mat_medic2_pairwise_use_eval)[mat_loc_i_j['j']])
      mat_loc_out <- intersect(dna_grp_opts_srt_by_aneu, mat_loc_out)
      return(c(mat_loc_out, co))
    })
    out <- t(out)
    return(out)
  }
})

df_comparisons <- as.data.frame(do.call('rbind', list_comparisons))
## FIXED terminology: 
## clone_w: clone with more distant to aneuploid
## clone_v: clone with less distant to aneuploid
colnames(df_comparisons) <- c('clone_w', 'clone_v', 'type')
df_comparisons <- tidyr::unite(
  df_comparisons, 'cp_code', c('clone_w', 'clone_v'), sep='_VS_', remove = T
)
df_comparisons <- left_join(df_comparisons_all, df_comparisons, by='cp_code')
df_comparisons$type[is.na(df_comparisons$type)] <- 'others'
df_comparisons <- tidyr::separate(
  df_comparisons, 'cp_code', into = c('clone_w', 'clone_v'), sep='_VS_', remove = F)
## remove the repeated cp
tmp <- table(df_comparisons$cp_code)
tmp <- names(tmp)[tmp>1]
tmp <- df_comparisons$cp_code %in% tmp & df_comparisons$type != 'abundant'
df_comparisons <- df_comparisons[!tmp, , drop=F]
df_comparisons$distanct_w_v <- sapply(1:nrow(df_comparisons), function(i) {
  mat_medic2_pairwise[df_comparisons$clone_w[i], df_comparisons$clone_v[i]]
})
rm(tmp)
df_comparisons$cp_code <- NULL
print(df_comparisons)
write_csv(df_comparisons, file.path(dir_res, 'df_comparisons.csv'))

# ~~~ 2. Run DAB and DAP for each comparision ~~~ --------------------------------------  
cat('# ~~~ 2. Run DAB and DAP for each comparision ~~~ --------------------------------------  \n')
n_clone_pairs <- nrow(df_comparisons)
cp_i <- 1

if (!'integer_fixed' %in% assayNames(objd)) {
  assay(objd, 'integer_fixed') <- assay(objd, 'integer')
}

for (cp_i in seq_len(n_clone_pairs)) {
  clone_w <- df_comparisons$clone_w[cp_i]  
  clone_v <- df_comparisons$clone_v[cp_i]  
  
  cp_code <- sprintf('%s.VS.%s', clone_w, clone_v)
  message(cp_code)
  dir_dab_cp <- file.path(dir_dab, cp_code); fs::dir_create(dir_dab_cp) 
  dir_dap_cp <- file.path(dir_dap, cp_code); fs::dir_create(dir_dap_cp)
  dir_res_cp <- file.path(dir_res, cp_code); fs::dir_create(dir_res_cp)
  
  cells_cp <- rownames(df_meta)[df_meta[[dna_grp_str]] %in% c(clone_w, clone_v)]
  print(head(cells_cp, 10))
  #------ create object ------
  if (!file.exists(file.path(dir_res_cp, 'metadata.df.rds'))) {
  # if (T) {
    cat('creating coda object for ', cp_code, '...\n')
    obja_cp <- subset(obja, cells = cells_cp)
    objd_cp <- objd[, cells_cp]
    if (!'integer_fixed' %in% assayNames(objd_cp)) {
      assay(objd_cp, 'integer_fixed') <- assay(objd_cp, 'integer')
    }
    df_meta_cp <- df_meta[cells_cp, ]
    objd_cp <- runDistMat(objd_cp, metric = 'manhattan', n_threads = 20)
    write_coda(dir_res_cp, df = df_meta_cp, A = obja_cp, D = objd_cp)
  } else {
    cat('loading coda object for ', cp_code, '...\n')
    obja_cp <- read_rds(file.path(dir_res_cp, 'obja.rds'))
    objd_cp <- read_rds(file.path(dir_res_cp, 'objd.rds'))
    df_meta_cp <- read_rds(file.path(dir_res_cp, 'metadata.df.rds'))
  }  
  
  #------ UMAP of clones on ATAC space ------
  pumap <- DimPlot(obja_cp, shuffle = T, group.by=ident_str, cols = pal_list[[ident_str]]) + 
    theme_void() + 
    theme(aspect.ratio = 1) + 
    ggpubr::border() + labs(caption = 'ATAC space')
  ggsave(file.path(dir_res_cp, sprintf('atac.umap.%s.pdf', ident_str)), pumap, width = 4, height = 3, useDingbats = F)
  #------ DAB ------
  cat('running DAB of', cp_code,  '...\n')
  f_dab0 <- file.path(dir_dab_cp, 'df_dab0.rds')
  if (! file.exists(f_dab0)) {
    df_dab0 <- findMarkerBins(
      objd, assay = 'integer_fixed', test_use = 'wilcox',
      cells.1 = colnames(objd)[objd@colData[[dna_grp_str]] == clone_w],
      cells.2 = colnames(objd)[objd@colData[[dna_grp_str]] == clone_v])
    # table(df_dab0$p_val_adj < 0.05, useNA='ifany')
    # colnames(df_dab0)
    varbin_is_dab <- with(
      df_dab0, p_val_adj < 0.05 & abs(integer.1 - integer.2) >= 0.1)
    varbin_is_dab[is.na(varbin_is_dab)] <- FALSE
    varbin_is_dab <- ifelse(varbin_is_dab, 'yes', 'no'); print(table(varbin_is_dab))
    df_dab0$varbin_is_dab <- varbin_is_dab
    write_rds(df_dab0, file.path(dir_dab_cp, 'df_dab0.rds'))
    write_csv(df_dab0, file.path(dir_dab_cp, 'df_dab0.csv'))
  } else {
    df_dab0 <- read_rds(file.path(dir_dab_cp, 'df_dab0.rds'))
    varbin_is_dab <- df_dab0$varbin_is_dab
    print(table(varbin_is_dab))
  }
  
  #------ DAP ------
  cat('running DAP of', cp_code,  '...\n')
  Idents(obja_cp) <- ident_str
  
  #------ remove less informative peaks to speed up computing ------
  atac_counts <- GetAssayData(obja_cp, slot='counts', assay='peaks')
  
  tmp <- rowSums2(atac_counts) > 0 & rowSums2(atac_counts>0) >= 3 & rowVars(atac_counts) > rowMeans2(atac_counts)
  table(tmp)
  peaks_intest_cp <- rownames(obja_cp)[which(tmp)]; str(peaks_intest_cp)
  cat('Testing', length(peaks_intest_cp), 'peaks...\n')
  
  if (!file.exists(file.path(dir_dap_cp, 'df_dap0.rds'))) {
    df_dap0 <- FindMarkers(
      obja_cp, 
      features = peaks_intest_cp, 
      ident.1 = clone_w, ident.2 = clone_v, 
      test.use='wilcox',
      # test.use='negbinom',
      # latent.vars = 'nFrags', 
      logfc.threshold = log2(1.01)
    )
    df_dap0 <- rownames_to_column(df_dap0, 'peak')
    write_rds(df_dap0, file.path(dir_dap_cp, 'df_dap0.rds'))
    write_csv(df_dap0, file.path(dir_dap_cp, 'df_dap0.csv'))
    
  } else {
    df_dap0 <- read_rds(file.path(dir_dap_cp, 'df_dap0.rds'))
  }

  
  table(df_dap0$p_val_adj < 0.05)
  peak_is_dap <- with(
    df_dap0, 
    p_val_adj < 0.05 & abs(avg_log2FC) > log(1.1)
  ); print(table(peak_is_dap, useNA='ifany'))
  peak_is_dap[is.na(peak_is_dap)] <- FALSE; print(table(peak_is_dap, useNA='ifany'))
  peak_is_dap <- ifelse(peak_is_dap, 'yes', 'no')
  df_dap0$peak_is_dap <- peak_is_dap
  
  write_rds(df_dap0, file.path(dir_dap_cp, 'df_dap0.rds'))
  write_csv(df_dap0, file.path(dir_dap_cp, 'df_dap0.csv'))
  
}  


# ~~~ 3. Collect DNA and ATAC information ~~~ --------------------------------------  
cat('# ~~~ 3. Collect DNA and ATAC information ~~~ --------------------------------------  \n')
for (cp_i in seq_len(nrow(df_comparisons))) {
  clone_w <- df_comparisons$clone_w[cp_i]  
  clone_v <- df_comparisons$clone_v[cp_i]  
  
  cp_code <- sprintf('%s.VS.%s', clone_w, clone_v)
  message(cp_code)

  dir_dab_cp <- file.path(dir_dab, cp_code)
  dir_dap_cp <- file.path(dir_dap, cp_code)
  dir_res_cp <- file.path(dir_res, cp_code)
  
  cat('collecting information', cp_code, '...\n')  
  obja_cp <- read_rds(file.path(dir_res_cp, 'obja.rds'))
  objd_cp <- read_rds(file.path(dir_res_cp, 'objd.rds'))
  df_meta_cp <- read_rds(file.path(dir_res_cp, 'metadata.df.rds'))
  df_dab0 <- read_rds(file.path(dir_dab_cp, 'df_dab0.rds'))
  df_dap0 <- read_rds(file.path(dir_dap_cp, 'df_dap0.rds'))
  # df_report_cp
  
  
  #------ DNA ------
  cat('range of integer copy number')
  objd_cp <- calcConsensus(objd_cp, assay = 'integer_fixed', 
                           consensus_by = dna_grp_str)
  print(range(consensus(objd_cp)))
  varbin_has_CNA <- is_CNA_varbin.int(as.matrix(consensus(objd_cp)))
  varbin_has_CNA <- ifelse(varbin_has_CNA, 'yes', 'no')
  try(varbin_is_dab <- df_dab0$varbin_is_dab)
  table(varbin_is_dab)
  #------ ATAC ------
  ## associate varbin with the peaks
  peaks_intest_cp <- rownames(obja_cp)
  head(df_dab0)
  dap_cp <- df_dap0$peak[df_dap0$peak_is_dap == 'yes']; str(dap_cp)
  if (length(dap_cp) != 0) {
    gr_dap_cp <- StringToGRanges(df_dap0$peak[df_dap0$peak_is_dap == 'yes'])
  } else {
    gr_dap_cp <- GRanges()
  }

  
  ovlp_dab_dap <- findOverlaps(
    query = rowRanges(objd_cp), 
    subject = gr_dap_cp, 
    ignore.strand = TRUE)
  
  head(ovlp_dab_dap)
  varbin_has_dap <- as.numeric(df_dab0$BIN) %in% queryHits(ovlp_dab_dap)
  varbin_has_dap <- ifelse(varbin_has_dap, 'yes', 'no')
  table(varbin_has_dap)
  
  ovlp_varbin_peaks <- findOverlaps(
    query = rowRanges(objd_cp), 
    subject = StringToGRanges(peaks_intest_cp), ignore.strand = TRUE)
  head(ovlp_varbin_peaks)
  varbin_has_peak <- as.numeric(df_dab0$BIN) %in% queryHits(ovlp_varbin_peaks)
  varbin_has_peak <- ifelse(varbin_has_peak, 'open', 'closed')
  
  df_report_cp <- data.frame(
    varbin_has_CNA,
    varbin_has_peak,
    varbin_is_dab, 
    varbin_has_dap,
    BIN_label = GRangesToString(rowRanges(objd_cp)),
    row.names = rownames(objd_cp)
  )
  
  write_rds(df_report_cp, file.path(dir_res_cp, 'df_report_cp.rds'))
  write_csv(df_report_cp, file.path(dir_res_cp, 'df_report_cp.csv'))
  
  print(table(varbin_has_CNA=df_report_cp$varbin_has_CNA))
  print(table(varbin_has_peak=df_report_cp$varbin_has_peak))
  print(table(varbin_is_dab=df_report_cp$varbin_is_dab))
  print(table(varbin_has_dap=df_report_cp$varbin_has_dap))
}  

# ~~~ 4. Viz DAB vs DAP for each comparison ~~~ --------------------------------------  
cat('# ~~~ 4. Viz DAB vs DAP for each comparison ~~~ --------------------------------------  \n')
library(UpSetR); library(VennDiagram); library(ComplexHeatmap)
df_comparisons <- read_csv(file.path(dir_res, 'df_comparisons.csv'))
for (cp_i in seq_len(nrow(df_comparisons))) {
  
  clone_w <- df_comparisons$clone_w[cp_i]  
  clone_v <- df_comparisons$clone_v[cp_i]  
  
  cp_code <- sprintf('%s.VS.%s', clone_w, clone_v)
  dir_dab_cp <- file.path(dir_dab, cp_code)
  dir_dap_cp <- file.path(dir_dap, cp_code)
  dir_res_cp <- file.path(dir_res, cp_code)
  
  message(cp_code)
  cat('visualizing ', cp_code, cp_i, '/', nrow(df_comparisons), '...\n')

  obja_cp <- read_rds(file.path(dir_res_cp, 'obja.rds'))
  objd_cp <- read_rds(file.path(dir_res_cp, 'objd.rds'))
  df_meta_cp <- read_rds(file.path(dir_res_cp, 'metadata.df.rds'))
  df_dab0 <- read_rds(file.path(dir_dab_cp, 'df_dab0.rds'))
  df_dap0 <- read_rds(file.path(dir_dap_cp, 'df_dap0.rds'))
  df_report_cp <- read_rds(file.path(dir_res_cp, 'df_report_cp.rds'))
  
  # UMAPPlot(obja_cp, group.by=ident_str)
  # plotUmap(objd_cp, label = ident_str)
  # ------ UpSetR / venndiagram ------
  cat('UpSetR / venndiagram...\n')

  colnames(df_report_cp)
  list_varbin_dab_dap_cp <- list(
    varbin_is_dab = df_report_cp$BIN_label[df_report_cp$varbin_is_dab=='yes'],
    varbin_has_dap = df_report_cp$BIN_label[df_report_cp$varbin_has_dap=='yes'])
  str(list_varbin_dab_dap_cp)
  p_upset <- try(upset(fromList(list_varbin_dab_dap_cp),
                       order.by = "degree", decreasing = T))
  pdf(file.path(dir_res_cp, 'plot.upset.pdf'),width = 6, height = 5, useDingbats = F)
  try(print(p_upset))
  dev.off()
  # ggVennDiagram(x = list_varbin_dab_dap_cp)
  p_venn <- try(VennDiagram::venn.diagram(
    list_varbin_dab_dap_cp, filename = NULL,
    col = c('goldenrod2', 'dodgerblue2'), cat.col = c('goldenrod2', 'dodgerblue2'),
    disable.logging=T, cat.default.pos = "text"))
  if (!'try-error' %in% class(p_venn)) { 
    pdf(file.path(dir_res_cp, 'plot.venn.pdf'), width = 5, height = 5, useDingbats = F)
    try(grid.draw(p_venn))
    dev.off()
  }  
  try(dev.off())
  #------ plotHeatmap ------
  cat('plotHeatmap...\n')
  objd_cp <- calcConsensus(objd_cp, assay = 'segment_ratios', 
                           consensus_by = dna_grp_str)
  objd_cp <- runConsensusPhylo(objd_cp)
  p1 <- plotHeatmap(
    objd_cp, assay = 'segment_ratios', 
    order_cells = 'consensus_tree',
    consensus = T, raster=T, col = pal_copykit_ratio_clip1, 
    label = dna_grp_str, label_colors = pal_list[dna_grp_str])
  
  if (!identical(assay(objd_cp, 'integer'), assay(objd_cp, 'integer_fixed'))) {
    message('Use integer_fixed instead!!!')
    assay(objd_cp, 'integer') <- assay(objd_cp, 'integer_fixed')
  }
  objd_cp <- calcConsensus(objd_cp, assay = 'integer', 
                           consensus_by = dna_grp_str)
  objd_cp <- runConsensusPhylo(objd_cp)
  p2 <- plotHeatmap(
    objd_cp, assay = 'integer', 
    order_cells = 'consensus_tree',
    consensus = T, raster=T, 
    label = dna_grp_str, label_colors = pal_list[dna_grp_str])
  tmp <- c('varbin_has_CNA', 'varbin_has_peak', 
           'varbin_is_dab', 'varbin_has_dap')
  ha_bin_class <- columnAnnotation(
    df = df_report_cp[, tmp, drop=F], 
    col = pal_list[tmp], 
    annotation_name_side = 'left')
  cat('draw...\n')
  pdf(file.path(dir_res_cp, 'cna.heatmap.consensus_DAB_DAP_birdview.pdf'), 
      width = 15, height = 8, useDingbats = F, onefile = TRUE)
  draw(p1 %v% p2 %v% ha_bin_class)
  dev.off()
  cat('[done]\n')
}

cat('Done', sample_name, '!\n')