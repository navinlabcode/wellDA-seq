#---------------------------
# Find DACH (differntial accessible chromatin hubs) between two groups              ----    
# 
# in a context of the dosage effect (geno-to-pheno) analysis. 
#---------------------------
setwd('/volumes/USR1/yyan/project/coda')
if (T) {
  library(Signac); library(copykit); library(Seurat); library(tidyverse)
  library(MatrixGenerics); library(Matrix)
  library(ComplexHeatmap); library(ggpubr)
  source('./rsrc/utils.R')
  source('./rsrc/utils.copykit.R')
  is_CNA_varbin.int <- function(csc, base=2) {
    csc[is.na(csc)] <- base
    o <- rowMin(csc) == base & rowMax(csc) == base
    return(!o)
  }
  library(cli)
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


dir_dab <- file.path(dir_coda, 'dab_two_clones') ## dab results [reusable]
dir_dap <- file.path(dir_coda, 'dap_two_clones') ## dap results [reusable]
dir_dach <- file.path(dir_coda, 'dach_two_clones') ## DACH results [reusable]
dir_dos <- file.path(dir_coda, 'dosage_effect_two_clones') ## dossage-effect analysis result
fs::dir_create(dir_dach)

#------ load copykit and signac co-assay objects ------
load_coda(dir_coda)
print(obja); print(objd); print(head(df_meta[, 1:3]))
#------ load CCAN object ------
objccan <- read_rds(file.path(dir_coda, 'cicero', 
                              'CCAN_assay.default.seurat.rds'))
print(objccan)
Idents(objccan) <- ident_str; print(table(Idents(objccan)))

#------ colors ------
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
# ~~~ Which two clones to compare ~~~ --------------------------------------  
# load the previously defined clones
df_comparisons <- read_csv(file.path(dir_dos, 'df_comparisons.csv'))
# ~~~ Run DACH ~~~ --------------------------------------  
# cp_i = 1
for (cp_i in seq_len(nrow(df_comparisons))) {
  clone_w <- df_comparisons$clone_w[cp_i]  
  clone_v <- df_comparisons$clone_v[cp_i]  
  
  cp_code <- sprintf('%s.VS.%s', clone_w, clone_v)
  message(cp_code)
  dir_dos_cp <- file.path(dir_dos, cp_code); fs::dir_create(dir_dos_cp)
  dir_dach_cp <- file.path(dir_dach, cp_code); fs::dir_create(dir_dach_cp)
  cells_cp <- rownames(df_meta)[df_meta[[dna_grp_str]] %in% c(clone_w, clone_v)]
  print(head(cells_cp, 10))
  
  #------ create object ------
  f_objccan_cp <- file.path(dir_dos_cp, 'CCAN_assay.default.seurat.rds')
  # if (!file.exists(f_objccan_cp)) {
  if (T) {
    objccan_cp <- subset(objccan, cells = cells_cp)
    write_rds(objccan_cp, f_objccan_cp)
  } else {
    objccan_cp <- read_rds(f_objccan_cp)
  }
  
  #------ DACH ------
  cat('running DACH', cp_code, '...\n')
  f_dach0 <- file.path(dir_dach_cp, 'df_dach0.rds')
  if (! file.exists(f_dach0)) {
  # if (T) {
    df_dach0 <- FindMarkers(
      objccan_cp, 
      ident.1 = clone_w, ident.2 = clone_v, 
      test.use='wilcox',
      # test.use='negbinom',
      # latent.vars = 'nFrags', 
      logfc.threshold = log2(1.01)
    )
    df_dach0 <- rownames_to_column(df_dach0, 'CH')
    write_rds(df_dach0, file.path(dir_dach_cp, 'df_dach0.rds'))
    write_csv(df_dach0, file.path(dir_dach_cp, 'df_dach0.csv'))
    
  } else {
    df_dach0 <- read_rds(f_dach0)
  }
  
  table(df_dach0$p_val_adj < 0.05)
  ch_is_dach <- with(
    df_dach0, 
    p_val_adj < 0.05 & abs(avg_log2FC) > log(1.1)
  ); ch_is_dach[is.na(ch_is_dach)] <- FALSE; print(table('chromatin_hub_is_DACH'=ch_is_dach))
  ch_is_dach <- ifelse(ch_is_dach, 'yes','no')
  
  df_dach0$CH_is_dach <- ch_is_dach
  
  write_rds(df_dach0, f_dach0)
  write_csv(df_dach0, file.path(dir_dach_cp, 'df_dach0.csv'))  
}

cat('Done DACH', sample_name, '\n')

# ~~~ 3. Collect DNA, ATAC (peak), ATAC (chromatin hub) information ~~~ --------------------------------------  
cat("# 3 ~~~ Collect DNA, ATAC (peak), ATAC (chromatin hub) information ~~~ --------------------------------------  ")
# curate these information:
# varbin_has_CNA,
# varbin_has_peak,
# varbin_is_dab, 
# varbin_has_dap,
# varbin_has_dach,
rm(cp_code)
for (cp_i in seq_len(nrow(df_comparisons))) {
  clone_w <- df_comparisons$clone_w[cp_i]  
  clone_v <- df_comparisons$clone_v[cp_i]  
  
  cp_code <- sprintf('%s.VS.%s', clone_w, clone_v)
  message(cp_code)
  
  dir_dab_cp <- file.path(dir_dab, cp_code)
  dir_dap_cp <- file.path(dir_dap, cp_code)
  dir_dach_cp <- file.path(dir_dach, cp_code)
  dir_dos_cp <- file.path(dir_dos, cp_code)
  f_report_cp <- file.path(dir_dos_cp, 'df_report_cp.rds')
  
  # if (file.exists(f_report_cp)) {
  #   df_report_cp <- read_rds(f_report_cp)
  #   if ('varbin_has_dach' %in% colnames(df_report_cp)) {
  #     cat('report has added varbin_has_dach.\n')
  #     next()
  #   } else {
  #     cat('adding varbin_has_dach...\n')
  #   }
  #   rm(df_report_cp)
  # }
  
  cat('collecting information', cp_code, '...\n')  
  obja_cp <- read_rds(file.path(dir_dos_cp, 'obja.rds'))
  objd_cp <- read_rds(file.path(dir_dos_cp, 'objd.rds'))
  objccan_cp <- read_rds(file.path(dir_dos_cp, 'CCAN_assay.default.seurat.rds'))
  df_meta_cp <- read_rds(file.path(dir_dos_cp, 'metadata.df.rds'))
  df_dab0 <- read_rds(file.path(dir_dab_cp, 'df_dab0.rds'))
  df_dap0 <- read_rds(file.path(dir_dap_cp, 'df_dap0.rds'))
  df_dach0 <- read_rds(file.path(dir_dach_cp, 'df_dach0.rds'))
  # df_report_cp
  
  
  #------ DNA ------
  cat('range of integer copy number')
  objd_cp <- calcConsensus(objd_cp, assay = 'integer_fixed', 
                           consensus_by = dna_grp_str)
  print(range(consensus(objd_cp)))
  varbin_has_CNA <- is_CNA_varbin.int(as.matrix(consensus(objd_cp)))
  varbin_has_CNA <- ifelse(varbin_has_CNA, 'yes', 'no')
  varbin_is_dab <- df_dab0$varbin_is_dab
  table(varbin_is_dab)
  
  #------ ATAC (chromatin hubs) ------
  ## DACH:{peaks} -%overlap%- varbin
  dach_cp <- df_dach0$CH[df_dach0$CH_is_dach == 'yes']; str(dach_cp)
  if (length(dach_cp) != 0) {
    dict_CH2peaks <- objccan_cp@assays$ccan@meta.features
    dict_CH2peaks <- ruok::deframe_to_list(dict_CH2peaks[, c('CH', 'peaks_set')])
    dach_peaks_cp <- lapply(dict_CH2peaks[dach_cp], function(s) {
      s <- unique(str_split_1(s, ','))
      s
    })
    dach_peaks_cp <- unlist(dach_peaks_cp)
    dach_peaks_cp <- unique(dach_peaks_cp); str(dach_peaks_cp)
    dach_peaks_cp <- StringToGRanges(dach_peaks_cp)
  } else {
    dach_peaks_cp <- GRanges()
  }
  ovlp_bin_dach_peaks <- findOverlaps(
    query = rowRanges(objd_cp), 
    subject = dach_peaks_cp, 
    ignore.strand = TRUE)
  varbin_has_dach <- as.numeric(df_dab0$BIN) %in% queryHits(ovlp_bin_dach_peaks)
  varbin_has_dach <- ifelse(varbin_has_dach, 'yes', 'no')
  print(table(varbin_has_dach))
  
  #------ ATAC (peaks) ------
  ## associate varbin with the peaks
  head(df_dab0)
  ## 1/3 if a varbin has 
  dap_cp <- df_dap0$peak[df_dap0$peak_is_dap == 'yes']; str(dap_cp)
  if (length(dap_cp) != 0) {
    gr_dap_cp <- StringToGRanges(df_dap0$peak[df_dap0$peak_is_dap == 'yes'])
  } else {
    gr_dap_cp <- GRanges()
  }
  ovlp_bin_dap <- findOverlaps(
    query = rowRanges(objd_cp), 
    subject = gr_dap_cp, 
    ignore.strand = TRUE)
  head(ovlp_bin_dap)
  varbin_has_dap <- as.numeric(df_dab0$BIN) %in% queryHits(ovlp_bin_dap)
  varbin_has_dap <- ifelse(varbin_has_dap, 'yes', 'no')
  print(table(varbin_has_dap))
  ## 2/3 if a varbin has open peaks
  peaks_intest_cp <- rownames(obja_cp)
  ovlp_varbin_peaks <- findOverlaps(
    query = rowRanges(objd_cp), 
    subject = StringToGRanges(peaks_intest_cp), ignore.strand = TRUE)
  head(ovlp_varbin_peaks)
  varbin_has_peak <- as.numeric(df_dab0$BIN) %in% queryHits(ovlp_varbin_peaks)
  varbin_has_peak <- ifelse(varbin_has_peak, 'open', 'closed')
  ## 3/3 sanity check
  n_top <- 1000
  n_top <- length(dach_peaks_cp)
  if (n_top == 0) {n_top <- 100}
  peak_top_fc <- head(
    df_dap0$peak[order(abs(df_dap0$avg_log2FC), decreasing = T)], 
    n_top)
  peak_code <- with(df_dap0, case_when(
    p_val_adj < 0.05 & abs(avg_log2FC) >= log2(1.1) ~ 2, 
    peak %in% peak_top_fc ~ 1, 
    .default = 0
  )); df_dap0$peak_code <- peak_code
  ovlp_varbin_peaks_code <- findOverlaps(
    query = rowRanges(objd_cp), 
    subject = StringToGRanges(df_dap0$peak, ignore.strand = TRUE)
  )
  head(ovlp_varbin_peaks_code)
  varbin_peak_code <- rep(-1, length(rowRanges(objd_cp)))
  ovlp_varbin_peaks_code <- as.data.frame(ovlp_varbin_peaks_code)
  ovlp_varbin_peaks_code$peak_code <- df_dap0$peak_code[ovlp_varbin_peaks_code$subjectHits]
  ovlp_varbin_peaks_code <- ovlp_varbin_peaks_code %>% dplyr::group_by(queryHits) %>%
    dplyr::summarise(peak_code = max(peak_code))
  idx <- match(ovlp_varbin_peaks_code$queryHits, 1:length(1:length(varbin_peak_code)))
  varbin_peak_code[idx] <- ovlp_varbin_peaks_code$peak_code; rm(idx)
  print(table(varbin_peak_code))
  varbin_peak_code <- as.character(varbin_peak_code)
  varbin_has_soft_dap <- ifelse(varbin_peak_code %in% c('-1', '0'), 'no', 'yes')
  
  
  
  df_report_cp <- data.frame(
    varbin_has_CNA,
    varbin_has_peak,
    varbin_is_dab, 
    varbin_has_dap,
    varbin_has_dach,
    varbin_peak_code,
    varbin_has_soft_dap, 
    BIN_label = GRangesToString(rowRanges(objd_cp)),
    row.names = rownames(objd_cp)
  )
  
  write_rds(df_report_cp, file.path(dir_dos_cp, 'df_report_cp.rds'))
  write_csv(df_report_cp, file.path(dir_dos_cp, 'df_report_cp.csv'))
  
  print(table(varbin_has_CNA=df_report_cp$varbin_has_CNA))
  print(table(varbin_has_peak=df_report_cp$varbin_has_peak))
  print(table(varbin_is_dab=df_report_cp$varbin_is_dab))
  print(table(varbin_has_dap=df_report_cp$varbin_has_dap))
  print(table(varbin_has_dach=df_report_cp$varbin_has_dach))
  
  rm(obja_cp); rm(objd_cp); rm(objccan_cp); rm(df_dab0); rm(df_dap0); rm(df_dach0); rm(df_report_cp)
}  

cat('Done report', sample_name, '\n')

# ~~~ 4. Viz geno vs pheno for each comparison ~~~ --------------------------------------  
cat('# ~~~ 4. Viz geno vs pheno for each comparison ~~~ --------------------------------------  \n')
library(UpSetR); library(VennDiagram); library(ComplexHeatmap)
df_comparisons <- read_csv(file.path(dir_dos, 'df_comparisons.csv'))

pal_list[['varbin_has_peak']] <- c('closed' = 'ghostwhite', 'open'='mistyrose4')
pal_list[['varbin_has_CNA']] <- c('no' = 'ghostwhite', 'yes'='mistyrose4')
pal_list[['varbin_is_dab']] <- c('no' = 'ghostwhite', 'yes'='#002fa7')
pal_list[['varbin_has_dap']] <- c('no' = 'ghostwhite', 'yes'='goldenrod2')
pal_list[['varbin_has_dach']] <- c('no' = 'ghostwhite', 'yes'='#800020')
pal_list[['varbin_peak_code']] <- c('-1' = 'ghostwhite', '0' = 'ghostwhite', 
                                    '1' = 'cyan', '2' = 'magenta', 
                                    '3' = 'red')
pal_list[['varbin_has_soft_dap']] <- c('no' = 'ghostwhite', 'yes'='cyan')


#   function(x) {
#   if (x == -1) return('ghostwhite')
#   if (x == 0)  return('grey')
#   if (x == 1)  return('cyan')
#   if (x == 2)  return('magenta')
#   if (x == 3)  return('red')
# }

for (cp_i in seq_len(nrow(df_comparisons))) {
  clone_w <- df_comparisons$clone_w[cp_i]  
  clone_v <- df_comparisons$clone_v[cp_i]  
  
  cp_code <- sprintf('%s.VS.%s', clone_w, clone_v)
  message(cp_code)
  
  dir_dab_cp <- file.path(dir_dab, cp_code)
  dir_dap_cp <- file.path(dir_dap, cp_code)
  dir_dach_cp <- file.path(dir_dach, cp_code)
  dir_dos_cp <- file.path(dir_dos, cp_code)
  
  obja_cp <- read_rds(file.path(dir_dos_cp, 'obja.rds'))
  objd_cp <- read_rds(file.path(dir_dos_cp, 'objd.rds'))
  objccan_cp <- read_rds(file.path(dir_dos_cp, 'CCAN_assay.default.seurat.rds'))
  df_meta_cp <- read_rds(file.path(dir_dos_cp, 'metadata.df.rds'))
  df_dab0 <- read_rds(file.path(dir_dab_cp, 'df_dab0.rds'))
  df_dap0 <- read_rds(file.path(dir_dap_cp, 'df_dap0.rds'))
  df_dach0 <- read_rds(file.path(dir_dach_cp, 'df_dach0.rds'))
  df_report_cp <- read_rds(file.path(dir_dos_cp, 'df_report_cp.rds'))
  
  cat('visualizing ', cp_code, '...\n')
  #------ venn diagram ------
  
  cat('UpSetR / venndiagram...\n')
  list_geno_pheno_cp <- list(
    varbin_is_dab = df_report_cp$BIN_label[df_report_cp$varbin_is_dab=='yes'],
    varbin_has_dap = df_report_cp$BIN_label[df_report_cp$varbin_has_dap=='yes'],
    varbin_has_dach = df_report_cp$BIN_label[df_report_cp$varbin_has_dach == 'yes'], 
    varbin_has_soft_dap = df_report_cp$BIN_label[df_report_cp$varbin_has_soft_dap == 'yes'])
  str(list_geno_pheno_cp)
  tmp <- combn(x=1:length(list_geno_pheno_cp), m=2)
  for (tmp_i in 1:ncol(tmp)){
    tmp_a <- names(list_geno_pheno_cp)[tmp[1, tmp_i]]
    tmp_b <- names(list_geno_pheno_cp)[tmp[2, tmp_i]]
    p_venn <- try(VennDiagram::venn.diagram(
      list_geno_pheno_cp[c(tmp_a, tmp_b)], filename = NULL,
      col = unname(c(pal_list[[tmp_a]]['yes'], pal_list[[tmp_b]]['yes'])), 
      cat.col = unname(c(pal_list[[tmp_a]]['yes'], pal_list[[tmp_b]]['yes'])), 
      disable.logging=T, cat.default.pos = "text"))
    if (!'try-error' %in% class(p_venn)) {
      pdf(file.path(dir_dos_cp, sprintf('plot.venn.%s.pdf', tmp_i)), 
          width = 5, height = 5, useDingbats = F); grid.draw(p_venn); dev.off()    
    } else {
      cat('did not create venndiagram...')
    }
  }; cat('finish venndigram code.\n')  
  #------ heatmap ------
  
  cat('Heatmap...\n')
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
           'varbin_is_dab', 'varbin_has_dap', 
           'varbin_has_dach')
  ha_bin_class <- columnAnnotation(
    df = df_report_cp[, tmp, drop=F], 
    col = pal_list[tmp], 
    annotation_name_side = 'left')
  cat('draw...\n')
  pdf(file.path(dir_dos_cp, 'cna.heatmap.consensus_DAB_DAP_DACH_birdview.pdf'), 
      width = 15, height = 8, useDingbats = F, onefile = TRUE)
  draw(p1 %v% p2 %v% ha_bin_class)
  dev.off()
  
  tmp <- c('varbin_has_CNA', 'varbin_has_peak', 
           'varbin_is_dab', 'varbin_has_dap', 'varbin_has_soft_dap', 
           'varbin_has_dach')
  ha_bin_class <- columnAnnotation(
    df = df_report_cp[, tmp, drop=F], 
    col = pal_list[tmp], 
    annotation_name_side = 'left')
  cat('draw...\n')
  pdf(file.path(dir_dos_cp, 'cna.heatmap.consensus_DAB_DAP_DACH_birdview.examine.pdf'), 
      width = 15, height = 7, useDingbats = F, onefile = TRUE)
  draw(p2 %v% ha_bin_class)
  dev.off()
  
  #------ trackplot ------
  cat('[done]\n')
  
  rm(obja_cp); rm(objd_cp); rm(objccan_cp); rm(df_dab0); rm(df_dap0); rm(df_dach0)
}


cat('Done all: ', sample_name, '\n')

