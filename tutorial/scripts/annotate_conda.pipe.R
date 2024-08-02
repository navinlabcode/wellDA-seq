## Written by: Yun Yan (https://github.com/Puriney)\n\n
#---------------------------
# Goal: Annotate conda raw object              ----    
# 
# Input: conda object without cell types and clones finalized
# Output: Finalize cell types and clones 
# 
# Progress: 
# 
# 1. Determine tentative cell types
# Pick a proper ATAC clustering
# Run DEGs to name cell types
#  
# 2. Determine clones
# Based on DNA's findAneuploid, split into two subsets: 
#   - diploid (white) => D + D_A (doublets)
#   - aneuploid (purple) => A + A_D (doublets)
# Merge the D and A subsets
#   - Compare similar subclones and merge into clones
# 
# 3. Examine if cell types and clones are consistent. Remove ATAC doublets. 
# 
# 4. Export to `coda_preready` folder because cell types and clones are finalized. 
# 
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
})
setwd("/volumes/USR1/yyan/project/coda")
source('/volumes/USR1/yyan/project/coda/rsrc/utils.R')

sample_name = 'DCIS66T_chip2'

dir_proj = file.path('./rds_coda/', sample_name)

f_a <- file.path(dir_proj, 'obja.rds')
f_d <- file.path(dir_proj, 'objd.rds')
f_meta <- file.path(dir_proj, 'metadata.df.rds')
obja <- read_rds(f_a)
objd <- read_rds(f_d)
df_meta <- read_rds(f_meta)

stopifnot( all.equal(Cells(obja), colnames(objd)) )
stopifnot( all.equal(Cells(obja), rownames(df_meta)) )
print(obja)
print(objd)

#------------------- ~~~ 1. cell typing ~~~ -------------------  
#------ 1/3. pick a proper ATAC clustering ------
ident_atac_clusters <- NULL
dir_res <- file.path(dir_proj, 'celltypes'); fs::dir_create(dir_res)
stop('beging manually celltpying...')
source("./rsrc/annotate_conda.s1_celltying_1_clustering.R")
message('Ident for the ATAC clusters: ', ident_atac_clusters)
#------ 2/3. iRNA DEG and manually assign cell types ------
ident_use <- ident_atac_clusters
f_atac_manual <- NULL
dir_snippet <- dir_res
source('./rsrc/snippet.ATAC_FindAllMarkers.R')
message('Todo: Fill in the ATAC cluster names in: ', f_atac_manual)
#------ 3/3. add `celltypes` column into the object ------

df_celltype <- read_csv(f_atac_manual, col_names = c('ident', 'celltypes'), 
                        col_types = c('c', 'c'))
dict_celltype2ident <- deframe(df_celltype[, c('celltypes', 'ident')])
obja$celltypes <- forcats::fct_recode(
  obja@meta.data[, ident_atac_clusters], 
  !!!dict_celltype2ident
)
obja$celltypes <- forcats::fct_relevel(
  obja$celltypes, 
  unique(names(dict_celltype2ident))
)
obja$seurat_clusters <- obja$celltypes
p <- UMAPPlot(obja, group.by='celltypes', label=T) + theme(aspect.ratio = 1)
print(p)
ggsave(file.path(dir_res, 'dr.celltypes.pdf'), p, width = 5, height = 4)

DefaultAssay(obja) <- 'peaks'
write_rds(obja, f_a) ## overwrite

update_coda(df_meta, obja, objd)

#------------------- ~~~ 2. Determine clones ~~~ -------------------  
if (!'run' %in% colnames(df_meta)) { df_meta$run <- sample_name; update_coda(df_meta, obja, objd)}
pal_ploidy_class = c('aneuploid'='#cd0bbc', 'diploid'='#61d04f', 'Unknown'='grey', `NA`='grey')
pal_ploidy_is_aneuploid = c(`TRUE`='#cd0bbc', `FALSE`='#61d04f', 'Unknown'='grey', `NA`='grey')
pal_cna_clones = new_palette_D(levels(df_meta$clones), pal = 'stallion')
pal_cna_superclones = new_palette_D(levels(df_meta$superclones))
pal_cna_subclones = new_palette_D(levels(df_meta$subclones), pal = 'ironMan')
pal_atac_clusters = ArchR::paletteDiscrete(df_meta$seurat_clusters, set = 'circus')
pal_run = new_palette_D(unique(df_meta$run), pal = 'scales_hue') # consistent with copykit
#------ split cells into green and purple cells ------
objd <- findAneuploidCells(objd, resolution = 'auto')
table(objd[['subclones']], useNA='always')
table(objd[['is_aneuploid']], useNA='always')
table(objd[['celltypes']], useNA='always')
try(objd[['ploidy_class']] <- ifelse(objd[['is_aneuploid']], 'aneuploid', 'diploid'))
table(objd[['ploidy_class']], useNA='always')
update_coda(df_meta, obja, objd)
## backup
objd0 <- objd; df_meta0 <- df_meta; obja0 <- obja
objd <- objd0; df_meta <- df_meta0; obja <- obja0
# objd@colData$sample <- colnames(objd)
## split
for (z in c('diploid', 'aneuploid')) {
# for (z in c('diploid')) {
  dir_snippet <- file.path(dir_proj, 'CNA_processing_ploidy_class', z)
  fs::dir_create(dir_snippet)
  tmp <- objd0[['ploidy_class']] == z; message(z, ' has ', sum(tmp), ' cells.')
  write_rds(objd0[, tmp], file.path(dir_snippet, 'objd.rds'))
  write_rds(obja0[, tmp], file.path(dir_snippet, 'obja.rds'))
  write_rds(df_meta0[tmp, ], file.path(dir_snippet, 'metadata.df.rds'))
  
}
for (cna_ploidy_class in c('diploid', 'aneuploid')) {
# for (cna_ploidy_class in c('diploid')) {
  message(cna_ploidy_class)
  dir_snippet <- file.path(dir_proj, 'CNA_processing_ploidy_class', cna_ploidy_class)
  objd    <- read_rds(file.path(dir_snippet, 'objd.rds'))
  objd <- runDistMat(objd, metric = 'manhattan', n_threads = 20)
  write_rds(objd, file.path(dir_snippet, 'objd.rds'))
}

#------ reclustering ------
for (cna_ploidy_class in c('diploid', 'aneuploid')) {
# for (cna_ploidy_class in c('diploid')) {
  message(cna_ploidy_class)
  dir_snippet <- file.path(dir_proj, 'CNA_processing_ploidy_class', cna_ploidy_class)
  load_coda(dir_snippet)
  print(ncol(objd))
  if (cna_ploidy_class == 'aneuploid') {
    source('./rsrc/snippet.CNA_clustering_in_aneuploid.R')
  } else {
    source('./rsrc/snippet.CNA_clustering_in_diploid.R')
  }
  update_coda(df_meta, obja, objd)
  write_coda(dir_snippet, df_meta, obja, objd)
}

#------ rerun scquantume ------
for (cna_ploidy_class in c('diploid', 'aneuploid')) {
# for (cna_ploidy_class in c('diploid')) {
  message(cna_ploidy_class)
  dir_snippet <- file.path(dir_proj, 'CNA_processing_ploidy_class', cna_ploidy_class)
  load_coda(dir_snippet)
  print(ncol(objd))
  
  objd <- calcInteger(objd, method = 'scquantum', assay = 'bincounts',
                      name = 'integer_scquantum')
  objd <- calcInteger(objd, method = 'fixed', assay = 'segment_ratios', 
                      ploidy_value = 2, 
                      name = 'integer_fixed')  
  
  update_coda(df_meta, obja, objd)
  write_coda(dir_snippet, df_meta, obja, objd)
}

#------ viz to prepare manual speculation [milestone] ------ 
for (cna_ploidy_class in c('diploid', 'aneuploid')) {
# for (cna_ploidy_class in c('aneuploid')) {
  message(cna_ploidy_class)
  dir_snippet <- file.path(dir_proj, 'CNA_processing_ploidy_class', 
                           cna_ploidy_class)
  load_coda(dir_snippet)
  z_tofig_opts <- c('subclones', 'celltypes')
  source('./rsrc/snippet.coda_viz.R')
}

# [milestone / drink water / grab coffee]

#------ Processing the green cells ------
dir_res <- file.path(dir_proj, 'CNA_processing_ploidy_class', 'diploid')
load_coda(o = dir_res)
stop('!manually examine and merge subclones!')
source('./rsrc/annotate_coda.s2_green.R')
print(table(objd[['clones']]))
update_coda(df_meta, obja, objd)
write_coda(dir_res, df_meta, obja, objd)
# load_coda(dir_res)

dir_snippet <- dir_res
z_tofig_opts <- c('clones')
source('./rsrc/snippet.coda_viz.R')

## subset the good dipoid cells if applicable
print(table(objd[['clones']]))
if ('X' %in% objd[['clones']]) {
  cat('removing noise diploid cells')
  cat(sum(objd[['clones']] %in% 'X'))
  write_lines(
    colnames(objd)[objd@colData$clones %in% 'X'], 
    file.path(dir_res, 'cell_names_noise.txt')
  )
  objd <- objd[, ! objd[['clones']] %in% 'X']
  objd[['clones']] <- fct_drop(objd[['clones']])
  subset_coda(df=df_meta, obja, objd, colnames(objd))
  update_coda(df_meta, obja, objd)
}
dir_res <- file.path(dir_proj, 'CNA_processing_ploidy_class', 
                     'diploid_rm_doublets')
fs::dir_create(dir_res)
write_coda(dir_res, df_meta, obja, objd)
dir_snippet <- dir_res
z_tofig_opts <- c('clones')
source('./rsrc/snippet.coda_figure.R')

# [milestone]

#------ Processing the purple cells ------
dir_res <- file.path(dir_proj, 'CNA_processing_ploidy_class', 'aneuploid')
load_coda(o = dir_res)

## re-fine clustering of purple cells if needed 
subclone_to_resubcluster <- NULL 
if (sample_name %in% c('DCIS28T') ) {subclone_to_resubcluster <- 'c4'}
if (!is.null(subclone_to_resubcluster)) {
  objd <- copykit_find_sub_clusters(
    objd, subcluster_by = 'subclones', cluster = subclone_to_resubcluster, 
    saveto_col_str = 'newsubclusters')
  table(objd[['newsubclusters']])
  plotUmap(objd, label = 'newsubclusters')
  objd[['subclones']] <- objd[['newsubclusters']]
  objd[['newsubclusters']] <- NULL
  plotUmap(objd, label = 'subclones')
  update_coda(df_meta, A = obja, D = objd)
  write_coda(o = dir_res, df = df_meta, A = obja, D = objd)
}
## Handling the low-quality cells (by scquantum) and doublets at the same time.
## Remove cells that are low-qual (by scquantum) doublets (lebeled as X)

## examine the scquantum strategy
source( './rsrc/annotate_coda.s2_purple_qc_by_scquantum.R' )
if (sample_name %in% c('DCIS28T', 'DCIS35T', 'DCIS67T', 'BCMDCIS69')) {
  objd[['pass_scquantum_QC']] <- objd[['pass_scquantum_sc']]
}
if (sample_name %in% 'BCMDCIS69') {
  objd[['pass_scquantum_QC']] <- TRUE
}
print(table(objd[['pass_scquantum_QC']], useNA='ifany'))
update_coda(df_meta, A = obja, D = objd)
ggsave(file.path(dir_res, 'coda_dimplot2.pass_scquantum_QC.pdf'), 
       coda_dimplot2(df_meta, 'pass_scquantum_QC', pal = 'Harmonic'), 
       width = 8, height = 5)
write_coda(o = dir_res, df = df_meta, A = obja, D = objd)
print(table(objd[['pass_scquantum_QC']]))

## examine doublets by merging subclones
source('./rsrc/annotate_coda.s2_purple_propose_doublets.R')
update_coda(df_meta, A = obja, D = objd)
write_coda(o = dir_res, df = df_meta, A = obja, D = objd)
print(table(clones=objd[['clones']]))

## subset the good aneuploid cells
source('./rsrc/annotate_coda.s2_purple_select_good_cells.R')
cname_good_aneu <- colnames(objd)[colData(objd)$is_good_aneu == TRUE]; str(cname_good_aneu)
subset_coda(df_meta, obja, objd, cname_good_aneu)
update_coda(df_meta, A = obja, D = objd)
objd <- runDistMat(objd, metric = 'manhattan', n_threads = 20)
update_coda(df_meta, A = obja, D = objd)
dir_res <- file.path(dir_proj, 'CNA_processing_ploidy_class', 
                     'aneuploid_pass_scquantum')
fs::dir_create(dir_res)
write_coda(o=dir_res, df=df_meta, A=obja, D=objd)

coda_dimplot2(df_meta, 'celltypes', pal = 'circus')

## re-clustering the good aneuploid cells of DNA [optional]
dir_res <- file.path(dir_proj, 'CNA_processing_ploidy_class', 
                     'aneuploid_pass_scquantum')
# load_coda(dir_res)
dir_snippet <- dir_res
source('./rsrc/snippet.CNA_clustering_in_aneuploid.R')
update_coda(df_meta, A = obja, D = objd)
z_tofig_opts <- c('subclones')
source('./rsrc/snippet.coda_viz.R')
objd <- findOutliers(objd, k = 2, resolution = 0.95)
pdf(file.path(dir_res, 'findOutliers.pdf'), width = 9, height = 6)
plotHeatmap(objd, label = 'outlier', row_split = 'outlier', n_threads = 40, 
            order_cells = 'consensus_tree', 
            col = pal_copykit_ratio_clip1)
dev.off()
table(objd@colData$outlier)
str(colnames(objd)[objd@colData$outlier == 'FALSE'])
subset_coda(df_meta, obja, objd, colnames(objd)[objd@colData$outlier == 'FALSE'])

## re-subclustering the clusters if needed
subclone_to_resubcluster <- NULL 
if (sample_name == c('DCIS66T_chip2') ) {subclone_to_resubcluster <- c('c0')}
if (!is.null(subclone_to_resubcluster)) {
  objd <- copykit_find_sub_clusters(
    objd, subcluster_by = 'subclones', cluster = subclone_to_resubcluster, 
    saveto_col_str = 'newsubclusters', k_nb = 2)
  table(objd[['newsubclusters']])
  plotUmap(objd, label = 'newsubclusters')
  objd[['subclones']] <- objd[['newsubclusters']]
  objd[['newsubclusters']] <- NULL
  objd[['subclones']] <- pretty_seamless_factor(
    objd[['subclones']], prefix='C'
  )
  plotUmap(objd, label = 'subclones')
  update_coda(df_meta, A = obja, D = objd)
  write_coda(o = dir_res, df = df_meta, A = obja, D = objd)
  z_tofig_opts <- c('subclones')
  dir_snippet <- dir_res
  source('./rsrc/snippet.coda_viz.R')
}

write_coda(o=dir_res, df=df_meta, A=obja, D=objd)

## merge similar subclones
dir_snippet <- dir_res <- file.path(
  dir_proj, 'CNA_processing_ploidy_class', 
  'aneuploid_pass_scquantum')
message('!manually examine and merge subclones!')
source('./rsrc/annotate_coda.s2_purple_merge_similar_subclones.R') 
print(table(colData(objd)[, 'clones']))
update_coda(df_meta, obja, objd)

z_tofig_opts <- c('clones')
source('./rsrc/snippet.coda_figure.R')
objd <- calcConsensus(objd, consensus_by = 'clones')
objd <- runConsensusPhylo(objd)
update_coda(df_meta, obja, objd)
write_coda(o=dir_res, df=df_meta, A=obja, D=objd)

#------ merge pure green cells and pure purple cells ------
## merge to './preready_coda/'
dir_green <- file.path(dir_proj, 'CNA_processing_ploidy_class', 'diploid')
dir_green <- file.path(dir_proj, 'CNA_processing_ploidy_class', 'diploid_rm_doublets')
dir_purple <- file.path(dir_proj, 'CNA_processing_ploidy_class', 'aneuploid_pass_scquantum')
dir_res <- file.path(dir_proj, 'preready_coda'); fs::dir_create(dir_res)
## merge CNA
objd_green <- read_rds(file.path(dir_green, 'objd.rds'))
objd_purple <- read_rds(file.path(dir_purple, 'objd.rds'))
print(table(objd_green[['clones']]))
print(table(objd_purple[['clones']]))
tmp <- setdiff( colnames(colData(objd_purple)), colnames(colData(objd_green)) )
for (tmp_x in tmp) {colData(objd_purple)[[tmp_x]] <- NULL}
tmp <- setdiff( colnames(colData(objd_green)), colnames(colData(objd_purple)) )
for (tmp_x in tmp) {colData(objd_green)[[tmp_x]] <- NULL}
rm(tmp); rm(tmp_x)

DietCopykit <- function(x){
  reducedDims(x) <- NULL
  x
}
objd_green <- DietCopykit(objd_green)
objd_purple <- DietCopykit(objd_purple)

objd <- cbind(objd_green, objd_purple)
objd <- runDistMat(objd, metric = 'manhattan', n_threads = 50) 
## merge ATAC
obja_green  <- read_rds(file.path(dir_green,  'obja.rds'))
obja_purple <- read_rds(file.path(dir_purple, 'obja.rds'))
obja <- merge(obja_green, obja_purple)

## merge df_meta
df_meta_green <- read_rds(file.path(dir_green, 'metadata.df.rds'))
df_meta_purple <- read_rds(file.path(dir_purple, 'metadata.df.rds'))
tmp <- intersect(colnames(df_meta_green), colnames(df_meta_purple))
df_meta <- rbind(df_meta_green[, tmp], df_meta_purple[, tmp])

## reclustering CNA
source('./rsrc/snippet.CNA_clustering_in_aneuploid.R')

## reclustering ATAC
source('./rsrc/snippet.ATAC_clustering.R')

## polish the metadata column 'clones' and 'celltypes'
objd[['clones']] <- fct_relevel(
  objd[['clones']], 
  gtools::mixedsort( levels(objd[['clones']]) ))
table(objd[['clones']])
obja$celltypes <- fct_relevel(
  obja$celltypes, 
  levels(df_meta$celltypes)
)
obja$seurat_clusters <- obja$celltypes

table(obja$seurat_clusters)
table(obja$clones)
# 
# coda_dimplot2(df_meta, 'ploidy_type_final', 'circus')
# coda_dimplot2(df_meta, 'clones', 'stallion',
#               # plot_order = levels(df_meta),
#               color_specific = c('D' = 'palegreen'))

## write coda
objd[['ploidy_type_final']] <- as.factor(ifelse(objd[['clones']]=='D', 'D', 'A'))
table(objd[['ploidy_type_final']])

update_coda(df=df_meta, A=obja, D=objd)

dir_res <- file.path(dir_proj, 'preready_coda'); fs::dir_create(dir_res)
write_coda(dir_res, df = df_meta, A = obja, D = objd)
load_coda(dir_res)
dir_snippet <- dir_res
z_tofig_opts <- c('clones', 'celltypes')
source('./rsrc/snippet.coda_figure.R') ##

load_coda(dir_res)
#------------------- ~~~ 3. Examine ATAC and CNA annotation~~~ -------------------  
#------ work on creating 'rds_coda_ready' object ------
dir_res <- file.path('./rds_coda_ready', sample_name); fs::dir_create(dir_res)
dir_snippet <- dir_res

#------ remove the noise cells ------
coda_dimplot2(df_meta, 'celltypes', pal = 'circus', pt.size = 0.5, color_specific = c('Noise' = 'magenta'))
Idents(obja) <- 'celltypes'
obja <- subset(obja, idents = setdiff(levels(Idents(obja)), 'Noise'))
table(Idents(obja))
obja$celltypes <- obja$seurat_clusters <- Idents(obja)
subset_coda(df_meta, obja, objd, cells = Cells(obja))
update_coda(df_meta, obja, objd)
write_coda(dir_res, df = df_meta, A = obja, D = objd)
#------ re-annotate cell types if needed ------
source('./rsrc/snippet.ATAC_clustering.R')
## manual work
ident_atac_clusters <- NULL
source('./rsrc/annotate_conda.s3_celltying_refine.R')
# idents absorb celltypes
ruok::qtable_scatter_hclust(
  pred  = obja@meta.data[, 'celltypes'],
  groundtruth = obja@meta.data[, ident_atac_clusters])
absorb_labels <- function(query, ref) {
  tab <- as.data.frame.matrix(table(query, ref))
  ptab <- prop.table(tab)
  ref_l <- colnames(ptab)
  query_l <- rownames(ptab)
  j <- apply(ptab, 1, nnet::which.is.max)
  res <- structure(names = query_l, ref_l[j])
}
dict_ident2celltypes <- absorb_labels(
  obja@meta.data[, ident_atac_clusters], 
  obja@meta.data[, 'celltypes']
)
if ('0' %in% names(dict_ident2celltypes)) {dict_ident2celltypes['0'] <- 'Noise'}
dict_ident2celltypes
obja$celltypes_before_refine <- obja$celltypes
obja$celltypes <- forcats::fct_recode(
  obja@meta.data[, ident_atac_clusters], 
  !!!structure(names = dict_ident2celltypes, as.character(names(dict_ident2celltypes)))
)
Idents(obja) <- 'celltypes'
obja <- subset(obja, idents=setdiff(levels(Idents(obja)), 'Noise'))
obja$celltypes <- Idents(obja)
obja$celltypes <- forcats::fct_relevel(
  obja$celltypes, levels(obja$celltypes_before_refine))
table(obja$celltypes, useNA='ifany')
obja$seurat_clusters <- obja$celltypes

subset_coda(df_meta, A = obja, D = objd, cells = Cells(obja))
update_coda(df_meta, obja, objd)
write_coda(dir_res, df = df_meta, A = obja, D = objd)

#------ reclustering DNA for a new UMAP ------
objd <- runDistMat(objd, metric = 'manhattan', n_threads = 50)
source('./rsrc/snippet.CNA_clustering_in_aneuploid.R')
update_coda(df_meta, obja, objd)
write_coda(dir_res, df = df_meta, A = obja, D = objd)
#------ re-annotate / merge similar clones if needed  [optional] ------
source('./rsrc/annotate_coda.s3_merge_similar_subclones.R') 
update_coda(df_meta, obja, objd)
write_coda(dir_res, df = df_meta, A = obja, D = objd)
#------ viz & export ------
dir_snippet <- dir_res
z_tofig_opts <- c('clones')#,  'seurat_clusters', 'celltypes')
source('./rsrc/snippet.coda_figure.R')

#------------------- ~~~ 4. Separate into A & D for viz ~~~ -------------------  
message('DONE annotate conda pipe.')
dir_res <- file.path('./rds_coda_ready', sample_name); fs::dir_create(dir_res)
load_coda(dir_res)
## backup
objd0 <- objd; df_meta0 <- df_meta; obja0 <- obja
types_split_to <- c('aneuploid_epi', 'aneuploid_nonepi', 'diploid')
# types_split_to <- c('aneuploid_nonepi', 'diploid')
# types_split_to <- c('diploid')
# types_split_to <- c('aneuploid_epi')
## split to the several halves: 
## aneuploid_epi = C* & LumSec|LumHR|Basal
## aneuploid_nonepi = C*|CD* & !LumSec|LumHR|Basal
## diploid = D
for (sp in types_split_to) {
  message('splitting into ', sp)
  dir_res <- file.path('./rds_coda_ready', sample_name, sp)
  fs::dir_create(dir_res)
  cells_wanted <- NULL
  epi_names <- c('LumHR', 'LumSec', 'Basal', 'TBD2')
  if (sp == 'aneuploid_epi') {
    idx <- str_detect(df_meta0$clones, 'C[0-9]+') & df_meta0$celltypes %in% epi_names
  } else if (sp == 'aneuploid_nonepi') {
    idx <- str_detect(df_meta0$clones, 'CD[0-9]+') & !df_meta0$celltypes %in% epi_names
  } else if (sp == 'diploid') {
    idx <- c(df_meta0$clones %in% 'D')
  }
  # load_coda(dir_res)
  # dir_snippet <- dir_res
  # source('./rsrc/snippet.coda_viz.R')
  # next()
  message(sp, ' has ', sum(idx), ' cells...')
  if (sum(idx) > 0) {
    subset_coda(df = df_meta0, A = obja0, D = objd0, 
                cells = rownames(df_meta0)[idx])
    objd <- runDistMat(objd, metric = 'manhattan', n_threads = 20) 
    write_coda(o = dir_res, df = df_meta, A = obja, D = objd)
    ## recluster CNA for a new UMAP
    ## recluster ATAC for a new UMAP
    ## viz each compartment 
    dir_snippet <- dir_res
    z_tofig_opts <- c('clones')#,  'seurat_clusters', 'celltypes')
    source('./rsrc/snippet.coda_figure.R')
  }
  
}



for (sp in c('aneuploid_epi')) {
  message('splitting into ', sp)
  dir_res <- file.path('./rds_coda_ready', sample_name, sp)
  load_coda(dir_res)
  dir_snippet <- dir_res
  z_tofig_opts <- c('clones')#,  'seurat_clusters', 'celltypes')
  source('./rsrc/snippet.coda_figure.R')
  
}


