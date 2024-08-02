## Written by: Yun Yan (https://github.com/Puriney)\n\n
# 
# Remove outlier cells
# Remove non-dispensed cells
# Common Clustering
# 
suppressPackageStartupMessages({
library(copykit)
library(BiocParallel)
library(RhpcBLASctl)
library(ComplexHeatmap); library(circlize)
library(tictoc)

library(tidyverse)
library(readr)
library(ggpubr); library(ggplot2)
setwd("/volumes/USR1/yyan/project/coda")

register(MulticoreParam(progressbar = T, workers = 20), default = T)
RhpcBLASctl::blas_set_num_threads(20)
source('/volumes/USR1/yyan/project/coda/rsrc/utils.R')
})


options <- commandArgs(trailingOnly = TRUE)


min_bincount <- 8
is_paired_end <- T

if (length(options) == 4) {
  indir <- options[[1]]#"/volumes/USR1/yyan/project/coda/rds/DCIS51T_chip1_low_dna/res_200_k_reg/final_result"
  dir_res <- options[[2]]#'/volumes/USR1/yyan/project/coda/rds/DCIS51T_chip1_low_dna'
  run_name <- options[[3]]#'DCIS51T_chip1_low_dna'
  sample_type <- options[[4]] # cell_line/DCIS/HBCA
}


#-------------------------- Read in --------------------------  

timestamp()
tumor <- read_rds( file.path(dir_res, 'ready.copykit.rds') )
colnames(colData(tumor))
table( outlier=colData(tumor)$outlier )

#------ Attach physical wells ------
source('/volumes/USR1/yyan/project/coda/rsrc/utils.wafargen_physical.R')
cnames <- colnames(tumor) %>% str_split(pattern = '_', n=3) %>% 
  lapply(., function(x) paste(x[1:2], collapse = '_')) %>%
  unlist(); str(cnames)
if(!all(cnames %in% lib_wafar$well_name))  {
  warning('[Wrong] Not all cell names are in the dictionary. Maybe wrong cell names.\n')
} else {
  message('attaching physical wells information')
  idx <- match(cnames, lib_wafar$well_name)
  colData(tumor) <- cbind(colData(tumor), lib_wafar[idx, ])
  write_rds(tumor,  file.path(dir_res, 'ready.copykit.rds') )
}
#------ Viz cell dispense on wafer chip ------
run_srch_key <- str_remove(run_name, '_low') %>% str_remove(., '_dna'); message(run_srch_key)
wafar_use <- read_tsv(Sys.glob(sprintf('/volumes/USR1/yyan/project/coda/wafar/%s/*.TXT', run_srch_key))); dim(wafar_use)
if (!purrr::is_empty(wafar_use)) {
  wafar_wid_use <- paste0(wafar_use$Row, '_', wafar_use$Col)
  
  wid_use <- lib_wafar$well_id[match(cnames, lib_wafar$well_name)]; str(wid_use)
  p <- VennDiagram::venn.diagram(
    list(DNA=wid_use, Dispense=wafar_wid_use),
    col = c('blue', 'gold'), cat.col = c('blue', 'gold'),
    cat.default.pos = "text",
    filename = NULL, disable.logging=T)
  pdf(file.path(dir_res, 'VennDiagram.well_id.pdf'), width = 5, height = 5, useDingbats = F)
  grid.draw(p)
  dev.off()
  wafar_df <- create_wafar_loc_dataframe(wid_use, wafar_wid_use)
  pa <- plot_wafar_capture(wafar_df) + labs(title=run_name); print(pa)
  ggsave(file.path(dir_res, 'wafargen_physical_dispense.pdf'), pa, width = 7, height = 7, useDingbats = F)
  
  str(setdiff(colData(tumor)$well_id, wafar_wid_use))
  
  colData(tumor)$is_dispensed <- c(colData(tumor)$well_id %in% wafar_wid_use); table(colData(tumor)$is_dispensed)
} else {
  message('assuming all cells are dispensed. ')
  colData(tumor)$is_dispensed <- TRUE  
}

write_rds(tumor,  file.path(dir_res, 'ready.copykit.rds') )
if (sum(!colData(tumor)$is_dispensed) >0){
  print(table(colData(tumor)$is_dispensed))
  cat('Some cells are false-positive!!.\n')
}

print(table(outlier=colData(tumor)$outlier, 
            is_dispensed=colData(tumor)$is_dispensed))
if ('is_dispensed' %in% colnames(colData(tumor))) {
  print(table(is_dispensed=colData(tumor)$is_dispensed))
  tumor2 <- tumor[, colData(tumor)$is_dispensed]
} else {
  tumor2 <- tumor
}

#------ for cell line only ------
if (sample_type %in% c('cell_line')) {
  tumor2 <- tumor[, colData(tumor)$outlier == 'FALSE']
  cat('--- removed', ncol(tumor) - ncol(tumor2), 'outlier cells.\n')
}
dim(tumor2)
write_rds( tumor2, file.path(dir_res, 'ready_clean.copykit.rds') )

#------ Run scquantum if absent ------
if (!any(c('ploidy_confidence', 'ploidy_score') %in% colnames(colData(tumor2))))  {
  tumor2 <- calcInteger(tumor2, method = 'scquantum')  
  write_rds( tumor2, file.path(dir_res, 'ready_clean.copykit.rds') )
} else {
  cat('scquantum already exist.\n')
}


tumor2 <- read_rds( file.path(dir_res, 'ready_clean.copykit.rds') )
try(colData(tumor2)$ploidy_confidence <- colData(tumor2)$confidence_ratio)

table(superclones=colData(tumor2)$superclones)
table(subclones=colData(tumor2)$subclones)
p <- ggscatter(as.data.frame(colData(tumor2)),
               x='ploidy', y='ploidy_confidence',
               facet.by= 'subclones') +
  geom_hline(yintercept = 1, lty='dashed', color='blue')  +
  geom_vline(xintercept = 2, lty='dashed', color='blue') 
print(p)
ggsave(file.path(dir_res, 'qc.scquantum.ploidy_by_subclones.pdf'), p, width = 8,height = 8, useDingbats = F)



#------ Examine doublets and viz on wafer chip [no big use] ------
pdf(file.path(dir_res, 'qc.wafergene_superclones.pdf'), width = 7, height = 7, useDingbats = F, onefile = T)
plot_wafar(df=as.data.frame(colData(tumor2)), z = 'superclones') %>% print()
for (z_hi in unique(colData(tumor2)[['superclones']])) {
  plot_wafar(df=as.data.frame(colData(tumor2)), z = 'superclones', z_hi = z_hi) %>% print()
}
dev.off()

#------ re-processing ------
source('/volumes/USR1/yyan/project/coda/rsrc/utils.R')
if (TRUE) {
  # Default
  k_superclones <- 20; k_subclones <- 10
  findClusters_embedding <- 'umap'
  findClusters_method <- 'hdbscan'
  findClusters_ncomponents <- 2
}
if (ncol(tumor) != ncol(tumor2)) {
  # k_superclones = 25
  cat(k_superclones, k_subclones)
  # cat('Reclustering because subseted cells.\n')
  # findClusters_embedding <- 'PCA'
  # findClusters_method <- 'leiden'
  # findClusters_ncomponents <- 50
  cat('Regenerate UMAP because of subseting cells\n')
  tumor2 <- copykit::runPca(tumor2)
  tumor2 <- runUmap(tumor2, n_neighbors=k_subclones, n_threads=20)
  tumor2  <- findClusters(
    tumor2,
    k_superclones = k_superclones,
    k_subclones = k_subclones,
    embedding = findClusters_embedding,
    method = findClusters_method,
    ncomponents = findClusters_ncomponents)
  print(table(colData(tumor2)$superclones))
  print(table(colData(tumor2)$subclones))
}

tumor2 <- calcConsensus(tumor2, consensus_by='subclones')
tumor2 <- runConsensusPhylo(tumor2)
tumor2 <- runPhylo(tumor2)

if (! 'clones' %in% colnames(colData(tumor2)) ) {
  # older objects
  colData(tumor2)$clones=colData(tumor2)$superclones  
}
if (! 'ploidy_class' %in% colnames(colData(tumor2)) ) {
  colData(tumor2)$ploidy_class=ifelse(colData(tumor2)$is_aneuploid, 'aneuploid', 'diploid')
}
write_rds(tumor2, file.path(dir_res, 'ready_clean.copykit.rds'))
tumor2 <- read_rds(file.path(dir_res, 'ready_clean.copykit.rds'))
metainfo <- as.data.frame(colData(tumor2)) %>% rownames_to_column('cname')
write_csv(metainfo, file.path(dir_res, 'ready_clean.copykit.metainfo.csv'))

#------ Viz: umap ------

for (z in c('subclones', 'superclones', 'pass_median_bincounts', 
            'is_aneuploid', 'clones', 'ploidy_class')) {
  pdf(file.path(dir_res, sprintf('clean.plotUmap.%s.pdf', z)), width = 4, height = 3)
  p <- plotUmap(tumor2, label = z) + theme(aspect.ratio = 1)
  try(print(p))
  dev.off()
}

#------ Viz: heatmap ------
pdf(file.path(dir_res, 'clean.plotHeatmap.clones.%03d.pdf'), 
    width = 12, height = 10, onefile = F, useDingbats = F)
try(plotHeatmap(
  tumor2, label = c('clones', 'subclones', 'superclones', 
                    'pass_median_bincounts', 'ploidy_class'), 
  row_split = 'ploidy_class',
  col = pal_copykit_ratio_clip1,
  order_cells = 'consensus_tree') %>% draw()
)
try(plotHeatmap(
  tumor2, label = c('clones', 'subclones', 'superclones', 
                    'pass_median_bincounts', 'ploidy_class'), 
  col = pal_copykit_ratio_clip1,
  row_split = 'subclones',
  order_cells = 'consensus_tree') %>% draw()
)
try(plotHeatmap(
  tumor2, label = c('clones', 'subclones', 'superclones', 
                    'pass_median_bincounts', 'ploidy_class'), 
  col = pal_copykit_ratio_clip1,
  row_split = 'superclones',
  order_cells = 'consensus_tree') %>% draw()
)
try(plotHeatmap(
  tumor2, label = c('clones', 'subclones', 'superclones', 
                    'pass_median_bincounts', 'ploidy_class'), 
  col = pal_copykit_ratio_clip1,
  row_split = 'clones',
  order_cells = 'consensus_tree') %>% draw()
)
dev.off()

pdf(file.path(dir_res, 'clean.plotHeatmap_integer.clones.%02d.pdf'), 
    width = 10, height = 10, onefile=F)
try(plotHeatmap(
  tumor2, label = c('clones', 'subclones', 'superclones', 
                    'pass_median_bincounts', 'ploidy_class'), 
  col = pal_copykit_int_from_ratio_clip1, 
  assay = 'integer',
  order_cells = 'consensus_tree') %>% draw()
)
try(plotHeatmap(
  tumor2, label = c('clones', 'subclones', 'superclones', 
                    'pass_median_bincounts', 'ploidy_class'), 
  col = pal_copykit_int_from_ratio_clip1, 
  assay = 'integer',
  row_split = 'clones',
  order_cells = 'consensus_tree') %>% draw()
)
dev.off()

#------ DONE ------

tumor2 <- read_rds( file.path(dir_res, 'ready_clean.copykit.rds') )
print(tumor2)
cat('DONE')
