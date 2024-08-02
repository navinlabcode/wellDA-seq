## Written by: Yun Yan (https://github.com/Puriney)\n\n
# Load library
# DNA object is created from uber files.
# 
library(copykit)
library(BiocParallel)
library(RhpcBLASctl)
library(ComplexHeatmap); library(circlize)
library(tictoc)
tictoc::tic()
register(MulticoreParam(progressbar = T, workers = 20), default = T)
RhpcBLASctl::blas_set_num_threads(20)
source('/volumes/USR1/yyan/project/coda/rsrc/utils.R')

options <- commandArgs(trailingOnly = TRUE)

if (length(options) == 3) {
  indir <- options[[1]]#"/volumes/USR1/yyan/project/coda/rds/DCIS51T_chip1_low_dna/res_200_k_reg/final_result"
  dir_res <- options[[2]]#'/volumes/USR1/yyan/project/coda/rds/DCIS51T_chip1_low_dna'
  run_name <- options[[3]]#'DCIS51T_chip1_low_dna'
}

library(tidyverse)
library(readr)
#-------------------------- START --------------------------  

#-------------------------- I --------------------------  
fs::dir_create(dir_res)
# Run pre-processing module
list.files(indir)
seg <- read.table(file.path(indir, 'uber.sample.seg.txt'), header = T)
rat <- read.table(file.path(indir, 'uber.sample.ratio.txt'), header = T)
bin <- read.table(file.path(indir, 'uber.sample.bin.txt'), header = T)
hg19_ranges <- makeGRangesFromDataFrame(copykit::hg19_rg, keep.extra.columns = TRUE)
dim(rat)
# remove Y chr
idx <- seg$chrom == '24'; table(idx); seg <- seg[!idx, ]
idx <- rat$chrom == '24'; table(idx); rat <- rat[!idx, ]
idx <- bin$chrom == '24'; table(idx); bin <- bin[!idx, ]
idx <- seqnames(hg19_ranges) == 'chrY'; table(idx); hg19_ranges <- hg19_ranges[!idx, ]
rm(idx)

tumor <- CopyKit(
  assay = list(bincounts = bin[,-c(1:3)], 
               segment_ratios = seg[,-c(1:3)], 
               ratios = rat[,-c(1:3)]),
  rowRanges = hg19_ranges)
print(tumor)
tumor <- runVst(tumor, transformation = 'ft')
tumor <- logNorm(tumor); print(tumor)
print(tumor)
metadata(tumor)$genome <- 'hg19'
metadata(tumor)$resolution <- "220kb"
colData(tumor)$sample <- colnames(tumor)
dim(tumor)
# idx <- str_detect(colnames(tumor), pattern = '.markdup'); table(idx)
# tumor <- tumor[, idx]
tumor <- runMetrics(tumor)

write_rds(tumor, file.path(dir_res, 'all.copykit.rds'))

#------ Remove the troublesome suffix ------
colnames(tumor) <- colnames(tumor) %>% str_split(pattern = '_', n=3) %>%
  lapply(., function(x) paste(x[1:2], collapse = '_')) %>%
  unlist()

#------ Basi processing ------
  
tumor <- read_rds(file.path(dir_res, 'all.copykit.rds'))
# Mark euploid cells if they exist
tumor <- findAneuploidCells(tumor)
tumor <- findOutliers(tumor, resolution = 0.8)


write_rds(tumor[, colData(tumor)$outlier == 'FALSE'], 
          file.path(dir_res, 'findOutliers0.8.copykit.rds'))

# findOutliers_res = 0.95
findOutliers_res = 0.9 ## all samples are using this. 
tumor <- findOutliers(tumor, resolution = findOutliers_res) ## high corr to make outliers
write_rds(tumor, file.path(dir_res, 'all.copykit.rds'))
# Mark low-quality cells for filtering
tumor <- runDistMat(tumor, metric = 'manhattan', n_threads = 20)
colData(tumor)$median_bincounts <- apply(bincounts(tumor), 2, median)

write_rds(tumor, file.path(dir_res, 'all.copykit.rds') )

#------ attach pipeline's metric ------
colnames(tumor)
colnames(colData(tumor))
meta_df <- read_tsv(file.path(dirname(indir), 'metrics', 'all_stat_metrics.txt'))
idx <- match(rownames(colData(tumor)), meta_df$`Sample Name`)
meta_df <- meta_df[idx, ]
meta_df$`Sample Name` <- NULL
colData(tumor) <- cbind(colData(tumor), meta_df)
range(colData(tumor)$MedianBinCount)
colnames(colData(tumor))
median( colData(tumor)$MedianBinCount )
table( colData(tumor)$MedianBinCount >10 )

for (iz in c('TotalReads', 'DupsRemoved', 'ReadsKept', 'MedianBinCount')) {
  pdf(file.path(dir_res, sprintf('qc.%s.pdf', iz)),
      width = 5, height = 5, useDingbats = F)
  print(plotMetrics(tumor, iz))
  rm(iz)
  dev.off()
}

colData(tumor)$pass_median_bincounts <- colData(tumor)$MedianBinCount>min_bincount
table(colData(tumor)$pass_median_bincounts)
write_rds(tumor, file.path(dir_res, 'all.copykit.rds') )

tumor <- read_rds(  file.path(dir_res, 'all.copykit.rds') )

# pal_copykit_ratio_clip2 <- circlize::colorRamp2(
#   seq(from=-2, to=2, length.out=11), 
#   hcl.colors(n=11, palette = 'Blue - Red 3'))
pal_copykit_ratio_clip1 <- circlize::colorRamp2(
  seq(from=-1, to=1, length.out=11), 
  hcl.colors(n=11, palette = 'RdBu', rev = T))

# Visualize cells labeled by filter and aneuploid status
pdf(file.path(dir_res, 'plotHeatmap.outlier.pdf'), width = 10, height = 10)
plotHeatmap(
  tumor, 
  label = intersect(c('outlier', 'pass_median_bincounts', 'is_aneuploid'), colnames(colData(tumor))), 
  row_split = 'outlier', n_threads = 20, 
  col = pal_copykit_ratio_clip1) %>% draw()
dev.off()
table( colData(tumor)$median_bincounts >10,  colData(tumor)$is_aneuploid)
# [optional] Remove cells marked as low-quality and/or aneuploid from the copykit object
# tumor <- tumor[,SummarizedExperiment::colData(tumor)$outlier == FALSE]
# tumor <- tumor[,SummarizedExperiment::colData(tumor)$is_aneuploid == TRUE]

#-------------------------- II --------------------------  

# Create a umap embedding 
tumor <- copykit::runPca(tumor)
## standardized
tumor <- runUmap(tumor, n_neighbors=50, n_threads=20)
plotUmap(tumor)

if (T) {
  # Default
  k_superclones <- 20; k_subclones <- 10
  findClusters_embedding <- 'umap'
  findClusters_method <- 'hdbscan'
  findClusters_ncomponents <- 2
}


tumor  <- findClusters(
  tumor, 
  k_superclones = k_superclones,
  k_subclones = k_subclones,
  embedding = findClusters_embedding, 
  method = findClusters_method, 
  ncomponents = findClusters_ncomponents)

pdf(file.path(dir_res, 'plotUmap.subclones.pdf'), width = 4, height = 3)
print(plotUmap(tumor, label = 'subclones') + theme(aspect.ratio = 1))
dev.off()

pdf(file.path(dir_res, 'plotUmap.superclones.pdf'), width = 4, height = 3)
print(plotUmap(tumor, label = 'superclones') + theme(aspect.ratio = 1))
dev.off()

pdf(file.path(dir_res, 'plotUmap.is_aneuploid.pdf'), width = 4, height = 3)
print(plotUmap(tumor, label = 'is_aneuploid') + theme(aspect.ratio = 1))
dev.off()

pdf(file.path(dir_res, 'plotUmap.pass_median_bincounts.pdf'), width = 4, height = 3)
print(plotUmap(tumor, label = 'pass_median_bincounts') + theme(aspect.ratio = 1))
dev.off()

# Calculate consensus profiles for each subclone, 
# and order cells by cluster for visualization with plotHeatmap
tumor <- calcConsensus(tumor)
tumor <- runConsensusPhylo(tumor)
tumor <- runPhylo(tumor)
# Run scquantum
try(tumor <- calcInteger(tumor, method = 'scquantum'))
try(colData(tumor2)$ploidy_confidence <- colData(tumor2)$confidence_ratio)
if (any(c('ploidy_confidence', 'ploidy_score', 'confidence_ratio') %in% colnames(colData(tumor)))) {
  cat('[good] scquantum done.\n')
}

write_rds(tumor, file.path(dir_res, 'ready.copykit.rds'))

# Plot a copy number heatmap with clustering annotation
pdf(file.path(dir_res, 'plotHeatmap.clones.%03d.pdf'), 
    width = 10, height = 10, onefile = F, useDingbats = F)
plotHeatmap(
  tumor, label = c('subclones', 'superclones', 'pass_median_bincounts', 'is_aneuploid'), 
  order_cells = 'consensus_tree', 
  col = pal_copykit_ratio_clip1 ) %>% draw()
plotHeatmap(
  tumor, label = c('subclones', 'superclones', 'pass_median_bincounts', 'is_aneuploid'), 
  row_split = 'subclones',
  col = pal_copykit_ratio_clip1,
  order_cells = 'consensus_tree') %>% draw()
plotHeatmap(
  tumor, label = c('subclones', 'superclones', 'pass_median_bincounts', 'is_aneuploid'), 
  row_split = 'superclones',
  col = pal_copykit_ratio_clip1,
  order_cells = 'consensus_tree') %>% draw()
dev.off()

