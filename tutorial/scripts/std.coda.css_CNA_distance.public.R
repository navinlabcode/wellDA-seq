#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# consensus CNA distance              ~~~~    
#~~~~~~~~~~~~~~~~~~~~~~~~~~~

suppressPackageStartupMessages({
  library(tidyverse); library(forcats)
  library(fs)
  library(copykit)
  library(VennDiagram)
  library(SummarizedExperiment)
  library(ggpubr); library(ruok); library(scales)
  library(ComplexHeatmap)
  library(ggtree)
  library(forcats)
  library(ape)
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

for (sample_name in c(names(pal_dcis_9samples))) {
# for (sample_name in c(sample_name)) {
  message(sample_name)
  dir_coda <- file.path(
    '/volumes/USR1/yyan/project/coda', 
    'rds_coda_ready', sample_name, 'aneuploid_epi')
  dir_res <- file.path(dir_coda, 'css_distance')
  fs::dir_create(dir_res)
  #------ manhaten distance ------
  
  ## using css mat from the MEDICC2 because I polished the CSS
  
  css_mat <- read_rds(file.path(dir_coda, 'medic2_polish', 'css_mat.rds'))
  print(dim(css_mat))
  print(colnames(css_mat))
  print(range(css_mat))
  
  css_dist_manhattan <- dist(t(css_mat), method = 'manhattan')
  
  tree <- ape::fastme.bal(css_dist_manhattan)
  # plot(tree)
  
  pdf(file = file.path(dir_res, 'plot.fastme_tree.pdf'), width = 5,height = 5, useDingbats = F)
  plot(tree %>% ladderize())
  dev.off()
  
  write_rds(tree, 
            file.path(dir_res, 'fastme_tree.rds'))
  
  css_dist_manhattan <- as.matrix(css_dist_manhattan)
  write_rds(
    css_dist_manhattan, 
    file.path(dir_res, 'css_integerCN_dist_manhattan.mat.rds'))
  rm(tree); rm(css_mat); rm(css_dist_manhattan)
}
