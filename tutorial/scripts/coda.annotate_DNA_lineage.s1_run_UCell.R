#---------------------------
# Annotate DNA lineage
# 1. quantify any signature for single-cell using UCELL for the ciceroRNA (preferred) or the iRNA assay
# 2. visualize signature on the DNA lineage tree              ----    

# Can use either iRNA or ciceroRNA
# Can use either UCell or Seurat's addmodulescore
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

options <- commandArgs(trailingOnly = TRUE)
if (length(options)>0) {
  cat('Reading user parameters:\n')
  sample_name <- options[[1]]
  rna_assay_use <- options[[2]]
  module_method <- options[[3]]
} else {
  sample_name = 'DCIS35T'
  rna_assay_use <- 'iRNA'
  rna_assay_use <- 'ciceroRNA'
  module_method <- 'ucell'
  module_method <- 'seurat'
}

#------ input / output ------
  
dir_coda <- file.path(
  '/volumes/USR1/yyan/project/coda', 
  'rds_coda_ready', sample_name, 'aneuploid_epi')


load_coda(dir_coda) # obja; objd; df_meta
# consensus_z <- 'clones'
print(obja)

#------ load ciceroRNA ------
f_crna <- file.path(dir_coda, 'cicero', 'ciceroRNA_assay.seurat.rds')
crna <- try(read_rds(f_crna))

#------ choose RNA assay ------
  
if ( 'try-error' %in% class(crna)) {
  DefaultAssay(obja) <- 'iRNA'
} else {
  stopifnot( identical(Cells(crna), Cells(obja)) )
  obja[['ciceroRNA']] <- crna[['ciceroRNA']]
  DefaultAssay(obja) <- 'ciceroRNA'
}

try(DefaultAssay(obja) <- rna_assay_use)

cat('ATAC uses RNA assay = ', DefaultAssay(obja), '\n')
print(obja)
#------ choose pathways ------
  
if (T) {
  library(msigdbr)
  db_msigdb_H <- msigdbr(species = "Homo sapiens", category = "H") %>% 
    dplyr::select(gs_name, gene_symbol)
  db_msigdb_H <- deframe_to_list(as.data.frame(db_msigdb_H))
  names(db_msigdb_H) <- make.names(names(db_msigdb_H))
}

db = db_msigdb_H
name_db <- 'msigdb_hallmark'

#------ setup dir_res ------

dir_res <- file.path(
  dir_coda, sprintf('%s_%s_%s', 
                    module_method,
                    name_db, 
                    DefaultAssay(obja)))
fs::dir_create(dir_res)

#------ run UCELL ------
if (module_method == 'ucell') {
  library(UCell)
  if (! file.exists(file.path(dir_res, 'addmodulescore.df.rds'))) {
    # if (T) { ## force to rerun
    message('Run AddModuleScore_UCell')
    obja <- AddModuleScore_UCell(obja, features = db, name="_UCell")
    names(db) <- paste0(names(db), '_UCell')  
    df_signature_res <- obja[[]][, names(db)]
    print(colnames(df_signature_res))
    write.csv(x=df_signature_res, file = file.path(dir_res, 'addmodulescore.df.csv'))
    write_rds(x=df_signature_res, file.path(dir_res, 'addmodulescore.df.rds'))
  } else {
    message('Load AddModuleScore_UCell')
    df_signature_res <- read_rds(file.path(dir_res, 'addmodulescore.df.rds'))
    stopifnot(identical(rownames(df_signature_res), Cells(obja)))
    obja <- AddMetaData(obja, metadata = df_signature_res)
    names(db) <- paste0(names(db), '_UCell')
  }
}  

#------ Run AddModulescore ------

if (module_method == 'seurat') {
  if (! file.exists(file.path(dir_res, 'addmodulescore.df.rds'))) {
    # if (T) { ## force to rerun
    message('Run AddModuleScore')
    obja <- AddModuleScore(obja, features = db, name='seuratmodule')
    for ( tmp in 1:length(db) ) {
      obja@meta.data[paste0(names(db)[tmp], '_seurat')] <- obja@meta.data[paste0('seuratmodule', tmp)]
      obja@meta.data[paste0('seuratmodule', tmp)] <- NULL
    }
    names(db) <- paste0(names(db), '_seurat')  
    df_signature_res <- obja[[]][, names(db)]
    print(colnames(df_signature_res))
    write.csv(
      x=df_signature_res, 
      file = file.path(dir_res, 'addmodulescore.df.csv'))
    write_rds(
      x=df_signature_res, 
      file.path(dir_res, 'addmodulescore.df.rds'))
  } else {
    message('Load AddModuleScore')
    df_signature_res <- read_rds(file.path(dir_res, 'addmodulescore.df.rds'))
    stopifnot(identical(rownames(df_signature_res), Cells(obja)))
    obja <- AddMetaData(obja, metadata = df_signature_res)
    names(db) <- paste0(names(db), '_seurat')
  }
}  



#------ viz Featureplot ------
update_coda(df_meta, obja, objd)

pdf(file.path(dir_res, sprintf('coda_featureplot.birdview_%s.pdf', name_db)), 
    width = 8, height = 5, onefile = T, useDingbats = F)
for (i in seq_along(names(db))) {
  print(coda_featureplot2(df = df_meta, feature = names(db)[i], pal='Plasma'))
}
dev.off()


#------ viz VlnPlot ------
pal_cna_clones = new_palette_D(levels(objd@colData$clones), pal = 'stallion')
n_clones <- length(unique(objd@colData$clones))
for (z in c('clones')) {
  message(z)
  pdf(file.path(dir_res, sprintf('vlnplot.birdview_%s.%s.pdf', name_db, z)), 
      width = n_clones*0.9, height = 4, useDingbats = F, onefile = T)
  for (i in seq_along(names(db))) {
    p <- VlnPlot(obja, feature = names(db)[i], 
                 group.by = z, pt.size = 0.5) &
      labs(x = z) &
      scale_fill_manual(values = pal_cna_clones) &
      stat_compare_means(label = "p.format") &
      stat_summary(fun = 'mean', colour = "cyan", size = 2, geom = "point", pch=4)
    print(p)

  }
  dev.off()
}
names(db)

cat('Done: ', sample_name, '-', rna_assay_use, '-', module_method, '\n')

