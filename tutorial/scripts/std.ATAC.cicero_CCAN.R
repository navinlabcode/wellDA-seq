#---------------------------
# Find cis-co-accessible networks (or chromatin hubs)              ----
# 
# - CCAN assay
# - dictionary: peaks ~ CCAN
#     
#---------------------------
setwd('/volumes/USR1/yyan/project/coda')
if (T) {
  library(Signac); library(copykit); library(Seurat)
  library(SeuratWrappers)
  library(tidyverse)
  library(MatrixGenerics); library(Matrix)
  library(ComplexHeatmap); library(ggpubr)
  library(pbapply)
  source('./rsrc/utils.R')
  source('./rsrc/utils.copykit.R')
  # Install Cicero
  # if (!requireNamespace("remotes", quietly = TRUE))
  #   install.packages("remotes")
  # devtools::install_github('cole-trapnell-lab/monocle3')
  # remotes::install_github("cole-trapnell-lab/cicero-release", ref = "monocle3")
  library(cicero)
  library(monocle3); library(SingleCellExperiment)
  
  ## Literally copied from https://rdrr.io/github/satijalab/seurat-wrappers/src/R/monocle3.R
  source("https://raw.githubusercontent.com/satijalab/seurat-wrappers/master/R/internal.R")
  as_cell_data_set <- function(
    x,
    assay = DefaultAssay(object = x),
    reductions = AssociatedDimReducs(object = x, assay = assay),
    default.reduction = DefaultDimReduc(object = x, assay = assay),
    graph = paste0(assay, '_snn'),
    group.by = NULL,
    ...
  ) {
    clusters.key <- 'monocle3_clusters'
    partitions.key <- 'monocle3_partitions'
    # Add assay data
    # Cheat and use as.SingleCellExperiment
    cds <- as(
      object = as.SingleCellExperiment(x = x, assay = assay),
      Class = 'cell_data_set'
    )
    # Ensure we have a counts assay
    if (is.null(x = SummarizedExperiment::assays(x = cds)$counts)) {
      SummarizedExperiment::assays(x = cds)$counts <- SummarizedExperiment::assays(x = cds)[[1]]
    }
    # Add Size_factor
    if (!"Size_Factor" %in% colnames(x = SummarizedExperiment::colData(x = cds))) {
      size.factor <- paste0('nCount_', assay)
      if (size.factor %in% colnames(x = x[[]])) {
        SummarizedExperiment::colData(x = cds)$Size_Factor <- x[[size.factor, drop = TRUE]]
      }
    }
    # Add DimReducs: Embeddings become a reduced dim, Loadings go to
    # reduce_dim_aux$gene_loadings, Stdev goes go reduce_dim_aux$prop_var_expl
    # First, reset the ones from as.SingleCellExperiment
    SingleCellExperiment::reducedDims(x = cds)[SingleCellExperiment::reducedDimNames(x = cds)] <- NULL
    reductions <- intersect(
      x = reductions,
      y = AssociatedDimReducs(object = x, assay = assay)
    )
    for (reduc in reductions) {
      SingleCellExperiment::reducedDims(x = cds)[[toupper(x = reduc)]] <- Embeddings(object = x[[reduc]])
      loadings <- Loadings(object = x[[reduc]])
      if (!IsMatrixEmpty(x = loadings)) {
        slot(object = cds, name = 'reduce_dim_aux')[['gene_loadings']] <- loadings
      }
      stdev <- Stdev(object = x[[reduc]])
      if (length(x = stdev)) {
        slot(object = cds, name = 'reduce_dim_aux')[['prop_var_expl']] <- stdev
      }
    }
    # Add clustering information
    # TODO: Figure out if I need to add relations, distMatrix, or clusters/partitions
    if (!is.null(x = group.by)) {
      Idents(object = x) <- group.by
    }
    # if (clusters.key %in% colnames(x = x[[]])) {
    clusters.list <- if (is.null(x = group.by) && all(c(clusters.key, partitions.key) %in% colnames(x = x[[]]))) {
      message("Using existing Monocle 3 cluster membership and partitions")
      list(
        partitions = factor(x = x[[partitions.key, drop = TRUE]]),
        clusters = factor(x = x[[clusters.key, drop = TRUE]])
      )
    } else if (graph %in% names(x = x)) {
      g <- igraph::graph_from_adjacency_matrix(
        adjmatrix = x[[graph]],
        weighted = TRUE
      )
      # TODO: figure out proper partitioning scheme
      # partitions <- igraph::components(graph = g)$membership[colnames(x = x)]
      warning(
        "Monocle 3 trajectories require cluster partitions, which Seurat does not calculate. Please run 'cluster_cells' on your cell_data_set object",
        call. = FALSE,
        immediate. = TRUE
      )
      partitions <- rep_len(x = 1, length.out = ncol(x = x))
      list(
        cluster_result = list(
          g = g,
          relations = NULL,
          distMatrix = 'matrix',
          coord = NULL,
          edge_links = NULL,
          optim_res = list(
            membership = as.integer(x = Idents(object = x)),
            modularity = NA_real_
          )
        ),
        partitions = factor(x = partitions),
        clusters = Idents(object = x)
      )
    } else {
      list()
    }
    if (length(x = clusters.list)) {
      slot(object = cds, name = 'clusters')[[toupper(x = default.reduction)]] <- clusters.list
    }
    # TODO: Add translated results from learn_graph
    return(cds)
  }
}
if (F) {
  library(future)
  plan("multisession", workers = 8)
  plan()
  options(future.globals.maxSize = 50 * 1024 ^ 3) # for 50 Gb RAM
}

coaccess_cutoff_override <- 0

# sample_name <- 'DCIS66T_chip2'
sample_name <- 'DCIS22T'
sample_name <- 'DCIS35T'
# sample_name = 'mda231wt_kaile'
args  <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) { 
  sample_name = args[[1]] 
  coaccess_cutoff_override = as.numeric(args[[2]]) # 0:use default; otherwise use coaccess_cutoff_override
}

message(sample_name)
dir_coda <- file.path("/volumes/USR1/yyan/project/coda/",
                      "rds_coda_ready", 
                      sample_name, "aneuploid_epi")

f_obja <- file.path(dir_coda, 'obja.rds')
# ident_str
# dir_res <- file.path(dirname(f_obja), 'cicero_demo'); fs::dir_create(dir_res)
dir_res <- file.path(dirname(f_obja), 'cicero'); fs::dir_create(dir_res)
obja <- read_rds(f_obja)
print(obja)
print(obja[['peaks']])

# ~~~ Create Cicero object ~~~ --------------------------------------  
#------ create the Monocle3 CDS object ------
f_monocle3 <- file.path(dir_res, 'obja.monocle3.rds')
if (! file.exists(f_monocle3)) {
# if ( T ) {
  cat('creating monocle3 object...\n')
  cds <- as_cell_data_set(obja)  
  cds <- cluster_cells(cds)
  write_rds(cds, f_monocle3)
} else {
  cat('loading monocle3 object...\n')
  cds <- read_rds(f_monocle3)
}

# plot_cells(cds)
# DimPlot(obja)

#------ preprocessing the CDS object ------
#------ export to cicero object ------
f_cicero <- file.path(dir_res, 'obja.cicero.rds')
if (! file.exists(f_cicero) ) {
# if ( T ) {
  cat('creating cicero object...\n')
  cds <- detect_genes(cds)
  ccr <- make_cicero_cds(cds, reduced_coordinates = reducedDims(cds)$UMAP)
  write_rds(ccr, f_cicero)
} else {
  cat('loading cicero object...\n')
  ccr <- read_rds(f_cicero)
}
# ~~~ Create Cicero connection ~~~ --------------------------------------  
f_conns <- file.path(dir_res, 'cicero_conns.df.rds')
if (! file.exists(f_conns)) {
  cat('creating cicero connection ...\n')
  data("human.hg19.genome")
  if (F) {
    cat('doing demo...\n')
    sample_genome <- subset(human.hg19.genome, V1 == "chr2")
    # sample_genome$V2[1] <- 10000000
    conns <- run_cicero(ccr, sample_genome, sample_num = 2) 
  } else {
    sample_genome <- human.hg19.genome
  }; head(sample_genome); dim(sample_genome)
  conns <- run_cicero(ccr, genomic_coords = sample_genome, 
                      window = 5e5, sample_num = 100)
  conns <- conns[!is.na(conns$coaccess), ]
  write_rds(conns, f_conns)
} else {
  cat('loading cicero connection ...\n')
  conns <- read_rds(f_conns)
}
conns <- conns[!is.na(conns$coaccess), ]
head(conns); nrow(conns); range(conns$coaccess, na.rm = T)
print(table('is_na'=is.na(conns$coaccess)))
print(table('great_than_zero'=conns$coaccess > 0))
hist(conns$coaccess[conns$coaccess>0], breaks = 500)

# sample  num of connection with coaccess > 0
# 1 BCMDCIS69     391,914       
# 2 DCIS22T       346,996     x
# 3 DCIS28T       484,472     
# 4 DCIS35T       776,622  
# 5 DCIS41T       1,319,154
# 6 DCIS51T       828,746     x
# 7 DCIS66T_chip2 446,604  
# 8 DCIS67T       710,242     x
# 9 DCIS68T       730,800     x

# ~~~ Create Cicero CCAN (chromatin hubs) ~~~ --------------------------------------  
f_ccan <- file.path(
  dir_res, 
  sprintf('cicero_CCAN.%s.df.rds', ifelse(coaccess_cutoff_override==0, 'default', coaccess_cutoff_override)))
if(coaccess_cutoff_override == 0) {
  
  xx = sort(conns$coaccess[conns$coaccess>0], decreasing = T)
  yy = cumsum(xx)
  df_xx <- data.frame(seq_along(xx), xx)
  df_yy <- data.frame(seq_along(yy), yy)
  
  library(pathviewr)
  #------ based on the elbowpoint of the sorted coaccess ------
  elbow_point_xx <- pathviewr::find_curve_elbow(df_xx)
  #------ based on the elbowpoint of the sorted accumulative coaccess [inuse]------
  elbow_point_yy <- pathviewr::find_curve_elbow(df_yy)
  
  coaccess_cutoff_override_xx <- xx[elbow_point_xx]
  pdf(file.path(dir_res, 'plot.decide_coaccess_cutoff_by_coaccess.pdf'), width = 6, height = 6, useDingbats = F)
  plot(seq_along(xx), xx, 'l', xlab='connection index', ylab='coaccess', 
       main = sprintf('%s / %s connections', sum(xx >= coaccess_cutoff_override_xx), length(xx)))
  abline(v=elbow_point_xx, col='red', lty='dashed')
  points(elbow_point_xx, coaccess_cutoff_override_xx, col='red')
  text(elbow_point_xx, coaccess_cutoff_override_xx, 
       signif(coaccess_cutoff_override_xx, 2), col='blue')
  dev.off()
  
  coaccess_cutoff_override_yy <- xx[elbow_point_yy]  
  pdf(file.path(dir_res, 'plot.decide_coaccess_cutoff_by_cumu_coaccess.pdf'), width = 6, height = 6, useDingbats = F)
  plot(seq_along(yy), yy, 'l', 
       main = sprintf('%s / %s connections', sum(xx>=coaccess_cutoff_override_yy), length(xx)), 
       xlab='connection index', ylab='cumulative coaccess')
  abline(v=elbow_point_yy, col='red')
  points(elbow_point_yy, yy[elbow_point_yy], col='red')
  text(elbow_point_yy, yy[elbow_point_yy], 
       signif(coaccess_cutoff_override_yy, 2), col='blue')
  dev.off()
  
  coaccess_cutoff_override_use <- coaccess_cutoff_override_xx
  coaccess_cutoff_override_use <- coaccess_cutoff_override_yy
  
} else {
  coaccess_cutoff_override_use <- coaccess_cutoff_override
}
cat('coaccess_cutoff_override_use=', signif(coaccess_cutoff_override_use, 3), '...\n')

# BCMDCIS69	0.17	47326 / 391914 connections  
# DCIS22T	0.22	55,088 / 346996 connections  
# DCIS28T	0.24	78220 / 484472 connections  
# DCIS35T	0.10	57690 / 776622 connections  
# DCIS41T	0.13	154572 / 1319154 connections
# DCIS51T	0.12	72,208 / 828746 connections  
# DCIS66T_chip2	0.19	70828 / 446604 connections  
# DCIS67T	0.13	65,224 / 710242 connections  
# DCIS68T	0.11	66,459 / 730800 connections  
write_rds(coaccess_cutoff_override_use, 
          file.path(dir_res, 'coaccess_cutoff_override_use.value.rds'))

if (! file.exists(f_ccan) ) {
# if (T) {
  cat('creating CCAN...\n')
  CCAN_assigns <- generate_ccans(
    conns, coaccess_cutoff_override = coaccess_cutoff_override_use)
  
  CCAN_assigns$CCAN <- paste0('CH', CCAN_assigns$CCAN)
  write_rds(CCAN_assigns, f_ccan)
} else {
  cat('loading CCAN...\n')
  CCAN_assigns <- read_rds(f_ccan)
} 
head(CCAN_assigns)
nrow(CCAN_assigns)
cat(sprintf('There are %s chromatin hubs.', length(unique(CCAN_assigns$CCAN))), '\n')
print(head(sort(table('CH with top size'=CCAN_assigns$CCAN), decreasing = T)))
cat('average CH size=',round(mean(table(CCAN_assigns$CCAN))), '\n')

try(dev.off())
pdf(file.path(dir_res, 'plot.hist.chromatin_hub_size.pdf'), width = 4, height = 3, useDingbats = F)
hist(table(CCAN_assigns$CCAN), breaks = 50, main = '', xlab='#peaks of a chromatin hub')
abline(v=mean(table(CCAN_assigns$CCAN)), col='red')
axis(side=1, at=round(mean(table(CCAN_assigns$CCAN))), label=round(mean(table(CCAN_assigns$CCAN))))
dev.off()
timestamp()
# cat('DONE\n')
# ~~~ Create assay of CCAN x cells ~~~ --------------------------------------  

head(CCAN_assigns)

CCAN_assigns_list <- ruok::deframe_to_list(CCAN_assigns[, c('CCAN', 'Peak')])
str(CCAN_assigns_list)
f_ccan_assay <- file.path(
  dir_res, 
  sprintf('CCAN_assay.%s.seurat.rds', ifelse(coaccess_cutoff_override==0, 'default', coaccess_cutoff_override)))
new_assay_name <- 'ccan'
if (! file.exists(f_ccan_assay) ) {
# if (T) {
  cat('creating CCAN assay for seurat...')
  
  ccan_assay_counts <- run_meta_feature(
    obja, CCAN_assigns_list, assay = 'peaks', slot='counts')
  # ccan_assay_data <- run_meta_feature(
  #   obja, CCAN_assigns_list, assay = 'peaks', slot='data')
  CCAN_assigns_str <- sapply(CCAN_assigns_list, function(x) {
    paste(x, collapse = ',')
  })
  CCAN_assigns_size <- sapply(CCAN_assigns_list, length)
  CCAN_assigns_df <- cbind(
    enframe(CCAN_assigns_str, name = 'CH', value = 'peaks_set'), 
    `peaks_num` = CCAN_assigns_size)
  ccan_assay <- CreateAssayObject(counts = ccan_assay_counts)
  ccan_assay <- CreateSeuratObject(
    counts = ccan_assay, assay = 'ccan', 
    meta.data = obja@meta.data[colnames(ccan_assay), ])
  
  ccan_assay <- RunTFIDF(ccan_assay)
  ccan_assay <- FindTopFeatures(ccan_assay, min.cutoff = 'q0')
  ccan_assay <- RunSVD(ccan_assay)  
  ccan_assay <- FindVariableFeatures(ccan_assay)
  ccan_assay <- Seurat::ScaleData(ccan_assay)
  ccan_assay <- RunPCA(ccan_assay)  
  print(range(ccan_assay@assays[[new_assay_name]]@data))
  if (!identical(rownames(ccan_assay@assays[[new_assay_name]]@meta.features), 
                 CCAN_assigns_df$CH)) {
    idx <- match(rownames(ccan_assay@assays[[new_assay_name]]@meta.features), 
                 CCAN_assigns_df$CH)
    CCAN_assigns_df <- CCAN_assigns_df[idx, ]
  }
  ccan_assay@assays[[new_assay_name]]@meta.features <- cbind(
    ccan_assay@assays[[new_assay_name]]@meta.features, 
    CCAN_assigns_df
  )
  write_rds(ccan_assay, f_ccan_assay)
  identical(rownames(ccan_assay@assays[[new_assay_name]]@meta.features), ccan_assay@assays[[new_assay_name]]@meta.features$CH)
  
  
  pdf(file.path(dir_res, 'plot.DepthCor.poor_lsi.pdf'), width = 7, height = 5, useDingbats = F)  
  p <- DepthCor(ccan_assay) + geom_hline(yintercept = c(.75, -.75), lty='dashed', col='red')
  print(p)
  dev.off()
  bad_lsi <- p$data$Component[abs(p$data$counts) > 0.75]
  dim_lsi <- head(setdiff(seq_len(ncol(Seurat::Reductions(ccan_assay, slot='lsi'))), 
                          bad_lsi), 30)
  # dim_lsi <- 1:30
  if (ncol(ccan_assay) < 100) {
    n.neighbors = 5
    k.param <- 5
  } else {
    n.neighbors = 20
    k.param <- 20
  }
  ccan_assay <- RunUMAP(
    ccan_assay, reduction = 'lsi', dims = dim_lsi, verbose = T,
    metric = 'correlation', n.neighbors=n.neighbors)
  umaplsi <- Reductions(ccan_assay, 'umap')
  umaplsi@key <- 'UMAPLSI_'
  colnames(umaplsi@cell.embeddings) <- str_replace_all(colnames(umaplsi@cell.embeddings), 'UMAP_', 'UMAPLSI_')
  ccan_assay[['umaplsi']] <- umaplsi
  
  pdf(file.path(dir_res, 'plot.ElbowPlot.poor_pca.pdf'), width = 7, height = 5, useDingbats = F)  
  print(ElbowPlot(ccan_assay))
  dev.off()
  ccan_assay <- RunUMAP(
    ccan_assay, reduction = 'pca', dims = 1:30, verbose = T,
    metric = 'correlation', n.neighbors=n.neighbors)
  umappca <- Reductions(ccan_assay, 'umap')
  umappca@key <- 'UMAPPCA_'
  colnames(umappca@cell.embeddings) <- str_replace_all(colnames(umappca@cell.embeddings), 'UMAP_', 'UMAPPCA_')
  ccan_assay[['umappca']] <- umappca
  
  write_rds(ccan_assay, f_ccan_assay)
  
  pdf(file.path(dir_res, 'plot.UMAP_clones.%01d.pdf'), width = 6, height = 5, useDingbats = F, onefile = F)
  try(print(DimPlot(ccan_assay, group.by = 'clones', reduction ='umaplsi') + theme(aspect.ratio = 1)))
  try(print(DimPlot(ccan_assay, group.by = 'clones', reduction ='umappca') + theme(aspect.ratio = 1)))
  # try(print(DimPlot(obja, group.by = 'clones')))
  dev.off()  
} else {
  cat('loading CCAN assay for seurat...')
  ccan_assay <- read_rds(f_ccan_assay)
}

cat('DONE', sample_name, '\n')