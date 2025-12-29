#---------------------------
# Develop functions for:
# - findMarkerBins 
# - findAllMarkerBins
# - findCNAEvents
#---------------------------
library(SummarizedExperiment)
library(SingleCellExperiment)
#-------------------------- misc --------------------------  
# Signac::GRangesToString
grange_to_string <- function(grange, sep = c("-", "-")) {
  regions <- paste0(
    as.character(seqnames(grange)), 
    sep[[1]], 
    start(grange), 
    sep[[2]], 
    end(grange))
  return(regions)
}

calc_fano_factor <- function(mat) {
  o <- matrixStats::rowVars(mat) / matrixStats::rowMeans2(mat)
  return(o)
}
calcFanoFactor <- function(x, bins = NULL, assay = 'bincounts') {
  if (is.null(bins)) { bins <- rownames(x) }
  mat <- assay(x, assay)
  mat <- as.matrix(mat[bins, ])
  
  ff <- calc_fano_factor(mat=mat)
  rowData(x)$fano_factor <- ff
  
  return(x)
}
# Signac::StringToGRanges
# string_to_grange <- function (regions, sep = c("-", "-"), ...) {
#   ranges.df <- data.frame(ranges = regions)
#   ranges.df <- tidyr::separate(
#     data = ranges.df, 
#     col = "ranges", 
#     sep = paste0(sep[[1]], "|", sep[[2]]), into = c("seqnames", "start", "end"))
#   granges <- makeGRangesFromDataFrame(df = ranges.df, ...)
#   return(granges)
# }

# calculate COV for bins
library(gtools)
findAllMarkerBins <- function(
    x, assay = NULL, group.by = 'subclones', 
    bins = NULL, 
    test_use = 'wilcox', 
    min_cov = 0.1, 
    latent_vars = NULL, 
    max_cells_per_ident = Inf, 
    random_seed = 42) {
  
  # require COV  todo
  
  # 
  if ('factor' %in% class(x[[group.by]])) {
    idents_opts <- levels(x[[group.by]])
  } else {
    idents_opts <- gtools::mixedsort(unique(x[[group.by]]))
  }
  
  # run findMarkerBins
  res <- list()
  for (xx in idents_opts) {
    message('Identifying DABs for ', group.by, '=', xx, '...')
    cells.1 <- colnames(x)[x[[group.by]] == xx]
    cells.2 <- setdiff(colnames(x), cells.1)
    df_xx <- try(findMarkerBins(
      x=x, assay = assay, bins = bins, 
      cells.1 = cells.1, cells.2 = cells.2, test_use = test_use, 
      latent_vars = latent_vars))
    
    if ('try-error' %in% class(df_xx)) {
      df_xx <- data.frame()
    } else {
      df_xx$cluster <- xx
    }
    res <- c(res, list(df_xx))
  }
  # merge results
  res <- do.call(rbind, res)
  res$cluster <- factor(as.character(res$cluster), levels = idents_opts)
  return(res)

}

library(Matrix.utils)
aggregate.Matrix.t <- function(x, groupings = NULL, fun = 'median') {
  # x: features x samples
  # o: features x groups
  x <- as.data.frame(t(x))
  if (fun %in% c('sum', 'mean', 'count')) {
    o <- t(Matrix.utils::aggregate.Matrix(
      x=x, groupings = groupings, fun = fun))
  } 
  if (fun %in% c('median')) {
    o <- lapply(split(x, groupings), function(data) {
      return(matrixStats::colMedians(as.matrix(data)))
    })
    o <- do.call(cbind, o)
  }
  return(o)
}

# x: copykit object
calc_consensus_assays <- function(
  x, assay = NULL, 
  bins = NULL, 
  cells.1 = NULL, 
  cells.2 = NULL) {
  
  if (is.null(bins)) {bins <- rownames(x)}
  
  assay_use <- intersect(c('bincounts', 'integer', 'logr', 'segment_ratios'), 
                         assayNames(x))
  
  groupings <- factor( rep(c('Group1', 'Group2'), 
                           c(length(cells.1), length(cells.2))) )
  bidx <- match(bins, rownames(x))
  cidx <- match(c(cells.1, cells.2), colnames(x))
  o <- lapply(assay_use, FUN = function(a) {
    # data: bins x cells
    data <- SummarizedExperiment::assay(x, a)
    data <- data[bidx, cidx, drop=F]
    df <- as.data.frame(aggregate.Matrix.t(data, groupings, fun = 'median'))
    colnames(df) <- paste0(a, c('.1', '.2'))
    
    return(df)
  })
  
  o <- do.call(cbind, o)
  
  avg_log2FC <- switch (assay,
    'bincounts' = log2(o$bincounts.1 + 1) - log2(o$bincounts.2 + 1), 
    'integer' = log2(o$integer.1 + 0.001) - log2(o$integer.2 + 0.001),
    'logr' = o$logr.1 - o$logr.2,
    'segment_ratios' = log2(o$segment_ratios.1 + 1e-3) - log2(o$segment_ratios.2 + 1e-3)
  )
  o$avg_log2FC <- avg_log2FC
  return(o)
}

library(pbapply)
# x: bins x cells
dab_test_Wilcox <- function(
    x, 
    cells.1, 
    cells.2,
    verbose = TRUE, ...
) {
  x <- x[, c(cells.1, cells.2), drop=F]
  groupings <- data.frame(row.names = c(cells.1, cells.2))
  groupings[cells.1, 'group'] <- 'Group1' # focus 
  groupings[cells.2, 'group'] <- 'Group2' # other
  # todo: BiocParallel::bplapply
  sapply_use <- pbapply::pbsapply
  pvals <- sapply_use(
    1:nrow(x), 
    function(b) { 
      wilcox.test(as.numeric(x[b, ]) ~ groupings[, 'group'])$p.value
    }
  )
  o <- data.frame(p_val = pvals, row.names = rownames(x))
  return(o)
}

# latent_vars: data frame
# adapted from https://github.com/satijalab/seurat/blob/763259d05991d40721dee99c9919ec6d4491d15e/R/differential_expression.R#L1696
dab_test_LR <- function(
    x, 
    cells.1, 
    cells.2,
    latent_vars, 
    verbose = TRUE, ...
) {
  group.info <- data.frame(row.names = c(cells.1, cells.2))
  group.info[cells.1, "group"] <- "Group1" # focus
  group.info[cells.2, "group"] <- "Group2" # other
  group.info[, "group"] <- factor(group.info[, "group"])
  
  x <- x[, rownames(group.info), drop = FALSE]
  
  latent_vars <- latent_vars[rownames(group.info), , drop = FALSE]
  
  sapply_use <- pbsapply
  p_val <- sapply_use(
    1:nrow(x),
    function(b) {
      if (is.null(latent_vars)) {
        model.data <- cbind(BIN = as.numeric(x[b, ]), group.info)
        fmla <- as.formula("group ~ BIN")
        fmla2 <- as.formula("group ~ 1")
      } else {
        model.data <- cbind(BIN = as.numeric(x[b, ]), group.info, latent_vars)
        fmla <- as.formula(object = paste(
          "group ~ BIN +",
          paste(colnames(latent_vars), collapse = "+")
        ))
        fmla2 <- as.formula(object = paste(
          "group ~",
          paste(colnames(latent_vars), collapse = "+")
        ))
      }
      
      model1 <- glm(formula = fmla, data = model.data, family = "binomial")
      model2 <- glm(formula = fmla2, data = model.data, family = "binomial")
      lrtest <- lmtest::lrtest(model1, model2)
      return(lrtest$Pr[2])
    }
  )
  o <- data.frame(p_val = p_val, row.names = rownames(x))
  return(o)
}

# x: matrix (bins x cells)
run_dab <- function(
    x,
    cells.1 = NULL, 
    cells.2 = NULL,
    test_use = 'wilcox', 
    latent_vars = NULL, 
    verbose = TRUE,
    ...) {

  o <- switch(
    test_use,
    'wilcox' = dab_test_Wilcox(
      x, cells.1 = cells.1, cells.2 = cells.2, 
      verbose = verbose, ...), 
    'LR' = dab_test_LR(
      x, cells.1 = cells.1, cells.2 = cells.2, 
      latent_vars = latent_vars,
      verbose = verbose, ...), 
    data.frame()
  )
  o <- cbind(BIN=rownames(o), o)
  return(o)
}

# x: object
findMarkerBins <- function(
    x, assay = NULL,
    bins = NULL, 
    cells.1 = NULL, 
    cells.2 = NULL,
    test_use = 'wilcox', 
    min_cov = 0.1, 
    latent_vars = NULL, 
    max_cells_per_ident = Inf, 
    random_seed = 42, 
    verbose = T, 
    ...) {
  
  # select bins based on COV
  if (is.null(bins)) { bins <- rownames(x) }
  
  # calculate FC and averages
  df_fc_avg <- calc_consensus_assays(
    x=x, assay = assay, 
    bins = bins, 
    cells.1 = cells.1, 
    cells.2 = cells.2)
  # run diff to detect DAB 
  mat <- assay(x, assay)
  mat <- mat[bins, , drop=F]
  latent_vars <- as.data.frame(colData(x)[, latent_vars, drop=F])
  df_dab <- run_dab(
    x = mat,
    cells.1 = cells.1, 
    cells.2 = cells.2,
    test_use = test_use, 
    latent_vars = latent_vars, 
    verbose = verbose, 
    ...)
  print(nrow(x))
  df_dab$p_val_adj <- p.adjust(df_dab$p_val, method = "BH", n = nrow(x))
  res <- cbind(df_fc_avg, df_dab)
  # fetch bin labels
  if (!'bin_labels' %in% colnames(rowData(x))) {
    rowData(x)$bin_labels <- grange_to_string(rowRanges(x))
  }
  bidx <- match(res$BIN, rownames(x))
  res <- cbind(BIN_label=rowData(x)$bin_labels[bidx], res)
  return(res)
}


#-------------------------- demo --------------------------  
if (F) {
  
  tmp = calc_consensus_assays(objd, assay = 'integer', 
                              cells.1 = colnames(objd)[objd[['subclones']]=='c1'], 
                              cells.2 = colnames(objd)[objd[['subclones']]!='c1'])
  tmp=aggregate.Matrix.t(assay(objd, 'bincounts'), objd[['subclones']], 
                         fun = 'median')
  head(tmp)
  objd <- calcConsensus(objd, assay = 'bincounts')
  objd <- calcConsensus(objd, assay = 'integer')
  copykit::consensus(objd)[1:6, ]
  # c1       c2       c3       c4
  # 1 50.24653 57.73852 38.14642 54.40918
  # 2 55.52665 59.66443 46.51447 62.28734
  # 3 62.11354 70.04213 56.21362 68.30558
  
  # c1 c2 c3 c4
  # 1  2  2  2  2
  # 2  2  2  2  2
  # 3  2  2  2  2
  # 4  2  2  2  2
  # 5  2  2  2  2
  # 6  2  2  2  2
  tmp = findMarkerBins(objd, assay = 'logr', 
                       cells.1 = colnames(objd)[objd[['subclones']]=='c1'], 
                       cells.2 = colnames(objd)[objd[['subclones']]!='c1'])
  tmp$qval <- p.adjust(tmp$p_val, method = 'BH', n=nrow(objd))
  all.equal(tmp$qval, tmp$p_val_adj)
  
  xx = matrix(1:12, nrow=4)
  aggregate(xx, c('A', 'B', 'A', 'B'), fun = 'median')
}