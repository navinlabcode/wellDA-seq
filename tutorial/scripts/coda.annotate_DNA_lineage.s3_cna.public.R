#---------------------------
# Annotate DNA lineage
# 
# - genotype and phenotype concordance based on the original 
# distance matrix, instead of trees [fixed]
# - plasticity score by Yun [fixed]
# - build a genetic tree
#   - medic2Orig: original MEDIC2 tree => delete diploid => reroot by the clone most similar to diploid => ladderize
#   - medic2ME: MEDIC2 distance => ME tree => delete diploid => ladderize
# - visualization of clone order depends on the tree
#---------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(tidyverse); library(forcats)
  library(stringr)
  library(fs)
  library(copykit)
  library(VennDiagram)
  library(SummarizedExperiment)
  library(ggpubr); library(ruok); library(scales)
  library(ComplexHeatmap)
  library(UCell)
  library(ggtree); library(phylogram); library(dendextend); library(phytools)
  library(ape)
  library(dendsort)
  library(rstatix)
  library(patchwork)
  theme_set(theme_pubr(legend = 'right', base_size = 6) %+replace% theme(axis.ticks.length=unit(0.1,"inch")))
  get_tree_tips_order <- function(phylo) {
    # Ref: https://github.com/Puriney/copykit/blob/faf47072f34130dda263037fda637520899ca013/R/plotHeatmap.R#L359-L362
    is_tip <- phylo$edge[, 2] <= length(phylo$tip.label)
    ordered_tips_index <- phylo$edge[is_tip, 2]
    phylo_tips_order <- phylo$tip.label[ordered_tips_index] %>% rev()
    return(phylo_tips_order)
  }
  
  calc_cophenetic_cor <- function(phylo, hc) { 
    a <- cophenetic(phylo) ## a matrix
    b <- cophenetic(hc) # a dist
    a <- a[labels(b), labels(b)]
    a <- as.dist(a)
    cor(a, b)
  }
  
  get_hclust_tip_order <- function(hc){ hc$labels[hc$order] }
  get_phylo_tip_order <- function(tree) {
    # ref: https://stackoverflow.com/a/34364914/1608734
    is_tip <- tree$edge[,2] <= length(tree$tip.label)
    ordered_tips <- tree$edge[is_tip, 2]
    return(tree$tip.label[ordered_tips])
  }
})
setwd("/volumes/USR1/yyan/project/coda")
source('/volumes/USR1/yyan/project/coda/rsrc/utils.R')
source('/volumes/USR1/yyan/project/coda/rsrc/color_pal.R')

options <- commandArgs(trailingOnly = TRUE)
if (length(options)>0) {
  cat('Reading user parameters:\n')
  sample_name <- options[[1]]
  rna_assay_use <- options[[2]]
  module_method <- options[[3]]
  tree_method <- options[[4]]
  cna_dist_assay <- options[[5]]
} else {
  sample_name = 'DCIS66T_chip2'
  rna_assay_use <- 'ciceroRNA'
  module_method <- 'ucell'
  tree_method <- 'medic2Orig'
  cna_dist_assay <- 'integerCNbased' # cna_dist_use <- 'eventCNbased'
}
message(sample_name)
ident_str <- 'clones'

#------------------ ~~~ Intro ~~~ --------------------

#------ input / output ------
  
dir_coda <- file.path(
  '/volumes/USR1/yyan/project/coda', 
  'rds_coda_ready', sample_name, 'aneuploid_epi')


load_coda(dir_coda) # obja; objd; df_meta

#------ load ciceroRNA ------
f_crna <- file.path(dir_coda, 'cicero', 'ciceroRNA_assay.seurat.rds')
crna <- try(read_rds(f_crna))

#------ choose RNA assay ------

if ( 'try-error' %in% class(crna)) {
  rna_assay_use = 'iRNA'
  DefaultAssay(obja) <- 'iRNA'
} else {
  stopifnot( identical(Cells(crna), Cells(obja)) )
  obja[['ciceroRNA']] <- crna[['ciceroRNA']]
  DefaultAssay(obja) <- 'ciceroRNA'
}

DefaultAssay(obja) <- rna_assay_use
cat('ATAC uses RNA assay = ', DefaultAssay(obja), '\n')
print(obja)

#------ load Modules scores ------
name_db <- 'msigdb_hallmark'

dir_modulescore <- file.path(
  dir_coda, 
  sprintf('%s_%s_%s', module_method, name_db, DefaultAssay(obja)))

message('Load AddModuleScore ', module_method)
df_signature_res <- read_rds(file.path(dir_modulescore, 'addmodulescore.df.rds'))
stopifnot(identical(rownames(df_signature_res), Cells(obja)))
obja <- AddMetaData(obja, metadata = df_signature_res)

update_coda(df_meta, obja, objd)

#------ Decide out dir ------
# medic2's original tree
dir_res <- file.path(
  dir_coda, 
  sprintf('%s_%s_%s_%s_%s', 
          cna_dist_assay,
          tree_method, 
          module_method,
          name_db, 
          DefaultAssay(obja)))
fs::dir_create(dir_res)
message(basename(dir_res))
#------------------ ~~~ Concordance ~~~ --------------------
ident_vec <- levels(forcats::fct_drop(df_meta[[ident_str]]))
my_sub <- factor(ident_vec, levels = ident_vec)
n_clones <- length(ident_vec)
#------ medic2 tree ------
dir_medic2 <- file.path(dir_coda, 'medic2')
medic_dist <- read_tsv(
  file.path(dir_medic2, 'medicc2_input_pairwise_distances.tsv'))
medic_dist <- medic_dist %>% as.data.frame() %>% 
  column_to_rownames(var="sample_id") %>% as.dist()
most_diploid_nb <- ident_vec[which.min(as.matrix(medic_dist)['diploid', ident_vec])]
print(most_diploid_nb)
#------ CNA css distance ------
dir_css_cna <- file.path(dir_coda, 'css_distance')
cna_dist_mat <- read_rds(file.path(dir_css_cna, 'css_integerCN_dist_manhattan.mat.rds'))
me_tree <- read_rds(file.path(dir_css_cna, 'fastme_tree.rds'))
me_tree <- phytools::reroot(
  me_tree %>% ladderize(), 
  node.number = which(me_tree$tip.label == most_diploid_nb))
me_tree_nodipoid <- ape::drop.tip(me_tree, tip = 'diploid')
cna_dist_val_range <- range(cna_dist_mat)
#------ atac distance ------
if (T) {
  ## based on LSI components
  DefaultAssay(obja) <- 'peaks'
  n_lsi <- ncol(Embeddings(obja, reduction = 'lsi'))
  p <- DepthCor(obja, n = n_lsi)
  lsi_idx <- head(setdiff(1:n_lsi, which(abs(p$data$counts)>0.75)), 30)
  atac_events <- Embeddings(obja, reduction = 'lsi')[, lsi_idx]
  DefaultAssay(obja) <- rna_assay_use
  dict_cell2ident <- obja@meta.data[[ident_str]]; names(dict_cell2ident) <- Cells(obja)
  dict_cell2ident <- forcats::fct_drop(dict_cell2ident)
  
  ## median is better than mean
  atac_events <- apply(atac_events, 2, function(vv) {
    tapply(vv, dict_cell2ident, median)
  })
  # atac_events <- apply(atac_events, 2, function(vv) {
  #   tapply(vv, dict_cell2ident, mean)
  # })
  atac_dist <- dist(atac_events)
  atac_dist <- as.matrix(atac_dist)
  atac_dist <- atac_dist[ident_vec, ident_vec]
}
print(atac_dist); print(range(atac_dist))
write_csv(as.data.frame(atac_dist), 
          file.path(dir_res, 'atac_pairwise_distance.csv'))
atac_dist <- scales::rescale(
  x = as.matrix(as.dist(atac_dist)), 
  to = cna_dist_val_range) ## good idea to keep the scale similar to get a better visualization
# atac_dist <- as.matrix(as.dist(atac_dist)) ## rescale and non-rescale have the same coph score
all(rownames(atac_dist) %in% ident_vec)
atac_dist <- atac_dist[intersect(ident_vec, rownames(atac_dist)),intersect(ident_vec, rownames(atac_dist))]

#------ concordance ------
stopifnot(all.equal(rownames(cna_dist_mat), rownames(atac_dist)))
stopifnot(all.equal(colnames(cna_dist_mat), colnames(atac_dist)))
score_GP_concordance_obj <- cor.test(
  as.dist(cna_dist_mat), 
  as.dist(atac_dist), method = 'pearson')
score_GP_concordance <- as.numeric(score_GP_concordance_obj$estimate)
cor.test(as.dist(cna_dist_mat), 
         as.dist(atac_dist), 
         method = 'pearson')
cat('score_GP_concordance = ', signif(score_GP_concordance, digits = 2), 'in sample', sample_name, '.\n')
write_lines(score_GP_concordance,
            file = file.path(dir_res, 'score_GP_concordance.txt'))
library(ade4)
run_mantel_rtest <- function(m1, m2){
  lbs <- rownames(m1)
  stopifnot(identical(intersect(lbs, rownames(m1)), lbs))
  stopifnot(identical(intersect(lbs, rownames(m2)), lbs))
  m1 <- m1[lbs, lbs]
  m2 <- m2[lbs, lbs]
  set.seed(42)
  mantel_res <- ade4::mantel.rtest(as.dist(m1), as.dist(m2))
  # print(mantel_res)
  return(c('pvalue'=mantel_res$pvalue, 'r'=mantel_res$obs))
}
score_GP_concordance_pval <- as.numeric(
  score_GP_concordance_obj$p.value)
#------ viz: concordance ------
if (tree_method == 'medic2ME') {
  genotree_tip_orders <- get_phylo_tip_order(me_tree_nodipoid)
} else if (tree_method == 'medic2Orig') {
  # always use the nondiploid tree  
  genotree_tip_orders <- get_phylo_tip_order(me_tree_nodipoid)
  # genotree_tip_orders <- get_phylo_tip_order(me_tree)
  genotree_tip_orders <- setdiff(genotree_tip_orders, 'diploid')
} else {
  warning('unsupported ', tree_method)
}
print(genotree_tip_orders)
write_lines(genotree_tip_orders, 
            file.path(dir_res, 'tree_tip_order.txt'))

if (T) {
  tmp_cna <- cna_dist_mat[rev(genotree_tip_orders), rev(genotree_tip_orders)]
  tmp_type <- 'CNA'
  pt_cna <- Heatmap(
    matrix = tmp_cna, 
    name = sprintf('%s dist', tmp_type), 
    cluster_rows = F, cluster_columns = F, 
    col = circlize::colorRamp2(
      quantile(c(c(cna_dist_mat), c(atac_dist)), c(0.01, 0.99)), 
      hcl.colors(n=2, 'Blues 3', rev = F)),
    show_row_dend = T, show_column_dend = T, 
    width = unit(n_clones*0.5, 'inch'),
    height = unit(n_clones*0.5 /ncol(tmp_cna)*nrow(tmp_cna), 'inch'),
    row_names_side  = 'left',
    rect_gp = gpar(type = "none"), 
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (i >= j) { 
        grid.rect(x, y, width, height, gp = gpar(fill = fill, col = fill))
        if (tmp_cna[i, j] > 0) {
          grid.text(sprintf("%d", tmp_cna[i, j]), x, y, gp = gpar(fontsize = 10, col='black'))
        }
      }
    },
    
    column_title = tmp_type, 
    column_title_gp = gpar(fontsize = 10), use_raster = F
  )
  
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  
  tmp_type <- 'ATAC'
  tmp_atac <- atac_dist[rev(genotree_tip_orders), rev(genotree_tip_orders)]
  pt_atac <- Heatmap(
    matrix = tmp_atac, 
    name = sprintf('%s dist', tmp_type), 
    cluster_rows = F, cluster_columns = F, 
    col = circlize::colorRamp2(
      quantile(c(c(cna_dist_mat), c(atac_dist)), c(0.01, 0.99)), 
      hcl.colors(n=2, 'Blues 3', rev = F)),
    show_row_dend = T, show_column_dend = T, 
    width = unit(n_clones*0.5, 'inch'),
    height = unit(n_clones*0.5 /ncol(tmp_atac)*nrow(tmp_atac), 'inch'),
    row_names_side  = 'left',
    rect_gp = gpar(type = "none"), 
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (i >= j) {
        grid.rect(x, y, width, height, gp = gpar(fill = fill, col = fill))
        if (tmp_atac[i, j] > 0) {
          if (is.wholenumber(tmp_atac[i, j])) {
            grid.text(sprintf("%d", as.integer(tmp_atac[i, j])), x, y, gp = gpar(fontsize = 10, col='black'))
          } else {
            grid.text(sprintf("%.2f", tmp_atac[i, j]), x, y, gp = gpar(fontsize = 10, col='black'))
          }
        }
      }
    },
    column_title = tmp_type, 
    column_title_gp = gpar(fontsize = 10), use_raster = F
  )
  # pt_cna
  # pt_atac
  tmp_diff <- abs(tmp_cna - tmp_atac)
  tmp_type <- 'abs diff'

  pt_diff <- Heatmap(
    matrix = tmp_diff, 
    name = sprintf('%s', tmp_type), 
    cluster_rows = F, cluster_columns = F, 
    col = circlize::colorRamp2(
      c(cna_dist_val_range[1], mean(cna_dist_val_range), cna_dist_val_range[2]), 
      hcl.colors(n=3, 'Cyan-Magenta', rev = F)),
    show_row_dend = T, show_column_dend = T, 
    width = unit(n_clones*0.5, 'inch'),
    height = unit(n_clones*0.5 /ncol(tmp_diff)*nrow(tmp_diff), 'inch'),
    row_names_side  = 'left',
    rect_gp = gpar(type = "none"), 
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (i >= j) {
        # only draw the low-triangle part
        grid.rect(x, y, width, height, gp = gpar(fill = fill, col = fill))
        if (tmp_diff[i, j] > 0) {
          if (is.wholenumber(tmp_diff[i, j])) {
            grid.text(sprintf("%d", as.integer(tmp_diff[i, j])), x, y, gp = gpar(fontsize = 10, col='black'))
          } else {
            grid.text(sprintf("%.2f", tmp_diff[i, j]), x, y, gp = gpar(fontsize = 10, col='black'))
          }
        }
      }      
    },
    column_title = tmp_type, 
    column_title_gp = gpar(fontsize = 10), use_raster = F
  )
  
  pdf(file.path(dir_res, 'heatmap.GP_concordance.pdf'), 
      width = n_clones * 0.8 *  3 + 1, 
      height = n_clones *0.8 + 2, useDingbats = F)
  draw(
    pt_cna + pt_atac + pt_diff, 
    column_title = sprintf('concordance = %.2f | pval=%.3g', 
                           score_GP_concordance, 
                           score_GP_concordance_pval), 
    column_title_gp = gpar(fontsize = 10))
  dev.off()
  
  pdf(file.path(dir_res, 'heatmap.GP_concordance.cor.pdf'), 
      width = n_clones * 0.5  + 1, 
      height = n_clones *0.5 + 1, useDingbats = F)
  (data.frame(CNA_dist=tmp_cna[lower.tri(tmp_cna)], 
              ATAC_dist=tmp_atac[lower.tri(tmp_atac)]) %>%
      ggplot(aes(x=CNA_dist, y=ATAC_dist)) + 
      geom_abline(slope = 1, intercept = 0, lty='dashed', color='grey') + 
      geom_point(pch=1, size=3) + 
      geom_smooth(method = 'lm', se=F) + 
      stat_cor() +
      coord_cartesian(xlim = cna_dist_val_range, ylim=cna_dist_val_range)
  ) %>% print()
  dev.off()
  
  # rm(tmp_atac); rm(tmp_cna); rm(tmp_type)
}


#------------------ ~~~ Plasticity ~~~ --------------------

#------ propose features ------
z = ident_str
phenotypes_opts <- colnames(df_signature_res)

phenotypes_sig_diverse_opts <- c()
pheno_ks_test_pval_vec <- c()
for (i in seq_along(phenotypes_opts)) {
  stat.test <- df_meta %>% kruskal_test( as.formula(sprintf('%s~%s', phenotypes_opts[i], z)))
  pheno_ks_test_pval_vec <- c(pheno_ks_test_pval_vec, stat.test$p)
}; names(pheno_ks_test_pval_vec) <- phenotypes_opts
pheno_ks_test_qval_vec <- p.adjust(pheno_ks_test_pval_vec, method = 'BH')
## pval < 0.05
phenotypes_sig_diverse_opts <- phenotypes_opts[which(pheno_ks_test_pval_vec < 0.05)]
if (length(phenotypes_sig_diverse_opts) <= 2) {
  cat('[warn] there are only ', length(phenotypes_sig_diverse_opts), 
      ' sig diff features in sample', sample_name, '.\n' )
}

#------ phenotype ordering ------
str(phenotypes_sig_diverse_opts)
df_pheno_viz <- df_meta %>% dplyr::select(all_of(c(phenotypes_sig_diverse_opts, z)))
df_pheno_viz <- df_pheno_viz %>% 
  dplyr::group_by_at(z) %>%
  dplyr::summarise_all(mean, na.rm=T)
df_pheno_viz$clones <- droplevels(df_pheno_viz$clones)
clones_opts <- levels(df_pheno_viz$clones)
mat_pheno_viz <- column_to_rownames(df_pheno_viz, var = z)
mat_pheno_viz <- as.matrix(mat_pheno_viz)
mat_pheno_viz_zscore <- scale(mat_pheno_viz)
df_pheno_viz_zscore <- as.data.frame(mat_pheno_viz_zscore) %>%
  rownames_to_column(var = z)
n_phenos <- ncol(mat_pheno_viz)

if (length(phenotypes_sig_diverse_opts)>=2) {
  pheno_hclust <- dendsort::dendsort(
    hclust(dist(t(mat_pheno_viz_zscore)), method = 'ward.D2'))
  pheno_viz_orders <- pheno_hclust$labels[pheno_hclust$order]
} else if (length(phenotypes_sig_diverse_opts) == 1) {
  pheno_viz_orders <- phenotypes_sig_diverse_opts
} else {
  warning('no sig diff features.')
  pheno_viz_orders <- phenotypes_sig_diverse_opts
}
pheno_viz_orders_pretty <- pheno_viz_orders %>% 
  str_replace_all('_', ' ') %>%
  str_remove_all('HALLMARK') %>% str_remove_all('UCell|seurat|ucell') %>% 
  str_squish() %>%
  tolower()
phenotree_tip_orders <- pheno_viz_orders

#------ plasticity scoring by Geary's C ------
calc_gearyc <- function(w, x) {
  x_mu <- mean(x, na.rm = T)
  mat_dist_x <- as.matrix(dist(x, method="euclidean"))
  shared_obs <- intersect(rownames(w), names(x))
  w <- w[shared_obs, shared_obs, drop=F]
  mat_dist_x <- mat_dist_x[shared_obs, shared_obs, drop=F]
  N <- length(x)
  W <- sum(w)
  C <- (N-1) * sum(mat_dist_x^2 * w) / (2 * W * sum((x - x_mu)^2))
  return(C)
}


calc_heritability_score_gearyc <- function(d, x) {
  ## d: distance between clones 
  ## x: values for the clones
  max_d <- apply(d, 1, max)
  w <- exp(-1 * d^2 / max_d^2)
  C <- calc_gearyc(w, x)
  return(1 - C)
}

library(pbmcapply)
calc_heritability_score_gearyc_hnull <- function(d, x, ntimes=1e3, ncores=8) {
  o <- pbmclapply(1:ntimes, function(r) {
    set.seed(r)
    d_names <- sample(x = rownames(d), size = nrow(d), replace = F)
    ## remake the dummy matrices
    d_new <- d
    rownames(d_new) <- colnames(d_new) <- d_names
    calc_heritability_score_gearyc(d_new, x)
  }, mc.cores = ncores, mc.set.seed = T)
  o <- sort(as.numeric(unlist(o)))
}

phenotypes_sig_diverse_heritability_gearyc <- c()
phenotypes_sig_diverse_heritability_gearyc_pval_up <- c()
phenotypes_sig_diverse_heritability_gearyc_pval_dn <- c()
pdf(file.path(dir_res, 'heatmap.heritability_gearyc_test.birdview.pdf'), 
    width = 5,
    height = 5, useDingbats = F, 
    onefile = T)
for (xx in phenotypes_sig_diverse_opts) {
  message(xx)
  xx_gearyc <- calc_heritability_score_gearyc(
    d = cna_dist_mat, 
    x=mat_pheno_viz[, xx, drop=T]
  )
  xx_null <- calc_heritability_score_gearyc_hnull(
    d = cna_dist_mat, 
    x = mat_pheno_viz[, xx, drop=T], 
    ntimes = 1e3, ncores = 8
  )
  xx_pval_up <- (sum(xx_null>xx_gearyc)+1)/(length(xx_null)+1)
  xx_pval_dn <- (sum(xx_null<xx_gearyc)+1)/(length(xx_null)+1)
  hist(xx_null, main = xx,
       xlab = sprintf('prob up = %s; prob dn = %s', 
                      signif(xx_pval_up, 3), 
                      signif(xx_pval_dn, 3)))
  abline(v = quantile(xx_null, c(0.05, 1-0.05)), col='red')
  abline(v = xx_gearyc, col='blue')
  
  
  phenotypes_sig_diverse_heritability_gearyc <- c(
    phenotypes_sig_diverse_heritability_gearyc, 
    xx_gearyc)
  phenotypes_sig_diverse_heritability_gearyc_pval_up <- c(
    phenotypes_sig_diverse_heritability_gearyc_pval_up, 
    xx_pval_up
  )
  phenotypes_sig_diverse_heritability_gearyc_pval_dn <- c(
    phenotypes_sig_diverse_heritability_gearyc_pval_dn, 
    xx_pval_dn
  )
  rm(xx_gearyc); rm(xx_pval_up); rm(xx_pval_dn)
}
names(phenotypes_sig_diverse_heritability_gearyc) <- phenotypes_sig_diverse_opts
names(phenotypes_sig_diverse_heritability_gearyc_pval_up) <- phenotypes_sig_diverse_opts
names(phenotypes_sig_diverse_heritability_gearyc_pval_dn) <- phenotypes_sig_diverse_opts
rm(xx)
dev.off()

#------ plasticity scoring by Yun ------

calc_heritability_score_orig <- function(d1, d2) {
  # d1, d2: two distance objects
  lb <- labels(d1)
  m2 <- as.matrix(d2)
  m2 <- m2[lb, lb]
  d2 <- as.dist(m2)
  # o <- cor.test(d1, d2, method = 'spearman')
  o <- cor.test(d1, d2, method = 'pearson')  # similar as spearman
  c('r' = as.numeric(o$estimate), 
    'pval' = as.numeric(o$p.value))
}

phenotypes_sig_diverse_heritability_orig_test <- sapply(
  phenotypes_sig_diverse_opts, function(xx) {
    calc_heritability_score_orig(
      as.dist(cna_dist_mat),
      dist(mat_pheno_viz_zscore[, xx, drop=F]) )
  }
)
phenotypes_sig_diverse_heritability_orig <- phenotypes_sig_diverse_heritability_orig_test['r', ]
phenotypes_sig_diverse_heritability_orig_pval <- phenotypes_sig_diverse_heritability_orig_test['pval', ]
names(phenotypes_sig_diverse_heritability_orig) <- names(phenotypes_sig_diverse_heritability_orig_pval) <- phenotypes_sig_diverse_opts

p_cna <- pheatmap(
  cna_dist_mat[rev(genotree_tip_orders), rev(genotree_tip_orders)], 
  color = circlize::colorRamp2(
    quantile(c(cna_dist_mat), c(0.01, 0.99)), 
    hcl.colors(n=2, 'Blues 2', rev = F)), 
  name = 'CNA dist', 
  display_numbers = T, number_color = 'black', 
  cluster_rows = F, cluster_cols = F, 
  border_color = 'black', main = 'CNA'
)
pdf(file.path(dir_res, 'heatmap.heritability_each_feature_vs_medic2.birdview.pdf'), 
    width = n_clones * 2 + 1, 
    height = n_clones + 1, useDingbats = F, 
    onefile = T)
for (xx in phenotypes_sig_diverse_opts) {
  message(xx)
  xx_dist_mat <- as.matrix(dist(mat_pheno_viz_zscore[, xx, drop=F], method = 'manhattan'))
  xx_dist_mat <- scales::rescale(xx_dist_mat, to=cna_dist_val_range)
  p_fea_xx <- pheatmap(
    xx_dist_mat[rev(genotree_tip_orders), rev(genotree_tip_orders)], 
    color = circlize::colorRamp2(
      quantile(c(xx_dist_mat), c(0.01, 0.99)), 
      hcl.colors(n=2, 'Blues 2', rev = F)), 
    name = 'feature dist', 
    display_numbers = T, number_color = 'black', 
    cluster_rows = F, cluster_cols = F, 
    border_color = 'black', 
    main = xx %>% str_remove_all('HALLMARK_|UCell|Seurat|seurat|ucell') %>%
      str_replace_all("_", " ") %>% str_squish() %>% tolower()
  )

  draw(p_cna + p_fea_xx, 
       column_title = sprintf(
         'heritability=%.2f pval=%.3g', 
         phenotypes_sig_diverse_heritability_orig[xx], 
         phenotypes_sig_diverse_heritability_orig_pval[xx]))
  
}
dev.off()


#------ plasticity scoring by PATH ------
library(PATH); library(castor); library(reshape2)
add_tip_edge_len <- function(tree, to=0.01) {
  tree$edge.length[sapply(seq_along(tree$tip.label),function(x,y) which(y==x),y=tree$edge[,2])] <- tree$edge.length[sapply(seq_along(tree$tip.label),function(x,y) which(y==x),y=tree$edge[,2])] + to
  return(tree)
}

# Winv <- inv.tree.dist(
#   add_tip_edge_len(me_tree_nodipoid, 0.01), 
#   node = TRUE, norm = FALSE)
Winv <- inv.tree.dist(
  me_tree_nodipoid, 
  node = TRUE, norm = FALSE)

rownames(Winv) <- colnames(Winv) <- me_tree_nodipoid$tip.label
if (module_method == 'ucell') {
  res_PATH <- xcor(
    data = mat_pheno_viz_zscore[rownames(Winv), , drop=F], 
    weight.matrix = Winv
  )
} 
if (module_method == 'seurat') {
  res_PATH <- xcor(
    data = mat_pheno_viz[rownames(Winv), , drop=F],
    weight.matrix = Winv
  )
}
if (n_clones == 3) {
  # https://github.com/landau-lab/PATH/issues/3#issuecomment-2023762205
  library(pbmcapply)

  permute_trees <- function(reps, z, w, cores) {
    pbmclapply(1:reps, function(r) xcor(z[sample(1:nrow(z), nrow(z)), , drop=F], w)$Moran, 
               mc.cores = cores) %>% 
      melt() %>% tibble
  }
  
  path_s0 <- permute_trees(
    reps = 10^3, 
    z = switch (module_method,
      ucell = mat_pheno_viz_zscore[rownames(Winv), , drop=F], 
      seurat = mat_pheno_viz[rownames(Winv), , drop=F]
    ), 
    w = Winv, 
    cores = 8) 
  
  
  dim(path_s0)  # reps x 4 # var1, var2, value, L1
  path_s <- path_s0 %>% group_by(Var1, Var2) %>% 
    summarize("s" = sd(value)) %>% 
    pivot_wider(id_cols = Var1, names_from = Var2, 
                values_from = s) %>% 
    ungroup() %>% 
    column_to_rownames("Var1") %>% as.matrix
  dim(path_s)
  path_z <- (res_PATH$Morans.I - res_PATH$Expected.I) / path_s
  Zdf <- reshape2::melt(path_z, value.name = 'Z')
} else {
  Zdf <- reshape2::melt(res_PATH$Z.score, value.name = "Z")
}  

Idf <- reshape2::melt(res_PATH$Morans.I, value.name = "I")


df_PATH <- full_join(Idf, Zdf, by=c("Var1", "Var2"))
df_PATH_use <- df_PATH %>% dplyr::filter(Var1 == Var2)
df_PATH_use$phenotype <- df_PATH_use$Var1

df_PATH_use %>%
  ggplot(aes_string(x='phenotype', y=factor(1), fill = 'Z')) +
  geom_tile() + 
  labs(fill='concordance') + 
  # theme(legend.position = 'bottom') + 
  scale_fill_continuous_diverging(palette='Tropic')

phenotypes_sig_diverse_xcor_scores <- df_PATH_use[, c('phenotype', 'Z')] %>% deframe()

report_heritability <- data.frame(
  phenotype = phenotypes_sig_diverse_opts, 
  heritability_cor = phenotypes_sig_diverse_heritability_orig[phenotypes_sig_diverse_opts], 
  heritability_cor_pval = phenotypes_sig_diverse_heritability_orig_pval[phenotypes_sig_diverse_opts], 
  heritability_path = phenotypes_sig_diverse_xcor_scores[phenotypes_sig_diverse_opts], 
  heritability_gearyc = phenotypes_sig_diverse_heritability_gearyc[phenotypes_sig_diverse_opts], 
  heritability_gearyc_pval_up = phenotypes_sig_diverse_heritability_gearyc_pval_up[phenotypes_sig_diverse_opts], 
  heritability_gearyc_pval_dn = phenotypes_sig_diverse_heritability_gearyc_pval_dn[phenotypes_sig_diverse_opts], 
  stringsAsFactors = F)
report_heritability$is_related <- ifelse(report_heritability$heritability_gearyc_pval_up < 0.05 | report_heritability$heritability_gearyc_pval_dn < 0.05, 'related', 'unrelated')
# table(report_heritability$heritability_cor_pval < 0.05)
# view(report_heritability)
write_rds(
  report_heritability, 
  file.path(dir_res, 'report_heritability.rds'))

#------ viz: phylo tree + module score matrix + plasticity ------

#------ 1/ phylo tree ------
pal_cna_clones = new_palette_D(levels(objd@colData$clones), pal = 'stallion')
n_clones <- length(unique(objd@colData$clones))
size_clones <- c(table(fct_drop(objd@colData$clones)))
prop_clones <- c(prop.table(table(fct_drop(objd@colData$clones))))
ggtree_misc_meta <- data.frame(
  Newick_label = c(names(size_clones)),
  n = c(size_clones),
  prop=c(prop_clones))

list_samples <- split(my_sub, my_sub)
if (tree_method == 'medic2Orig') {
  # tree <- ggtree::groupOTU(me_tree, list_samples)
  ## always use the diploid-deleted tree because here this
  ## is highlighting the clone-wise difference
  tree <- ggtree::groupOTU(me_tree_nodipoid, list_samples)
} else if (tree_method == 'medic2ME') {
  tree <- ggtree::groupOTU(me_tree_nodipoid, list_samples)
} else {
  warning('unsupported ', tree_method)
}
treeplt <- ggtree::ggtree(tree, ladderize=F, size = 0.2)
treeplt <- treeplt %<+% ggtree_misc_meta
# view(ggtree_misc_meta)


treeplt <- treeplt + 
  ggtree::geom_tiplab(size=3, aes(color=group),hjust = -0.4, alpha=1)+
  ggtree::geom_tippoint(aes(color=group, size = prop), alpha=1)+
  scale_colour_manual(values = pal_cna_clones) +
  geom_text(aes(x=branch,
                label=plyr::mapvalues(round(branch.length, 1), from = "0",to = ""),
                vjust=-.5), size = 2) +
  theme(legend.position = "none") +
  ggtree::geom_rootpoint() +
  ggtree::theme_tree2() + 
  guides(color = guide_legend(override.aes = list(size = 3)))

ggsave(file.path(dir_res, 'tree.genotree_medic2_viz.pdf'), 
       treeplt, width = 5, height = 4, useDingbats = F)
ggsave(file.path(dir_res, 'paper.tree.genotree_medic2_viz.pdf'), 
       treeplt + scale_colour_manual(values = pal_cna_clones_paper),
       width = 5, height = 4, useDingbats = F)

pL <- treeplt + theme(legend.position = 'left')
pL
#------ 2/ plot of pheno matrix ------

if (n_phenos > 0) {
  pmat_pheno_zscore <- df_pheno_viz_zscore %>% 
    tidyr::pivot_longer(!z, names_to = 'phenotype', values_to = 'val') %>%
    ggplot(aes_string(y=z, x='phenotype', fill='val')) + 
    geom_tile() + 
    scale_x_discrete(
      position = "top", 
      limits = pheno_viz_orders, 
      labels = str_wrap(pheno_viz_orders_pretty, width = 35),
      guide = guide_axis(angle = 45)) + 
    scale_fill_gradient2(  low = muted("blue"),
                           mid = "white",
                           high = muted("red")) + 
    rremove('x.title') + rremove('y.title') + labs(fill='zscore')
  pmat_pheno_score <- df_pheno_viz %>% 
    tidyr::pivot_longer(!z, names_to = 'phenotype', values_to = 'val') %>%
    ggplot(aes_string(y=z, x='phenotype', fill='val')) + 
    geom_tile() + 
    scale_x_discrete(
      position = "top", 
      limits = pheno_viz_orders, 
      labels = str_wrap(pheno_viz_orders_pretty, width = 35),
      guide = guide_axis(angle = 45)) + 
    # scale_fill_viridis_c(option = 'B') +
    scale_fill_gradient2(  low = muted("blue"),
                           mid = "white",
                           high = muted("red")) + 
    rremove('x.title') + rremove('y.title') + labs(fill='score')
}  else {
  pmat_pheno_zscore <- pmat_pheno_score <- NULL
}

#------ 3/ plot of feature plasticity ------
cutoff_herit_corr <- 0.18
cutoff_herit_path <- 0
cutoff_herit_gearyc <- 0.2
viz_herit_corr_min <- -0.5
viz_herit_corr_max <- 1
viz_herit_path_min <- -2
viz_herit_path_max <- 4
viz_herit_gearyc_min <- 0
viz_herit_gearyc_max <- 0.6

pal_heritability_is_related <- c('related' = 'cyan4', 
                                 'unrelated'='orchid')

pTileHerit_orig <- enframe(
  report_heritability %>% dplyr::select(phenotype, heritability_cor) %>% deframe(), 
  name = 'phenotype', value = 'val') %>%
  ggplot(aes_string(x = 'phenotype', y=factor(1), fill='val')) + 
  geom_tile(col='black') + 
  labs(fill='heritability') + 
  scale_fill_continuous_diverging(
    palette='Tropic', 
    limits = c(viz_herit_corr_min, viz_herit_corr_max),
    oob = scales::squish, 
    mid=cutoff_herit_corr, rev=T, na.value='black')
pTileHerit_pval_bi <- enframe(
  report_heritability %>% dplyr::select(phenotype, is_related) %>% deframe(), 
  name = 'phenotype', value = 'val') %>%
  ggplot(aes_string(x = 'phenotype', y=factor(1), fill='val')) + 
  geom_tile(col='black') + 
  labs(fill='is_related') + 
  scale_fill_manual(values = pal_heritability_is_related, na.value='black')
pTileHerit_xcor <- enframe(
  report_heritability %>% dplyr::select(phenotype, heritability_path) %>% deframe(), 
  name = 'phenotype', value = 'val') %>%
  ggplot(aes_string(x = 'phenotype', y=factor(1), fill='val')) + 
  geom_tile(col='black') + 
  labs(fill='heritability (PATH)') + 
  scale_fill_continuous_diverging(
    palette='Tropic', 
    limits = c(viz_herit_path_min, viz_herit_path_max),
    oob = scales::squish, 
    mid=cutoff_herit_path, rev=T, na.value='black')
pTileHerit_gearyc <- report_heritability %>%
  dplyr::mutate(star = ifelse(heritability_gearyc_pval_up < 0.05 | heritability_gearyc_pval_dn < 0.05, '*', '')) %>%
  ggplot(aes_string(x = 'phenotype', y=factor(1), fill='heritability_gearyc')) + 
  geom_tile(col='black') + 
  geom_text(aes(label=star))+
  labs(fill='heritability (gearyC)') + 
  scale_fill_continuous_diverging(
    palette='Tropic', 
    limits = c(viz_herit_gearyc_min, viz_herit_gearyc_max),
    oob = scales::squish, 
    mid=cutoff_herit_gearyc, rev=T, na.value='black')


pTileHerit <- wrap_plots(
  pTileHerit_orig, 
  # pTileHerit_pval_bi,
  pTileHerit_gearyc,
  pTileHerit_xcor, 
  ncol=1)
pTileHerit <- pTileHerit & scale_x_discrete(
  position = "bottom", 
  limits = pheno_viz_orders, 
  labels = str_wrap(pheno_viz_orders_pretty, width = 35),
  guide = guide_axis(angle = 45)) &
  rremove('y.text') & rremove('y.title') & rremove('x.title')
pTileHerit[[1]] <- pTileHerit[[1]] + rremove('x.text')
pTileHerit[[2]] <- pTileHerit[[2]] + rremove('x.text')
pTileHerit <- pTileHerit &
  guides(fill = guide_colourbar(
    frame.colour = 'black', ticks.colour = "black"))

pout_design <- "ZLM
                ##P"
pout1 <- patchwork::wrap_plots(
  Z = as_ggplot(get_legend(pL)), 
  L = pL + rremove('legend'), 
  M = pmat_pheno_zscore + scale_y_discrete(limits = genotree_tip_orders), 
  P = pTileHerit, 
  design = pout_design, 
  heights = c(n_clones, 3), guides = 'collect',
  widths = c(0.3, 1, 3)
)
pout2 <- patchwork::wrap_plots(
  Z = as_ggplot(get_legend(pL)), 
  L = pL + rremove('legend'), 
  M = pmat_pheno_score + scale_y_discrete(limits = genotree_tip_orders), 
  P = pTileHerit, 
  design = pout_design, 
  heights = c(n_clones, 3), guides = 'collect',
  widths = c(0.3, 1, 3)
)

pout1 <- pout1 + plot_annotation(
  title = sprintf('concordance=%.2f', score_GP_concordance)
)
pout2 <- pout2 + plot_annotation(
  title = sprintf('concordance=%.2f', score_GP_concordance)
)

pdfmat_h <- pmax((n_clones+3)*0.35 + 1, 2.5)
pdfmat_w <- pdfmat_h/(n_clones+3) * (n_phenos+1) + 1
pdf(file.path(dir_res, sprintf('concordance.heritability.medic2.%s.zscore.birdview_SigOnly.combo.pdf', name_db)), 
    height = pdfmat_h + 2, width = pdfmat_w + 2, 
    useDingbats = F, onefile = T)
print(pout1)
dev.off()

pdf(file.path(dir_res, sprintf('concordance.heritability.medic2.%s.score.birdview_SigOnly.combo.pdf', name_db)), 
    height = pdfmat_h + 2, width = pdfmat_w + 2, 
    useDingbats = F, onefile = T)
print(pout2)
dev.off()

pdfmat_h <- pmax((n_clones+1)*0.35 + 1, 2.5)
pdfmat_w <- pdfmat_h/(n_clones+1) * (n_phenos+1) + 1
pout3 <- patchwork::wrap_plots(
  Z = as_ggplot(get_legend(pL)), 
  L = pL + rremove('legend'), 
  M = pmat_pheno_zscore + scale_y_discrete(limits = genotree_tip_orders), 
  P = pTileHerit[[3]], 
  design = pout_design, 
  heights = c(n_clones, 1), guides = 'collect',
  widths = c(0.3, 1, 3)
)
pout3 <- pout3 + plot_annotation(
  title = sprintf('concordance=%.2f', score_GP_concordance)
)
pdf(file.path(dir_res, sprintf('concordance.heritability.medic2.%s.zscore.birdview_SigOnly.PATH.pdf', name_db)), 
    height = pdfmat_h + 2, width = pdfmat_w + 2, 
    useDingbats = F, onefile = T)
print(pout3)
dev.off()

pdfmat_h <- pmax((n_clones+2)*0.35 + 1, 2.5)
pdfmat_w <- pdfmat_h/(n_clones+2) * (n_phenos+1) + 1
pTileHerit <- wrap_plots(
  pTileHerit_orig, 
  pTileHerit_xcor, 
  ncol=1)
pTileHerit <- pTileHerit & scale_x_discrete(
  position = "bottom", 
  limits = pheno_viz_orders, 
  labels = str_wrap(pheno_viz_orders_pretty, width = 35),
  guide = guide_axis(angle = 45)) &
  rremove('y.text') & rremove('y.title') & rremove('x.title')
pTileHerit[[1]] <- pTileHerit[[1]] + rremove('x.text')
pTileHerit
pout4 <- patchwork::wrap_plots(
  Z = as_ggplot(get_legend(pL)), 
  L = pL + rremove('legend'), 
  M = pmat_pheno_zscore + scale_y_discrete(limits = genotree_tip_orders), 
  P = pTileHerit, 
  design = pout_design, 
  heights = c(n_clones, 2), guides = 'collect',
  widths = c(0.3, 1, 3)
)
pout4 <- pout4 + plot_annotation(
  title = sprintf('concordance=%.2f', score_GP_concordance)
)
pdf(file.path(dir_res, sprintf('concordance.heritability.medic2.%s.zscore.birdview_SigOnly.combo2.pdf', name_db)), 
    height = pdfmat_h + 2, width = pdfmat_w + 2, 
    useDingbats = F, onefile = T)
print(pout4)
dev.off()

#------ replot the consensus heatmap ------
## following the clone orders
css_mat <- read_rds(file.path(dir_medic2, 'css_mat.rds'))
objd@consensus <- css_mat

para.plotHeatmap.group <- 'run'
if (length(unique(objd@colData[['run']])) == 1) {
  para.plotHeatmap.group = NULL
} 
objd@consensusPhylo <- me_tree_nodipoid
pdf(file.path(dir_res, sprintf('cna.heatmap_consensus_int.combo_coda.%s.pdf', z)), 
    width = 15, height = 5, onefile = F, useDingbats = F)
p <- try( draw( plotHeatmap(
  objd, label = z, 
  assay = 'integer',
  order_cells = 'consensus_tree',
  # col = pal_copykit_int_clip_0_2_6,
  label_colors = list(clones=pal_cna_clones),
  group = para.plotHeatmap.group, consensus=T) ) )
dev.off()

pdf(file.path(dir_res, sprintf('paper.cna.heatmap_consensus_int.combo_coda.%s.pdf', z)), 
    width = 15, height = 5, onefile = F, useDingbats = F)
p <- try( draw( plotHeatmap(
  objd, label = z, 
  assay = 'integer',
  order_cells = 'consensus_tree',
  # col = pal_copykit_int_clip_0_2_6,
  label_colors = list(clones=pal_cna_clones_paper),
  group = para.plotHeatmap.group, consensus=T) ) )
dev.off()
message(paste('Done: ', sample_name, rna_assay_use, module_method, name_db, tree_method))
