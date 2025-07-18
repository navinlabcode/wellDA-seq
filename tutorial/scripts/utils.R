library(patchwork); library(ggpubr)
library(tidyverse)
#------------------ ~~~ Misc ~~~ --------------------
write_gmx2 <- function(L, path) {
  column_names <- names(L)
  if (purrr::is_empty(column_names)) {
    column_names <- as.character(seq_along(L))
  }
  n_row <- max(sapply(L, length))
  n_col <- length(L)
  df <- as.data.frame(matrix(NA, nrow = n_row, ncol = n_col))
  for (j in seq_len(n_col)) {
    i_to <- length(L[[j]])
    if (i_to == 0) {
      next()
    }
    df[1:i_to, j] <- L[[j]]
  }
  colnames(df) <- column_names
  readr::write_csv(x = df, file = path, col_names = TRUE, na = "")
}
#-------------------------- coda --------------------------  
library(tidyverse)
library(ggrastr)
library(scales)
library(ggplot2)

if (F) {
  library(ggrastr)
  library(scales)
  library(ggplot2)
  data("iris")
  head(iris)

  emb_a_x <- 'Sepal.Length' 
  emb_a_y <- 'Sepal.Width'
  emb_b_x <- 'Petal.Length'
  emb_b_y <- 'Petal.Width'
  
  crossdimplot <- function(
    df, emb_a_x, emb_a_y, emb_b_x, emb_b_y, group.by, 
    assay_a_label='A', assay_b_label='B',
    pt.size=2, 
    seg_alpha = 0.5, seg_lwd = 0.5,
    byrow=NULL, bycol=NULL,
    do_shuffle = NULL) {
    
    
    emb_a <- as_tibble(df[, c(emb_a_x, emb_a_y)])
    emb_b <- as_tibble(df[, c(emb_b_x, emb_b_y)])
    z <- df[, group.by, drop=T]
    if (all(is.null(c(byrow, bycol)))) {
      byrow <- T; bycol <- F
    } else {
      if (is.null(bycol)) {bycol <- !byrow}
      if (is.null(byrow)) {byrow <- !bycol}
    }
    spacing_h <- 0.3; spacing_v <- 0
    if (byrow & !bycol) {spacing_h <- 0.3; spacing_v <- 0}
    if (!byrow & bycol) {spacing_h <- 0; spacing_v <- 0.3}
    
    emb_a <- emb_a %>% mutate(across(where(is.numeric), function(x) {scales::rescale(x, to=c(0.1, 1.1))} ))
    emb_b <- emb_b %>% mutate(across(where(is.numeric), function(x) {scales::rescale(x, to=c(0.1, 1.1))} ))
    if (spacing_h > 0) {
      emb_b[, emb_b_x] <- emb_b[, emb_b_x] + diff(range(emb_a[, emb_a_x])) + spacing_h
    }
    if (spacing_v > 0) {
      emb_b[, emb_b_y] <- emb_b[, emb_b_y] + diff(range(emb_a[, emb_a_y])) + spacing_v
    }
    
    df <- cbind(emb_a, emb_b, z)
    
    if (!is.null(do_shuffle)) {
      set.seed(42)
      df <- df[sample(seq_len(nrow(df)), size = nrow(df), replace = F), ]
    }
    p <- ggplot() +
      ggrastr::rasterize(
        geom_segment(
          data = df, 
          aes_string(x=emb_a_x, y=emb_a_y,
                     xend=emb_b_x, yend=emb_b_y,
                     color='z'), 
          alpha=seg_alpha, 
          lwd=seg_lwd), dpi=300)+
      
      geom_point_rast(
        data=df, aes_string(x=emb_a_x, y=emb_a_y, fill='z'),
        size = pt.size, pch=21, alpha=1, stroke=0.1, color='black') +
      
      geom_point_rast(
        data=df, aes_string(x=emb_b_x, y=emb_b_y, fill='z'), 
        size = pt.size, pch=21, alpha=1, stroke=0.1, color='black') + 
      coord_equal() + 
      theme_void() +
      geom_rect(aes(xmin=min(emb_a[, emb_a_x])-0.1,
                    xmax=max(emb_a[, emb_a_x])+0.1,
                    ymin=min(emb_a[, emb_a_y])-0.1,
                    ymax=max(emb_a[, emb_a_y])+0.1), 
                fill=NA, color='black'
      ) +
      geom_rect(aes(xmin=min(emb_b[, emb_b_x])-0.1,
                    xmax=max(emb_b[, emb_b_x])+0.1,
                    ymin=min(emb_b[, emb_b_y])-0.1,
                    ymax=max(emb_b[, emb_b_y])+0.1), 
                fill=NA, color='black'
      )
    p <- p + 
      annotate('text', 
               x=min(emb_a[, emb_a_x])-0.1,
               y=max(emb_a[, emb_a_y])+0.1,
               label=assay_a_label) +
      annotate('text',
               x=min(emb_b[, emb_b_x])-0.1,
               y=max(emb_b[, emb_b_y])+0.1,
               label=assay_b_label)
    p <- p + labs(fill=group.by, color=group.by)
    return(p)
  }
  
  crossdimplot(
    df=iris, 
    emb_a_x, emb_a_y, emb_b_x, emb_b_y, 
    group.by = 'Species', byrow = T, 
    pt.size = 2, seg_lwd = 1, seg_alpha = .1) +
    scale_color_viridis_d() + 
    scale_fill_viridis_d()
  
  crossdimplot(
    df=iris, 
    emb_a_x, emb_a_y, emb_b_x, emb_b_y, 
    group.by = 'Species', byrow = T, 
    pt.size = 1, seg_lwd = 2, seg_alpha = .1) +
    scale_color_viridis_d() + 
    scale_fill_viridis_d()
  
  crossdimplot(
    df=iris, 
    emb_a_x, emb_a_y, emb_b_x, emb_b_y, 
    group.by = 'Species', byrow = F, 
    pt.size = 1, seg_lwd = 2, seg_alpha = .1) +
    scale_color_viridis_d() + 
    scale_fill_viridis_d()
}

coda_crossdimplot2 <- function(
    df, group.by, 
    pal = c('stallion'), 
    pt.size=2, 
    seg_alpha = 0.5, seg_lwd = 0.5,
    plot_order = NULL,
    color_specific = NULL, 
    nrow = NULL, ncol = NULL, do_shuffle = NULL) {

  emb_a_x <- 'UMAP_1_CNA'
  emb_a_y <- 'UMAP_2_CNA'  
  emb_b_x <- 'UMAP_1_ATAC' 
  emb_b_y <- 'UMAP_2_ATAC'

  
  emb_a <- as_tibble(df[, c(emb_a_x, emb_a_y)])
  emb_b <- as_tibble(df[, c(emb_b_x, emb_b_y)])
  z <- df[, group.by, drop=T]
  
  spacing_h <- 0.3
  spacing_v <- 0
  emb_a <- emb_a %>% mutate(across(where(is.numeric), function(x) {scales::rescale(x, to=c(0.1, 1.1))} ))
  emb_b <- emb_b %>% mutate(across(where(is.numeric), function(x) {scales::rescale(x, to=c(0.1, 1.1))} ))
  if (spacing_h > 0) {
    emb_b[, emb_b_x] <- emb_b[, emb_b_x] + diff(range(emb_a[, emb_a_x])) + spacing_h
  }
  if (spacing_v > 0) {
    emb_b[, emb_b_y] <- emb_b[, emb_b_y] + diff(range(emb_a[, emb_a_y])) + spacing_v
  }
  
  df <- cbind(emb_a, emb_b, z)
  
  if (is.null(plot_order) & !is.null(do_shuffle)) {
    set.seed(42)
    df <- df[sample(seq_len(nrow(df)), size = nrow(df), replace = F), ]
  }
  p <- ggplot() +
    ggrastr::rasterize(
      geom_segment(
        data = df, 
        aes_string(x=emb_a_x, y=emb_a_y,
                   xend=emb_b_x, yend=emb_b_y,
                   color='z'), 
        alpha=seg_alpha, 
        lwd=seg_lwd), dpi=300)+
    
    geom_point_rast(
      data=df, aes_string(x=emb_a_x, y=emb_a_y, fill='z'),
      size = pt.size, pch=21, alpha=1, stroke=0.1, color='black') +
    
    geom_point_rast(
      data=df, aes_string(x=emb_b_x, y=emb_b_y, fill='z'), 
      size = pt.size, pch=21, alpha=1, stroke=0.1, color='black') + 
    coord_equal() + 
    theme_void() +
    geom_rect(aes(xmin=min(emb_a[, emb_a_x])-0.1,
                  xmax=max(emb_a[, emb_a_x])+0.1,
                  ymin=min(emb_a[, emb_a_y])-0.1,
                  ymax=max(emb_a[, emb_a_y])+0.1), 
              fill=NA, color='black'
    ) +
    geom_rect(aes(xmin=min(emb_b[, emb_b_x])-0.1,
                  xmax=max(emb_b[, emb_b_x])+0.1,
                  ymin=min(emb_b[, emb_b_y])-0.1,
                  ymax=max(emb_b[, emb_b_y])+0.1), 
              fill=NA, color='black'
    )
  p <- p + 
    annotate('text', 
             x=min(emb_a[, emb_a_x])-0.1,
             y=max(emb_a[, emb_a_y])+0.1,
             label='DNA') +
    annotate('text',
             x=min(emb_b[, emb_b_x])-0.1,
             y=max(emb_b[, emb_b_y])+0.1,
             label='ATAC')
  p <- p + labs(fill=group.by)
  return(p)
}





coda_dimplot2 <- function(df, group.by, 
                          pal = c('stallion'), 
                          pt.size=1.5, 
                          plot_order = NULL,
                          color_specific = NULL, 
                          nrow = NULL, ncol = NULL, do_shuffle = NULL) {
  umap1_a <- 'UMAP_1_ATAC'
  umap2_a <- 'UMAP_2_ATAC'
  umap1_d <- 'UMAP_1_CNA'
  umap2_d <- 'UMAP_2_CNA'
  
  stopifnot(all(c(umap1_a, umap2_a, umap1_d, umap2_d) %in% colnames(df)))
  
  ident_lv <- levels(df[, group.by])
  if (any(is.null(ident_lv))) {ident_lv <- gtools::mixedsort(unique(df[, group.by]))}
  pal <- new_palette_D(ident_lv, pal=pal)
  
  if (!is.null(color_specific)) {
    color_specific <- color_specific[intersect(names(color_specific), names(pal))]
    for (i in 1:length(color_specific)) {
      pal[names(color_specific)[i]] <- color_specific[i]
    }
  }
  
  if (!is.null(plot_order) & is.null(do_shuffle)) {
    df[, group.by] <- fct_relevel(df[, group.by], plot_order)
    df <- df[order(df[, group.by], decreasing = T), ]
  }
  if (is.null(plot_order) & !is.null(do_shuffle)) {
    set.seed(42)
    df <- df[sample(seq_len(nrow(df)), size = nrow(df), replace = F), ]
  }
  pa <- ggplot(df, aes_string(x=umap1_a, y=umap2_a, color=group.by)) + 
    geom_point_rast(pch=16, size=pt.size) + 
    labs(title = '', caption = 'ATAC space', fill='', color='')
  pd <- ggplot(df, aes_string(x=umap1_d, y=umap2_d, color=group.by)) + 
    geom_point_rast(pch=16, size=pt.size) + 
    labs(title = '', caption = 'CNA space', fill='', color='')
  
  pa <- pa + theme_void() + border() + theme(aspect.ratio = 1) +
    guides(color = guide_legend(override.aes = list(size = 3)))
  pd <- pd + theme_void() + border() + theme(aspect.ratio = 1) +
    guides(color = guide_legend(override.aes = list(size = 3)))
  
  ident_tab <- ruok::pretty_table2str(table(df[, group.by]))
  
  po <-  wrap_plots(
    pa + scale_color_manual(values = pal, labels=ident_tab) + 
      theme(plot.margin = unit(c(0,10,0,0), "pt")) + 
      labs(title = sprintf('%s (n=%s)', group.by, nrow(df))), 
    pd + scale_color_manual(values = pal, labels=ident_tab) + 
      theme(plot.margin = unit(c(0,0,0,10), "pt")), 
    guides = 'collect', nrow=nrow, ncol=ncol
  )
  return(po)
}

coda_dimplot <- function(A, D, group.by, pal = c('stallion', 'circus', 'ironMan')) {
  pal <- match.arg(pal)
  pal <- ArchR::paletteDiscrete(unique(A@meta.data[, group.by]), set = pal)
  
  p1 = DimPlot(A, reduction = 'umap', group.by = group.by, pt.size = 1.5) + 
    labs(title = '', caption = 'ATAC', fill='', color='')
  p1 <- p1 + theme_void() + border() + theme(aspect.ratio = 1) 
  
  p2 = plotUmap(D, label= group.by) + 
    labs(caption = 'DNA', fill='') + theme_void() + border() + theme(aspect.ratio = 1)
  
  po = wrap_plots(
    p1 + scale_color_manual(values = pal) + theme(plot.margin = unit(c(0,10,0,0), "pt")), 
    p2 + scale_fill_manual(values = pal) + theme(plot.margin = unit(c(0,0,0,10), "pt")), 
    guides = 'collect')
  po + labs(title = group.by)
}

library(colorspace)
coda_featureplot <- function(A, D, feature, pal='Blues') {

  #ArchR::ArchRPalettes
  # [1] "stallion"     "stallion2"    "calm"        
  # [4] "kelly"        "bear"         "ironMan"     
  # [7] "circus"       "paired"       "grove"       
  # [10] "summerNight"  "zissou"       "darjeeling"  
  # [13] "rushmore"     "captain"      
  # "horizon"     
  # [16] "horizonExtra" "blueYellow"   "sambaNight"  
  # [19] "solarExtra"   "whitePurple"  "whiteBlue"   
  # [22] "whiteRed"     "comet"        "greenBlue"   
  # [25] "beach"        "coolwarm"     "fireworks"   
  # [28] "greyMagma"    "fireworks2"   "purpleOrange"
  
  # pal <- match.arg(pal)
  
  p1 = FeaturePlot(A, reduction = 'umap', features = feature, pt.size = 1.5) + 
    labs(title = '', caption = 'ATAC', fill='', color='')
  p1 <- p1 + theme_void() + border() + theme(aspect.ratio = 1) 
  
  p2 = plotUmap(D, label= feature)+labs(caption = 'DNA', fill='') 
  p2 <-  p2 + theme_void() + border() + theme(aspect.ratio = 1)
  
  po = wrap_plots(
    p1 + scale_color_continuous_sequential(palette = pal) + theme(plot.margin = unit(c(0,10,0,0), "pt")),
    p2 + scale_fill_continuous_sequential(palette = pal) + theme(plot.margin = unit(c(0,0,0,10), "pt")), 
    guides = 'collect')
  
  po + labs(title = feature)
}


coda_featureplot2 <- function(df, feature, 
                              pal = c('Plasma'), 
                              pt.size=1.5, 
                              plot_order = c('rand', 'sort'),
                              nrow = NULL, ncol = NULL) {
  umap1_a <- 'UMAP_1_ATAC'
  umap2_a <- 'UMAP_2_ATAC'
  umap1_d <- 'UMAP_1_CNA'
  umap2_d <- 'UMAP_2_CNA'
  plot_order <- match.arg(plot_order)
  stopifnot(all(c(umap1_a, umap2_a, umap1_d, umap2_d) %in% colnames(df)))
  
  if (plot_order == 'rand') {
    idx <- sample(1:nrow(df), size = nrow(df), replace = F)
  } else if (plot_order == 'sort') {
    idx <- order(df[, feature], decreasing = F)
  } else {
    idx <- 1:nrow(df)
  }
  df <- df[idx, ]
  
  pa <- ggplot(df, aes_string(x=umap1_a, y=umap2_a, color=feature)) + 
    geom_point_rast(pch=16, size=pt.size) + 
    labs(title = '', caption = 'ATAC space', fill='', color='')
  pd <- ggplot(df, aes_string(x=umap1_d, y=umap2_d, color=feature)) + 
    geom_point_rast(pch=16, size=pt.size) + 
    labs(title = '', caption = 'CNA space', fill='', color='')
  
  pa <- pa + theme_void() + border() + theme(aspect.ratio = 1) 
  pd <- pd + theme_void() + border() + theme(aspect.ratio = 1) 
  
  po = wrap_plots(
    pa + scale_color_continuous_sequential(palette = pal) + 
      theme(plot.margin = unit(c(0,10,0,0), "pt")) + 
      labs(title = feature),
    pd + scale_color_continuous_sequential(palette = pal) + 
      theme(plot.margin = unit(c(0,0,0,10), "pt")), 
    guides = 'collect', nrow=nrow, ncol=ncol)
  

  return(po)
}

left_update <- function(x, y) {
  ## x: old data frame
  ## y: new data frame to borrow information
  stopifnot(all.equal(rownames(x), rownames(y)))
  cols_x <- colnames(x)
  cols_y <- colnames(y)
  
  ## new columns in y
  cols_y_uniq <- setdiff(cols_y, cols_x)
  cols_xy <- intersect(cols_y, cols_x)
  if (length(cols_y_uniq)!=0) {
    cat('adding', cols_y_uniq, '...\n')
    x <- cbind(x, y[, cols_y_uniq, drop=F])
  }
  ## updated columns by y
  if (length(cols_xy) > 0) {
    for (j in cols_xy) {
      if (identical(x[, j], y[, j])) { next() }
      cat('updating', j, '...')
      x[,j] <- y[,j]
    }
  }
  cat('\n')
  return(x)
}

left_update2 <- function(x, y1, y2) {
  # x: old
  # y1: new
  # y2: new
  stopifnot(all.equal(rownames(x), rownames(y1)))
  stopifnot(all.equal(rownames(x), rownames(y2)))
  o <- x
  cols_x <- colnames(x)
  cols_y1 <- colnames(y1)
  cols_y2 <- colnames(y2)
  uniq_cols_xy1 <- setdiff(cols_y1, cols_x)
  uniq_cols_xy2 <- setdiff(cols_y2, cols_x)
  # sh_cols_xy1 <- intersect(cols_y1, cols_x)
  # sh_cols_xy2 <- intersect(cols_y2, cols_x)
  sh_cols_all <- intersect(cols_x, intersect(cols_y1, cols_y2))
  ## add new columns
  if (length(uniq_cols_xy1)!=0) {
    cat('adding', uniq_cols_xy1, ' from y1 ...\n')
    o <- cbind(o, y1[, uniq_cols_xy1, drop=F])
  }
  if (length(uniq_cols_xy2)!=0) {
    cat('adding', uniq_cols_xy2, ' from y2 ...\n')
    o <- cbind(o, y2[, uniq_cols_xy2, drop=F])
  }
  ## update columns
  if (length(sh_cols_all) > 0) {
    for (j in sh_cols_all) {
      
      if (!identical(x[, j], y1[, j]) & identical(x[, j], y2[, j]) ) {
        cat('updating', j, ' from y1 ...\n')
        o[,j] <- y1[,j]
      }
      
      if (!identical(x[, j], y2[, j]) & identical(x[, j], y1[, j]) ) {
        cat('updating', j, ' from y2 ...\n')
        o[,j] <- y2[,j]
      }
    
      if (identical(y1[, j], y2[, j])) { 
        if (!identical(x[, j], y1[, j])) {
          cat('updating', j, ' that is same in y1 and y2...\n')
          o[,j] <- y1[,j]
        }
      }
      
      
    }
  }
  cat('\n')
  ## remove duplicated column names
  cat('remove columns with duplicated names:', 
      colnames(o)[duplicated(colnames(o))], '...\n')
  o <- o[, !duplicated(colnames(o)), drop=F]

  return(o)
}

coda_trackplot <- function(A, D, feature, group.by, pal) {
  DefaultAssay(A) <- 'peaks'
  p <- CoveragePlot(
    A, region = feature, 
    # features = 'SYNGR3',
    group.by=group.by,
    # expression.assay='iRNA', 
    # tile = T, tile.cells = 10,
    extend.upstream = 2e3, extend.downstream=2e3)
  # track plot
  p[[1]][[1]] <- p[[1]][[1]] + scale_fill_manual(values = pal)
  # tile plot
  p[[1]][[2]] <- p[[1]][[2]] + rremove('legend')
  # gene anno
  p[[1]][[3]] <- p[[1]][[3]] + scale_color_manual(values=c('slategray', 'slategray2'))
  # peak anno
  p[[1]][[4]] <- p[[1]][[4]] + scale_color_manual(values=c('black'))
  # p
  
  p2 <- copykit::plotGeneCopy(D, genes = feature, label = group.by, geom='swarm') +
    scale_fill_manual(values = pal)
  
  DefaultAssay(A) <- 'iRNA'
  p3 <- VlnPlot(A, feature, group.by = group.by, pt.size = 0.1) + 
    scale_fill_manual(values = pal) + 
    rremove('x.title')

  wrap_plots(p, 
             p3 + rremove('legend'), 
             p2 + rremove('legend'), widths = c(7, 3.5, 1))

}


update_coda <- function(df, A, D) {
  df_a <- A@meta.data
  df_d <- as.data.frame(colData(D))
  umap_a <- as.data.frame(A@reductions$umap@cell.embeddings)
  colnames(umap_a) <- c('UMAP_1_ATAC', 'UMAP_2_ATAC')
  umap_d <- as.data.frame(reducedDim(D, 'umap'))
  colnames(umap_d) <- c('UMAP_1_CNA', 'UMAP_2_CNA')
  stopifnot( all.equal( rownames(umap_a), rownames(umap_d) ) )
  
  coda_metadf <- cbind(umap_a, umap_d, df)
  coda_metadf <- left_update2(coda_metadf, df_a, df_d) #// resolve conflicts and update info
  
  A@meta.data <- subset(coda_metadf, select = -c(UMAP_1_ATAC, UMAP_2_ATAC, UMAP_1_CNA, UMAP_2_CNA))
  colData(D) <- as(subset(coda_metadf, select = -c(UMAP_1_ATAC, UMAP_2_ATAC, UMAP_1_CNA, UMAP_2_CNA)), 'DFrame')
  
  # return(list(df_meta = coda_metadf, obja=A, objd=D))
  assign('obja', A, envir = .GlobalEnv); message('See ATAC result as obja.')
  assign('objd', D, envir = .GlobalEnv); message('See CNA result in objd.')
  assign('df_meta', coda_metadf, envir = .GlobalEnv); message('See meta info in df_meta.')
  return(invisible(df_meta))
}

subset_coda <- function(df, A, D, cells) {
  idx <- match(cells, rownames(df))
  df <- df[idx, ]
  A <- subset(A, cells = cells)
  D <- D[, cells]
  stopifnot(identical(rownames(df), as.character(Cells(A))))
  stopifnot(identical(rownames(df), as.character(colnames(D))))
  # D <- runDistMat(D, metric = 'manhatten', n_threads = 20)
  assign('obja', A, envir = .GlobalEnv); message('See ATAC result as obja.')
  assign('objd', D, envir = .GlobalEnv); message('See CNA result in objd.')
  assign('df_meta', df, envir = .GlobalEnv); message('See meta info in df_meta.')
  return(invisible(0))
}

write_coda <- function(o, df, A, D) {
  write_rds(df, file.path(o, 'metadata.df.rds'))
  write_csv(df, file.path(o, 'metadata.df.csv'))
  write_rds(A,  file.path(o, 'obja.rds'))
  write_rds(D,  file.path(o, 'objd.rds'))
  return(invisible(0))
}

load_coda <- function(o) {
  objd    <- read_rds(file.path(o, 'objd.rds'))
  obja    <- read_rds(file.path(o, 'obja.rds'))
  df_meta <- read_rds(file.path(o, 'metadata.df.rds')) 
  assign('obja', obja, envir = .GlobalEnv); message('See ATAC result as obja.')
  assign('objd', objd, envir = .GlobalEnv); message('See CNA result in objd.')
  assign('df_meta', df_meta, envir = .GlobalEnv); message('See meta info in df_meta.')
  return(invisible(0))
}


CNALinePlot <- function(objd, region, group.by, ident) {
  region <- StringToGRanges(region)
  varbins_intersect <- subsetByOverlaps(x = rowRanges(objd), ranges = region, ignore.strand=TRUE)
  if (length(varbins_intersect)==0) {
    cat('cannot find overlap of varbins and the requested region')
  }
  varbins_df <- as.data.frame(x = varbins_intersect)
  varbins_df$BIN_str <- GRangesToString(varbins_intersect)
  start.pos <- start(region)
  end.pos <- end(region)
  chromosome <- seqnames(region)
  
  
  lineplot <- ggplot()
  if (nrow(varbins_df) > 0) {
    
    csc_mat <- consensus(objd)[, ident, drop=F]
    mdx <- match(varbins_df$BIN_str, GRangesToString(rowRanges(objd)))
    csc_mat <- csc_mat[mdx, , drop=F]
    csc_mat <- as.data.frame(csc_mat)
    varbins_df <- cbind(varbins_df, csc_mat)
    
    varbins_df$start[varbins_df$start < start.pos] <- start.pos
    varbins_df$end[varbins_df$end > end.pos] <- end.pos
    varbins_df <- tidyr::pivot_longer(
      varbins_df, cols = ident, names_to = group.by, values_to = 'CNA')
    varbins_df[[group.by]] <- factor(
      as.character(varbins_df[[group.by]]), 
      levels = ident
    )
    lineplot <- ggplot(varbins_df, aes_string(color=group.by)) + 
      geom_hline(yintercept = c(0, 2)) +
      geom_segment(aes(x = start, y = CNA, xend = end, yend = CNA),
                   size = 2,
                   data = varbins_df) +
      xlim(c(start.pos, end.pos)) +
      xlab(label = paste0(chromosome, " position (bp)")) 
  }
  return(lineplot)
}
#-------------------------- ATAC --------------------------  

bgzip_tabix_fragments <- function(f, bgzip, tabix, run_bgzip=T, run_tabix=T){
    fo <- paste0(f, '.gz')
    cmd <- sprintf('%s -f -c %s > %s', bgzip, f, fo)
    if (run_bgzip) { system(cmd) }
    cmd <- sprintf('%s -f -p bed %s', tabix, fo)
    if (run_tabix) { system(cmd) }
    return(fo)    
}

# }
# bgzip_tabix_fragments(
#     f=file.path(dir_tosignac, 'fragments.tsv'),
#     bgzip='/volumes/lab/users/yyan/apps/anaconda3/bin/bgzip',
#     run_bgzip = T,
#     tabix='/volumes/lab/users/yyan/apps/anaconda3/bin/tabix', 
#     run_tabix = T)
# fpath_frag <- file.path(dir_tosignac, 'fragments.tsv.gz')

get_density2 <- function(x = NULL, y = NULL, n = 100, sample = NULL, densityMax = 0.95){
    #modified from http://slowkow.com/notes/ggplot2-color-by-density/
    #https://rdrr.io/github/GreenleafLab/ArchR/src/R/GgplotUtils.R
    df <- data.frame(x=x,y=y)
    dens <- MASS::kde2d(x = x, y = y, n = n)
    ix <- findInterval(x, dens$x)
    iy <- findInterval(y, dens$y)
    ii <- cbind(ix, iy)
    df$density <- dens$z[ii]
    df$density[df$density > quantile(unique(df$density),densityMax)] <- quantile(unique(df$density),densityMax) #make sure the higher end doesnt bias colors
    if(!is.null(sample)){
        df <- df[sample(nrow(df), min(sample,nrow(df))),]
    }
    return(df)
}

## ArchR to Signac object
library(Signac); library(Seurat); library(tibble); library(readr)

#' ConvertArchRtoSignac
#' @param archr Archr object
#' @param dir_tosignac folder path to save all the results
#' @param do_binarization TRUE if the orignal ArchR's count matrix is binarized. 
#' @param bgzip path to the bgzip command
#' @param tabix path to the tabix command
#' @return A signac object
#' @examples
#' archr <- loadArchRProject(dir_archr)
#' signac <- ConvertArchRtoSignac(
#'   archr, 'path/to/out', do_binarization=TRUE,
#'   bgzip='/volumes/USR1/yyan/anaconda3/bin/bgzip', 
#'   tabix='/volumes/USR1/yyan/anaconda3/bin/tabix')
#' @author yun yan (yun.yan at uth.tmc.edu)
ConvertArchRtoSignac <- function(archr, dir_tosignac, do_binarization=TRUE, 
                                 bgzip = 'bgzip', tabix = 'tabix') {
  ## 1/5. peaks.bed
  message('## 1/5. peaks.bed')
  gr_peaks <- getPeakSet(archr)
  length(gr_peaks)
  data.frame(
    chr=seqnames(gr_peaks),
    start=start(gr_peaks)-1,
    end = end(gr_peaks)) %>%
    readr::write_tsv(., path = file.path(dir_tosignac, 'peaks.bed'), col_names = F)
  
  getAvailableMatrices(archr)
  
  ## 2/5. singlecell.csv
  message('## 2/5. singlecell.csv')
  write.csv(x = as.data.frame(getCellColData(archr)), 
            file = file.path(dir_tosignac, 'singlecell.csv'))
  ## 3/5. peak x cells
  message('## 3/5. peak x cells')
  mat_peaks_cells <- getMatrixFromProject(archr, useMatrix='PeakMatrix', 
                                          binarize = do_binarization)
  dim(mat_peaks_cells)
  write_rds(mat_peaks_cells, 
            file.path(dir_tosignac, 'filtered_peak_bc_matrix.se.rds'))
  ## option. Genome Bin matrix. 
  # mat_tile_cells <- getMatrixFromProject(archr, useMatrix = 'TileMatrix',
  #                                        binarize = do_binarization)
  # write_rds(mat_tile_cells,
  #           file.path(dir_tosignac, 'filtered_tile_bc_matrix.se.rds'))
  ## 4/5. gene_score x cells
  message('## 4/5. gene_score x cells')
  mat_irna_cells <- getMatrixFromProject(archr, useMatrix = 'GeneScoreMatrix')
  write_rds(mat_irna_cells, 
            file.path(dir_tosignac, 'irna_bc_matrix.se.rds'))
  ## 5/5. fragments... takes time...
  message('## 5/5. fragments... takes time...')
  fragment_gr <- getFragmentsFromArrow(getArrowFiles(archr))
  seqlevelsStyle(fragment_gr) <- 'UCSC'
  fragment_gr <- sort(sortSeqlevels(fragment_gr))
  tibble(seqnames = as.character(seqnames(fragment_gr)),
         start = start(fragment_gr)-1,
         end = end(fragment_gr),
         RG= as.character(fragment_gr$RG),
         N=1) %>% 
    write_tsv(file.path(dir_tosignac, 'fragments.tsv'), col_names = F)
  
  bgzip_tabix_fragments(
    f=file.path(dir_tosignac, 'fragments.tsv'),
    bgzip=bgzip,
    run_bgzip = T,
    tabix=tabix,
    run_tabix = T)
  fpath_frag <- file.path(dir_tosignac, 'fragments.tsv.gz')
  ## 6/ optional. iLSI matrix
  message('## 6/ optional. iLSI matrix')
  mat_ilsi <- getReducedDims(archr, reducedDims = 'IterativeLSI', 
                             corCutOff=0.75)
  colnames(mat_ilsi) <- paste0('iLSI_', 1:ncol(mat_ilsi))
  
  message('Creating Signac Object')
  signac_peaks <- read_tsv(file.path(dir_tosignac, 'peaks.bed'), col_names = F)
  colnames(signac_peaks) <- c('seqname', 'start', 'end')
  signac_peaks <- makeGRangesFromDataFrame(signac_peaks, starts.in.df.are.0based = TRUE)
  signac_counts <- assays(mat_peaks_cells)$PeakMatrix
  rownames(signac_counts) <- GRangesToString(signac_peaks, sep=c(':', '-'))
  signac_cellmeta <- read.csv(file.path(dir_tosignac, 'singlecell.csv'), row.names = 1, stringsAsFactors = F)
  signac_genemat <- assays(mat_irna_cells)$GeneScoreMatrix
  rownames(signac_genemat) <- make.unique(rowData(mat_irna_cells)$name)
  
  
  chrom_assay <- CreateChromatinAssay(
    counts = signac_counts,
    sep = c(":", "-"),
    fragments = fpath_frag,
    min.cells = 1,
    min.features = 0
  )
  idx <- match(Cells(chrom_assay), rownames(signac_cellmeta))
  sr <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks",
    meta.data = signac_cellmeta[idx, ]
  )
  sr[['iRNA']] <- CreateAssayObject(data = signac_genemat) ## This is a normalized-to-1e4 mat
  sr[['iLSI']] <- Seurat::CreateDimReducObject(
    embeddings = mat_ilsi, key='iLSI_', assay='peaks')
  
  message('Saving Signac Object')
  write_rds(sr, file.path(dir_tosignac, 'ready.signac_default.rds'))
  
  return(sr)
  
}


if (F) {
  # https://github.com/fzhaouf/scaDA
estParams <- function(object, group.1=NULL, group.2=NULL){
  message("start initial parameter estiamte")
  
  count <- object@count
  # normlaization factor
  sfs <- apply(count,2,mean); sfs=sfs/median(sfs)
  
  if(is.null(group.2)){
    cells_loc <- sample.loc(object=object, group.1 = group.1, method="balanced")
    group.1.loc <- cells_loc$group.1.loc  # assuming group 1 is case
    group.2.loc <- cells_loc$group.2.loc  # assuming group 2 is reference
  } else { # if both group.1 and group.2 are specified
    group.1.loc <- which(object@colData==group.1)  # assuming c1 is case
    group.2.loc <- which(object@colData==group.2)  # assuming c2 is reference
  }
  # subset dat
  dat <- count[,c(group.1.loc, group.2.loc)]
  npeak <- dim(dat)[1]
  nsam <- dim(dat)[2]
  nsam1 <- length(group.1.loc)
  nsam2 <- length(group.2.loc)
  poolCol <- c(1:nsam)
  cond1Col <- c(1:nsam1)
  cond2Col <- c((nsam1+1):nsam)
  
  counts_pooled <- dat[,poolCol]
  counts_cell1 <- dat[,cond1Col]
  counts_cell2 <- dat[,cond2Col]
  # seq depth adj factors
  sfs_pooled <- sfs[c(group.1.loc, group.2.loc)]
  sfs_cell1 <- sfs[group.1.loc]
  sfs_cell2 <- sfs[group.2.loc]
  # df store estimates
  est_params_pooled <- matrix(NA,dim(counts_pooled)[1],ncol=3)
  est_params_cell1 <- matrix(NA,dim(counts_cell1)[1],ncol=3)
  est_params_cell2 <- matrix(NA,dim(counts_cell2)[1],ncol=3)
  colnames(est_params_pooled) <- c("mu","phi","p0")
  colnames(est_params_cell1) <- c("mu","phi","p0")
  colnames(est_params_cell2) <- c("mu","phi","p0")
  message("ready run pb")
  pval_zinb3p <- NULL
  tstats <- NULL
  pb <- progress::progress_bar$new(
    format = "  [:bar] :percent :elapsedfull",
    total = npeak, clear = FALSE, width = 60
  )
  suppressWarnings({
    for (i in 1:npeak){
      cat(i, ',')
      counts.peak <- as.numeric(counts_pooled[i, ])
      dat.whole.peak <- data.frame(counts = counts.peak)
      
      ctrl <- pscl::zeroinfl.control(method = "L-BFGS-B")
      ctrl$reltol <- NULL
      ctrl$factr <- 1e-3/.Machine$double.eps
      #fit zero infl ng to whole sample
      m1 <- pscl::zeroinfl(counts ~ 1+offset(log(sfs_pooled)) | 1+offset(log(sfs_pooled)),
                           data = dat.whole.peak, dist = "negbin" ,control = ctrl)
      mu <- exp(m1$coefficients$count)
      if (m1$theta>150){
        theta <- 150
      } else {
        theta <- m1$theta
      }
      p0 <- plogis(m1$coefficients$zero)
      est_params_pooled[i,] <- c(mu,theta,p0)
      logL_null <- zinb.loglink(counts=counts.peak,p=p0,u=mu,k=1/theta)
      #fit zinb to cond1
      counts.peak <- as.numeric(counts_cell1[i, ])
      dat.cell1.peak <- data.frame(counts = counts.peak)
      m2 <- pscl::zeroinfl(counts ~ 1+offset(log(sfs_cell1)) | 1+offset(log(sfs_cell1)),
                           data = dat.cell1.peak, dist = "negbin" ,control = ctrl)
      mu <- exp(m2$coefficients$count)
      if (m2$theta>150){
        theta <- 150
      } else {
        theta <- m2$theta
      }
      p0 <- plogis(m2$coefficients$zero)
      est_params_cell1[i,] <- c(mu,theta,p0)
      logL_alter_1 <- zinb.loglink(counts=counts.peak,p=p0,u=mu,k=1/theta)
      #fit zinb to cond2
      counts.peak <- as.numeric(counts_cell2[i, ])
      dat.cell2.peak <- data.frame(counts = counts.peak)
      m3 <- pscl::zeroinfl(counts ~ 1 + offset(log(sfs_cell2)) | 1+offset(log(sfs_cell2)),
                           data = dat.cell2.peak, dist = "negbin" ,control = ctrl)
      mu <- exp(m3$coefficients$count)
      if (m3$theta>150){
        theta <- 150
      } else {
        theta <- m3$theta
      }
      p0 <- plogis(m3$coefficients$zero)
      est_params_cell2[i,] <- c(mu,theta,p0)
      logL_alter_2 <- zinb.loglink(counts=counts.peak,p=p0,u=mu,k=1/theta)
      logL_alter <- logL_alter_1+logL_alter_2
      
      #LRT test follow x 3 degress of freedom
      test.stats <- -2*(logL_null - logL_alter)
      pvl <- pchisq(test.stats, df=3, lower.tail = FALSE)
      pval_zinb3p <- c(pval_zinb3p,pvl)
      tstats <- c(tstats, test.stats)
      pb$tick()
    }
    
    est_params_pooled <- data.frame(est_params_pooled)
    est_params_cell1 <- data.frame(est_params_cell1)
    est_params_cell2 <- data.frame(est_params_cell2)
    
    object@params <- list(g1 = group.1.loc,
                          g2 = group.2.loc,
                          param_pooled = est_params_pooled,
                          param_g1 = est_params_cell1,
                          param_g2 = est_params_cell2)
    return(object)
  })
}
library(pscl)
zinb.loglink <- function(counts,p,u,k){
  counts <- as.numeric(counts)
  dens <- numeric(length(counts))
  for (i in 1:length(counts)){
    if (counts[i]==0){
      dens[i] <- p+(1-p)*(1/(1+k*u))^(1/k)
    } else {
      g <- gamma(counts[i]+1/k)/(gamma(1/k)*gamma(counts[i]+1))
      dens[i] <- (1-p)*g*(1/(1+k*u))^(1/k)*((k*u)/(1+k*u))^counts[i]
    }
    dens <- unlist(dens)
  }
  loglink <- sum(log(dens))
  return(loglink)
}
}

#-------------------------- Viz and misc --------------------------  

new_palette_D <- function(x, pal='stallion') {
  if (pal %in% c('stallion', 'ironMan', 'circus', 'calm', 'kelly', 'bear')) {
    o <- structure(as.character(ArchR::paletteDiscrete(x, set=pal)), names = x)
  } else if (pal %in% c('scales_hue')) {
    o <- structure(scales::hue_pal()(length(x)), names=x)
  } else {
    o <- structure(hcl.colors(n=length(x), palette = pal), names=x)
  }
  return(o)
}

library(gtools)
library(forcats)

fct_anon_simple <- function(f, prefix = "") {
  ## https://github.com/tidyverse/forcats/blob/7cade9e7e3527c1dc85ca7b48af0f09e3a8c50a1/R/anon.R
  digits <- function(x) nchar(max(x, na.rm = TRUE))
  
  zero_pad <- function(x) {
    sprintf(paste0("%0", digits(x), "d"), x)
  }
  
  # levels <- paste0(prefix, zero_pad(seq_len(nlevels(f))))
  levels <- paste0(prefix, seq_len(nlevels(f)))
  
  f <- lvls_revalue(f, levels)
  lvls_reorder(f, match(levels, levels(f)))
}

pretty_seamless_factor <- function(x, drop_empty_levels = T, prefix=''){
  ## pretty_seamless_factor(factor(c('s1', 's10', 's8', 's4'))) => 1 2 4 3
  ## pretty_seamless_factor(c('s1', 's10', 's8', 's4')) => 1 4 3 2
  if (!class(x) %in% 'factor') {
    x_lvs <- gtools::mixedsort(unique(x))
    x <- factor(x, levels = x_lvs)
  }
  if (drop_empty_levels) { x <- forcats::fct_drop(x) }
  return( fct_anon_simple(x, prefix=prefix) )
}

quick_summary_stats <- function(v, name='metric', value='value') {
  stderror <- function(x) sd(x)/sqrt(length(x))
  qt <- quantile(v, c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975))
  names(qt) <- str_remove_all(names(qt), '%')
  names(qt) <- paste0('qt', names(qt))
  enframe(
    c(c(mean = mean(v, na.rm=T),
        sd = sd(v), 
        sem = stderror(v), 
        min = min(v), 
        max = max(v)), 
      qt
    ), name = name, value = value
  )
}

#-------------------------- Seurat/Signac --------------------------  
library(pbapply); library(Matrix)
library(dbscan)
run_hdbscan_sr <- function(sr, reduction='umap', dims=1:2, n_neighbors=30) {
  emb <- Embeddings(sr, reduction)[, dims]
  tmp <- hdbscan(emb, minPts = n_neighbors, verbose = T)
  o <- tmp$cluster; o <- as.factor(o); names(o) <- Cells(sr)
  o_str <- sprintf('%s_hdbscan_%s_minPts.%s', 
                   DefaultAssay(sr), reduction, n_neighbors)
  sr <- AddMetaData(
    sr, metadata = o, 
    col.name=o_str)
  Idents(sr) <- o_str
  message(length(unique(Idents(sr))), ' clusters found.')
  return(sr)
}
FindSubCluster_hdbscan <- function(object,
                                   cluster,
                                   reduction = 'umap', dims = 1:2, 
                                   subcluster.name = "sub.cluster",
                                   n_neighbors = 20) {
  o_sub <- subset(object, idents = cluster)
  o_sub <- run_hdbscan_sr(o_sub, reduction = reduction, dims = dims, n_neighbors = n_neighbors)
  o_str <- sprintf('%s_hdbscan_%s_minPts.%s', 
                   DefaultAssay(o_sub), reduction, n_neighbors)
  idx <- match(Cells(o_sub), Cells(object))
  object@meta.data[, subcluster.name] <- as.character(Idents(object))
  object@meta.data[idx, subcluster.name] <- paste0(cluster, '_', as.character(o_sub@meta.data[, o_str]))
  object@meta.data[, subcluster.name] <- as.factor(object@meta.data[, subcluster.name])
  return(object)
  
}

sr.mean.fxn.CNA.bincounts <- function(x) {
  return(log(x = rowMeans(x = x) + 1, base = 2))
}
sr.mean.fxn.CNA.integer <- function(x) {
  return(log(x = rowMeans(x = x) + 0.001, base = 2))
}
sr.mean.fxn.CNA.logr <- rowMeans

PlotGeneExpression <- function(x, genes, label, geom='swarm', dodge.width=0) {
  ## Equivalent to copykit::plotGeneCopy
  data <- FetchData(x, vars = genes) # cells x genes
  data$cellname <- rownames(data)
  data <- pivot_longer(data, cols = setdiff(colnames(data), 'cellname'))
  grp <- x@meta.data[data$cellname, label]
  data$group <- grp
  p <- ggplot(data, aes(x=name, y=value, fill=group)) + 
    ggbeeswarm::geom_quasirandom(shape = 21, size = 2.2, stroke=0.2, dodge.width=dodge.width)
  return(p)
}

GetGeneSeqname <- function(anno, genes) {
  genes <- intersect(genes, mcols(anno)$gene_name)
  idx <- mcols(anno)$gene_name %in% genes
  o <- data.frame(seqname = as.character(seqnames(anno[idx, ])), 
                  gene_name = anno[idx, ]$gene_name)
  o <- o[!duplicated(o$gene_name), ]
  o <- o[match(genes, o$gene_name), ]
  return(structure(o$seqname, names = o$gene_name))
}


# Average multiple features into a single feature
# Ref: motivated by Seurat's MetaFeatures but the calculation is different here.
run_meta_feature <- function(object, features, cells = NULL, assay = NULL, slot = 'data') {
  if (all(is.null(cells))) {cells <- Cells(object)}
  if (is.null(assay)) {assay <- DefaultAssay(object) }
  data_mat <- GetAssayData(object = object, assay = assay, slot = slot)
  res_mat <- pbapply::pbsapply(features, function(r) {
    newmat <- data_mat[r, , drop=F]
    if (slot == 'scale.data') {
      newdata <- Matrix::colMeans(newmat)
    } else if (slot == 'counts') {
      newdata <- Matrix::colSums(newmat) ## or colMean? TBD...
    } else {
      newdata <- Matrix::colMeans(newmat)
      # rowtotals <- Matrix::rowSums(newmat)
      # newmat <- newmat / rowtotals
      # newdata <- Matrix::colMeans(newmat)
    }
    return(newdata)
  }, simplify = T)
  res_mat <- t(res_mat)
  res_mat <- as(res_mat, "sparseMatrix")    
  return(res_mat)
}


#------------------- ~~~ Copykit ~~~ -------------------  

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


calc_consensus_ident <- function(x, assay = NULL, group.by = 'subclones', fun='median') {
  data <- assay(x, assay)
  o <- aggregate.Matrix.t(data, factor(x[[group.by]]), fun = fun)
  return(o)
}

ploidy_scale <- function(ploidy_VAL, df, round = TRUE) {
  ## from Darlan
  # population segmenter returns log ratios, transforming back to ratios values
  popseg_long_exp <- as.data.frame(2^df)
  
  # correcting to ploidy
  popseg_ploidy_cor <- as.data.frame(popseg_long_exp * ploidy_VAL)
  
  if (round == TRUE) {
    # rounding to the nearest integer values
    popseg_round <- as.data.frame(round(popseg_ploidy_cor, 0))
    return(popseg_round)
  } else return(popseg_ploidy_cor)
  
}

#' For every bin, returns the inferred genomic classification class (cCNA, sCNA, uCNA)
consensus_genomic_classes <- function(consensus_int, ploidy_VAL = 2) {
  # consensus_int: features x cells
  consensus_int <- t(consensus_int) # cells x features
  ## from Darlan
  percent_clonal <- 1
  percent_extant <- 1 / nrow(consensus_int)
  
  # for every bin
  ps_percents_list <- pbapply::pbapply(consensus_int, 2, function(x) {
    return(prop.table(table(x)))
  }, simplify=F)
  
  bin_classes <- pbapply::pbsapply(ps_percents_list, function(x) {
    if (any(x == percent_extant) & all(names(x) != as.character(ploidy_VAL)) ) { return("uCNA") } 
    if (any(x == percent_clonal) & all(names(x) == as.character(ploidy_VAL)) ) { return("noCNA") } 
    if (any(x == percent_clonal)) { return("cCNA") } 
    return("sCNA")
  })
  
  # is_neutral <- apply(consensus_int, 2, function(x) all(x == ploidy_VAL))
  # bin_classes[is_neutral] <- 'noCNA'
  
  bin_classes <- unname(bin_classes)
  names(bin_classes) <- colnames(consensus_int)
  bin_classes <- factor(bin_classes, levels = paste0(c('no', 'c', 's', 'u'), 'CNA'))
  return(bin_classes)
}

copykit_find_intra_outliers <- function(x, group.by='subclones', idents=NULL, ...) {
  # x: copykit object
  # group.by: groups to search
  # idents: particular idents to search
  # ...: parameters to run copykit::findOutliers such as k and resolution
  # res: `intra_outlier`
  if (is.null(idents)) { idents <- unique(colData(x)[[group.by]]) }
  
  colData(x)$intra_outlier <- 'FALSE'
  
  for (s in idents) {
    idx <- colData(x)[[group.by]] == s
    xx <- x[, idx]
    xx <- findOutliers(scCNA=xx, ...)
    colData(x)$intra_outlier[idx] <- colData(xx)$outlier
  }
  colData(x)$intra_outlier <- as.logical( colData(x)$intra_outlier )
  message('see intra_outlier')
  return(x)
}

copykit_find_intra_aneuploid_cells <- function(x, group.by='subclones', idents=NULL, ...) {
  # x: copykit object
  # group.by: groups to search
  # idents: particular idents to search
  # ...: parameters to run copykit::findOutliers such as k and resolution
  # res: `intra_aneuploid`
  if (is.null(idents)) { idents <- unique(colData(x)[[group.by]]) }
  
  colData(x)$intra_aneuploid <- FALSE
  
  for (s in idents) {
    idx <- colData(x)[[group.by]] == s
    xx <- x[, idx]
    xx <- findAneuploidCells(scCNA=xx, ...)
    colData(x)$intra_aneuploid[idx] <- colData(xx)$is_aneuploid
  }
  message('see intra_aneuploid')
  return(x)
}



copykit_plotUmap_hi <- function(x, cells.highlight, cols.highlight='red') {
  mycol <- c(`TRUE` = cols.highlight, `FALSE` = 'grey')
  colData(x)$tohighlight <- colnames(x) %in% cells.highlight
  p <- plotUmap(x, label = 'tohighlight') + 
    scale_fill_manual(values = mycol) +
    rremove('legend')
  return(p)
}


copykit_depth_cor <- function(x, reduction = 'PCA', dims = NULL, 
                              latent_var = 'ReadsKept', do_log = F) {
  embed_cell <- reducedDim(x, reduction)
  if (is.null(dims)) { dims <- head(seq_len(ncol(embed_cell)), 10) }
  dims <- sort(as.numeric(dims))
  embed_cell <- embed_cell[, dims, drop=F]
  
  latent_var_val <- x[[latent_var]]
  if (do_log) { latent_var_val <- log1p(latent_var_val) }
  df <- as.data.frame(cor(embed_cell, latent_var_val))
  df$value <- df[, 1]
  df$comp <- dims #factor(dims, labels = rownames(df), levels = dims)
  p <- ggplot(df, aes(x=comp, y=value)) +
    geom_point() + 
    scale_y_continuous(limits = c(-1, 1)) + 
    # scale_x_continuous(n.breaks = length(dims), limits = c(1, length(dims))) +
    scale_x_continuous(breaks = dims, limits = c(1, max(dims))) + 
    labs(y='correlation', x=reduction, 
         title = sprintf('correlation between %s and %s', reduction, latent_var)) + 
    theme_light()
  # p
  return(p)  
}


# copykit_propose_merge_consensus <- function(x) { 
#   
#   data <- consensus(x)
#   clone_names <- colnames(data)
# 
# }

copykit_correct_naughty_leiden <- copykit_correct_naughty_louvain <- function(x, a_col = 'subclones', b_col = 'superclones', do_what = c('rm')) {
  ## some cells of the leiden/louvain clusters are running away to other superclones 
  ## Pick out these cells
  a <- colData(x)[[a_col]]
  b <- colData(x)[[b_col]]
  mat <- as.data.frame.matrix(table(a, b))
  print(mat)
  best_a2b <- apply(mat, 1, function(qx) {
    best_r <- which.max(qx == max(qx))
    return(colnames(mat)[best_r])
  })
  print(best_a2b)
  
  best_y <- best_a2b[as.character(colData(x)[[a_col]])]
  return(best_y != b)
}

copykit_find_sub_clusters <- function(
    x, subcluster_by, cluster, k_nb=3, saveto_col_str='newsubclusters') {
  subcells <- colnames(x)[colData(x)[, subcluster_by] %in% cluster]
  obj <- x[, subcells]
  obj[['subclones_orig']] <- obj[['subclones']]
  obj[['superclones_orig']] <- obj[['superclones']]
  obj <- findClusters(obj, k_superclones = k_nb, k_subclones = k_nb)
  
  subcells_clusters_new <- obj[['subclones']]
  subcells_clusters_orig <- obj[['subclones_orig']]
  subcells_clusters_new <- structure(
    paste0(subcells_clusters_orig, '.', subcells_clusters_new),
    names = subcells)
  print(table(subcells_clusters_new), useNA='ifany')
  x[[saveto_col_str]] <- as.character(x[['subclones']])
  i <- match(subcells, colnames(x))
  colData(x)[i, saveto_col_str] <- subcells_clusters_new[subcells]
  x[[saveto_col_str]] <- factor(x[[saveto_col_str]], levels = gtools::mixedsort(unique(x[[saveto_col_str]])))
  return(x)  
}




pal_copykit_ratio_clip2 <- circlize::colorRamp2(
  seq(from=-2, to=2, length.out=11),
  hcl.colors(n=11, palette = 'RdBu', rev = T))
pal_copykit_ratio_clip1 <- circlize::colorRamp2(
  seq(from=-1, to=1, length.out=11), 
  hcl.colors(n=11, palette = 'RdBu', rev = T))
# pal_copykit_int_from_ratio_clip1 <- circlize::colorRamp2(
#   c(0, 2* 2^seq(from=-1, to=1, length.out=11)), 
#   c('black', hcl.colors(n=11, palette = 'RdBu', rev = T)))
# pal_copykit_int_from_ratio_clip1 <- structure(
#   pal_copykit_int_from_ratio_clip1(c(0, 2* 2^seq(from=-1, to=1, length.out=11))), 
#   names=as.character(c(0, 2* 2^seq(from=-1, to=1, length.out=11)))
# )
pal_copykit_int_from_ratio_clip1 <- structure(
  pals::ocean.balance(5), 
  names = 0:4
)
# pal_copykit_int_clip_0_2_6 <- circlize::colorRamp2(
#   0:6,
#   hcl.colors(n=9, palette = 'RdBu', rev = T)[c(1,3,5:9)])
pal_copykit_int_clip_0_2_6 <- structure(
  hcl.colors(n=9, palette = 'RdBu', rev = T)[c(1,3,5:9)],
  names = 0:6
)
pal_heatmap_auto_c <- function(x, pal='Plasma', rev=F) {
  circlize::colorRamp2(
    seq(from=min(x), to=max(x), length.out=9),
    hcl.colors(n=9, palette =pal, rev = rev))
}

#------ calcualte MPD ------
quick_scTree <- function(df,
                         method = "nj",
                         metric = "manhattan",
                         assay = "ratio",
                         add_fake_dipoid = TRUE, 
                         dipoid_val = 1, 
                         n_threads = parallel::detectCores() / 4) {
  # df: bins x cells
  
  # Originally written by Kris Wang
  # Add the option to include the fake diploid
  
  # cores check
  if (n_threads < 1) {n_threads <- 1}

  seg_data <- df
  
  if (add_fake_dipoid) {
    seg_data_diploid <- rep(dipoid_val, times = nrow(df))
    seg_data$FAKE_DIPLOID <- seg_data_diploid
  } 
  
  if (assay == "integer") {
    ## with integers
    message("Using integer data...")
    seg_data <- t(seg_data) %>% as.matrix()
    ## recommend using hamming distance for integer profiles
    distMat <- as.matrix(parallelDist::parDist(seg_data, method= "hamming", diag=T, upper=T,n_threads=n_threads))
    
    if (metric != "hamming") {
      stop("Recommend only using hamming distance for integer profiles")
    }
    
    
  } else {
    # with ratios
    # message("Using ratio data...")
    seg_data <- t(seg_data) %>% as.data.frame()
    # calculating distance matrix
    # message("Calculating distance matrix")
    distMat <- amap::Dist(seg_data,
                          method = metric,
                          nbproc = n_threads)
  }
  # ordering cells
  if (method %in% c("nj", "me")) {
    if (method == "nj") {
      # message("Creating neighbor-joining tree.")
      tree <- ape::nj(distMat)
    }
    
    if (method == "me") {
      # message("Creating minimum evolution tree.")
      tree <- ape::fastme.bal(distMat)
    }
    
  } else {
    stop("Currently only nj and me trees are supported.")
  }
  
  return(tree)
}

calc_mpd_phylo_tree <- function(tree) {
  # tree: a phylo tree
  # originally developed by Kris Wang
  n<-length(tree$tip.label)
  ## removing end node
  tree$edge.length[sapply(1:n,function(x,y) which(y==x),y=tree$edge[,2])] <- 0
  return(mean(ape::cophenetic.phylo(tree)))  
}

calc_mean_distance_to_ref <- function(tree, ref_tip='FAKE_DIPLOID') {
  return(mean(ape::cophenetic.phylo(tree)[, ref_tip]))
}

#------------------ ~~~ Stats ~~~ --------------------
calc_histogram <- function(array, nbins = 100) {
  # Ensure input is a numeric vector
  array <- as.numeric(array)
  
  # Compute histogram
  hist_data <- hist(array, breaks = nbins, plot = FALSE)
  
  # Compute bin centers
  bin_centers <- (head(hist_data$breaks, -1) + tail(hist_data$breaks, -1)) / 2
  
  return(list(hist = hist_data$counts, bin_centers = bin_centers))
}


threshold_otsu <- function(array, nbins = 100, min_value = 100) {
  ## Ref: https://github.com/aertslab/scATAC-seq_benchmark/blob/36dd41912b56460fc546daca12b797d11bebd713/0_resources/scripts/qc_plots_public.py#L34
  # Filter the array based on min_value
  cat("# Filter the array based on min_value")
  array <- array[array >= min_value]
  cat("# Get histogram and bin centers")
  # Get histogram and bin centers
  hist_info <- calc_histogram(array, nbins)
  hist <- as.numeric(hist_info$hist)
  bin_centers <- hist_info$bin_centers
  
  # Convert to probabilities
  weight1 <- cumsum(hist)
  weight2 <- rev(cumsum(rev(hist)))
  
  # Avoid division by zero
  safe_divide <- function(a, b) {
    out <- rep(0, length(a))
    nonzero <- b != 0
    out[nonzero] <- a[nonzero] / b[nonzero]
    return(out)
  }
  
  mean1 <- safe_divide(cumsum(hist * bin_centers), weight1)
  mean2 <- rev(safe_divide(cumsum(rev(hist * bin_centers)), rev(weight2)))
  
  # Compute inter-class variance
  variance12 <- weight1[-length(weight1)] * weight2[-1] * 
    (mean1[-length(mean1)] - mean2[-1])^2
  
  # Find the index of the maximum variance
  idx <- which.max(variance12)
  
  # Return the threshold
  threshold <- bin_centers[idx]
  return(threshold)
}

if (F) {
  set.seed(42)
  data <- c(rnorm(1000, 400, sd=200), rpois(1000, 1000))
  hist(data)
  abline(v=threshold_otsu(data, nbins = 100), col='red')
}


# #------ viz CODA ggtree ------
# library(ggnewscale)
# ## Copied from https://dmnfarrell.github.io/r/ggtree-heatmaps
# gettreedata <- function(tree, meta){
#   #get treedata object
#   d<-meta[row.names(meta) %in% tree$tip.label,]
#   d$label <- row.names(d)
#   y <- full_join(as_tibble(tree), d, by='label')
#   y <- as.treedata(y)
#   return(y)
# }
# 
# get_color_mapping <- function(data, col, cmap){
#   labels <- (data[[col]])   
#   names <- levels(as.factor(labels))
#   n <- length(names)
#   if (n<10){      
#     colors <- suppressWarnings(c(brewer.pal(n, cmap)))[1:n]
#   }
#   else {
#     colors <- colorRampPalette(brewer.pal(8, cmap))(n)
#   }
#   names(colors) = names
#   return (colors)
# }
# 
# ggplottree <- function(tree, meta, cols=NULL, cmaps=NULL, layout="rectangular",
#                        offset=10, tiplabel=FALSE, tipsize=3) {
#   
#   y <- gettreedata(tree, meta)
#   p <- ggtree(y, layout=layout)   
#   if (is.null(cols)){
#     return (p)
#   }
#   
#   col <- cols[1]
#   cmap <- cmaps[1]
#   df<-meta[tree$tip.label,][col]
#   colors <- get_color_mapping(df, col, cmap)
#   
#   #tip formatting    
#   p1 <- p + new_scale_fill() +    
#     geom_tippoint(mapping=aes(fill=.data[[col]]),size=tipsize,shape=21,stroke=0) +
#     scale_fill_manual(values=colors, na.value="white")
#   
#   p2 <- p1
#   if (length(cols)>1){
#     for (i in 2:length(cols)){
#       col <- cols[i]
#       cmap <- cmaps[i]
#       df <- meta[tree$tip.label,][col]
#       type <- class(df[col,])            
#       p2 <- p2 + new_scale_fill()
#       p2 <- gheatmap(p2, df, offset=i*offset, width=.08,
#                      colnames_angle=0, colnames_offset_y = .05)  
#       #deal with continuous values
#       if (type == 'numeric'){               
#         p2 <- p2 + scale_color_brewer(type="div", palette=cmap)
#       }
#       else {
#         colors <- get_color_mapping(df, col, cmap)
#         p2 <- p2 + scale_fill_manual(values=colors, name=col)
#       }          
#     }
#   }
#   
#   p2 <- p2 + theme_tree2(legend.text = element_text(size=20), legend.key.size = unit(1, 'cm'),
#                          legend.position="left", plot.title = element_text(size=40))     
#   guides(color = guide_legend(override.aes = list(size=10)))
#   
#   return(p2)
# }
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


cat('!!! source util.R -- ')
cat('The following functions are created:\n')
cat('- bgzip_tabix_fragments\n')
cat('- get_density2 \n')
cat('- run_hdbscan_sr\n')
cat('- new_palette_D\n')
cat('- pretty_seamless_factor\n')
cat('- copykit_find_outliers\n')
cat('The following color functions: \n')
cat('- pal_copykit_int_from_ratio_clip1\n')
cat('- pal_copykit_ratio_clip2\n')
cat('- pal_copykit_ratio_clip1\n')