#--------------------------
# Goal: Create ready-signac object purely by Signac
# The only input is the fragment file
## Written by: Yun Yan (https://github.com/Puriney)\n\n
#--------------------------
library(future)
options(future.globals.maxSize = 50 * 1024 ^ 3) # for 50 Gb RAM
plan("multisession", workers = 4)
library(Signac); library(readr); library(tidyverse); library(fs); library(ruok)
library(Seurat); library(ggrastr); library(ggplot2); library(ggpubr)
library(cli); library(tictoc)
theme_set(theme_pubr(base_family = 'Helvetica', legend = 'right'))
suppressPackageStartupMessages(library(MASS))
get_density2 <- function(x = NULL, y = NULL, n = 100, sample = NULL, densityMax = 0.95, fillna=0){
    #modified from http://slowkow.com/notes/ggplot2-color-by-density/
    #https://rdrr.io/github/GreenleafLab/ArchR/src/R/GgplotUtils.R
    x[is.na(x)] <- fillna
    y[is.na(y)] <- fillna
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
#
#
#------ Param (mostly no need to change) ------
options <- commandArgs(trailingOnly = TRUE)
## suggested for the common pair-end seq
tss_score_min_init <- 0
fragment_num_min_init <- 10
tss_score_min_use <- 8 # 6-8 for hg19 Read: https://www.encodeproject.org/atac-seq/
# tss_score_min_use <- 2 # Signac is different from ArchR
fragment_num_min_use <- 1e3
fragment_num_max_use <- 1e5
## suggested for the single-end seq
tss_score_min_init <- 0
fragment_num_min_init <- 10
tss_score_min_use <- 3
fragment_num_min_use <- 1e4
fragment_num_max_use <- 1e6

do_binarization <- FALSE
qc_by_archr <- TRUE
#------ inputs ------

if (length(options)>0) {
    cat('Reading user parameters:\n')
    run_name             <- options[[1]]
    fpath_fragment       <- options[[2]]
    tss_score_min_use    <- as.numeric(options[[3]]) # Signac is different from ArchR
    fragment_num_min_use <- as.numeric(options[[4]])
    fragment_num_max_use <- as.numeric(options[[5]])
    do_binarization      <- as.numeric(options[[6]]) == 1 # 1:do binarize
    qc_by_archr          <- as.numeric(options[[7]]) == 1 # 1:do qc by archr
}

message(run_name)
message(fpath_fragment)
message(sprintf('TSS min = %s', tss_score_min_use))
message(sprintf('nFrags min = %s', fragment_num_min_use))
message(sprintf('nFrags max = %s', fragment_num_max_use))
message(sprintf('do_binarization = %s', do_binarization))

stopifnot(file.exists(fpath_fragment))

#-------------------------- Outputs --------------------------  
## Outputs
dir_rawsignac <- file.path(
    "/volumes/USR1/yyan/project/coda/data",
    run_name, "scATAC", 'signac')
dir_tosignac <- file.path(
    "/volumes/USR1/yyan/project/coda/rds",
    run_name, "signac")
fs::dir_create(dir_rawsignac)
fs::dir_create(dir_tosignac)

#------ Env and macs2 path Load txdb ------
macs2_cmd <- '/volumes/USR1/yyan/anaconda3/envs/macs2/bin/macs2'
if (species_genome == 'mm10') {
    library(BSgenome.Mmusculus.UCSC.mm10)
    bsgenome <- BSgenome.Mmusculus.UCSC.mm10
    library(TxDb.Mmusculus.UCSC.mm10.knownGene)
    txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
    library(EnsDb.Mmusculus.v79)
    endb <- EnsDb.Mmusculus.v79
    library(org.Mm.eg.db)
    orggenome <- org.Mm.eg.db
    species_full_name <- "Mus musculus"
} 
if (species_genome == 'hg19') {
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    library(EnsDb.Hsapiens.v75)
    endb <- EnsDb.Hsapiens.v75
    species_full_name = "Homo sapiens"
}
if (species_genome == 'hg38') {
    # library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    # txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    library(EnsDb.Hsapiens.v86)
    endb <- EnsDb.Hsapiens.v86
    species_full_name = "Homo sapiens"
}



#-------------------------- part-1 --------------------------    

# fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
# cells <- colnames(x = atac_small)
# frags <- CreateFragmentObject(path = fpath, cells = cells, validate.fragments = T)
# frags <- CreateFragmentObject(path = fpath, validate.fragments = F)

#------ Raw output ------
if (file.exists(file = file.path(dir_rawsignac, 'signac.rds'))) {
  message('load existing objects')
    sr3 <- read_rds(file = file.path(dir_rawsignac, 'signac.rds'))
} else {
    fragments <- CreateFragmentObject(
        path = fpath_fragment, validate.fragments = F
    )
    fragments
    tmpdir <- file.path(dir_rawsignac, 'tmp_callpeaks'); fs::dir_create(tmpdir)
    peaks <- CallPeaks(
        object = fragments,
        macs2.path = macs2_cmd, 
        outdir = tmpdir
    )
    fs::dir_delete(tmpdir)
    length(peaks)
    counts <- FeatureMatrix(
        fragments = fragments,
        features = peaks,
        sep = c("-", "-")
    )
    print(dim(counts))
    # nfrags
    # TSS
    # FRIP
    chrom_assay <- CreateChromatinAssay(
        counts = counts,
        sep = c("-", "-"),
        genome = species_genome,
        fragments = fpath_fragment, 
        min.cells = 1,
        min.features = 1
    )
    
    total_fragments <- CountFragments(fragments = fpath_fragment)
    total_fragments$nFrags <- total_fragments$frequency_count
    total_fragments$frequency_count <- NULL
    colnames(total_fragments)
    rownames(total_fragments) <- total_fragments$CB
    stopifnot(all(Cells(chrom_assay) %in% rownames(total_fragments)))
    sr3 <- CreateSeuratObject(
        counts = chrom_assay,
        assay = "peaks",
        meta.data = total_fragments[Cells(chrom_assay), ]
    )
    annotations <- GetGRangesFromEnsDb(ensdb = endb)
    seqlevelsStyle(annotations) <- 'UCSC'
    Annotation(sr3) <- annotations
    sr3 <- NucleosomeSignal(object = sr3)
    sr3 <- TSSEnrichment(object = sr3, fast = FALSE)
    
    sr3 <- FRiP(
        object = sr3,
        assay = 'peaks',
        total.fragments = 'nFrags'
    )
    sr3$log10nFrags <- log10(sr3$nFrags)
    colnames(sr3@meta.data)
    # export 
    # - raw object
    write_rds(sr3, file = file.path(dir_rawsignac, 'signac.rds'))
    # - matrix.rds
    # - matrix.mtx
    # - barcodes.txt
    # - features.txt
    # - cellmeta.csv
}

#------ Attach physical wells ------
source('/volumes/USR1/yyan/project/coda/rsrc/utils.wafargen_physical.R')
cnames <- Cells(sr3) ; str(cnames)
if(!all(cnames %in% lib_wafar$well_name))  {
  warning('[Wrong] Not all cell names are in the dictionary. Maybe wrong cell names.\n')
} else {
  message('attaching the available physical well locations')  
  idx <- match(cnames, lib_wafar$well_name)
  sr3 <- AddMetaData(sr3, data.frame(lib_wafar[idx, ], row.names = cnames))
  write_rds(sr3, file = file.path(dir_rawsignac, 'signac.rds'))
}
sr3 <- read_rds( file.path(dir_rawsignac, 'signac.rds') )
#------ Add dispense information ------
run_srch_key <- str_remove(run_name, '_low') %>% str_remove(., '_atac') %>% str_remove(., '_chise'); message(run_srch_key)
wafar_use <- read_tsv(Sys.glob(sprintf('/volumes/USR1/yyan/project/coda/wafar/%s/*.TXT', run_srch_key))); dim(wafar_use)
if (!purrr::is_empty(wafar_use)) {
  wafar_wid_use <- paste0(wafar_use$Row, '_', wafar_use$Col)
  
  wid_use <- lib_wafar$well_id[match(cnames, lib_wafar$well_name)]; str(wid_use)

  str(setdiff(sr3$well_id, wafar_wid_use))
  
  sr3$is_dispensed <- c(sr3$well_id %in% wafar_wid_use); print(table(sr3$is_dispensed))
  
  if (sum(!sr3$is_dispensed) >0){
    print(table(is_dispensed=sr3$is_dispensed))
    cat('Some cells are false-positive!!.\n')
  }
} else {
  message('assuming all cells are dispensed. ')
  sr3$is_dispensed <- TRUE
}

write_rds(sr3, file = file.path(dir_rawsignac, 'signac.rds'))
#------ QC cells ------
df <- sr3@meta.data
f_archr_qc <- file.path(dirname(dir_tosignac), 'archr_project', 'overviewQC.csv')
qc_method <- 'signac'
try(df$TSS.enrichment.signac <- df$TSS.enrichment)
try(sr3$TSS.enrichment.signac <- sr3$TSS.enrichment)
if ( qc_by_archr & file.exists(f_archr_qc) ) {
  qc_archr <- read.csv(f_archr_qc, row.names = 1, stringsAsFactors = F)
  qc_archr$cellname <- stringr::str_remove(rownames(qc_archr), pattern = paste0(run_name, '#') )
  str(df$CB); str(qc_archr$cellname); 
  str(intersect(df$CB, qc_archr$cellname))
  idx <- match(df$CB, qc_archr$cellname)
  df$TSS.enrichment.archr <- qc_archr$TSSEnrichment[idx]
  df$TSS.enrichment <- df$TSS.enrichment.archr ## override the signac TSS
  sr3$TSS.enrichment <- df$TSS.enrichment      ## override the signac TSS
  # df <- df[!is.na(df$TSS.enrichment.use), ]
  qc_method <- 'archr'  
  
} else {
  # do nothing
}
is_pass <- with(df, TSS.enrichment >= tss_score_min_use & nFrags >= fragment_num_min_use & nFrags <= fragment_num_max_use); print(table(is_pass))
names(is_pass) <- df$CB
tmp <- get_density2(df$log10nFrags, df$TSS.enrichment, n=100)
df$density <- tmp$density
df$density <- df$density / max(df$density) ### !!! Normalized to max to unify visualization

p <- ggplot(df, aes(x=log10nFrags, y=TSS.enrichment, color=density)) +
  ggrastr::geom_point_rast(size=1) +
  # border() +
  scale_color_viridis_c(option = 'E')+
  labs(x="Log10 Unique Fragments", 
       y= sprintf("TSS Enrichment (%s)", qc_method), 
       title=sprintf('%s (%d/%d cells)', 
                     run_name, sum(is_pass, na.rm = T), nrow(df))) +
  geom_hline(yintercept = tss_score_min_use, lty = "dashed") + 
  geom_vline(xintercept = c(log10(fragment_num_min_use), log10(fragment_num_max_use)),
             lty = "dashed") +
  theme(aspect.ratio = 1)
print(p)
ggsave2(file.path(dir_rawsignac, 'overviewQC'), p, width = 6.5, height = 6)

#------ Export QC for all possible cells ------
df$is_pass <- ifelse(is_pass, 'good', 'bad')
sr3$is_pass <- ifelse(is_pass, 'good', 'bad')
table(is_pass=df$is_pass, is_dispensed=df$is_dispensed, useNA='ifany')
print(table(is_pass=sr3$is_pass, is_dispensed=sr3$is_dispensed, useNA='ifany'))

write_rds(sr3, file = file.path(dir_rawsignac, 'signac.rds'))
write.csv(df, file.path(dir_rawsignac, 'overviewQC.csv'))

#-------------------------- part-2 --------------------------    

#------ Select good & dispensed cells ------
# examine fragment content
# sr3 <- read_rds( file.path(dir_rawsignac, 'signac.rds') )
# p <- VennDiagram::venn.diagram(
#   list(ATAC=sr3$is_pass, Dispense=ifelse(sr3$is_dispensed, 'good', 'bad')),
#   col = c('blue', 'gold'), cat.col = c('blue', 'gold'),
#   cat.default.pos = "text",
#   filename = NULL, disable.logging=T)

wid_use <- df[which(df$is_pass == 'good'), 'well_id']
wafar_wid_use <- df[which(df$is_dispensed), 'well_id']
p <- ggvenn::ggvenn(list(ATAC=wid_use, Dispense=wafar_wid_use))
pdf(file.path(dir_rawsignac, 'VennDiagram.well_id.pdf'), width = 5, height = 5, useDingbats = F)
print(p)
dev.off()
wafar_df <- create_wafar_loc_dataframe(wid_use, wafar_wid_use)
pa <- plot_wafar_capture(wafar_df) + labs(title=run_name); print(pa)
ggsave(file.path(dir_rawsignac, 'wafargen_physical_dispense.pdf'), pa, width = 7, height = 7, useDingbats = F)

sr3$is_pass <- ifelse(is_pass, 'good', 'bad')
sr3 <- subset(sr3, subset = is_pass == 'good')
if (sum(!sr3$is_dispensed)>0){
  sr3 <- subset(sr3, subset = is_dispensed == TRUE)
  print(sr3)
  sr3$is_dispensed <- NULL # because since now all cells are dispensed
}
sr3$is_pass <- NULL # because since now all cells are good


old_frag_list <- Fragments(sr3)
new_frag_list <- vector(mode='list', length = length(old_frag_list))
## actually it is only 1 frag file
for (i in seq_along(old_frag_list)) {
    old_frag_fpath <- GetFragmentData(old_frag_list[[i]], slot = 'path')
    new_frag_fpath <- file.path(dir_tosignac, 'fragments.tsv.gz')
    if (file.exists(new_frag_fpath)) {fs::file_delete(new_frag_fpath)}
    FilterCells(old_frag_fpath, cells=Cells(sr3), outfile = new_frag_fpath)
    new_frag_list[[i]] <- new_frag_fpath
}
fragments <- CreateFragmentObject(
    path = new_frag_list[[1]], cells = Cells(sr3), validate.fragments = T
)
Fragments(sr3) <- NULL
Fragments(sr3) <- fragments

#------ export ------
# - object.rds
write_rds(sr3, file.path(dir_tosignac, 'ready.signac_default.rds'))
# - matrix.rds
# - matrix.mtx
# - barcodes.txt
# - features.txt
# - cellmeta.csv

#-------------------------- part-3 --------------------------    
#------ Analysis-ready signac object ------
# sr3 <- read_rds(file.path(dir_tosignac, 'ready.signac_default.rds'))
if (!do_binarization){
    cli_rule('Analysis-ready signac without binarization')
} else {
    cli_rule('Analysis-ready signac with binarization')
    sr3 <- BinarizeCounts(sr3)
}
sr3 <- RunTFIDF(sr3)
sr3 <- FindTopFeatures(sr3, min.cutoff = 'q5')
sr3 <- RunSVD(sr3, assay='peaks', 
              n = ifelse(ncol(sr3) < 100, ncol(sr3)-1, 100), 
              reduction.key = 'LSI_', reduction.name = 'lsi')

p <- DepthCor(sr3)+geom_hline(yintercept = c(0.75, -0.75) , lty='dashed', color='red')
ggsave2(file.path(dir_tosignac, 'decide_poor_lsi'), p, width = 5, height = 5)
if ( 'iLSI' %in% names(sr3@reductions) ){
    n_ilsi <- ncol(mat_ilsi)
    sr3 <- RunUMAP(sr3, reduction = 'iLSI', dims = 1:n_ilsi, verbose = T)
    sr3 <- RunTSNE(sr3, reduction = 'iLSI', dims = 1:n_ilsi, verbose = T)
    sr3 <- FindNeighbors(sr3, reduction = "iLSI", dims = 1:n_ilsi, verbose = T)
} else {
  # FeatureScatter(sr3, 'nFrags', 'LSI_1') ## LSI-1 is the proxy of nFragments
  bad_lsi <- p$data$Component[abs(p$data$counts) > 0.75]
  dim_lsi <- head(setdiff(c(1:100), bad_lsi), 50)
  sr3 <- RunUMAP(sr3, reduction = 'lsi', dims = dim_lsi, verbose = F, 
                 n.neighbors=ifelse(ncol(sr3) < 30, 5, 20))
  sr3 <- RunTSNE(sr3, reduction = 'lsi', dims = dim_lsi, verbose = F)
  sr3 <- FindNeighbors(sr3, reduction = "lsi", dims = dim_lsi, verbose = F, 
                       k.param=ifelse(ncol(sr3) < 30, 5, 20))
}



for (snn_res in rev(seq(0.1, 0.6, by=0.1))){
    sr3 <- FindClusters(object = sr3, resolution = snn_res, 
                        verbose = F, n.start = 10, algorithm = 1)
    snn_str <- sprintf('peaks_snn_res.%s', snn_res)
}

library(clustree)
p <- clustree::clustree(sr3, prefix='peaks_snn_res.')
p
ggsave2(file.path(dir_tosignac, 'clustree'), p, width = 10, height = 10)


if (!do_binarization){
    write_rds(sr3, file.path(dir_tosignac, 'ready.signac_default.rds'))
} else {
    write_rds(sr3, file.path(dir_tosignac, 'ready.signac_binary.rds'))
}

# p1 <- UMAPPlot(sr3, group.by='peaks_snn_res.0.3') + coord_fixed()
p1 <- UMAPPlot(sr3) + coord_fixed()
pdf(file.path(dir_tosignac, 'overview_umap_snn.pdf'), width = 7.5, useDingbats = F)
print(p1)
dev.off()

cat('[DONE]')


sr3 <- read_rds( file.path(dir_tosignac, 'ready.signac_default.rds') )
plan('sequential')
