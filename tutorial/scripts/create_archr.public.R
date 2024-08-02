#/volumes/USR1/yyan/anaconda3/envs/r_env/bin/Rscript
suppressPackageStartupMessages(library(magrittr))
"
#--------------------------
# Goal: Given a fragment.tsv.gz file, create an ArchR arrow file and export to a Signac-ready object
# 1. Arrow file
# 2. Cell QC plot
# 3. Good cell names
# 4. peak.bed
# History:
# v2: set maxnFrags to 1e6 because of the single-end fragments
#--------------------------
# Author: Yun Yan (yun.yan@uth.tmc.edu)
#--------------------------
" -> doc_help
timestamp()
suppressPackageStartupMessages({
    
    library(tidyverse); library(ggplot2); library(ggpubr); library(patchwork)
    theme_set(theme_pubr(base_size = 20, legend='right'))
    library(cli); library(tictoc); library(glue); library(scales); library(tools); library(ruok)
    library(GenomeInfoDb); library(GenomicRanges)
    library(ArchR)
    library(parallel); library(parallelly)
    library(ggrastr)
    my_scatter_theme <- theme_pubr(base_size = 12, legend = 'right') %+replace% theme(axis.text=element_text(size=rel(.2)), axis.title=element_text(size=rel(.3)), axis.ticks=element_line(size = rel(.3)), panel.border = element_rect(fill = NA, size=rel(.5)), axis.line = element_blank())  
    my_scatter_themevoid <- theme_pubr(base_size = 12, legend = 'bottom') %+replace% theme(axis.text=element_blank(), axis.title=element_blank(), axis.ticks=element_blank(), panel.border = element_blank(), axis.line = element_blank())
    
})
options <- commandArgs(trailingOnly = TRUE)
## Inputs

if (length(options)==2) {
    run_name       <- options[[1]]
    fpath_fragment <- options[[2]]
}
message(run_name)
message(fpath_fragment)
stopifnot(file.exists(fpath_fragment))
#-------------------------- Outputs --------------------------  
## Outputs
dir_arrow <- file.path(
    "/volumes/USR1/yyan/project/coda/data",
    run_name, "scATAC", 'archr')
dir_archr <- file.path(
    "/volumes/USR1/yyan/project/coda/rds",
    run_name, "archr_project")
dir_tosignac <- file.path(
    "/volumes/USR1/yyan/project/coda/rds",
    run_name, "tosignac")
fs::dir_create(dir_arrow)
fs::dir_create(dir_archr)
fs::dir_create(dir_tosignac)
#--------------------------
# Parameters (mostly no need to change)
#--------------------------
tss_score_min_init <- 0
fragment_num_min_init <- 1
fragment_num_max_init <- 1e7

tss_score_min_use <- 8#8 # Read: https://www.encodeproject.org/atac-seq/
fragment_num_min_use <- 1000
fragment_num_max_use <- 1e5
#--------------------------
# Envs (no need to change)
#--------------------------
n_threads <- 20
macs2_cmd <- '/volumes/USR1/yyan/anaconda3/envs/macs2/bin/macs2'
addArchRGenome(species_genome)
addArchRThreads(threads = n_threads)


#--------------------------
# Helpers
#--------------------------
suppressPackageStartupMessages(library(MASS))
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
#-------------------------- START --------------------------  
fs::dir_create(dir_archr)
fs::dir_create(dir_arrow)
#-------------------------- Create Arrow files --------------------------  
# arrow files can always be re-used.
if (T) {
    setwd(dir_arrow)
    ArrowFiles <- createArrowFiles(
        inputFiles = fpath_fragment,
        sampleNames = run_name,
        filterTSS = tss_score_min_init,
        minFrags = fragment_num_min_init, 
        maxFrags = fragment_num_max_init,
        addTileMat = TRUE,
        TileMatParams = list(binarize = do_binarization),
        # removeFilteredCells = FALSE, 
        nChunk = 20, # higher number reduces memeory usage.
        force = TRUE,
        addGeneScoreMat = TRUE)
    stopifnot(file.exists((ArrowFiles)))
    getwd()
    
    doubScores <- addDoubletScores(
        input = ArrowFiles,
        k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
        knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
        LSIMethod = 1
    )
}
#-------------------------- Prepare Analysis --------------------------  
archr <- ArchRProject(
    ArrowFiles = file.path(dir_arrow, sprintf('%s.arrow', run_name)), 
    outputDirectory = dir_archr,
    copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)

setwd(dir_archr)

#-------------------------- START --------------------------  
# archr <- read_rds(file.path(dir_archr, 'Save-ArchR-Project.rds'))
str(getCellNames(archr))
colnames(getCellColData(archr))
min(archr$TSSEnrichment)
#-------------------------- QC cells -------------------------- 
if (T) {
    df <- getCellColData(archr, select = c("log10(nFrags)", "TSSEnrichment", 
                                           "nFrags", "Sample"))
    df <- as.data.frame(df)
    
    is_pass <- with(df, TSSEnrichment >= tss_score_min_use & nFrags >= fragment_num_min_use & nFrags <= fragment_num_max_use)
    df$log10nFrags <- log10(df$nFrags)
    write.csv(df, file.path(dir_archr, 'overviewQC.csv'))
    # df <- read.csv(file.path(dir_archr, 'overviewQC.csv'))
    # df <- df[df$TSSEnrichment>3,]
    tmp <- get_density2(df$log10nFrags, df$TSSEnrichment, 
                        n=100)
    df$density <- tmp$density
    df$density <- df$density / max(df$density) ### !!! Normalized to max to unify visualization
    
    p <- ggplot(df, aes(x=log10nFrags, y=TSSEnrichment, color=density)) +
        geom_point_rast(size=1) +
        border() +
        scale_color_viridis_c(option = 'E')+
        labs(x="Log10 Unique Fragments", y= "TSS Enrichment",
             title=sprintf('%s (%d/%d cells)', 
                           unique(df$Sample), sum(is_pass, na.rm = T), nrow(df))) +
        geom_hline(yintercept = tss_score_min_use, lty = "dashed") + 
        geom_vline(xintercept = log10(fragment_num_min_use), lty = "dashed") +
        theme(aspect.ratio = 1)
    ggsave2(file.path(dir_archr, 'overviewQC'), p, width = 6.5, height = 6)
    
    archr <- ArchR::subsetCells(
        archr, 
        cellNames = getCellNames(archr)[is_pass])
    
    str(getCellNames(archr))
}
n_cells = length(getCellNames(archr))

archr <- saveArchRProject(ArchRProj = archr, 
                          outputDirectory = dir_archr, 
                          load = TRUE)

#-------------------------- Doublets --------------------------  
if (T) {
    archr <- addIterativeLSI(
        ArchRProj = archr,
        useMatrix = "TileMatrix", 
        name = "IterativeLSI", 
        force=T,
        iterations = 2, 
        #See Seurat::FindClusters to set clusterParams
        clusterParams = list( 
            resolution = c(0.2), 
            sampleCells = 10e3, 
            n.start = 10), 
        varFeatures = 25000, 
        scaleTo = 1e4,
        binarize = do_binarization, 
        dimsToUse = 1:20)
    archr <- addTSNE(archr, force=TRUE)
    archr <- addDoubletScores(archr, force=TRUE)
    archr_tmp <- filterDoublets(archr, filterRatio = doublet_filter_ratio)
    
    df_meta <- getEmbedding(archr, 'TSNE')
    colnames(df_meta) <- c(sprintf('Dim_%d', 1:2))
    cnames_pass_doublet <- getCellNames(archr_tmp); str(cnames_pass_doublet)
    df_meta$doublet_type <- ifelse(rownames(df_meta) %in% cnames_pass_doublet, 'singlet', 'doublet')
    try(dev.off())
    write.csv(df_meta, file.path(dir_archr, 'doublet_detection.csv'))
    df_meta <- read.csv(file.path(dir_archr, 'doublet_detection.csv'))
    p1 <- ggplot(df_meta, aes_string(x='Dim_1', y='Dim_2')) +
        geom_point_rast(aes(color = doublet_type), size=1) +
        scale_color_discrete(labels=pretty_table2str(table(df_meta$doublet_type)))+
        coord_equal()+ guides(color = guide_legend(override.aes=list(size=3))) +
        labs(color='')
    
    ggsave2(file.path(dir_archr, 'doublet_detection'),
            p1 + my_scatter_themevoid + border(), width = 6, height = 6.5)
    
    archr <- archr_tmp; rm(archr_tmp)
    
}


archr <- saveArchRProject(ArchRProj = archr, 
                          outputDirectory = dir_archr, 
                          load = TRUE)

cat('[DONE]')
