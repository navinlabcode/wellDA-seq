source('/volumes/USR1/yyan/project/coda/rsrc/r_libs.R')
source('/volumes/USR1/yyan/project/coda/rsrc/utils.R')
if (T) {
  library(cicero); library(monocle3)
  library(data.table);
  setDTthreads(threads = 4)
}

cat('load input...\n')
# sample_name <- 'DCIS66T_chip2'
# sample_name <- 'DCIS22T'
# sample_name <- 'DCIS51T'
# sample_name <- 'DCIS68T'
sample_name <- 'DCIS35T'

args  <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) { 
  sample_name = args[[1]] 
}

dir_coda <- file.path(
  '/volumes/USR1/yyan/project/coda/rds_coda_ready', 
  sample_name, 'aneuploid_epi')
dir_cicero <- file.path(
  '/volumes/USR1/yyan/project/coda/rds_coda_ready', 
  sample_name, 'aneuploid_epi', 'cicero')


load_coda(dir_coda)
cicero <- read_rds(file.path(dir_cicero, 'obja.cicero.rds'))
monocle <- read_rds(file.path(dir_cicero, 'obja.monocle3.rds'))
conns <- read_rds(file.path(dir_cicero, 'cicero_conns.df.rds'))
## REF: https://cole-trapnell-lab.github.io/cicero-release/docs_m3/#cicero-gene-activity-scores

fpath_hg19_gtf <- 'https://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz'
fpath_hg19_gtf <- '/volumes/USR1/yyan/public/hg19/Homo_sapiens.GRCh37.87.gtf.gz'
fpath_hg19_gtf <- '/volumes/USR1/yyan/public/hg19/Homo_sapiens.GRCh37.87.rtracklayer.dataframe.rds'
if (F) {
  gene_anno <- rtracklayer::readGFF(fpath_hg19_gtf); class(gene_anno)
  write_rds(gene_anno, fpath_hg19_gtf)
} else {
  gene_anno <- read_rds( fpath_hg19_gtf )
}
colnames(gene_anno)
gene_anno$chromosome <- paste0("chr", gene_anno$seqid)
gene_anno$gene <- gene_anno$gene_id
gene_anno$transcript <- gene_anno$transcript_id
gene_anno$symbol <- gene_anno$gene_name

if (T) {
  pos <- subset(gene_anno, strand == "+")
  pos <- pos[order(pos$start),] 
  # remove all but the first exons per transcript
  pos <- pos[!duplicated(pos$transcript),] 
  # make a 1 base pair marker of the TSS
  pos$end <- pos$start + 1 
  neg <- subset(gene_anno, strand == "-")
  neg <- neg[order(neg$start, decreasing = TRUE),] 
  # remove all but the first exons per transcript
  neg <- neg[!duplicated(neg$transcript),] 
  neg$start <- neg$end - 1
  gene_annotation_sub <- rbind(pos, neg)
  # Make a subset of the TSS annotation columns containing just the coordinates 
  # and the gene name
  gene_annotation_sub <- gene_annotation_sub[,c("chromosome", "start", "end", "symbol")]
  # Rename the gene symbol column to "gene"
  names(gene_annotation_sub)[4] <- "gene"
  gene_annotation_sub$chromosome[gene_annotation_sub$chromosome=='chrMT'] <- 'chrM'
}

cat('annotate_cds_by_site...\n')
monocle <- annotate_cds_by_site(monocle, gene_annotation_sub, verbose = TRUE)
tail(fData(monocle))


coaccess_cutoff_override_use <- readRDS(
  file.path(dir_cicero, 'coaccess_cutoff_override_use.value.rds'))
cat('using coaccesss cutoff=', coaccess_cutoff_override_use, '...\n')

cat('build_gene_activity_matrix...\n')
unnorm_ga <- build_gene_activity_matrix(
  monocle, conns, 
  dist_thresh=100e3,
  coaccess_cutoff = coaccess_cutoff_override_use)
class(unnorm_ga) # dgCMatrix
print(dim(unnorm_ga) )
print(range(unnorm_ga)) # count matrix
print(unnorm_ga[1:3, 1:3])
write_rds(unnorm_ga, file.path(dir_cicero, 'ciceroRNA_counts.mat.rds'))
unnorm_ga <- readRDS( file.path(dir_cicero, 'ciceroRNA_counts.mat.rds') )
#------ create seurat object assay ------
cat('CreateSeuratObject...\n')
f_ciceroRNA_assay <- file.path(dir_cicero, 'ciceroRNA_assay.seurat.rds')
ciceroRNA_assay_name <- 'ciceroRNA'
cRNA_assay <- CreateAssayObject(counts = unnorm_ga)
stopifnot(identical(colnames(cRNA_assay), Cells(obja)))
cRNA_assay <- CreateSeuratObject(
  counts = cRNA_assay, assay = 'ciceroRNA', 
  meta.data = obja@meta.data[colnames(cRNA_assay), ])
DefaultAssay(cRNA_assay) <- ciceroRNA_assay_name
cRNA_assay <- NormalizeData(cRNA_assay)
write_rds(cRNA_assay, f_ciceroRNA_assay)

cRNA_assay <- FindVariableFeatures(cRNA_assay)
cRNA_assay <- ScaleData(cRNA_assay)
cRNA_assay <- RunPCA(cRNA_assay)
pdf(file.path(dir_cicero, 'plot.cRNA_assay.ElbowPlot.pdf'), width = 5, height = 4)
print(ElbowPlot(cRNA_assay))
dev.off()
write_rds(cRNA_assay, f_ciceroRNA_assay)

cRNA_assay <- RunUMAP(cRNA_assay, dims = 1:30)
write_rds(cRNA_assay, f_ciceroRNA_assay)
DimPlot(cRNA_assay, group.by = 'clones', label=T)

cat('Done:', sample_name, '\n')
