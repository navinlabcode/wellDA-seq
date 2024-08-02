#---------------------------
# Examine the overlaps of the cells              ----    
# -- This is different from the coda result. 
# -- The cells used here are the cells passing the loose qc criteria.

#---------------------------
library(Seurat)
library(Signac)
library(tidyverse)
library(fs)
library(copykit)
library(VennDiagram)
library(SummarizedExperiment)

setwd("/volumes/USR1/yyan/project/coda")

run_name <- 's7'
f_obja <- './rds/s7_atac/signac/ready.signac_default.rds'
f_objd <- './rds/s7_dna/ready.copykit.rds'

run_name <- 's8'
f_obja <- './rds/s8_atac/signac/ready.signac_default.rds'
f_objd <- './rds/s8_dna/ready.copykit.rds' 



dir_res <- sprintf('./rds_coda/%s', run_name)


dir_create(dir_res)



# if (T) {
message('Run...')
obja <- read_rds(f_obja)
objd <- read_rds(f_objd)

str(Cells(obja))
str(colnames(objd))

cnames_a <- Cells(obja)
cnames_d <- colnames(objd)

if (length(intersect(cnames_a, cnames_d)) == 0) {
  head(colnames(objd))
  cnames_d <- colnames(objd) %>% str_split(pattern = '_', n=3) %>% 
    lapply(., function(x) paste(x[1:2], collapse = '_')) %>%
    unlist()
  head(cnames_d)
  colnames(objd) <- cnames_d
  # objd <- copykit::findClusters(objd, k_superclones=20)
  # objd <- objd[,colData(objd)$subclones != 'c0']  
}
table(colData(objd)$is_dispensed)
objd <- objd[, colData(objd)$is_dispensed]

cnames_a <- Cells(obja); str(cnames_a)
cnames_d <- colnames(objd); str(cnames_d)

library(grDevices)
p <- VennDiagram::venn.diagram(
  list(ATAC=cnames_a, DNA=cnames_d), filename = NULL,
  # col = c('blue', 'gold'), cat.col = c('blue', 'gold'),
  disable.logging=T, cat.default.pos = "text")
pdf(file.path(dir_res, 'vennplot.raw_cellnames.beforeinteresection.pdf'), width = 5, height = 5, useDingbats = F)
grid.draw(p)
dev.off()


setwd("/volumes/USR1/yyan/project/coda")
source('./rsrc/utils.wafargen_physical.R')
run_srch_key <- run_name
if (run_name == 'DCIS22T_chip2_low') {run_srch_key='DCIS22T_chip2'}
if (str_detect(run_name, 'DCIS41T_chip1') ) {run_srch_key='DCIS41T_chip1'}
run_srch_key <- str_remove(run_name, '_low')
message(run_srch_key)

wafar_use <- read_tsv(Sys.glob(sprintf('./wafar/%s/*.TXT', run_srch_key))); dim(wafar_use)
wafar_wid_use <- paste0(wafar_use$Row, '_', wafar_use$Col)
stopifnot( all(wafar_wid_use %in% lib_wafar$well_id) )

# cnames <- rownames(coda_metadf); str(cnames)
cnames <- union(cnames_a, cnames_d)
wid_use <- lib_wafar$well_id[match(cnames, lib_wafar$well_name)]
wid_a <- lib_wafar$well_id[match(cnames_a, lib_wafar$well_name)]
wid_d <- lib_wafar$well_id[match(cnames_d, lib_wafar$well_name)]
library(VennDiagram)
str(setdiff(wid_use, wafar_wid_use))
ggvenn::ggvenn(list(CODA=wid_use, Dispense=wafar_wid_use))
ggvenn::ggvenn(list(DNA=wid_d, ATAC=wid_a, Dispense=wafar_wid_use))

p <- VennDiagram::venn.diagram(
  list(CODA=wid_use, Dispense=wafar_wid_use),
  col = c('blue', 'gold'), cat.col = c('blue', 'gold'),
  cat.default.pos = "text",
  filename = NULL, disable.logging=T)
pdf(file.path(dir_res, 'VennDiagram.well_id.coda.pdf'), width = 5, height = 5, useDingbats = F)
grid.draw(p)
dev.off()

wafar_df <- create_wafar_loc_dataframe(wid_use, wafar_wid_use)
table(wafar_df$well_type)
pa <- plot_wafar_capture(wafar_df) + labs(title='coda'); print(pa)
ggsave(file.path(dir_res, 'wafargen_physical_dispense.pdf'), pa, width = 7, height = 7, useDingbats = F)

