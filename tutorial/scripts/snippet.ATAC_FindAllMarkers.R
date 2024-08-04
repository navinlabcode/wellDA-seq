#---------------------------
# Snippet of lines to run FindAllMarkers              ----    
# No need to change
# 
# Variables in use: 
# obja; dir_snippet; ident_use
#---------------------------
message('Run FindAllMarkers on ident: ', ident_use)

DefaultAssay(obja) <- 'iRNA'
Idents(obja) <- ident_use
print(table(Idents(obja)))
message(ident_use)
f_deg <- file.path(dir_snippet, sprintf('deg0.%s.df.rds', ident_use))

if (! file.exists(f_deg)) {
  df_deg0 <- FindAllMarkers(
    obja, min.pct=1/100,
    test.use = 'wilcox',
    assay = 'iRNA', 
    slot = 'data',
    verbose = T
  )
  write_rds(df_deg0, f_deg)
} else {
  df_deg0 <- read_rds(f_deg)
}
df_deg <- df_deg0 %>% dplyr::filter(
  p_val_adj < 0.05, avg_log2FC > log2(1)) %>%
  dplyr::arrange(desc(avg_log2FC))


## manual rename clusters
f_atac_manual <- file.path(dir_snippet, sprintf('deg.iRNA.%s.manual_celltype.csv', ident_use))
write_lines(levels(Idents(obja)), file = f_atac_manual)

message('Todo: Fill in the ATAC cluster names in: ', f_atac_manual)

DefaultAssay(obja) <- 'peaks'
