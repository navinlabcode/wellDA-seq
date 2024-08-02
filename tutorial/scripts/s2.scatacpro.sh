#!/usr/bin/bash

## local version
set -e

sample_name=$1
assay=$2

# sample_name='hbca_s3_c'
# assay='atac'



proj='/volumes/USR1/yyan/project/coda'
indir=${proj}/fq_aggre/${sample_name}/${assay}
outdir=${proj}/scatac_pro/${sample_name}/${assay}
outdir=${proj}/scatac_pro/${sample_name}/${assay}
if [[ ! -e $dir ]]; then
mkdir -p $outdir
fi
config_path=$proj/rsrc/config_scatac_pro/${sample_name}.config.txt 

## BEGIN
cd $outdir ## scATAC_pro uses working_dir (not programmer friendly)
echo $outdir

STATIC_FOLDER='output'
##
## process_no_dex will do everything
##
if [[ ! -e _done_1_demplx_fastq.txt ]]; then
echo 'pre-processing'
scATAC-pro -s trimming \
-i ${indir}/R1.fastq.gz,${indir}/R2.fastq.gz \
-c ${config_path}
touch ${outdir}/_done_1_demplx_fastq.txt
fi

echo 'done'

if [[ ! -e _done_2_mapping.txt ]]; then
echo 'mapping'
scATAC-pro -s mapping \
-i ${outdir}/${STATIC_FOLDER}/trimmed_fastq/${sample_name}.trimmed.demplxed.PE1.fastq.gz,${outdir}/${STATIC_FOLDER}/trimmed_fastq/${sample_name}.trimmed.demplxed.PE2.fastq.gz \
-c ${config_path}
touch ${outdir}/_done_2_mapping.txt
# -i ${outdir}/${STATIC_FOLDER}/trimmed_fastq/R1_val_1.fq.gz,${outdir}/${STATIC_FOLDER}/trimmed_fastq/R2_val_2.fq.gz \
fi

if [[ ! -e _done_3_peakcalling.txt ]]; then
echo 'call_peak'
scATAC-pro -s call_peak \
-i ${outdir}/${STATIC_FOLDER}/mapping_result/${sample_name}.positionsort.MAPQ30.bam \
-c ${config_path}

# Mannually change the source code to run the PE mode
# ~/apps/scATAC-pro.installed/scATAC-pro_1.4.0/scripts/call_peak.sh
# macs2 callpeak \
# -t ${outdir}/${STATIC_FOLDER}/mapping_result/${sample_name}.positionsort.MAPQ30.bam \
# --outdir ${outdir}/${STATIC_FOLDER}/peaks/MACS2 -n s1 -f BAMPE -q 0.05 -g hs
touch ${outdir}/_done_3_peakcalling.txt
fi

if [[ ! -e _done_4_aggr_signal.txt ]]; then
## Time consuming
echo 'aggr_signal'
scATAC-pro -s aggr_signal \
-i ${outdir}/${STATIC_FOLDER}/mapping_result/${sample_name}.positionsort.MAPQ30.bam \
-c ${config_path}
touch ${outdir}/_done_4_aggr_signal.txt
fi

if [[ ! -e _done_5_get_mtx.txt ]]; then
echo 'get_mtx'
scATAC-pro -s get_mtx \
-i ${outdir}/${STATIC_FOLDER}//summary/${sample_name}.fragments.tsv.gz,${outdir}/${STATIC_FOLDER}//peaks/MACS2/${sample_name}_features_BlacklistRemoved.bed \
-c ${config_path}
touch ${outdir}/_done_5_get_mtx.txt
fi  

if [[ ! -e _done_6_qc_per_barcode.txt ]]; then
echo qc_per_barcode
scATAC-pro -s qc_per_barcode \
-i ${outdir}/${STATIC_FOLDER}//summary/${sample_name}.fragments.tsv.gz,${outdir}/${STATIC_FOLDER}//peaks/MACS2/${sample_name}_features_BlacklistRemoved.bed \
-c ${config_path}
touch ${outdir}/_done_6_qc_per_barcode.txt
fi


# Different PEAK_CALLER will save different matrix.rds. Here: MACS2
if [[ ! -e _done_7_call_cell.txt ]]; then
echo call_cell
scATAC-pro -s call_cell \
-i ${outdir}/${STATIC_FOLDER}//raw_matrix/MACS2/matrix.rds \
-c ${config_path} 
touch ${outdir}/_done_7_call_cell.txt
fi

# Different CELL_CALLER will save different results. Here: FILTER
## Create a new BAM files for the good cells and save inside the mapping_result folder 
if [[ ! -e _done_8_get_bam4Cells.txt ]]; then
echo get_bam4Cells
scATAC-pro -s get_bam4Cells \
-i ${outdir}/${STATIC_FOLDER}//mapping_result/${sample_name}.positionsort.bam,${outdir}/${STATIC_FOLDER}//filtered_matrix/MACS2/FILTER/barcodes.txt \
-c ${config_path}
touch ${outdir}/_done_8_get_bam4Cells.txt
fi


# if [[ ! -e _done_.txt ]]; then
#   echo
# 
#   touch ${outdir}/_done_.txt
# fi
#     ## after running the above module, you can run module report (list below)
#     ## to generate first page of the summary report
#     $ scATAC-pro -s rmDoublets
#                  -i ${outdir}/${STATIC_FOLDER}//filtered_matrix/PEAK_CALLER/CELL_CALLER/matrix.rds,0.03 (0.03 is the expected fraction of doublets ) 
#                  -c ${config_path}


# Then below I can use Signac
if [[ ! -e _done_clustering.txt ]]; then
echo clustering
scATAC-pro -s clustering \
-i ${outdir}/${STATIC_FOLDER}//filtered_matrix/MACS2/FILTER/matrix.rds \
-c ${config_path}
touch ${outdir}/_done_clustering.txt
fi

# 
#     $ scATAC-pro -s motif_analysis
#                  -i ${outdir}/${STATIC_FOLDER}//filtered_matrix/PEAK_CALLER/CELL_CALLER/matrix.rds (or matrix.mtx) 
#                  -c ${config_path}
#                  
#     $ scATAC-pro -s split_bam
#                  -i ${outdir}/${STATIC_FOLDER}//downstream_analysis/PEAK_CALLER/CELL_CALLER/cell_cluster_table.tsv
#                  -c ${config_path}
# 
#     $ scATAC-pro -s footprint ## supporting comparison two groups of cell clusters, and one-vs-rest
#                  -i 0,1  ## or '0:3,1:2' (group1 consist of cluster0,3, and group2 for cluster1,2)) or 'one,rest' (all one-vs-rest comparison)
#                  -c ${config_path}
#                  
#     $ scATAC-pro -s runCicero
#                  -i ${outdir}/${STATIC_FOLDER}//downstream_analysis/PEAK_CALLER/CELL_CALLER/seurat_obj.rds
#                  -c ${config_path}
# 
#     $ scATAC-pro -s runDA
#                  -i ${outdir}/${STATIC_FOLDER}//downstream_analysis/PEAK_CALLER/CELL_CALLER/seurat_obj.rds,0:1:3,2  ## group1 consist of cluster 0,1,and 3; group2 cluster2 
#                  -c ${config_path}
#                  
#     $ scATAC-pro -s runGO
#                  -i ${outdir}/${STATIC_FOLDER}//filtered_matrix/PEAK_CALLER/CELL_CALLER/differential_accessible_features_0:1:3_vs_2.tsv,  
#                  -c ${config_path}
#                  
#     $ scATAC-pro -s report
#                  -i ${outdir}/${STATIC_FOLDER}//summary
#                  -c ${config_path}
