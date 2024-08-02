#!/usr/bin/bash

## local version
set -e

sample_name=$1
assay=$2

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

