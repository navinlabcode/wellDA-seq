<!-- Written by: Yun Yan (https://github.com/Puriney) -->

# A minimal, reproducible example of analyzing wellDA-seq data

In this minimal, reproducible example of wellDA-seq, we will address the following questions/tasks: 
1) how to preprocess the wellDA-seq data to create analysis-ready data objects?
2) how to identify subclones and cell types (or states)?
3) how to investigate the interplay of CNA events and open/closed chromatin regions (e.g., the GtoE and EbyG scores, and the plasticity/heritability score)?

This tutorial contains all the analysis that we performed on each of the samples in the manuscript. 

## Software requirement

In brief, the following key tools are required: 
- preprocessing the DNA data by single-cell [CNV pipeline](https://github.com/navinlabcode/CNV_pipeline).
- preprocessing the ATAC data by [scATAC-pro](https://github.com/Puriney/scATAC-pro) and [ArchR](https://github.com/GreenleafLab/ArchR).
- creating an analysis-ready DNA object by [copykit](https://github.com/navinlabcode/copykit).
- creating an analysis-ready ATAC object by [Signac](https://stuartlab.org/signac/).

See the folder [wellDA-seq/install](https://github.com/navinlabcode/wellDA-seq/tree/main/install) for detailed instruction of installing software requirement. 

## Outline

Here, we use the real sample P8 (DCIS66T_chip2) in the manuscript for demonstration.

### 1. Proprocessing

See the detailed instruction in [01.preprocessing.md](https://github.com/navinlabcode/wellDA-seq/blob/main/tutorial/01.preprocessing.md)

Input: FASTQ files of the DNA and ATAC modality data. 

*Expected result*:
- a preliminary Signac object for ATAC
- a preliminary Copykit object for DNA


| item                  | ATAC                                                                                                                                                | DNA                                                                                                                                                           |
| --------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Cell Dispensing       | <img src="https://github.com/navinlabcode/wellDA-seq/raw/main/website_resource/tutorial/atac.wafargen_physical_dispense.png?raw=true" width="250">  | <img src="https://github.com/navinlabcode/wellDA-seq/blob/main/website_resource/tutorial/dna_dispense.png?raw=true" alt="txt" width="250">                    |
| Basic Quality Control | <img src="https://github.com/navinlabcode/wellDA-seq/raw/main/website_resource/tutorial/atac_QC_basic.png?raw=true" alt="txt" width="250">          | <img src="https://github.com/navinlabcode/wellDA-seq/raw/main/website_resource/tutorial/dna_snapshot_CNApipeline.png?raw=true" alt="txt" width="250">         |
| Basic Clustering      | <img src="https://github.com/navinlabcode/wellDA-seq/raw/main/website_resource/tutorial/atac_overview_umap_snn.png?raw=true" alt="txt" width="250"> | <img src="https://github.com/navinlabcode/wellDA-seq/raw/main/website_resource/tutorial/dna.clean.plotHeatmap.clones.004.png?raw=true" alt="txt" width="250"> |




### 2. Initiating the wellDA data

See the detailed instruction in [02.wellDA_initiation.md](https://github.com/navinlabcode/wellDA-seq/blob/main/tutorial/02.wellDA_initiation.md)

*Expected result*:
- a folder of wellDA data

```
├── metadata.csv      <--- data frame of single-cell metadata
├── metadata.df.rds   <--- data frame of single-cell metadata
├── obja.rds          <--- Signac object
├── objd.rds          <--- Copykit object
```

<img src="https://github.com/navinlabcode/wellDA-seq/blob/main/website_resource/tutorial/02.coda.vennplot.cellnames.beforeinteresection.png?raw=true" alt="txt" width="350">

### 3. Further analyze the ATAC part

See the detailed instruction in [03.wellDA_scATAC_annotation.md](https://github.com/navinlabcode/wellDA-seq/blob/main/tutorial/03.wellDA_scATAC_annotation.md)

*Expected result*:
- a Signac object with cell types / states annotated. 

<img src="https://github.com/navinlabcode/wellDA-seq/raw/main/website_resource/tutorial/03.coda_celltype.png?raw=true" alt="txt" width="350">


### 4. Further analyze the CNA part

See the detailed instruction in [04.wellDA_scCNA_annotation.md](https://github.com/navinlabcode/wellDA-seq/blob/main/tutorial/04.wellDA_scCNA_annotation.md)

*Expected result*:
- a Copykit object with the diploid and low-quality cells removed and the tentative subclones determined. 

| Before data cleaning the aneuploid cells                                                                                                                                | After data cleaning the aneuploid cells                                                                                                                                |
| ----------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| <img src="https://github.com/navinlabcode/wellDA-seq/blob/main/website_resource/tutorial/04.aneu.before.cna.heatmap.combo_coda.subclones.001.png?raw=true" width="400"> | <img src="https://github.com/navinlabcode/wellDA-seq/blob/main/website_resource/tutorial/04.aneu.after.cna.heatmap.combo_coda.subclones.001.png?raw=true" width="400"> |

### 5. Refine the wellDA data

See the detailed instruction in [05.wellDA_refine.md](https://github.com/navinlabcode/wellDA-seq/blob/main/tutorial/05.wellDA_refine.md)

*Expected result*:
- an analysis-ready folder of wellDA-seq data for each sample

<img src="https://github.com/navinlabcode/wellDA-seq/blob/main/website_resource/tutorial/05.cna.heatmap.combo_coda.clones.001.png?raw=true" alt="txt" width="700">

### 6. GtoE and EbyG score
See the detailed instruction in [06.GtoE_EbyG.md](https://github.com/navinlabcode/wellDA-seq/blob/main/tutorial/06.GtoE_EbyG.md)

*Expected result*:
- GtoE and EbyG scores for each sample

<img src="https://github.com/navinlabcode/wellDA-seq/blob/main/website_resource/tutorial/06.cna.heatmap.consensus_DAB_DAP_DACH_birdview.png?raw=true" width="800">


### 7. Global concordance score between genotypes and chromatin accessibility profiles
See the detailed instruction in [07.global_concordance.md](https://github.com/navinlabcode/wellDA-seq/blob/main/tutorial/07.global_concordance.md)


*Expected result*:
- a global concordance score for each sample


| UMAP of CNA and ATAC of single-cell colored by clones       | <img src="https://github.com/navinlabcode/wellDA-seq/blob/main/website_resource/tutorial/07.coda_crossdimplot.clones.png?raw=true" width="600"> |
| ----------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------- |
| Global concordance score between genotypes and epigenotypes | <img src="https://github.com/navinlabcode/wellDA-seq/blob/main/website_resource/tutorial/07.heatmap.GP_concordance.png?raw=true" width="600">   |



### 8. Heritable and plastic tumor phenotypes

See the detailed instruction in [08.plasticity_heritability.md](https://github.com/navinlabcode/wellDA-seq/blob/main/tutorial/08.plasticity_heritability.md)

*Expected result*:

<img src="" width='600'>

