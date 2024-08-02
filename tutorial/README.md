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

Expected result:
- preliminary Signac object for ATAC
- preliminary Copykit object for DNA


### 2. Creating an analysis-ready scCNA data object

See the detailed instruction in [02.scATAC.md](xxx)

Expected result:
(xxx)

### 3. Creating an analysis-ready scATAC data object

See the detailed instruction in [03.scCNA.md](xxx)

Expected result:
(xxx)

### 4. Creating an analysis-ready co-assay data object

See the detailed instruction in [04.wellDA.md](xxx)

Expected result:
(xxx)

### 5. GtoE and EbyG score

### 6. Global concordance score between genotypes and chromatin accessibility profiles

### 7. Heritable and plastic tumor phenotypes

