# A minimal, producible example of wellDA-seq data

Here, we will show you how to 1) preprocess a wellDA-seq data to create analysis-ready objects, 2) perform analysis of identifying subclones and cell types, and 3) visualizing and analyzing the interplay of CNA and ATAC data.

This tutorial form the foundation of all analysis that we performed on all other samples in the manuscript

## Software requirement

In brief, we will start proprocessing by using the single-cell [CNV pipeline](https://github.com/navinlabcode/CNV_pipeline) to preprocess the DNA data and using [scATAC-pro](https://github.com/Puriney/scATAC-pro) to proprocess the ATAC data. Then we will create analysis-ready data object for DNA using [copykit](https://github.com/navinlabcode/copykit) and for ATAC using [ArchR](https://github.com/GreenleafLab/ArchR) and [Signac](https://stuartlab.org/signac/). 

See the folder "[wellDA-seq/install](https://github.com/navinlabcode/wellDA-seq/tree/main/install)" for detailed instruction. 

## Proprocessing

## Creating an analysis-ready scCNA data object

## Creating an analysis-ready scATAC data object