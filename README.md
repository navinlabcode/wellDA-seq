# wellDA-seq

This repository includes the codes for the manuscript: '*Single cell genome and epigenome co-profiling reveals hardwiring and plasticity in breast cancer*'. 

It includes: 1) codes for processing and analyzing wellDA-seq data, and 2) codes for repeating the analysis and figures in the manuscript. Specifically, it is organized as following folders: 

1. [<kbd>install</kbd>](https://github.com/navinlabcode/wellDA-seq/tree/main/install). Installing the required the command line tools and packages of R and Python languages. 
2. <kbd>MRE</kbd> ([link](https://github.com/navinlabcode/wellDA-seq/tree/main/MRE)). Tutorial of showing a minimal, reproducible example (https://stackoverflow.com/help/minimal-reproducible-example) of preprocessing and analyzing wellDA-seq data. 
3. <kbd>preprocessing</kbd>. Scripts to take the BCL files of wellDA-seq to create the single-cell CNA and ATAC matrices files. This is for the generic wellDA-seq users to pre-process the wellDA-seq data.
4. <kbd>analysis</kbd>. Codes of computational analysis used in the manuscript. This folder is organized in the order of applications. 
5. <kbd>figure</kbd>. Codes of data visualization shown in the manuscript. This folder is organized in the order of figures. 
- <kbd>data_portal</kbd>. Data objects of input and output. 
