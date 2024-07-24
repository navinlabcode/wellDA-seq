# Preprocessing


## Preprocessing the scATAC modality data

### Generating scATAC fragment file by scATAC-pro

With FASTQ files as input, scATAC-pro maps reads to genome and generates the ATAC fragment file of Tn5 insertions. The ATAC fragment file is the common format of inputs of other downstream analysis tools, including ArchR and Signac. 

To run scATAC-pro, it first needs to make scATAC-pro recognize reads for each cell. We attach the cell barcode names to the reads in the FASTQ files. 

For example, the cell barcode name 'A13_A10' was inserted into the specific read. This is the format scATAC-pro require. 

```
@VH00219:49:AACYKCYM5:1:1101:20674:7759 1:N:0:TAAGGCGAATC+CGATACACTCA
CTAAAGGGTTGTAGCTGTGAAGGAATAATGGGCGGAGTGGTGGCACACACGGCATT
+
C;C-CC-CC-CCC-;CCCCCC;CCCCCCCCC-CC;CC;CCCC;;CCCCCCCCCCCC
@VH00219:49:AACYKCYM5:1:1101:23420:27239 1:N:0:TAAGGCGAATC+CGATACACTCA
AGGCTGGGAAGCCCAAGACCAAAGTACCAGCAGATTTGGTGTCTGGTGAGGGGCTG
+
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
```

```
@A13_A10:VH00219:49:AACYKCYM5:1:1101:20674:7759 1:N:0:TAAGGCGAATC+CGATACACTCA
CTAAAGGGTTGTAGCTGTGAAGGAATAATGGGCGGAGTGGTGGCACACACGGCATT
+
C;C-CC-CC-CCC-;CCCCCC;CCCCCCCCC-CC;CC;CCCC;;CCCCCCCCCCCC
@A13_A10:VH00219:49:AACYKCYM5:1:1101:23420:27239 1:N:0:TAAGGCGAATC+CGATACACTCA
AGGCTGGGAAGCCCAAGACCAAAGTACCAGCAGATTTGGTGTCTGGTGAGGGGCTG
+
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
```


### Quality-control by ArchR

### Creating Signac object for analysis

## Preprocessing the scDNA modality data

### Generating scCNA matrix and quality control by CNA-pipeline

### Creating copykit object for analysis