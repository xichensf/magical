# MAGICAL (Multiome Accessible Gene Integration Calling And Looping)

MAGICAL analyzes paired scRNA-seq and scATAC-seq datasets from different conditions using a hierarchical Bayesian framework that improves model robustness by leveraging prior transcription factor motif and chromatin domain information. MAGICAL explicitly models errors both in chromatin accessibility and in gene expression measures, which improves the accuracy of the condition-specific 3D triads comprising distal chromatin sites, regulators, and their interactions with genes that are identified. 

Please check our paper "Mapping cell type regulatory triads modulated by disease from single-cell multiomics data" for more details. Any questions about the technical details of the paper or about this MAGICAL package can be emailed to: Xi Chen, xchen@flatironinstute.org

MATLAB and R scripts are provided seperately for ***two condition analysis*** (inferring disease modulated regulatory triads) or ***one condition analysis*** (inferring active regulatory triad). 


# MAGICAL input

Users need to prepare following input files (which can be easily obtained from popular single cell data processing packages like Seurat, ArchR (or Signac))

**scRNA-seq Input files**:

*Read count table*: A three-column matrix with *gene index*, *cell index*, and *RNA read count*
*Gene list*: The order of *gene names* should be corresponding to the order in the read count table
*Cell meta*: A four-column matrix with *cell barcode*, *cell type label*, *sample/subject ID*, and *condition* (can be more than two conditions but only up to two conditions will be analyzed in each run)

**scATAC-seq Input files**:

*Read count table*: A three-column matrix with *peak index*, *cell index*, and *ATAC read count*
*Peak list*: A three-column matrix with *chr*, *peak_point1*, *peak_point2*. The order of peaks (chr, point1, point2) should be corresponding to the order in the read count table
*Cell meta*: A four-column matrix with *cell barcode* (can be different from the scRNA-seq cell barcodes), *cell type label* (must be the same as scRNA-seq cell type label), *sample/subject ID* (must be the same as scRNA-seq sample ID), and *condition* (must be the same as scRNA-seq condition)


**Transcriotion factor motif mapping file (prior)**
A binary matrix (can also be continuous but will be converted to binary) with peaks as rows and transcription factor motifs as columns. (This can be easily obtained using ArchR or Signac during scATAC-seq data processing) 


**Topological assciated domain file (prior)**:
A six column matrix with genome coordinates (left_chr, left_point1, left_point2, right_chr, right_point1, right_point2) for the two boundaries of each domain. A no proper TAD information or HiC profile is available for the context being studied. We also provide an option to use 200kb to TSS or 500kb to TSS as prior to initally pair peaks and genes. 


**Candidate gene and peak input files**:

These two files must be prepared in a similar way (i.e., DEG and DAS differential at level, DEG and DAS at both single cell and psuedobulk levels, active genes and active peaks in the selected cell type). The file names should be (Cell type name) candidate genes.txt and (Cell type name) candidate peaks.txt, where the cell type name must be the same as the cell type label used in the scRNA-seq and scATAX-seq datasets. Note, MAGICAL integrates data and infer triads for each cell type. Thus, it is fine to only provide candidate genes and peaks selection for the cell type to be analyzed. 




# MAGICAL output
