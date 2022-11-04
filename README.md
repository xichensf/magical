# MAGICAL analysis

MAGICAL (Multiome Accessible Gene Integration Calling And Looping) analyzes scRNA-seq and scATAC-seq datasets from different conditions. It features a hierarchical Bayesian framework that improves model robustness by leveraging prior transcription factor motif and chromatin domain information. MAGICAL explicitly models signal and data noise in both chromatin accessibility and gene expression measures, which improves the accuracy of regulatory circuits comprising chromatin sites, transcription factors, and target genes that are identified. 

![alt text](https://github.com/xichensf/magical/blob/main/MAGICAL.png)

Please check our paper "*Xi Chen et al. **Mapping disease regulatory circuits at cell-type resolution from single-cell multiomics data**. 2022*" for more details. Questions can be emailed to: Xi Chen, xchen@flatironinstute.org.


## Input files

Tools for scRNA-seq and scATAC-seq data processing are widely available, e.g. Seurat, ArchR. Researchers may have different preference on data processing especially when there are multiple conditions, batches and samples involved. MAGICAL only needs gene symbols, peak coordinates, read count and cell meta information including cell type, sample/subject ID and sample group/condition. These information are very fundamental and can be easily obtained from any single cell multioimc dataset. We provide a [R script](https://github.com/xichensf/magical/blob/main/Multiomics_input_for_MAGICAL.R) to demo how to extra the necessary input files from the single cell multiomics data for use with MAGICAL. [(download demo input files)](https://drive.google.com/file/d/1CerwMHMnS1PNFNMy00OoHQjn6T30M1j4/view?usp=sharing)


#### **Cell type**

MAGICAL infers regulatory circuits for each cell type. Therefore, the input scRNA-seq and scATAC-seq data should be cell type labelled. Users would need to select one cell type and use the provided R script to prepare the following input files for MAGICAL.


#### **Candidate genes (DEG) and chromatin sites (DAS)**

As differential calling is usually done seperately during the scRNA-seq and scATAC-seq data processing, we highly recommand preparing these two files using similar differential statistics cutoffs.  

  * *Candidate gene file*: a list of ``` gene symbols ```
  * *Candidate chromatin site file*: a three-column matrix of ```chr```, ```point1```, and ```point2``` 

#### **scRNA-seq read count**
The scRNA-seq read count information from cells labelled to the selected cell type.   

  * *scRNA read count file*: a three-column matrix with ```gene index```, ```cell index```, and ```RNA read count```  
  * *scRNA gene name file*: a two-column matrix with ```gene index``` and ```gene name```.
  * *scRNA cell meta file*: a five-column matrix with ```cell index```, ```cell barcode```, ```cell type label```, ```sample/subject ID```, and ```condition```

Note, each sample must have a unique name and this name should be the same in the scATAC-seq data to allow MAGICAL to pair the two modalities together. 


#### **scATAC-seq read count**
The scATAC-seq read count information from cells labelled to the selected cell type. 

  * *scATAC read count file*: a three-column matrix with ```peak index```, ```cell index```, and ```ATAC read count```
  * *scATAC peak coordinate file*: a four-column matrix with ```peak index```, ```chr```, ```peak_point1```, ```peak_point2```.
  * *scATAC cell meta file*: a five-column matrix with ```cell index```, ```cell barcode```, ```cell type label```, ```sample/subject ID```, and ```condition```


#### **TF motif mapping (prior)**
We map 870 motifs from the [chromVARmotifs](https://github.com/GreenleafLab/chromVARmotifs) library to all peaks and get the binary mapping relationship. 

  * *Motif mapping file*: a three-column matrix with ```peak index```, ```motif index```, and ```binary binding```.
  * *Motif name file*: a two-column matrix with ```motif index``` and ```motif names```.

Note, the motif names of the same TF can be very different if they are collected from different databases. To avoid duplicated motifs and minimize redundance, we required using the corresponding gene name for each TF motif. 

#### **TAD boundary (prior)**
TAD boundary information can be easily obatined from HiC profiles or similar experiments. We include TAD file (~6000 domains with median size 400kb) in our demo for blood context analysis. 
  * *TAD file*: a three column matrix with ```chr```, ```left_boundary```, and ```right_boundary``` 

In case no proper TAD information or HiC profile is available for the context being studied, another option to use relative distance to TSS (e.g. 500kb) as prior to initally pair peaks and genes. Hg38 RefSeq file is provided for TSS reference.  


A demo (~10mins run):

```
loading all input data ...

We detected 2 conditions from the meta file.

The input scRNAseq data includes 36601 genes, 13248 cells from 12 samples/subecjts.

The input scATACseq data includes 144387 peaks, 13248 cells from 12 samples/subecjts.

There are paired data for 12 samples/subecjts. (check sample IDs if this number is lower than expected)

870 motifs, 945 candidate chromatin sites and 1143 candidate genes are provided.
```


## Regulatory circuit inference

#### **1. Build candidate circuits**  
To identify candidate disease-modulated triads, candidate chromatin sites are associated with TFs by motif sequence matching. These sites are then linked to the candidate genes by requiring them to be within the same TAD or within a user controlled distance. 
```
Candidate regulatory circuits construction ...

MAGICAL model initialization ...
```
#### **2. Infer circuit linkages** 
MAGICAL uses a Bayesian framework to iteratively model chromatin accessibility and gene expression variation across cells and conditions and estimate the strength of the circuit TF-peak-gene linkages. The TF-peak binding confidence and the hidden TF activity are optimized to fit to the chromatin accessibility data. The estimated TF binding strength, TF activity and the input gene expression data are used to infer the peak-gene looping. The updated states of TF-peak-gene linkages based on the estimated strength are used to initialize the next round of estimations. 
```
MAGICAL work starts ...

MAGICAL finished 10 percent

MAGICAL finished 20 percent

MAGICAL finished 30 percent

MAGICAL finished 40 percent

MAGICAL finished 50 percent

MAGICAL finished 60 percent

MAGICAL finished 70 percent

MAGICAL finished 80 percent

MAGICAL finished 90 percent

MAGICAL finished 100 percent
```
#### **3. Output disease modulated circuits** 
Finally, optimized circuits fitting the variation in both modalities are selected. Circuit genes, assocaited chromatin sites and the regulatory TFs will be written to "MAGICAL_selected_regulatory_circuits.txt". MAGICAL uses default thresholds (posterior probabilities on TF-peak binding and peak-gene looping) for circuit selection. Users can adjust these thresholds in the provided scripts to allow more or fewer outputs.  

```
MAGICAL selected regulatory circuits with 93 TFs, 390 peaks and 310 genes.
```

