
# PBMC single cell multiomics data generated in the paper

We provide download links for datasets generated in the paper "***Mapping disease regulatory circuits at cell-type resolution from single-cell multiomics data***". The sample-paired PBMC scRNA-seq and scATAC-seq data of *S. aureus* infected subjects and uninfected controls were respectively processed using Seurat and ArchR. The PBMC scATAC-seq data of mild COVID-19 subjects and uninfected controls were processed using ArchR. 

  * The Seurat R object of the intergated ***S. aureus*** scRNA-seq data can be downloaded [here](https://wisp.princeton.edu/media/magical/MRSA-MSSA-CTRL-all-combine-20210908.RData.gz) (15GB). 
  * The ArchR R project (including arrow files) of the integrated ***S. aureus*** scATAC-seq data can be downloaded [here](https://wisp.princeton.edu/media/magical/Staph_scATAC_integration.tar.gz) (34GB).
  * The ArchR R project (including arrow files) of the integrated **COVID-19** scATAC-seq data can be downloaded [here](https://wisp.princeton.edu/media/magical/COVID19_scATAC_integration.tar.gz) (7GB).


![alt text](https://github.com/xichensf/magical/blob/main/UMAP.png)



# MAGICAL analysis

MAGICAL (Multiome Accessible Gene Integration Calling And Looping) analyzes scRNA-seq and scATAC-seq datasets from different conditions. It features a hierarchical Bayesian framework that improves model robustness by leveraging prior transcription factor motif and chromatin domain information. MAGICAL explicitly models signal and data noise in both chromatin accessibility and gene expression measures, which improves the accuracy of regulatory circuits comprising chromatin sites, transcription factors, and target genes that are identified. 

![alt text](https://github.com/xichensf/magical/blob/main/MAGICAL.png)







## MAGICAL input

MAGICAL only requires gene symbols, peak coordinates, read count and cell meta information including cell type, sample/subject ID and sample group/condition. These information are very fundamental and can be easily obtained from any single cell multioimc dataset. We provide [Multiomics_input_for_MAGICAL.R](https://github.com/xichensf/magical/blob/main/Multiomics_input_for_MAGICAL.R) to demo how to extra the necessary input files from the single cell multiomics data for use with MAGICAL. The script includes code to extract needed information from Seurat processed scRNA-seq data and ArchR or Signac processed scATAC-seq data. Since the single cell data files are usually very large, for testing MAGICAL, we provide demo input files that users can easily download and run MAGICAL on their local machine. [(download demo input files)](https://drive.google.com/file/d/1CerwMHMnS1PNFNMy00OoHQjn6T30M1j4/view?usp=sharing)


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
TAD boundary information can be easily obatined from HiC profiles or similar experiments. We include in our demo a GM12878 HiC-based TAD file including ~6000 domains with median size 400kb for blood context analysis. 
  * *TAD file*: a three column matrix with ```chr```, ```left_boundary```, and ```right_boundary``` 

In case no proper TAD information or HiC profile is available for the context being studied, another option to use relative distance to TSS (e.g. 500kb) as prior to initally pair peaks and genes. A hg38 RefSeq file is included in the demo for TSS reference.  


A demo (~10mins run):

```
loading all input data ...

We detected 2 conditions from the meta file.

The input scRNAseq data includes 36601 genes, 13248 cells from 12 samples/subecjts.

The input scATACseq data includes 144387 peaks, 13248 cells from 12 samples/subecjts.

There are paired data for 12 samples/subecjts. (check sample IDs if this number is lower than expected)

870 motifs, 945 candidate chromatin sites and 1143 candidate genes are provided.
```


## MAGICAL circuit inference

#### **Build candidate circuits**  
To identify candidate disease-modulated triads, candidate chromatin sites are associated with TFs by motif sequence matching. These sites are then linked to the candidate genes by requiring them to be within the same TAD or within a user controlled distance. 
```
Candidate regulatory circuits construction ...

MAGICAL model initialization ...
```
#### **Infer circuit linkages** 
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
## MAGICAL output
Finally, optimized circuits fitting the variation in both modalities are selected. Circuit genes, assocaited chromatin sites and the regulatory TFs will be written to [MAGICAL_selected_regulatory_circuits.txt](https://github.com/xichensf/magical/blob/main/MAGICAL_selected_regulatory_circuits.txt). MAGICAL uses default thresholds (posterior probabilities on TF-peak binding and peak-gene looping) for circuit selection. Users can adjust these thresholds in the provided scripts to allow more or fewer outputs.  

```
MAGICAL selected regulatory circuits with 93 TFs, 390 peaks and 310 genes.
```


## Contact
Questions about scATAC-seq data processing and MAGICAL can be emailed to Xi Chen (<xchen@flatironinstitute.org>).
Questions about scRNA-seq data processing can be emailed to Yuan Wang (<yuanwang@cs.princeton.edu>).
Other questions about our work should be emailed to Olga G. Troyanskaya (<ogt@genomics.princeton.edu>) and Stuart C. Sealfon (<stuart.sealfon@mssm.edu>).
