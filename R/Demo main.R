rm(list = ls())
cat("\f")

setwd('~/Desktop/MAGICAL demo/R')
source('~/Desktop/MAGICAL demo/R/MAGICAL_functions.R')


library(Matrix)
library(dplyr)


#****************************MAGICAL input*******************************
# pre-selected candidate genes and peaks for the cell type
Candidate_gene_file_path = '~/Desktop/MAGICAL demo/Demo input files/Cell type candidate genes.txt'
Candidate_peak_file_path = '~/Desktop/MAGICAL demo/Demo input files/Cell type candidate peaks.txt'


# filtered scRNA data of the cell type
scRNA_readcount_file_path = '~/Desktop/MAGICAL demo/Demo input files/Cell type scRNA read count.txt'
scRNA_gene_file_path = '~/Desktop/MAGICAL demo/Demo input files/scRNA genes.txt'
scRNA_cellmeta_file_path = '~/Desktop/MAGICAL demo/Demo input files/Cell type scRNA cell meta.txt'


# filtered scATAC data of the cell type
scATAC_readcount_file_path = '~/Desktop/MAGICAL demo/Demo input files/Cell type scATAC read count.txt'
scATAC_peak_file_path = '~/Desktop/MAGICAL demo/Demo input files/scATAC peaks.txt'
scATAC_cellmeta_file_path = '~/Desktop/MAGICAL demo/Demo input files/Cell type scATAC cell meta.txt'


# TF motif prior on all ATAC peaks
Motif_mapping_file_path = '~/Desktop/MAGICAL demo/Demo input files/Motif mapping prior.txt'
Motif_name_file_path = '~/Desktop/MAGICAL demo/Demo input files/Motifs.txt'


# TAD prior
TAD_file_path = '~/Desktop/MAGICAL demo/Demo input files/RaoGM12878_40kb_TopDomTADs_filtered_hg38.txt'
distance_control=5e5


# Refseq file for transcription starting site extraction
Ref_seq_file_path = '~/Desktop/MAGICAL demo/Demo input files/hg38_Refseq.txt'



loaded_data <- Data_loading(Candidate_gene_file_path, Candidate_peak_file_path,
                            scRNA_readcount_file_path, scRNA_gene_file_path, scRNA_cellmeta_file_path,
                            scATAC_readcount_file_path, scATAC_peak_file_path, scATAC_cellmeta_file_path,
                            Motif_mapping_file_path, Motif_name_file_path, Ref_seq_file_path)



#****************************MAGICAL analysis*******************************

#Candidate circuits construction
Candidate_circuits <- Candidate_circuits_construction_with_TAD(loaded_data, TAD_file_path)

#Candidate_circuits <- Candidate_circuits_construction_without_TAD(loaded_data, distance_control)


#Model initialization
Initial_model<-MAGICAL_initialization(loaded_data, Candidate_circuits)


#Model parameter estimation
source('~/Desktop/MAGICAL demo/R/MAGICAL_functions.R')

Circuits_linkage_posterior<-MAGICAL_estimation(loaded_data, Candidate_circuits, Initial_model, iteration_num = 1000)



#****************************MAGICAL output******************************* 
#source('~/Desktop/MAGICAL demo/R/MAGICAL_functions.R')
source('~/Desktop/MAGICAL demo/R/MAGICAL_functions.R')

MAGICAL_circuits_output(Output_file_path = 'MAGICAL_selected_regulatory_circuits.txt', 
                        Candidate_circuits, Circuits_linkage_posterior)

save(Candidate_circuits, Circuits_linkage_posterior, file = "MAGICAL_results.RData")

