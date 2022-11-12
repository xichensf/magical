rm(list = ls())
setwd('~/Desktop/MAGICAL demo/R')
source('~/Desktop/MAGICAL demo/R/MAGICAL_functions.R')

library(Matrix)

#***************set input file path***************************************
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
TAD_flag = 1 #if no TAD provided, simply set the path as empty and set the flag to 0
TAD_file_path = '~/Desktop/MAGICAL demo/Demo input files/RaoGM12878_40kb_TopDomTADs_filtered_hg38.txt'

#TAD_flag=0 #if no TAD provided, simply set the path as empty and set the flag to 0
distance_control=5e5

# Refseq file for transcription starting site extraction
Ref_seq_file_path = '~/Desktop/MAGICAL demo/Demo input files/hg38_Refseq.txt'

#Output file 
Output_file_path = 'MAGICAL_selected_regulatory_circuits.txt'
prob_threshold_TF_peak_binding=0.8
prob_threshold_peak_gene_looping=0.8


#**************************load all data into the R environment**********

loaded_data <- Data_loading(Candidate_gene_file_path, Candidate_peak_file_path,
                            scRNA_readcount_file_path, scRNA_gene_file_path, scRNA_cellmeta_file_path,
                            scATAC_readcount_file_path, scATAC_peak_file_path, scATAC_cellmeta_file_path,
                            Motif_mapping_file_path, Motif_name_file_path, TAD_flag, TAD_file_path, Ref_seq_file_path,
                            Output_file_path, prob_threshold_TF_peak_binding, prob_threshold_peak_gene_looping)

 
#********************** candidate circuits construction ****************

if (TAD_flag==1){
  Candidate_circuits <- Candidate_circuits_construction_with_TAD(loaded_data)
}else{
  Candidate_circuits <- Candidate_circuits_construction_without_TAD(loaded_data, distance_control)
}


#************************* MAGICAL initialization *********************
  
Initial_model<-MAGICAL_initialization(loaded_data, Candidate_circuits)


#************************* MAGICAL sampling ***************************
  
Circuits_linkage_posterior<-MAGICAL_estimation(loaded_data, Candidate_circuits, Initial_model)


#************************* MAGICAL output *****************************
#*
#*
#*