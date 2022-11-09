rm(list = ls())
setwd('~/Desktop/MAGICAL demo/R')
#source('~/Desktop/MAGICAL demo/R/MAGICAL_functions.R')

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
TAD_flag=0 #if no TAD provided, simply set the path as empty and set the flag to 0
TAD_file_path = '~/Desktop/MAGICAL demo/Demo input files/RaoGM12878_40kb_TopDomTADs_filtered_hg38.txt'


# Refseq file for transcription starting site extraction
Ref_seq_file_path = '~/Desktop/MAGICAL demo/Demo input files/hg38_Refseq.txt'


#Output file 
Output_file_path = 'MAGICAL_selected_regulatory_circuits.txt'
prob_threshold_TF_peak_binding=0.8
prob_threshold_peak_gene_looping=0.8




#*********************load all input data**********************
print('loading all input data ...')

#**********load candidate genes**********
Candidate_genes <- read.table(Candidate_gene_file_path, sep='\t')
colnames(Candidate_genes) = c("gene_symbols") 


#**********load scRNAseq data**********
scRNA_genes<-read.table(scRNA_gene_file_path, sep='\t')
colnames(scRNA_genes) = c("gene_index", "gene_symbols")

scRNA_cells<-read.table(scRNA_cellmeta_file_path, sep='\t')
colnames(scRNA_cells) = c("cell_index", "cell_barcode", "cell_type", "subject_ID", "condition")

scRNA_read_count_table<-read.table(scRNA_readcount_file_path, sep='\t')
scRNA_read_count_matrix <- sparseMatrix(i = scRNA_read_count_table[,1], j=scRNA_read_count_table[,2], x=scRNA_read_count_table[,3],
             dimnames=list(levels(scRNA_read_count_table[,1]), levels(scRNA_read_count_table[,2])), repr = "T")
rm(scRNA_read_count_table)

Group_conditions = unique(scRNA_cells$condition)
print(paste('We detected', length(Group_conditions), 'conditions from the meta file', sep = ' '))

scRNA_samples = unique(scRNA_cells$subject_ID)
print(paste('The input scRNAseq data includes', length(scRNA_genes$gene_symbols),'genes and', length(scRNA_cells$cell_barcode),'cells from',length(scRNA_samples),'samples', sep = ' '))



#**********load candidate peaks**********
Candidate_peaks <- read.table(Candidate_peak_file_path, sep='\t')
colnames(Candidate_peaks) = c("chr", "point1", 'point2')


#**********load scATACseq data**********
scATAC_peaks <- read.table(scATAC_peak_file_path, sep='\t')
colnames(scATAC_peaks) = c("peak_index", "chr", "point1", 'point2')

scATAC_cells <- read.table(scATAC_cellmeta_file_path, sep ="\t")
colnames(scATAC_cells) = c("cell_index", "cell_barcode", "cell_type", 'subject_ID')

scATAC_read_count_table <- read.table(scATAC_readcount_file_path, sep ="\t")
scATAC_read_count_matrix <- sparseMatrix(i = scATAC_read_count_table[,1], j=scATAC_read_count_table[,2], x=scATAC_read_count_table[,3],
                                        dimnames=list(levels(scATAC_read_count_table[,1]), levels(scATAC_read_count_table[,2])), repr = "T")
rm(scATAC_read_count_table)

scATAC_samples = unique(scATAC_cells$subject_ID)
print(paste('The input scATACseq data includes',length(scATAC_peaks$peak_index), 'peaks and', length(scATAC_cells$cell_barcode), 'cells from', length(scATAC_samples), 'samples', sep = ' '))

Common_samples=intersect(scRNA_samples, scATAC_samples);
print(paste('There are paired data for', length(Common_samples), 'samples. (Please check sample IDs if this number is lower than expected).', sep = ' '))



#**********load TF motif prior**********
Motifs<-read.table(Motif_name_file_path, sep ="\t")
colnames(Motifs) = c("motif_index", "name")
Motif_mapping_table <- read.table(Motif_mapping_file_path)
TF_peak_binding_matrix <- sparseMatrix(i = Motif_mapping_table[,1], j=Motif_mapping_table[,2], x=Motif_mapping_table[,3],
                                         dimnames=list(levels(Motif_mapping_table[,1]), levels(Motif_mapping_table[,2])), repr = "T")
rm(Motif_mapping_table)


#**********load TAD prior**********
#TAD_flag=1 if TAD prior is provided, TAD_flag=0 if not
if (TAD_flag==1) {
  TAD <- read.table(TAD_file_path, sep='\t')
  colnames(TAD) = c("chr", "left_boundary", 'right_boundary')  
} else {
  Distance_control=5e5 #500kb to TSS
}

Refseq <- read.table(Ref_seq_file_path, header = TRUE, sep='\t')


print(paste(length(Motifs$name), 'motifs,', length(Candidate_peaks$chr), 'candidate chromatin sites and', length(Candidate_genes$gene_symbols), 'candidate genes are provided.', sep = ' '))


