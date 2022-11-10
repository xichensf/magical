#**************************************************************************
#                         MAGICAL functions
#                            11/10/2022
#
#***************************************************************************




#***************************load all input files***************************
Data_loading <- function(Candidate_gene_file_path, Candidate_peak_file_path,
                         scRNA_readcount_file_path, scRNA_gene_file_path, scRNA_cellmeta_file_path,
                         scATAC_readcount_file_path, scATAC_peak_file_path, scATAC_cellmeta_file_path,
                         Motif_mapping_file_path, Motif_name_file_path, TAD_flag, TAD_file_path, Ref_seq_file_path,
                         Output_file_path, prob_threshold_TF_peak_binding, prob_threshold_peak_gene_looping)
  {
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
  scATAC_peaks <- scATAC_peaks[,c(1,2,3,4)]
  colnames(scATAC_peaks) = c("peak_index", "chr", "point1", 'point2')
  
  scATAC_cells <- read.table(scATAC_cellmeta_file_path, sep ="\t")
  colnames(scATAC_cells) = c("cell_index", "cell_barcode", "cell_type", 'subject_ID')
  
  scATAC_read_count_table <- read.table(scATAC_readcount_file_path, sep ="\t")
  scATAC_read_count_matrix <- sparseMatrix(i = scATAC_read_count_table[,1], j=scATAC_read_count_table[,2], x=scATAC_read_count_table[,3],
                                           dimnames=list(levels(scATAC_read_count_table[,1]), levels(scATAC_read_count_table[,2])), repr = "T")
  rm(scATAC_read_count_table)
  
  scATAC_samples = unique(scATAC_cells$subject_ID)
  print(paste('The input scATACseq data includes',length(scATAC_peaks$peak_index), 'peaks and', length(scATAC_cells$cell_barcode), 'cells from', length(scATAC_samples), 'samples', sep = ' '))
  
  Common_samples=intersect(scRNA_samples, scATAC_samples)
  print(paste('There are sample-paired single cell multiomics data for', length(Common_samples), 'samples. (Please check sample IDs if this number is lower than expected).', sep = ' '))
  
  
  
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
  
  loaded_data = list('Common_samples'=Common_samples, 
                     'Candidate_genes'=Candidate_genes, 
                     'Candidate_peaks'=Candidate_peaks,
                     'scRNA_genes'=scRNA_genes,
                     'scRNA_cells'=scRNA_cells,
                     'scRNA_read_count_matrix'=scRNA_read_count_matrix,
                     'scATAC_peaks'=scATAC_peaks, 
                     'scATAC_cells'=scATAC_cells,
                     'scATAC_read_count_matrix'=scATAC_read_count_matrix,
                     'Motifs'=Motifs, 
                     'TF_peak_binding_matrix'=TF_peak_binding_matrix,
                     'Refseq'=Refseq, 
                     'TAD'=TAD)
  return(loaded_data)
}






#***************************build candidate circuits***************************
Candidate_circuits_construction_with_TAD <- function(loaded_data){
  print('Candidate regulatory circuits constrcution ...')
  
  Common_samples=loaded_data$Common_samples
  Candidate_genes=loaded_data$Candidate_genes
  Candidate_peaks=loaded_data$Candidate_peaks
  scRNA_genes=loaded_data$scRNA_genes 
  scRNA_cells=loaded_data$scRNA_cells 
  scRNA_read_count_matrix=loaded_data$scRNA_read_count_matrix
  scATAC_peaks=loaded_data$scATAC_peaks 
  scATAC_cells=loaded_data$scATAC_cells 
  scATAC_read_count_matrix=loaded_data$scATAC_read_count_matrix
  Motifs=loaded_data$Motifs 
  TF_peak_binding_matrix=loaded_data$TF_peak_binding_matrix
  Refseq=loaded_data$Refseq
  TAD=loaded_data$TAD
  
  
  # TF-peak binding
  library(dplyr)
  Candidate_peaks<-inner_join(scATAC_peaks, Candidate_peaks)
  index<-match(Candidate_peaks$peak_index, scATAC_peaks$peak_index)
  Candidate_TF_Peak_Binding<-TF_peak_binding_matrix[index,]
  TF_num<-colSums(Candidate_TF_Peak_Binding)
  TF_pct_in_candidate_peaks<-colSums(Candidate_TF_Peak_Binding)/nrow(Candidate_peaks)
  TF_pct_in_all_peaks<-colSums(TF_peak_binding_matrix)/nrow(scATAC_peaks)
  TF_enrichment_FC<-TF_pct_in_candidate_peaks/TF_pct_in_all_peaks
  TF_index<-which(TF_pct_in_candidate_peaks>0.05 & TF_num>30 & TF_enrichment_FC>1.2)
  if (length(TF_index) == 0){
    print('Too few peaks with TF binding sites.')
    print('MAGICAL not applicable to this cell type!')
  }else{
    Candidate_TFs <- Motifs[TF_index,]
    Candidate_TF_Peak_Binding <- Candidate_TF_Peak_Binding[,TF_index]
  }
  peak_index<-which(rowSums(Candidate_TF_Peak_Binding)>0)
  Candidate_peaks <- Candidate_peaks[peak_index,]
  Candidate_TF_Peak_Binding <- Candidate_TF_Peak_Binding[peak_index,]
  
  
  
  # extract gene TSS
  gene_symbols<-intersect(intersect(Candidate_genes$gene_symbols, Refseq$gene_name), scRNA_genes$gene_symbols)
  gene_chr=character(length(gene_symbols))
  gene_TSS=integer(length(gene_symbols))
  Candidate_genes=data.frame(gene_symbols, gene_chr, gene_TSS)
  colnames(Candidate_genes) = c("gene_symbols", "chr", "TSS")
  for (g in 1:length(gene_symbols)){
    index = which(Refseq$gene_name==Candidate_genes$gene_symbols[g])
    Candidate_genes$chr[g] <- Refseq$chr[index[1]]
    if (Refseq$strand[index[1]]=='+'){
      Candidate_genes$TSS[g] <- min(Refseq$start[index])
    }else{
      Candidate_genes$TSS[g] <- max(Refseq$end[index])
    }
  }
  
  
  
  # Peak-gene looping
  Candidate_Peak_Gene_looping_TAD=sparseMatrix(i=1,j=1,x=0, dims=c(nrow(Candidate_peaks), nrow(Candidate_genes)), repr = "T")
  for (t in 1:nrow(TAD)){
#    if (TAD$right_boundary[t]-TAD$left_boundary[t]<2e6){
      peak_index=which(Candidate_peaks$chr==TAD$chr[t] & Candidate_peaks$point1>TAD$left_boundary[t] & Candidate_peaks$point2<TAD$right_boundary[t])
      gene_index=which(Candidate_genes$chr==TAD$chr[t] & Candidate_genes$TSS>TAD$left_boundary[t] & Candidate_genes$TSS<TAD$right_boundary[t])
      if (length(peak_index)>0 && length(gene_index)>0){
        Candidate_Peak_Gene_looping_TAD[peak_index,gene_index]=1
      }
#    }
  }
  #some TAD could be extra wide, based on HiC interaction scale, we limit the peak-gene distance upto 1M bps
  Candidate_Peak_Gene_looping_distance_control=sparseMatrix(i=1,j=1,x=0, dims=c(nrow(Candidate_peaks), nrow(Candidate_genes)), repr = "T")
  for (g in 1:length(Candidate_genes$gene_symbols)){
    index=which(Candidate_peaks$chr==Candidate_genes$chr[g] & abs((Candidate_peaks$point1+Candidate_peaks$point2)/2-Candidate_genes$TSS[g])<1e6)
    if (length(index)>0){
      Candidate_Peak_Gene_looping_distance_control[index,g]=1
    }
  }
  Candidate_Peak_Gene_looping<-Candidate_Peak_Gene_looping_TAD*Candidate_Peak_Gene_looping_distance_control
  peak_index<-which(rowSums(Candidate_Peak_Gene_looping)>0)
  gene_index<-which(colSums(Candidate_Peak_Gene_looping)>0)
  Candidate_Peak_Gene_looping<-Candidate_Peak_Gene_looping[peak_index,gene_index]
  Candidate_peaks<-Candidate_peaks[peak_index,]
  Candidate_genes<-Candidate_genes[gene_index,]
  TF_index=which(colSums(Candidate_TF_Peak_Binding[peak_index,])>0)
  Candidate_TFs=Candidate_TFs[TF_index,]
  Candidate_TF_Peak_Binding=Candidate_TF_Peak_Binding[peak_index, TF_index]
  
  
  
  
  #pseudobulk data is calculated for model initialization
  ATAC_count=matrix(0,nrow=nrow(scATAC_read_count_matrix), ncol=length(Common_samples))
  for (s in 1:length(Common_samples)){
    ATAC_count[,s]=rowSums(scATAC_read_count_matrix[,which(scATAC_cells$subject_ID==Common_samples[s])])
  }
  colnames(ATAC_count)<-Common_samples
  rownames(ATAC_count)<-scATAC_peaks$peak_index
  
  RNA_count=matrix(0,nrow=nrow(scRNA_read_count_matrix), ncol=length(Common_samples))
  for (s in 1:length(Common_samples)){
    RNA_count[,s]=rowSums(scRNA_read_count_matrix[,which(scRNA_cells$subject_ID==Common_samples[s])])
  }
  colnames(RNA_count)<-Common_samples
  rownames(RNA_count)<-scRNA_genes$gene_symbols
  
  #normalization
  total_ATAC_reads=5e6
  ATAC_scale_facter=total_ATAC_reads/(colSums(ATAC_count)+1)
  ATAC_count=sweep(ATAC_count, 2, ATAC_scale_facter, "*")
  ATAC_log2=log2(ATAC_count+1)
    
  total_RNA_reads=5e6
  RNA_scale_facter=total_RNA_reads/(colSums(RNA_count)+1)
  RNA_count=sweep(RNA_count, 2, RNA_scale_facter, "*")
  RNA_log2=log2(RNA_count+1)
  
  
  # select actively accessbile peaks
  index=match(Candidate_peaks$peak_index, scATAC_peaks$peak_index)
  Candidate_Peak_log2Count=ATAC_log2[index,]
  
  peak_index=which(rowSums(Candidate_Peak_log2Count)>0)
  Candidate_peaks=Candidate_peaks[peak_index,]
  Candidate_TF_Peak_Binding=Candidate_TF_Peak_Binding[peak_index,]
  Candidate_Peak_Gene_looping=Candidate_Peak_Gene_looping[peak_index,]
  Candidate_Peak_log2Count=Candidate_Peak_log2Count[peak_index,]
  Candidate_Peak_log2Count=sweep(Candidate_Peak_log2Count, 1, rowMeans(Candidate_Peak_log2Count), '-')
  
  
  # select actively expressed genes
  index=match(Candidate_genes$gene_symbols, scRNA_genes$gene_symbols)
  Candidate_Gene_log2Count=RNA_log2[index,]
  
  gene_index=which(rowSums(Candidate_Gene_log2Count)>0)
  Candidate_genes=Candidate_genes[gene_index,]
  Candidate_Peak_Gene_looping=Candidate_Peak_Gene_looping[,gene_index]
  Candidate_Gene_log2Count=Candidate_Gene_log2Count[gene_index,]
  Candidate_Gene_log2Count=sweep(Candidate_Gene_log2Count, 1, rowMeans(Candidate_Gene_log2Count), '-')
  
  
  # select actively expressed TFs
  selected_TFs=intersect(Candidate_TFs$name, scRNA_genes$gene_symbols)
  index=match(selected_TFs, Candidate_TFs$name)
  Candidate_TFs=Candidate_TFs[index,]
  
  index=match(selected_TFs, scRNA_genes$gene_symbols)
  Candidate_TF_Peak_Binding=Candidate_TF_Peak_Binding[,index]
  Candidate_TF_log2Count=RNA_log2[index,]
  
  TF_index=which(rowSums(Candidate_TF_log2Count)>0)
  Candidate_TFs=Candidate_TFs[TF_index,]
  Candidate_TF_Peak_Binding=Candidate_TF_Peak_Binding[,TF_index]
  Candidate_TF_log2Count=Candidate_TF_log2Count[TF_index,]
  Candidate_TF_log2Count=sweep(Candidate_TF_log2Count, 1, rowMeans(Candidate_TF_log2Count), '-')
  
  
  
  # final filtering for candidate circuits
  peak_index=which(rowSums(Candidate_Peak_Gene_looping)>0 & rowSums(Candidate_TF_Peak_Binding)>0)
  gene_index=which(colSums(Candidate_Peak_Gene_looping[peak_index,])>0)
  TF_index=which(colSums(Candidate_TF_Peak_Binding[peak_index,])>0)
  
  Candidate_TF_Peak_Binding=Candidate_TF_Peak_Binding[peak_index,TF_index]
  Candidate_Peak_Gene_looping=Candidate_Peak_Gene_looping[peak_index,gene_index]
  
  Candidate_TFs=Candidate_TFs[TF_index,]
  Candidate_TF_log2Count=Candidate_TF_log2Count[TF_index,]
  
  Candidate_peaks=Candidate_peaks[peak_index,]
  Candidate_Peak_log2Count=Candidate_Peak_log2Count[peak_index,]  
  
  Candidate_genes=Candidate_genes[gene_index,]
  Candidate_Gene_log2Count=Candidate_Gene_log2Count[gene_index,]
  
  Candidate_circuits=list('TFs'=Candidate_TFs, 
                          'TF_log2Count'=Candidate_TF_log2Count,
                          'Peaks'=Candidate_peaks, 
                          'Peak_log2Count'=Candidate_Peak_log2Count,
                          'Genes'=Candidate_genes, 
                          'Gene_log2Count'=Candidate_Gene_log2Count,
                          'TF_Peak_Binding'=Candidate_TF_Peak_Binding, 
                          'Peak_Gene_looping'=Candidate_Peak_Gene_looping)
  return(Candidate_circuits)
}



#***************************build candidate circuits***************************
Candidate_circuits_construction_without_TAD <- function(loaded_data, distance_control){
  print('Candidate regulatory circuits constrcution ...')
  
  Common_samples=loaded_data$Common_samples
  Candidate_genes=loaded_data$Candidate_genes
  Candidate_peaks=loaded_data$Candidate_peaks
  scRNA_genes=loaded_data$scRNA_genes 
  scRNA_cells=loaded_data$scRNA_cells 
  scRNA_read_count_matrix=loaded_data$scRNA_read_count_matrix
  scATAC_peaks=loaded_data$scATAC_peaks 
  scATAC_cells=loaded_data$scATAC_cells 
  scATAC_read_count_matrix=loaded_data$scATAC_read_count_matrix
  Motifs=loaded_data$Motifs 
  TF_peak_binding_matrix=loaded_data$TF_peak_binding_matrix
  Refseq=loaded_data$Refseq
  
  
  # TF-peak binding
  library(dplyr)
  Candidate_peaks<-inner_join(scATAC_peaks, Candidate_peaks)
  index<-match(Candidate_peaks$peak_index, scATAC_peaks$peak_index)
  Candidate_TF_Peak_Binding<-TF_peak_binding_matrix[index,]
  TF_num<-colSums(Candidate_TF_Peak_Binding)
  TF_pct_in_candidate_peaks<-colSums(Candidate_TF_Peak_Binding)/nrow(Candidate_peaks)
  TF_pct_in_all_peaks<-colSums(TF_peak_binding_matrix)/nrow(scATAC_peaks)
  TF_enrichment_FC<-TF_pct_in_candidate_peaks/TF_pct_in_all_peaks
  TF_index<-which(TF_pct_in_candidate_peaks>0.05 & TF_num>30 & TF_enrichment_FC>1.2)
  if (length(TF_index) == 0){
    print('Too few peaks with TF binding sites.')
    print('MAGICAL not applicable to this cell type!')
  }else{
    Candidate_TFs <- Motifs[TF_index,]
    Candidate_TF_Peak_Binding <- Candidate_TF_Peak_Binding[,TF_index]
  }
  peak_index<-which(rowSums(Candidate_TF_Peak_Binding)>0)
  Candidate_peaks <- Candidate_peaks[peak_index,]
  Candidate_TF_Peak_Binding <- Candidate_TF_Peak_Binding[peak_index,]
  
  
  
  # extract gene TSS
  gene_symbols<-intersect(intersect(Candidate_genes$gene_symbols, Refseq$gene_name), scRNA_genes$gene_symbols)
  gene_chr=character(length(gene_symbols))
  gene_TSS=integer(length(gene_symbols))
  Candidate_genes=data.frame(gene_symbols, gene_chr, gene_TSS)
  colnames(Candidate_genes) = c("gene_symbols", "chr", "TSS")
  for (g in 1:length(gene_symbols)){
    index = which(Refseq$gene_name==Candidate_genes$gene_symbols[g])
    Candidate_genes$chr[g] <- Refseq$chr[index[1]]
    if (Refseq$strand[index[1]]=='+'){
      Candidate_genes$TSS[g] <- min(Refseq$start[index])
    }else{
      Candidate_genes$TSS[g] <- max(Refseq$end[index])
    }
  }
  
  
  
  # Peak-gene looping
  Candidate_Peak_Gene_looping_distance_control=sparseMatrix(i=1,j=1,x=0, dims=c(nrow(Candidate_peaks), nrow(Candidate_genes)), repr = "T")
  for (g in 1:length(Candidate_genes$gene_symbols)){
    index=which(Candidate_peaks$chr==Candidate_genes$chr[g] & abs((Candidate_peaks$point1+Candidate_peaks$point2)/2-Candidate_genes$TSS[g])<distance_control)
    if (length(index)>0){
      Candidate_Peak_Gene_looping_distance_control[index,g]=1
    }
  }
  Candidate_Peak_Gene_looping<-Candidate_Peak_Gene_looping_distance_control
  peak_index<-which(rowSums(Candidate_Peak_Gene_looping)>0)
  gene_index<-which(colSums(Candidate_Peak_Gene_looping)>0)
  Candidate_Peak_Gene_looping<-Candidate_Peak_Gene_looping[peak_index,gene_index]
  Candidate_peaks<-Candidate_peaks[peak_index,]
  Candidate_genes<-Candidate_genes[gene_index,]
  TF_index=which(colSums(Candidate_TF_Peak_Binding[peak_index,])>0)
  Candidate_TFs=Candidate_TFs[TF_index,]
  Candidate_TF_Peak_Binding=Candidate_TF_Peak_Binding[peak_index, TF_index]
  
  
  
  
  #pseudobulk data is calculated for model initialization
  ATAC_count=matrix(0,nrow=nrow(scATAC_read_count_matrix), ncol=length(Common_samples))
  for (s in 1:length(Common_samples)){
    ATAC_count[,s]=rowSums(scATAC_read_count_matrix[,which(scATAC_cells$subject_ID==Common_samples[s])])
  }
  colnames(ATAC_count)<-Common_samples
  rownames(ATAC_count)<-scATAC_peaks$peak_index
  
  RNA_count=matrix(0,nrow=nrow(scRNA_read_count_matrix), ncol=length(Common_samples))
  for (s in 1:length(Common_samples)){
    RNA_count[,s]=rowSums(scRNA_read_count_matrix[,which(scRNA_cells$subject_ID==Common_samples[s])])
  }
  colnames(RNA_count)<-Common_samples
  rownames(RNA_count)<-scRNA_genes$gene_symbols
  
  #normalization
  total_ATAC_reads=5e6
  ATAC_scale_facter=total_ATAC_reads/(colSums(ATAC_count)+1)
  ATAC_count=sweep(ATAC_count, 2, ATAC_scale_facter, "*")
  ATAC_log2=log2(ATAC_count+1)
  
  total_RNA_reads=5e6
  RNA_scale_facter=total_RNA_reads/(colSums(RNA_count)+1)
  RNA_count=sweep(RNA_count, 2, RNA_scale_facter, "*")
  RNA_log2=log2(RNA_count+1)
  
  
  # select actively accessbile peaks
  index=match(Candidate_peaks$peak_index, scATAC_peaks$peak_index)
  Candidate_Peak_log2Count=ATAC_log2[index,]
  
  peak_index=which(rowSums(Candidate_Peak_log2Count)>0)
  Candidate_peaks=Candidate_peaks[peak_index,]
  Candidate_TF_Peak_Binding=Candidate_TF_Peak_Binding[peak_index,]
  Candidate_Peak_Gene_looping=Candidate_Peak_Gene_looping[peak_index,]
  Candidate_Peak_log2Count=Candidate_Peak_log2Count[peak_index,]
  Candidate_Peak_log2Count=sweep(Candidate_Peak_log2Count, 1, rowMeans(Candidate_Peak_log2Count), '-')
  
  
  # select actively expressed genes
  index=match(Candidate_genes$gene_symbols, scRNA_genes$gene_symbols)
  Candidate_Gene_log2Count=RNA_log2[index,]
  
  gene_index=which(rowSums(Candidate_Gene_log2Count)>0)
  Candidate_genes=Candidate_genes[gene_index,]
  Candidate_Peak_Gene_looping=Candidate_Peak_Gene_looping[,gene_index]
  Candidate_Gene_log2Count=Candidate_Gene_log2Count[gene_index,]
  Candidate_Gene_log2Count=sweep(Candidate_Gene_log2Count, 1, rowMeans(Candidate_Gene_log2Count), '-')
  
  
  # select actively expressed TFs
  selected_TFs=intersect(Candidate_TFs$name, scRNA_genes$gene_symbols)
  index=match(selected_TFs, Candidate_TFs$name)
  Candidate_TFs=Candidate_TFs[index,]
  
  index=match(selected_TFs, scRNA_genes$gene_symbols)
  Candidate_TF_Peak_Binding=Candidate_TF_Peak_Binding[,index]
  Candidate_TF_log2Count=RNA_log2[index,]
  
  TF_index=which(rowSums(Candidate_TF_log2Count)>0)
  Candidate_TFs=Candidate_TFs[TF_index,]
  Candidate_TF_Peak_Binding=Candidate_TF_Peak_Binding[,TF_index]
  Candidate_TF_log2Count=Candidate_TF_log2Count[TF_index,]
  Candidate_TF_log2Count=sweep(Candidate_TF_log2Count, 1, rowMeans(Candidate_TF_log2Count), '-')
  
  
  
  # final filtering for candidate circuits
  peak_index=which(rowSums(Candidate_Peak_Gene_looping)>0 & rowSums(Candidate_TF_Peak_Binding)>0)
  gene_index=which(colSums(Candidate_Peak_Gene_looping[peak_index,])>0)
  TF_index=which(colSums(Candidate_TF_Peak_Binding[peak_index,])>0)
  
  Candidate_TF_Peak_Binding=Candidate_TF_Peak_Binding[peak_index,TF_index]
  Candidate_Peak_Gene_looping=Candidate_Peak_Gene_looping[peak_index,gene_index]
  
  Candidate_TFs=Candidate_TFs[TF_index,]
  Candidate_TF_log2Count=Candidate_TF_log2Count[TF_index,]
  
  Candidate_peaks=Candidate_peaks[peak_index,]
  Candidate_Peak_log2Count=Candidate_Peak_log2Count[peak_index,]  
  
  Candidate_genes=Candidate_genes[gene_index,]
  Candidate_Gene_log2Count=Candidate_Gene_log2Count[gene_index,]
  
  Candidate_circuits=list('TFs'=Candidate_TFs, 
                          'TF_log2Count'=Candidate_TF_log2Count,
                          'Peaks'=Candidate_peaks, 
                          'Peak_log2Count'=Candidate_Peak_log2Count,
                          'Genes'=Candidate_genes, 
                          'Gene_log2Count'=Candidate_Gene_log2Count,
                          'TF_Peak_Binding'=Candidate_TF_Peak_Binding, 
                          'Peak_Gene_looping'=Candidate_Peak_Gene_looping)
  return(Candidate_circuits)
}

