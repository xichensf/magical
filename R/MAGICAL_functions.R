#**************************************************************************
#                         MAGICAL functions
#***************************************************************************



#***************************load all input files***************************
Data_loading <- function(Candidate_Gene_file_path, Candidate_Peak_file_path,
                         scRNA_readcount_file_path, scRNA_Gene_file_path, scRNA_cellmeta_file_path,
                         scATAC_readcount_file_path, scATAC_Peak_file_path, scATAC_cellmeta_file_path,
                         Motif_mapping_file_path, Motif_name_file_path, Ref_seq_file_path){
  
  #load candidate Genes
  Candidate_Genes <- read.table(Candidate_Gene_file_path, sep='\t')
  colnames(Candidate_Genes) = c("Gene_symbols") 
  print('loading candidate genes ...')
  cat('\n')
  print(paste(length(Candidate_Genes$Gene_symbols), 'candidate genes are provided.', sep = ' '))
  cat('\n')
  
  
  #load scRNAseq data
  print('loading scRNAseq data ...')
  cat('\n')
  
  scRNA_Genes<-read.table(scRNA_Gene_file_path, sep='\t')
  colnames(scRNA_Genes) = c("Gene_index", "Gene_symbols")
  
  scRNA_cells<-read.table(scRNA_cellmeta_file_path, sep='\t')
  colnames(scRNA_cells) = c("cell_index", "cell_barcode", "cell_type", "subject_ID", "condition")
  
  scRNA_read_count_table<-read.table(scRNA_readcount_file_path, sep='\t')
  scRNA_read_count_matrix <- sparseMatrix(i = scRNA_read_count_table[,1], j=scRNA_read_count_table[,2], x=scRNA_read_count_table[,3],
                                          dimnames=list(levels(scRNA_read_count_table[,1]), levels(scRNA_read_count_table[,2])), repr = "T")
  rm(scRNA_read_count_table)
  
  Group_conditions = unique(scRNA_cells$condition)
  print(paste('We detected', length(Group_conditions), 'conditions from the meta file', sep = ' '))
  cat('\n')
  
  scRNA_samples = unique(scRNA_cells$subject_ID)
  print(paste('The input scRNAseq data includes', length(scRNA_Genes$Gene_symbols),'Genes and', length(scRNA_cells$cell_barcode),'cells from',length(scRNA_samples),'samples', sep = ' '))
  cat('\n')
  
  #print('RNA cells of each sample')
  #table(scRNA_cells$subject_ID)
  
  
  #load candidate Peaks
  print('loading candidate peaks ...')
  cat('\n')
  Candidate_Peaks <- read.table(Candidate_Peak_file_path, sep='\t')
  colnames(Candidate_Peaks) = c("chr", "point1", 'point2')
  print(paste(nrow(Candidate_Peaks), 'candidate peaks are provided.', sep = ' '))
  cat('\n')
  
  
  #load scATACseq data
  print('loading scATACseq data ...')
  cat('\n')
  scATAC_Peaks <- read.table(scATAC_Peak_file_path, sep='\t')
  scATAC_Peaks <- scATAC_Peaks[,c(1,2,3,4)]
  colnames(scATAC_Peaks) = c("Peak_index", "chr", "point1", 'point2')
  
  scATAC_cells <- read.table(scATAC_cellmeta_file_path, sep ="\t")
  colnames(scATAC_cells) = c("cell_index", "cell_barcode", "cell_type", 'subject_ID')
  
  scATAC_read_count_table <- read.table(scATAC_readcount_file_path, sep ="\t")
  scATAC_read_count_matrix <- sparseMatrix(i = scATAC_read_count_table[,1], j=scATAC_read_count_table[,2], x=scATAC_read_count_table[,3],
                                           dimnames=list(levels(scATAC_read_count_table[,1]), levels(scATAC_read_count_table[,2])), repr = "T")
  rm(scATAC_read_count_table)
  
  scATAC_samples = unique(scATAC_cells$subject_ID)
  print(paste('The input scATACseq data includes',length(scATAC_Peaks$Peak_index), 'Peaks and', length(scATAC_cells$cell_barcode), 'cells from', length(scATAC_samples), 'samples', sep = ' '))
  cat('\n')

  #print('ATAC cells of each sample')
  #table(scATAC_cells$subject_ID)
  
  
  Common_samples=intersect(scRNA_samples, scATAC_samples)
  print(paste('There are sample-paired single cell multiomics data for', length(Common_samples), 'samples.', sep = ' '))
  print('Please check sample IDs if this number is lower than expected.')
  cat('\n')

  #load TF motif prior
  print('loading TF motif match data ...')
  cat('\n')
  Motifs<-read.table(Motif_name_file_path, sep ="\t")
  colnames(Motifs) = c("motif_index", "name")
  Motif_mapping_table <- read.table(Motif_mapping_file_path)
  TF_Peak_binding_matrix <- sparseMatrix(i = Motif_mapping_table[,1], j=Motif_mapping_table[,2], x=Motif_mapping_table[,3],
                                         dimnames=list(levels(Motif_mapping_table[,1]), levels(Motif_mapping_table[,2])), repr = "T")
  rm(Motif_mapping_table)
  print(paste(nrow(Motifs), 'motifs are screened.', sep = ' '))
  cat('\n')
  
  
  #load RefSeq
  Refseq <- read.table(Ref_seq_file_path, header = TRUE, sep='\t')
  colnames(Refseq) = c("chr", "strand", "start", 'end', 'Gene_symbols')
  
  
  loaded_data = list('Common_samples'=Common_samples, 
                     'Candidate_Genes'=Candidate_Genes, 
                     'Candidate_Peaks'=Candidate_Peaks,
                     'scRNA_Genes'=scRNA_Genes,
                     'scRNA_cells'=scRNA_cells,
                     'scRNA_read_count_matrix'=scRNA_read_count_matrix,
                     'scATAC_Peaks'=scATAC_Peaks, 
                     'scATAC_cells'=scATAC_cells,
                     'scATAC_read_count_matrix'=scATAC_read_count_matrix,
                     'Motifs'=Motifs, 
                     'TF_Peak_binding_matrix'=TF_Peak_binding_matrix,
                     'Refseq'=Refseq)
  
  return(loaded_data)
}




#***************************Build candidate circuits***************************
Candidate_circuits_construction_with_TAD <- function(loaded_data, TAD_file_path){
  
  Common_samples=loaded_data$Common_samples
  Candidate_Genes=loaded_data$Candidate_Genes
  Candidate_Peaks=loaded_data$Candidate_Peaks
  scRNA_Genes=loaded_data$scRNA_Genes 
  scRNA_cells=loaded_data$scRNA_cells 
  scRNA_read_count_matrix=loaded_data$scRNA_read_count_matrix
  scATAC_Peaks=loaded_data$scATAC_Peaks 
  scATAC_cells=loaded_data$scATAC_cells 
  scATAC_read_count_matrix=loaded_data$scATAC_read_count_matrix
  Motifs=loaded_data$Motifs 
  TF_Peak_binding_matrix=loaded_data$TF_Peak_binding_matrix
  Refseq=loaded_data$Refseq
  
  #**********load TAD prior**********
  print('loading TAD data ...')
  cat('\n')
  TAD <- read.table(TAD_file_path, sep='\t')
  colnames(TAD) = c("chr", "left_boundary", 'right_boundary')  
  print(paste(nrow(TAD), 'TAD segments provided', sep = ' '))
  cat('\n')
  print('MAGICAL model initialization ...')
  cat('\n')
  
  # TF-Peak binding
  library(dplyr)
  Candidate_Peaks<-inner_join(scATAC_Peaks, Candidate_Peaks, c("chr", "point1", "point2"))
  index<-match(Candidate_Peaks$Peak_index, scATAC_Peaks$Peak_index)
  Candidate_TF_Peak_Binding<-TF_Peak_binding_matrix[index,]
  TF_num<-colSums(Candidate_TF_Peak_Binding)
  TF_pct_in_candidate_Peaks<-colSums(Candidate_TF_Peak_Binding)/nrow(Candidate_Peaks)
  TF_pct_in_all_Peaks<-colSums(TF_Peak_binding_matrix)/nrow(scATAC_Peaks)
  TF_enrichment_FC<-TF_pct_in_candidate_Peaks/TF_pct_in_all_Peaks
  TF_index<-which(TF_pct_in_candidate_Peaks>0.05 & TF_num>30 & TF_enrichment_FC>1.2)
  if (length(TF_index) == 0){
    print('Too few Peaks with TF binding sites.')
    print('MAGICAL not applicable to this cell type!')
  }else{
    Candidate_TFs <- Motifs[TF_index,]
    Candidate_TF_Peak_Binding <- Candidate_TF_Peak_Binding[,TF_index]
  }
  Peak_index<-which(rowSums(Candidate_TF_Peak_Binding)>0)
  Candidate_Peaks <- Candidate_Peaks[Peak_index,]
  Candidate_TF_Peak_Binding <- Candidate_TF_Peak_Binding[Peak_index,]
  
  
  # extract Gene TSS
  Gene_symbols<-intersect(intersect(Candidate_Genes$Gene_symbols, Refseq$Gene_symbols), scRNA_Genes$Gene_symbols)
  Gene_chr=character(length(Gene_symbols))
  Gene_TSS=integer(length(Gene_symbols))
  Candidate_Genes=data.frame(Gene_symbols, Gene_chr, Gene_TSS)
  colnames(Candidate_Genes) = c("Gene_symbols", "chr", "TSS")
  for (g in 1:length(Gene_symbols)){
    index = which(Refseq$Gene_symbols==Candidate_Genes$Gene_symbols[g])
    Candidate_Genes$chr[g] <- Refseq$chr[index[1]]
    if (Refseq$strand[index[1]]=='+'){
      Candidate_Genes$TSS[g] <- min(Refseq$start[index])
    }else{
      Candidate_Genes$TSS[g] <- max(Refseq$end[index])
    }
  }

  
  # Peak-Gene looping
  Candidate_Peak_Gene_looping_TAD=sparseMatrix(i=1,j=1,x=0, dims=c(nrow(Candidate_Peaks), nrow(Candidate_Genes)), repr = "T")
  for (t in 1:nrow(TAD)){
      Peak_index=which(Candidate_Peaks$chr==TAD$chr[t] & Candidate_Peaks$point1>TAD$left_boundary[t] & Candidate_Peaks$point2<TAD$right_boundary[t])
      Gene_index=which(Candidate_Genes$chr==TAD$chr[t] & Candidate_Genes$TSS>TAD$left_boundary[t] & Candidate_Genes$TSS<TAD$right_boundary[t])
      if (length(Peak_index)>0 && length(Gene_index)>0){
        Candidate_Peak_Gene_looping_TAD[Peak_index,Gene_index]=1
      }
  }
  
  
  #some TAD could be extra wide, based on HiC interaction scale, we limit the Peak-Gene distance upto 1M bps
  Candidate_Peak_Gene_looping_distance_control=sparseMatrix(i=1,j=1,x=0, dims=c(nrow(Candidate_Peaks), nrow(Candidate_Genes)), repr = "T")
  for (g in 1:length(Candidate_Genes$Gene_symbols)){
    index=which(Candidate_Peaks$chr==Candidate_Genes$chr[g] & abs((Candidate_Peaks$point1+Candidate_Peaks$point2)/2-Candidate_Genes$TSS[g])<1e6)
    if (length(index)>0){
      Candidate_Peak_Gene_looping_distance_control[index,g]=1
    }
  }
  
  Candidate_Peak_Gene_looping<-as(Candidate_Peak_Gene_looping_TAD*Candidate_Peak_Gene_looping_distance_control, "dgTMatrix")
  
  
  Peak_index<-which(rowSums(Candidate_Peak_Gene_looping)>0)
  Gene_index<-which(colSums(Candidate_Peak_Gene_looping)>0)
  Candidate_Peak_Gene_looping<-Candidate_Peak_Gene_looping[Peak_index,Gene_index]
  Candidate_Peaks<-Candidate_Peaks[Peak_index,]
  Candidate_Genes<-Candidate_Genes[Gene_index,]
  TF_index=which(colSums(Candidate_TF_Peak_Binding[Peak_index,])>0)
  Candidate_TFs=Candidate_TFs[TF_index,]
  Candidate_TF_Peak_Binding=Candidate_TF_Peak_Binding[Peak_index, TF_index]
  
   
  #pseudobulk ATAC data is calculated for model initialization
  ATAC_count=matrix(0,nrow=nrow(scATAC_read_count_matrix), ncol=length(Common_samples))
  for (s in 1:length(Common_samples)){
    ATAC_count[,s]=rowSums(scATAC_read_count_matrix[,which(scATAC_cells$subject_ID==Common_samples[s])])
  }
  colnames(ATAC_count)<-Common_samples
  rownames(ATAC_count)<-scATAC_Peaks$Peak_index
 
  total_ATAC_reads=5e6
  ATAC_scale_facter=total_ATAC_reads/(colSums(ATAC_count)+1)
  ATAC_count=sweep(ATAC_count, 2, ATAC_scale_facter, "*")
  ATAC_log2=log2(ATAC_count+1)
  
  index=match(Candidate_Peaks$Peak_index, scATAC_Peaks$Peak_index)
  Candidate_Peak_log2Count=ATAC_log2[index,]
  Candidate_Peak_log2Count=sweep(Candidate_Peak_log2Count, 1, rowMeans(Candidate_Peak_log2Count), '-')
  
  
  #pseudobulk RNA data is calculated for model initialization
  RNA_count=matrix(0,nrow=nrow(scRNA_read_count_matrix), ncol=length(Common_samples))
  for (s in 1:length(Common_samples)){
    RNA_count[,s]=rowSums(scRNA_read_count_matrix[,which(scRNA_cells$subject_ID==Common_samples[s])])
  }
  colnames(RNA_count)<-Common_samples
  rownames(RNA_count)<-scRNA_Genes$Gene_symbols
  
  total_RNA_reads=5e6
  RNA_scale_facter=total_RNA_reads/(colSums(RNA_count)+1)
  RNA_count=sweep(RNA_count, 2, RNA_scale_facter, "*")
  RNA_log2=log2(RNA_count+1)
  
  index=match(Candidate_Genes$Gene_symbols, scRNA_Genes$Gene_symbols)
  Candidate_Gene_log2Count=RNA_log2[index,]
  Candidate_Gene_log2Count=sweep(Candidate_Gene_log2Count, 1, rowMeans(Candidate_Gene_log2Count), '-')
  
  
  # TFs that are not expressed in the current data are not considered
  selected_TFs=intersect(Candidate_TFs$name, scRNA_Genes$Gene_symbols)
  
  index=match(selected_TFs, Candidate_TFs$name)
  Candidate_TFs=Candidate_TFs[index,]
  Candidate_TF_Peak_Binding=Candidate_TF_Peak_Binding[,index]
  
  index=match(selected_TFs, scRNA_Genes$Gene_symbols)
  Candidate_TF_log2Count=RNA_log2[index,]
  
  TF_index=which(rowSums(Candidate_TF_log2Count)>0)
  Candidate_TFs=Candidate_TFs[TF_index,]
  Candidate_TF_Peak_Binding=Candidate_TF_Peak_Binding[,TF_index]
  Candidate_TF_log2Count=Candidate_TF_log2Count[TF_index,]
  Candidate_TF_log2Count=sweep(Candidate_TF_log2Count, 1, rowMeans(Candidate_TF_log2Count), '-')
  
  
  # final filtering for candidate circuits
  Peak_index=which(rowSums(Candidate_Peak_Gene_looping)>0 & rowSums(Candidate_TF_Peak_Binding)>0)
  Gene_index=which(colSums(Candidate_Peak_Gene_looping[Peak_index,])>0)
  TF_index=which(colSums(Candidate_TF_Peak_Binding[Peak_index,])>0)
  
  Candidate_TF_Peak_Binding=Candidate_TF_Peak_Binding[Peak_index,TF_index]
  Candidate_Peak_Gene_looping=Candidate_Peak_Gene_looping[Peak_index,Gene_index]
  
  Candidate_TFs=Candidate_TFs[TF_index,]
  Candidate_TF_log2Count=Candidate_TF_log2Count[TF_index,]
  
  Candidate_Peaks=Candidate_Peaks[Peak_index,]
  Candidate_Peak_log2Count=Candidate_Peak_log2Count[Peak_index,]  
  
  Candidate_Genes=Candidate_Genes[Gene_index,]
  Candidate_Gene_log2Count=Candidate_Gene_log2Count[Gene_index,]
  
  Candidate_circuits=list('TFs'=Candidate_TFs, 
                          'TF_log2Count'=Candidate_TF_log2Count,
                          'Peaks'=Candidate_Peaks, 
                          'Peak_log2Count'=Candidate_Peak_log2Count,
                          'Genes'=Candidate_Genes, 
                          'Gene_log2Count'=Candidate_Gene_log2Count,
                          'TF_Peak_Binding'=Candidate_TF_Peak_Binding, 
                          'Peak_Gene_looping'=Candidate_Peak_Gene_looping)
  
  return(Candidate_circuits)
  
}




#***************************Build candidate circuits***************************
Candidate_circuits_construction_without_TAD <- function(loaded_data, distance_control){

  Common_samples=loaded_data$Common_samples
  Candidate_Genes=loaded_data$Candidate_Genes
  Candidate_Peaks=loaded_data$Candidate_Peaks
  scRNA_Genes=loaded_data$scRNA_Genes 
  scRNA_cells=loaded_data$scRNA_cells 
  scRNA_read_count_matrix=loaded_data$scRNA_read_count_matrix
  scATAC_Peaks=loaded_data$scATAC_Peaks 
  scATAC_cells=loaded_data$scATAC_cells 
  scATAC_read_count_matrix=loaded_data$scATAC_read_count_matrix
  Motifs=loaded_data$Motifs 
  TF_Peak_binding_matrix=loaded_data$TF_Peak_binding_matrix
  Refseq=loaded_data$Refseq
  
  #**********load TAD prior**********
  print(paste('Naive pairing of peaks and genes if they are within', distance_control,'bps', sep = ' '))
  cat('\n')
  print('MAGICAL model initialization ...')
  cat('\n')
  
  # TF-Peak binding
  Candidate_Peaks<-inner_join(scATAC_Peaks, Candidate_Peaks)
  index<-match(Candidate_Peaks$Peak_index, scATAC_Peaks$Peak_index)
  Candidate_TF_Peak_Binding<-TF_Peak_binding_matrix[index,]
  TF_num<-colSums(Candidate_TF_Peak_Binding)
  TF_pct_in_candidate_Peaks<-colSums(Candidate_TF_Peak_Binding)/nrow(Candidate_Peaks)
  TF_pct_in_all_Peaks<-colSums(TF_Peak_binding_matrix)/nrow(scATAC_Peaks)
  TF_enrichment_FC<-TF_pct_in_candidate_Peaks/TF_pct_in_all_Peaks
  TF_index<-which(TF_pct_in_candidate_Peaks>0.05 & TF_num>30 & TF_enrichment_FC>1.2)
  if (length(TF_index) == 0){
    print('Too few Peaks with TF binding sites.')
    print('MAGICAL not applicable to this cell type!')
  }else{
    Candidate_TFs <- Motifs[TF_index,]
    Candidate_TF_Peak_Binding <- Candidate_TF_Peak_Binding[,TF_index]
  }
  Peak_index<-which(rowSums(Candidate_TF_Peak_Binding)>0)
  Candidate_Peaks <- Candidate_Peaks[Peak_index,]
  Candidate_TF_Peak_Binding <- Candidate_TF_Peak_Binding[Peak_index,]
  
  
  # extract Gene TSS
  Gene_symbols<-intersect(intersect(Candidate_Genes$Gene_symbols, Refseq$Gene_symbols), scRNA_Genes$Gene_symbols)
  Gene_chr=character(length(Gene_symbols))
  Gene_TSS=integer(length(Gene_symbols))
  Candidate_Genes=data.frame(Gene_symbols, Gene_chr, Gene_TSS)
  colnames(Candidate_Genes) = c("Gene_symbols", "chr", "TSS")
  for (g in 1:length(Gene_symbols)){
    index = which(Refseq$Gene_symbols==Candidate_Genes$Gene_symbols[g])
    Candidate_Genes$chr[g] <- Refseq$chr[index[1]]
    if (Refseq$strand[index[1]]=='+'){
      Candidate_Genes$TSS[g] <- min(Refseq$start[index])
    }else{
      Candidate_Genes$TSS[g] <- max(Refseq$end[index])
    }
  }
  
  
  # Peak-Gene looping
  Candidate_Peak_Gene_looping_distance_control=sparseMatrix(i=1,j=1,x=0, dims=c(nrow(Candidate_Peaks), nrow(Candidate_Genes)), repr = "T")
  for (g in 1:length(Candidate_Genes$Gene_symbols)){
    index=which(Candidate_Peaks$chr==Candidate_Genes$chr[g] & abs((Candidate_Peaks$point1+Candidate_Peaks$point2)/2-Candidate_Genes$TSS[g])<distance_control)
    if (length(index)>0){
      Candidate_Peak_Gene_looping_distance_control[index,g]=1
    }
  }
  Candidate_Peak_Gene_looping<-Candidate_Peak_Gene_looping_distance_control
  Peak_index<-which(rowSums(Candidate_Peak_Gene_looping)>0)
  Gene_index<-which(colSums(Candidate_Peak_Gene_looping)>0)
  Candidate_Peak_Gene_looping<-Candidate_Peak_Gene_looping[Peak_index,Gene_index]
  Candidate_Peaks<-Candidate_Peaks[Peak_index,]
  Candidate_Genes<-Candidate_Genes[Gene_index,]
  TF_index=which(colSums(Candidate_TF_Peak_Binding[Peak_index,])>0)
  Candidate_TFs=Candidate_TFs[TF_index,]
  Candidate_TF_Peak_Binding=Candidate_TF_Peak_Binding[Peak_index, TF_index]
  
  
  #pseudobulk ATAC data is calculated for model initialization
  ATAC_count=matrix(0,nrow=nrow(scATAC_read_count_matrix), ncol=length(Common_samples))
  for (s in 1:length(Common_samples)){
    ATAC_count[,s]=rowSums(scATAC_read_count_matrix[,which(scATAC_cells$subject_ID==Common_samples[s])])
  }
  colnames(ATAC_count)<-Common_samples
  rownames(ATAC_count)<-scATAC_Peaks$Peak_index
  
  total_ATAC_reads=5e6
  ATAC_scale_facter=total_ATAC_reads/(colSums(ATAC_count)+1)
  ATAC_count=sweep(ATAC_count, 2, ATAC_scale_facter, "*")
  ATAC_log2=log2(ATAC_count+1)
  
  index=match(Candidate_Peaks$Peak_index, scATAC_Peaks$Peak_index)
  Candidate_Peak_log2Count=ATAC_log2[index,]
  Candidate_Peak_log2Count=sweep(Candidate_Peak_log2Count, 1, rowMeans(Candidate_Peak_log2Count), '-')
  
  
  #pseudobulk RNA data is calculated for model initialization
  RNA_count=matrix(0,nrow=nrow(scRNA_read_count_matrix), ncol=length(Common_samples))
  for (s in 1:length(Common_samples)){
    RNA_count[,s]=rowSums(scRNA_read_count_matrix[,which(scRNA_cells$subject_ID==Common_samples[s])])
  }
  colnames(RNA_count)<-Common_samples
  rownames(RNA_count)<-scRNA_Genes$Gene_symbols
  
  total_RNA_reads=5e6
  RNA_scale_facter=total_RNA_reads/(colSums(RNA_count)+1)
  RNA_count=sweep(RNA_count, 2, RNA_scale_facter, "*")
  RNA_log2=log2(RNA_count+1)
  
  index=match(Candidate_Genes$Gene_symbols, scRNA_Genes$Gene_symbols)
  Candidate_Gene_log2Count=RNA_log2[index,]
  Candidate_Gene_log2Count=sweep(Candidate_Gene_log2Count, 1, rowMeans(Candidate_Gene_log2Count), '-')
  
  
  # TFs that are not expressed in the current data are not considered
  selected_TFs=intersect(Candidate_TFs$name, scRNA_Genes$Gene_symbols)
  index=match(selected_TFs, Candidate_TFs$name)
  Candidate_TFs=Candidate_TFs[index,]
  Candidate_TF_Peak_Binding=Candidate_TF_Peak_Binding[,index]
  
  index=match(selected_TFs, scRNA_Genes$Gene_symbols)
  Candidate_TF_log2Count=RNA_log2[index,]
  
  TF_index=which(rowSums(Candidate_TF_log2Count)>0)
  Candidate_TFs=Candidate_TFs[TF_index,]
  Candidate_TF_Peak_Binding=Candidate_TF_Peak_Binding[,TF_index]
  Candidate_TF_log2Count=Candidate_TF_log2Count[TF_index,]
  Candidate_TF_log2Count=sweep(Candidate_TF_log2Count, 1, rowMeans(Candidate_TF_log2Count), '-')
  
  
  # final filtering for candidate circuits
  Peak_index=which(rowSums(Candidate_Peak_Gene_looping)>0 & rowSums(Candidate_TF_Peak_Binding)>0)
  Gene_index=which(colSums(Candidate_Peak_Gene_looping[Peak_index,])>0)
  TF_index=which(colSums(Candidate_TF_Peak_Binding[Peak_index,])>0)
  
  Candidate_TF_Peak_Binding=Candidate_TF_Peak_Binding[Peak_index,TF_index]
  Candidate_Peak_Gene_looping=Candidate_Peak_Gene_looping[Peak_index,Gene_index]
  
  Candidate_TFs=Candidate_TFs[TF_index,]
  Candidate_TF_log2Count=Candidate_TF_log2Count[TF_index,]
  
  Candidate_Peaks=Candidate_Peaks[Peak_index,]
  Candidate_Peak_log2Count=Candidate_Peak_log2Count[Peak_index,]  
  
  Candidate_Genes=Candidate_Genes[Gene_index,]
  Candidate_Gene_log2Count=Candidate_Gene_log2Count[Gene_index,]
  
  Candidate_circuits=list('TFs'=Candidate_TFs, 
                          'TF_log2Count'=Candidate_TF_log2Count,
                          'Peaks'=Candidate_Peaks, 
                          'Peak_log2Count'=Candidate_Peak_log2Count,
                          'Genes'=Candidate_Genes, 
                          'Gene_log2Count'=Candidate_Gene_log2Count,
                          'TF_Peak_Binding'=Candidate_TF_Peak_Binding, 
                          'Peak_Gene_looping'=Candidate_Peak_Gene_looping)
  
  return(Candidate_circuits)
}




#*************************** MCMC model initialization ***************************
#for MCMC method, the initial values of parameters are not deterministic for the results 
#use pseudo bulk data to initialize the model as only at this dimension all data are matched

MAGICAL_initialization <- function(loaded_data, Candidate_circuits){
  
  Candidate_TFs <- Candidate_circuits$TFs
  Candidate_TF_log2Count <- Candidate_circuits$TF_log2Count 
  
  Candidate_Peaks <- Candidate_circuits$Peaks
  Candidate_Peak_log2Count <- Candidate_circuits$Peak_log2Count
  
  Candidate_Genes <- Candidate_circuits$Genes
  Candidate_Gene_log2Count <- Candidate_circuits$Gene_log2Count 
  
  Candidate_TF_Peak_Binding <- Candidate_circuits$TF_Peak_Binding 
  Candidate_Peak_Gene_looping <- Candidate_circuits$Peak_Gene_looping
  
  
  S <- length(loaded_data$Common_samples)
  M <- nrow(Candidate_TFs)#total number of candidate TFs
  P <- nrow(Candidate_Peaks)#total number of candidate Peaks
  G <- nrow(Candidate_Genes)#total number of candidate Genes

  
  # TF activity prior: TF RNA expression
  T_prior <- Candidate_TF_log2Count
  T_mean <- Candidate_TF_log2Count
  T_var = matrix(0, nrow=M, ncol=S)
  for (m in 1:M){
    T_var[m,] <- var(Candidate_TF_log2Count[m,])
  }
  
  
  # TF-Peak binding prior: regression weight between TF activity (expression) and Peak ATAC activity
  B_prior = as.matrix(Candidate_TF_Peak_Binding)
  B_prob = as.matrix(Candidate_TF_Peak_Binding)
  TF_Peak_index <- which(as.matrix(Candidate_TF_Peak_Binding)>0, arr.ind = TRUE)
  
  for (b in 1:length(Candidate_TF_Peak_Binding@x)){
    mdl = summary(lm(Candidate_Peak_log2Count[TF_Peak_index[b,1],]~Candidate_TF_log2Count[TF_Peak_index[b,2],]))
    B_prior[TF_Peak_index[b,1], TF_Peak_index[b,2]] = mdl$coefficients[2,1]
    B_prob[TF_Peak_index[b,1], TF_Peak_index[b,2]] = 1-mdl$coefficients[2,4]
  }
  B_mean = B_prior
  
  B_var = matrix(0.5, nrow=1, ncol=M)
  
  for (m in 1:M){
    if(sum(Candidate_TF_Peak_Binding[,m])>1){
      B_var[m]=var(B_prior[which(Candidate_TF_Peak_Binding[,m]>0),m])#variance of binding events weights of TF t
    }else{
      B_var[m]=0.5
    }
   }
  
  
  # Peak-Gene looping prior
  L_prior=as.matrix(Candidate_Peak_Gene_looping)
  L_prob=as.matrix(Candidate_Peak_Gene_looping)
  Peak_Gene_index <- which(as.matrix(Candidate_Peak_Gene_looping)>0, arr.ind = TRUE)
  
  for (l in 1:nrow(Peak_Gene_index)){
    mdl=summary(lm(Candidate_Gene_log2Count[Peak_Gene_index[l,2],]~Candidate_Peak_log2Count[Peak_Gene_index[l,1],]))
    L_prior[Peak_Gene_index[l,1], Peak_Gene_index[l,2]] = mdl$coefficients[2,1]
    L_prob[Peak_Gene_index[l,1], Peak_Gene_index[l,2]] = 1-mdl$coefficients[2,4]
  }
  
  L_mean=L_prior
  L_var=var(L_prior[which(L_prob>0)])
  
  Initial_model=list('T_prior'=T_prior, 
                     'T_mean'=T_mean, 
                     'T_var'=T_var, 
                     'B_prior'=B_prior, 
                     'B_mean'=B_mean, 
                     'B_var'=B_var, 
                     'B_prob'=B_prob, 
                     'L_prior'=L_prior, 
                     'L_mean'=L_mean, 
                     'L_var'=L_var, 
                     'L_prob'=L_prob)
  
  return(Initial_model)
}




#*************************** TFA estimation ***************************
TF_activity_T_sampling <- function(A, A_sample, ATAC_Cell_Sample_vector, R, R_sample, RNA_Cell_Sample_vector, TFA, T_prior_mean, T_prior_var, B, B_state, sigma_A_noise, P, G, M, S){
  
  TF_index=sample(M)
  for (i in 1:M){
    #during the looping, TFA is iteratively updated for all M TFs and the update order is random in each round of sampling
    m=TF_index[i]
    if (sum(B_state[,m])>0){
      
      temp_var=sum(B[,m]^2)*T_prior_var[m,]/sum(B_state[,m])+sigma_A_noise
      mean_T=(t(B[,m])%*%(A_sample-B%*%TFA$T_sample+B[,m]%*%t(TFA$T_sample[m,]))/sum(B_state[,m])+T_prior_mean[m,]*sigma_A_noise)/temp_var
      variance_T=T_prior_var[m,]*sigma_A_noise/temp_var
      
      for (s in 1:S){
        aa=rnorm(1)
        if (aa-3>0){
          aa=3
        }
        if (aa+3<0){
          aa=-3
        }
        
        TFA$T_sample[m,s]=aa*sqrt(abs(variance_T[s]))+mean_T[s]
        
        ATAC_cell_index=which(ATAC_Cell_Sample_vector==s)
        TFA$T_A[m,ATAC_cell_index]=rnorm(length(ATAC_cell_index))*sqrt(abs(variance_T[s]))+TFA$T_sample[m,s]
        
        RNA_cell_index=which(RNA_Cell_Sample_vector==s)
        TFA$T_R[m,RNA_cell_index]=rnorm(length(RNA_cell_index))*sqrt(abs(variance_T[s]))+TFA$T_sample[m,s]
      }
    }
  }
  return(TFA)
  
}





#*************************** TF-peak binding confidence estimation ***************************
TF_peak_binding_B_sampling <- function(A, A_sample, ATAC_Cell_Sample_vector, TFA, B, B_state, B_prior_mean, B_prior_var, sigma_A_noise, P, G, M, S){

  T_sample=matrix(0, nrow=M, ncol=S)
  for (s in 1:S){
    T_sample[,s]=rowMeans(TFA$T_A[,which(ATAC_Cell_Sample_vector==s)])
  }
  
  TF_index=sample(M)
  
  for (i in 1:M){
    
    m=TF_index[i]
    temp_var=sum(TFA$T_sample[m,]^2)*B_prior_var[m]/S+sigma_A_noise
    mean_B=((A_sample-B%*%TFA$T_sample+B[,m]%*%t(TFA$T_sample[m,]))%*%TFA$T_sample[m,]*B_prior_var[m]/S+B_prior_mean[,m]*sigma_A_noise)/temp_var
    vairance_B=B_prior_var[m]*sigma_A_noise/temp_var
    
    bb = rnorm(P)
    bb[which(bb-3 > 0)] = 3
    bb[which(bb+3 < 0)] = -3
    
    B[,m]=(bb*sqrt(abs(vairance_B))+mean_B)*B_state[,m]
  }
  
  return(B)
}




#*************************** Peak-gene looping confidence estimation ***************************
Peak_gene_looping_L_samping <- function(R, R_sample, RNA_Cell_Sample_vector, TFA, B, L, L_state, L_prior_mean, L_prior_var, sigma_R_noise, P, G, M, S){
  
  T_sample=matrix(0, nrow=M, ncol=S)
  for (s in 1:S){
    T_sample[,s]=rowMeans(TFA$T_R[,which(RNA_Cell_Sample_vector==s)])
  }
  
  A_estimate=B%*%TFA$T_sample
  
  Peak_index=sample(P)
  for (i in 1:P){
    f=Peak_index[i]
    temp_var=sum(A_estimate[f,]^2)*L_prior_var/S+sigma_R_noise
    mean_L=(A_estimate[f,]%*%t(R_sample-t(L)%*%A_estimate+L[f,]%*%t(A_estimate[f,]))*L_prior_var/S+L_prior_mean[f,]*sigma_R_noise)/temp_var
    vairance_L=L_prior_var*sigma_R_noise/temp_var
    
    ll=rnorm(G)
    ll[which(ll-3>0)]=3
    ll[which(ll+3<0)]=-3
    
    L[f,]=(ll*sqrt(vairance_L)+mean_L)*L_state[f,]
  }
  return(L)
}




#*************************** TF-peak binding state update ***************************
TF_peak_binary_binding_B_state_sampling <-function(A, A_sample, ATAC_Cell_Sample_vector, TFA,
                                                   B, B_state, B_prior_mean, B_prior_var, B_prior_prob, sigma_A_noise, P, G, M, S){
  T_sample=matrix(0, nrow=M, ncol=S)
  for (s in 1:S){
    T_sample[,s]=rowMeans(TFA$T_A[,which(ATAC_Cell_Sample_vector==s)])
  }
  
  TF_index=sample(M)
  Peak_index=sample(P)
  
  for (i in 1:P){
    
    f=Peak_index[i]
    temp=t(B[f,])%*%TFA$T_sample
    
    for (j in 1:M){
      
      m=TF_index[j]
      
      temp_var=sum(TFA$T_sample[m,]^2)*B_prior_var[m]/S+sigma_A_noise
      
      if (B_state[f,m]>0){
        
        mean_B=(sum((A_sample[f,]-temp+B[f,m]*TFA$T_sample[m,])*TFA$T_sample[m,])*B_prior_var[m]/S+B_prior_mean[f,m]*sigma_A_noise)/temp_var
       
        vairance_B=B_prior_var[m]*sigma_A_noise/temp_var
        
        post_b1 = exp(-(B[f,m]*1-mean_B)^2/(2*vairance_B))*(B_prior_prob[f,m]+0.25)+1e-6
        
        post_b0 = exp(-(B[f,m]*0-mean_B)^2/(2*vairance_B))*(1-B_prior_prob[f,m]+0.25)+1e-6
        
        P1=post_b1/(post_b1+post_b0)
        
        threshold_c = runif(1)
        
        if (is.na(P1) || is.nan(P1)){
          P1=0.5
        }
        
        if (P1 < threshold_c){
          
          B[f,m] = 0
          B_state[f,m] = 0
          
        }else{
          
          B[f,m] = B[f,m]
          B_state[f,m] = 1
          
        }
      } 
      
      if (B_state[f,m]==0 && B_prior_prob[f,m]>0){
        mean_B = (sum((A_sample[f,]-temp)*TFA$T_sample[m,])*B_prior_var[m]/S+B_prior_mean[f,m]*sigma_A_noise)/temp_var
        vairance_B = B_prior_var[m]*sigma_A_noise/temp_var
        
        bb=rnorm(1)
        
        if (bb-3 > 0){
          bb=3
        }
        
        if (bb+3 < 0){
          bb=-3
        }
        
        B_temp = bb*sqrt(abs(vairance_B))+mean_B
        
        post_b1 = exp(-(bb^2)/2)*(B_prior_prob[f,m]+0.25)+1e-6
        
        post_b0 = exp(-(mean_B)^2/(2*vairance_B))*(1-B_prior_prob[f,m]+0.25)+1e-6
        
        P1 = post_b1/(post_b1+post_b0)
        
        threshold_c = runif(1)
        
        if (is.na(P1) || is.nan(P1)){
          P1=0.5
        }
        
        if (P1 < threshold_c){
          
          B[f,m] = 0
          B_state[f,m] = 0
          
        } else {
          
          B[f,m] = B_temp
          B_state[f,m] = 1          
  
        }
      }
    }
  }
  Updated_B_and_State = list('B'=B, 'B_state'=B_state)
  return(Updated_B_and_State)
}




#*************************** Peak-gene looping state update ***************************
Peak_gene_binary_looping_L_state_samping <- function(R, R_sample, RNA_Cell_Sample_vector, TFA, B, L, L_state, L_prior_mean, L_prior_var, L_prior_prob, sigma_R_noise, P, G, M, S){
  
  Peak_index=sample(P)
  Gene_index=sample(G)
  
  T_sample=matrix(0, nrow=M, ncol=S)
  for (s in 1:S){
    T_sample[,s]=rowMeans(TFA$T_R[,which(RNA_Cell_Sample_vector==s)])
  }

  
  A_estimate=B%*%TFA$T_sample
  
  
  for (i in 1:G){
    g=Gene_index[i]
    temp=L[,g]%*%A_estimate
    
    for (j in 1:P){
      f=Peak_index[j]
      temp_var=sum(A_estimate[f,]^2)*L_prior_var/S+sigma_R_noise
      
      if (L_state[f,g]>0){
        
        mean_L=(sum((R_sample[g,]-temp+L[f,g]*A_estimate[f,])*A_estimate[f,])*L_prior_var/S+L_prior_mean[f,g]*sigma_R_noise)/temp_var
        vairance_L=L_prior_var*sigma_R_noise/temp_var
        
        post_l1 = exp(-(L[f,g]*1-mean_L)^2/(2*vairance_L))*(L_prior_prob[f,g]+0.25)+1e-6
        
        post_l0 = exp(-(L[f,g]*0-mean_L)^2/(2*vairance_L))*(1-L_prior_prob[f,g]+0.25)+1e-6
        
        P1=post_l1/(post_l1+post_l0)
        
        threshold_c=runif(1)
        
        if (is.na(P1) || is.nan(P1)){
          P1=0.5
        }
        
        if (P1 < threshold_c){
          
          L[f,g] = 0
          L_state[f,g] = 0
          
        }else{
          
          L[f,g] = L[f,g]
          L_state[f,g] = 1
          
        }
      }
      
      if (L_state[f,g]==0 && L_prior_prob[f,g]>0){
        
        mean_L=(sum((R_sample[g,]-temp)*A_estimate[f,])*L_prior_var/S+L_prior_mean[f,g]*sigma_R_noise)/temp_var
        vairance_L=L_prior_var*sigma_R_noise/temp_var
        
        ll = rnorm(1)
        
        if (ll-3 > 0){
          ll = 3
        }
        if (ll+3 < 0){
          ll = -3
        }
        
        L_temp = ll*sqrt(abs(vairance_L))+mean_L
        
        post_l1 = exp(-ll^2/2)*(L_prior_prob[f,g]+0.1)+1e-6
        
        post_l0 = exp(-(mean_L)^2/(2*vairance_L))*(1-L_prior_prob[f,g]+0.1)+1e-6
        
        P1 = post_l1/(post_l1+post_l0)
        
        threshold_c=runif(1)
        
        if (is.na(P1) || is.nan(P1)){
          P1=0.5
        }
        
        if (P1<threshold_c){
          
          L[f,g] = 0
          L_state[f,g] = 0          
          
        }else{
          
          L[f,g] = L_temp
          L_state[f,g] = 1
          
        }
      }
    }
  }
  Updated_L_and_State = list('L'=L, 'L_state'=L_state)
  return(Updated_L_and_State)
}




#*************************** MCMC sampling***************************

MAGICAL_estimation <- function(loaded_data, Candidate_circuits, Initial_model, iteration_num){
  

  Common_samples=loaded_data$Common_samples
  S=length(Common_samples)#number of samples
  
  # ATAC data
  Peak_index=match(Candidate_circuits$Peaks$Peak_index, loaded_data$scATAC_Peaks$Peak_index)
  A=loaded_data$scATAC_read_count_matrix[Peak_index,]
  ATAC_Cell_Sample_vector=matrix(0, nrow=1, ncol=nrow(loaded_data$scATAC_cells))
  for (s in 1:S){
    ATAC_Cell_Sample_vector[which(loaded_data$scATAC_cells$subject_ID==Common_samples[s])]=s
  }#ATAC cells of each sample
  A_sample=Candidate_circuits$Peak_log2Count
  P=length(Peak_index)#number of peaks
  
  
  
  # RNA data
  Gene_index=match(Candidate_circuits$Genes$Gene_symbols, loaded_data$scRNA_Genes$Gene_symbols)
  R=loaded_data$scRNA_read_count_matrix[Gene_index,]
  RNA_Cell_Sample_vector=matrix(0, nrow=1, ncol=nrow(loaded_data$scRNA_cells))
  for (s in 1:S){
    RNA_Cell_Sample_vector[1,which(loaded_data$scRNA_cells$subject_ID==Common_samples[s])]=s
  }#RNA cells of each sample
  R_sample=Candidate_circuits$Gene_log2Count
  G=length(Gene_index)#number of genes
  
  
  # TF-peak binding initial values
  B_prior_mean=Initial_model$B_mean 
  B_prior_var=Initial_model$B_var 
  B_prior_prob=Initial_model$B_prob
  B=Initial_model$B_prior
  B_state=as.matrix(Candidate_circuits$TF_Peak_Binding) 
  
  
  
  # Peak-gene looping initial values
  L_prior_mean=Initial_model$L_mean
  L_prior_var=Initial_model$L_var 
  L_prior_prob=Initial_model$L_prob
  L=Initial_model$L_prior
  L_state=as.matrix(Candidate_circuits$Peak_Gene_looping) 
  
  
  
  # TFA initial values
  T_prior_mean=Initial_model$T_mean 
  T_prior_var=Initial_model$T_var
  T_sample=T_prior_mean
  M=nrow(T_sample)#number of TFs
  
  T_A=matrix(0, nrow=nrow(Candidate_circuits$TFs), ncol=nrow(loaded_data$scATAC_cells))
  for (s in 1:S){
    for (m in 1:M){
      index=which(ATAC_Cell_Sample_vector==s)
      T_A[m,index]=rnorm(length(index), mean=T_prior_mean[m,s], sd=sqrt(abs(T_prior_var[m,s])))
    }
  }#initial values of the hidden TF activity of ATAC cells for each sample
  
  T_R=matrix(0, nrow=nrow(Candidate_circuits$TFs), ncol=nrow(loaded_data$scRNA_cells))
  for (s in 1:S){
    for (m in 1:M){
      index=which(RNA_Cell_Sample_vector==s)
      T_R[m,index]=rnorm(length(index), mean=T_prior_mean[m,s], sd=sqrt(abs(T_prior_var[m,s])))
    }
  }#initial values of the hidden TF activity of RNA cells for each sample
  
  TFA=list('T_A'=T_A, 'T_R'=T_R, 'T_sample'=T_sample)
  
  
  
  # ATAC fitting residual variance initial values
  alpha_A=1
  beta_A=1
  sigma_A_noise=1/rgamma(1, shape=alpha_A+1/2, rate=(beta_A+sum((A_sample-B%*%TFA$T_sample)^2))/(2*P*S))
  
  
  
  # RNA fitting residual variance initial values
  alpha_R=1
  beta_R=1
  sigma_R_noise=1/rgamma(1, shape=alpha_R+1/2, scale=1/((beta_R+sum((R_sample-t(L)%*%(B%*%TFA$T_sample))^2))/(2*G*S)))
  
  
  
  #MCMC sampling initial state
  iteration_seg=round(iteration_num/10)
  B_state_frq=B_state
  L_state_frq=L_state
  
  print('MAGICAL integration starts ...')
  cat('\n')
  
  
  Noise_parameters=matrix(0, nrow=iteration_num, ncol=2)
  for (i in 1:iteration_num){
    
    #Step 1: TF activity sampling
    TFA <- TF_activity_T_sampling(A, A_sample, ATAC_Cell_Sample_vector, R, R_sample, RNA_Cell_Sample_vector,
                                  TFA, T_prior_mean, T_prior_var, B, B_state, sigma_A_noise, P, G, M, S)
    
    
    #Step 2: TF-peak binding weight sampling
    B <- TF_peak_binding_B_sampling(A, A_sample, ATAC_Cell_Sample_vector, TFA,
                                    B, B_state, B_prior_mean, B_prior_var, sigma_A_noise, P, G, M, S)
    
    
    #Step 3: TF-peak binding state update
    Updated_B_and_State <- TF_peak_binary_binding_B_state_sampling(A, A_sample, ATAC_Cell_Sample_vector,
                                                                   TFA, B, B_state, B_prior_mean, B_prior_var, B_prior_prob,
                                                                   sigma_A_noise, P, G, M, S)
    B <- Updated_B_and_State$B
    B_state <- Updated_B_and_State$B_state
    B_state_frq = B_state_frq+B_state
    
    
    #Step 4: Peak-Gene looping weight sampling
    L <- Peak_gene_looping_L_samping(R, R_sample, RNA_Cell_Sample_vector, TFA, B, L, L_state, L_prior_mean, L_prior_var, sigma_R_noise, P, G, M, S)
    
    
    #Step 5: Peak-Gene looping state update
    Updated_L_and_State <- Peak_gene_binary_looping_L_state_samping(R, R_sample, RNA_Cell_Sample_vector, TFA, B, L, L_state, L_prior_mean, L_prior_var, L_prior_prob, sigma_R_noise, P, G, M, S)
    L <- Updated_L_and_State$L
    L_state <- Updated_L_and_State$L_state    
    L_state_frq=L_state_frq+L_state
    
    
    #Step 6: fitting residue variance control
    sigma_A_noise = 1/rgamma(1, shape=alpha_A+1/2, rate=(beta_A+sum((A_sample-B%*%TFA$T_sample)^2))/(2*P*S))
    
    sigma_R_noise = 1/rgamma(1, shape=alpha_R+1/2, rate=(beta_R+sum((R_sample-t(L)%*%(B%*%TFA$T_sample))^2))/(2*G*S))
    
    Noise_parameters[i,1]=sigma_A_noise
    Noise_parameters[i,2]=sigma_R_noise
    
    if (i%%iteration_seg==0){
      print(paste('MAGICAL finished', 10*i/iteration_seg,  'percent', sep=' '))
      cat('\n')
    }
    
    
  }
  
  #Sample Summary
  Candidate_TF_Peak_Binding_prob = B_state_frq/(iteration_num+1)
  Candidate_Peak_Gene_Looping_prob = L_state_frq/(iteration_num+1)
  
  Circuits_linkage_posterior<-list('TF_Peak_Binding_prob'=Candidate_TF_Peak_Binding_prob,
                                   'Peak_Gene_Looping_prob'=Candidate_Peak_Gene_Looping_prob,
                                   'Noise_parameters'=Noise_parameters)
  
  return(Circuits_linkage_posterior)
}





#*************************** MAGICAL output ***************************

MAGICAL_circuits_output<- function(Output_file_path, Candidate_circuits, Circuits_linkage_posterior,
                                   prob_threshold_TF_peak_binding=0.8, prob_threshold_peak_gene_looping=0.95){
  
  Peak_Gene_index <- which(Circuits_linkage_posterior$Peak_Gene_Looping_prob>prob_threshold_peak_gene_looping, arr.ind = TRUE)
  circuit_flag=matrix(0, nrow=nrow(Peak_Gene_index), ncol=1)
  writeLines('Gene_symbol\tGene_chr\tGene_TSS\tPeak_chr\tPeak_start\tPeak_end\tLooping_prob\tTFs(binding prob)\n', Output_file_path)
  
  TF_vector=matrix(0, nrow(Candidate_circuits$TFs), ncol=1)
  for (i in 1:nrow(Peak_Gene_index)){
    
    TF_prob <- sort.int(Circuits_linkage_posterior$TF_Peak_Binding_prob[Peak_Gene_index[i,1],], decreasing = TRUE, index.return = TRUE)
    index=which(TF_prob$x>prob_threshold_TF_peak_binding)
    
    if (length(index)>0){
      
      circuit_flag[i]=1
      
      cat(Candidate_circuits$Genes$Gene_symbols[Peak_Gene_index[i,2]],
          Candidate_circuits$Genes$chr[Peak_Gene_index[i,2]],
          Candidate_circuits$Genes$TSS[Peak_Gene_index[i,2]],
          Candidate_circuits$Peaks$chr[Peak_Gene_index[i,1]],
          Candidate_circuits$Peaks$point1[Peak_Gene_index[i,1]],
          Candidate_circuits$Peaks$point2[Peak_Gene_index[i,1]],
          Circuits_linkage_posterior$Peak_Gene_Looping_prob[Peak_Gene_index[i,1],Peak_Gene_index[i,2]],
          file=Output_file_path,sep="\t",append=TRUE)
      
      cat('\t', file=Output_file_path,append=TRUE)
      
      for (j in 1:length(index)){
        cat(paste(Candidate_circuits$TFs[TF_prob$ix[j],2], ',(', TF_prob$x[j],'),', 'sep'= ''), file=Output_file_path, sep=" ", append=TRUE)
        
        TF_vector[TF_prob$ix[j]]=TF_vector[TF_prob$ix[j]]+1
      }
    }
    cat('\n', file=Output_file_path, sep="", append=TRUE)
  }   
  print(paste('MAGICAL selected regulatory circuits with', length(which(TF_vector>1)), 'TFs,', length(unique(Peak_Gene_index[circuit_flag>0,1])), 'peaks, and',  length(unique(Peak_Gene_index[circuit_flag>0,2])), 'genes', sep=' '))
}

