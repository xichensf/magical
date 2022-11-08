setwd("~/Desktop/Muscle/")

library(Seurat)
library(Matrix)
library(stringr)
library(Signac)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(chromVARmotifs)
library(dplyr)
library(tidyr)
library(data.table)

#*****************load single cell multiomics data*************************
scRNA <- readRDS(file = "scRNA.rds")
scATAC <- readRDS(file = "scATAC.rds")

#***********************select cell type **********************************
Cell_type_selected = 'xxxx'


#******extract and output scRNA-seq info of the selected cell type**********
cell_index=which(scRNA$refined_annotations==Cell_type_selected)
#gene_index=which(match(scRNA@assays$RNA@counts@Dimnames[[1]], Candidate_genes)>0)

Cell_type_scRNA_counts = as(scRNA@assays$RNA@counts[, cell_index], "dgTMatrix")
Cell_type_scRNA_meta = scRNA@meta.data[cell_index,]

write.table(summary(Cell_type_scRNA_counts), file=paste(Cell_type_selected, 'scRNA read count.txt',sep=' '),
            quote = FALSE, row.names = FALSE, col.names = FALSE,  sep = '\t')

write.table(Cell_type_scRNA_counts@Dimnames[[1]], file = 'scRNA gene name.txt',
            quote = FALSE, row.names = TRUE, col.names = FALSE, sep = '\t')

write.table(data.frame(rownames(Cell_type_scRNA_meta), Cell_type_scRNA_meta$refined_annotations, Cell_type_scRNA_meta$Sample, Cell_type_scRNA_meta$Group),
            file = paste(Cell_type_selected, 'scRNA cell meta.txt', sep=' '),
            quote = FALSE, row.names = TRUE, col.names = FALSE,  sep = '\t')


#******extract and output scATAC-seq info of the selected cell type**********
DefaultAssay(scATAC) <- "ATAC"
main.chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
keep.peaks <- as.logical(seqnames(granges(scATAC)) %in% main.chroms)
scATAC <- scATAC[keep.peaks, ]

cell_index=which(scATAC$refined_annotations==Cell_type_selected)
Cell_type_scATAC_counts = as(scATAC@assays$ATAC@counts[, cell_index], "dgTMatrix")
Cell_type_scATAC_meta = scATAC@meta.data[cell_index,]

write.table(summary(Cell_type_scATAC_counts), file=paste(Cell_type_selected, 'scATAC read count.txt',sep=''),
            quote = FALSE, row.names = FALSE, col.names = FALSE,  sep = '\t')

write.table(data.frame(rownames(Cell_type_scATAC_meta), Cell_type_scATAC_meta$refined_annotations, Cell_type_scATAC_meta$Sample),
            file = paste(Cell_type_selected, 'scATAC cell meta.txt', sep='_'),
            quote = FALSE, row.names = TRUE, col.names = FALSE,  sep = '\t')


#******************************TF motif mapping**********************************
scATAC <- AddMotifs(object = scATAC, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = human_pwms_v2)

Peak_motif_mapping  = as(scATAC@assays$ATAC@motifs@data, "dgTMatrix")

write.table(summary(Peak_motif_mapping), file='Motif mapping prior.txt',
            quote = FALSE, row.names = FALSE, col.names = FALSE,  sep = '\t')
write.table(Peak_motif_mapping@Dimnames[[2]], file = 'Motif names.txt', quote = FALSE, col.names = FALSE,  sep = '\t')

write.table(Peak_motif_mapping@Dimnames[[1]], file = 'scATAC peaks.txt', quote = FALSE, col.names = FALSE,  sep = '\t')

