#*****************load Seurat integrated scRNA-seq data*******************
library(Seurat)
library(Matrix)
library(stringr)
library(dplyr)
library(tidyr)
library(data.table)

scRNA <- readRDS(file = 'path to scRNA.rds')

Cell_type_selected = 'xxxx'

#RNA assay genes
write.table(scRNA@assays$RNA@counts@Dimnames[[1]], file = 'scRNA genes.txt',
            quote = FALSE, row.names = TRUE, col.names = FALSE, sep = '\t')


#RNA assay cell read counts
cell_index=which(scRNA$refined_annotations==Cell_type_selected)
Cell_type_scRNA_counts = as(scRNA@assays$RNA@counts[, cell_index], "dgTMatrix")
write.table(summary(Cell_type_scRNA_counts), file='Cell type scRNA read count.txt',
            quote = FALSE, row.names = FALSE, col.names = FALSE,  sep = '\t')


#RNA assay cell meta including barcode, cell type annotation, sample names/IDs and group/condition
Cell_type_scRNA_meta = scRNA@meta.data[cell_index,]
write.table(data.frame(rownames(Cell_type_scRNA_meta), Cell_type_scRNA_meta$refined_annotations, Cell_type_scRNA_meta$Sample, Cell_type_scRNA_meta$Group),
            file = 'Cell type scRNA cell meta.txt',
            quote = FALSE, row.names = TRUE, col.names = FALSE,  sep = '\t')





#*****************load Signac integrated scATAC-seq data*******************

library(Signac)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(chromVARmotifs)

scATAC <- readRDS(file = 'path to scATAC.rds')

Cell_type_selected = 'xxxx'

DefaultAssay(scATAC) <- "ATAC"
main.chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
keep.peaks <- as.logical(seqnames(granges(scATAC)) %in% main.chroms)
scATAC <- scATAC[keep.peaks, ]


#ATAC assay peaks
write.table(scATAC@assays$ATAC@ranges, file = 'scATAC peaks.txt', quote = FALSE, col.names = FALSE,  sep = '\t')


#ATAC assay cell read counts
cell_index=which(scATAC$refined_annotations==Cell_type_selected)
Cell_type_scATAC_counts = as(scATAC@assays$ATAC@counts[, cell_index], "dgTMatrix")
write.table(summary(Cell_type_scATAC_counts), file='Cell type scATAC read count.txt',
            quote = FALSE, row.names = FALSE, col.names = FALSE,  sep = '\t')


#ATAC assay cell meta including barcode, cell type annotation, and sample names/IDs
Cell_type_scATAC_meta = scATAC@meta.data[cell_index,]
write.table(data.frame(rownames(Cell_type_scATAC_meta), Cell_type_scATAC_meta$refined_annotations, Cell_type_scATAC_meta$Sample),
            file = 'Cell type scATAC cell meta.txt',
            quote = FALSE, row.names = TRUE, col.names = FALSE,  sep = '\t')


#TF motif mapping
scATAC <- AddMotifs(object = scATAC, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = human_pwms_v2)
Peak_motif_mapping  = as(scATAC@assays$ATAC@motifs@data, "dgTMatrix")
write.table(summary(Peak_motif_mapping), file='Motif mapping prior.txt',
            quote = FALSE, row.names = FALSE, col.names = FALSE,  sep = '\t')
write.table(Peak_motif_mapping@Dimnames[[2]], file = 'Motifs.txt', quote = FALSE, col.names = FALSE,  sep = '\t')
#make sure to check the motif name file to clean up the names and use the standard gene symbol for each motif






#*****************load ArchR integrated scATAC-seq data*******************

library(ArchR)
library(chromVARmotifs)

proj<-loadArchRProject('path to scATAC proj', force = FALSE, showLogo = TRUE)

Cell_type_selected = 'xxxx'

idxPass <- which(str_detect(proj$Cell_type_voting,Cell_type_selected))
proj_temp<-proj[idxPass, ]

#ATAC assay peaks
Peaks <- getPeakSet(proj_temp)
write.table(Peaks[,1], file = 'scATAC peaks.txt', quote = FALSE, col.names = FALSE,  sep = '\t')


#ATAC assay cell count
Cell_type_scATAC_counts<-getMatrixFromProject(ArchRProj = proj_temp, useMatrix = "PeakMatrix", useSeqnames = NULL,
                                verbose = TRUE,binarize = FALSE,threads = getArchRThreads(),
                                logFile = createLogFile("getMatrixFromProject"))
write.table(summary(Cell_type_scATAC_counts@assays@data@listData$PeakMatrix), 
            file='Cell type scATAC read count.txt',
            quote = FALSE, row.names = FALSE, col.names = FALSE,  sep = '\t')


#ATAC assay cell meta
write.table(data.frame(proj_temp$cellNames, proj_temp$Cell_type_voting, proj_temp$Sample),
            file = 'Cell type scATAC cell meta.txt',
            quote = FALSE, row.names = TRUE, col.names = FALSE,  sep = '\t')


#TF motif mapping
proj_temp <- addMotifAnnotations(ArchRProj = proj_temp, motifSet = "cisbp", name = "Motif")
#Different from Signac, motif match in ArchR is not stored in the data space, 
#instead, a Motif-Matches-In-Peaks.rds file will be created under the Annotations folder
Peak_motif_mapping <- readRDS(file = "~/Desktop/ECHO/Staph/Staph_scATAC_integration/Annotations/Motif-Matches-In-Peaks.rds")
write.table(lapply(summary(Peak_motif_mapping@assays@data@listData$matches), as.numeric), file='Motif mapping prior.txt',
            quote = FALSE, row.names = FALSE, col.names = FALSE,  sep = '\t')
write.table(Peak_motif_mapping@assays@data@listData$matches@Dimnames[[2]], file = 'Motifs.txt', quote = FALSE, col.names = FALSE,  sep = '\t')
#make sure to check the motif name file to clean up the names and use the standard gene symbol for each motif
