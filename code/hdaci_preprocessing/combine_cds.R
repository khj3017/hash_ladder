suppressPackageStartupMessages({
    library(monocle3)
    library(MASS)
    library(reshape2)
    library(VGAM)
    library(dplyr)
    library(ggplot2)
    library(grid)
    library(gridExtra)
    library(Matrix)
    library(IRdisplay)
    library(tibble)
    library(ggridges)
    library(RColorBrewer)
})

wd = getwd()
wd

### load.cds from sci-RNA-seq output
load.cds = function(mat.path, gene.annotation.path, cell.annotation.path, 
                    hashTable.path, Experiment = 1, Batch = 1,
                    min_expr_cutoff = 0, num_cells = 0) {
  df = read.table(
    mat.path,
    col.names = c("gene.idx", "cell.idx", "count"),
    colClasses = c("integer", "integer", "integer"))
  
  gene.annotations = read.table(
    gene.annotation.path,
    col.names = c("id", "gene_short_name"),
    colClasses = c("character", "character"))
  
  cell.annotations = read.table(
    cell.annotation.path,
    col.names = c("Cell", "Sample"),
    colClasses = c("character", "factor"))
  
  hashTable <- read.table(hashTable.path, header=T, stringsAsFactors = F)
  
  
  rownames(gene.annotations) = gene.annotations$id
  rownames(cell.annotations) = cell.annotations$Cell
    
  
  # add a dummy cell to ensure that all genes are included in the matrix
  # even if a gene isn't expressed in any cell
  df = rbind(df, data.frame(
    gene.idx = c(1, nrow(gene.annotations)),
    cell.idx = rep(nrow(cell.annotations)+1, 2),
    count = c(1, 1)))
  

  mat = sparseMatrix(i = df$gene.idx, j = df$cell.idx, x = df$count)
  mat = mat[, 1:(ncol(mat)-1)]
    
  rownames(mat) = gene.annotations$id
  colnames(mat) = cell.annotations$Cell
  
  matchNames <- match(hashTable$Cell, rownames(cell.annotations))
  mat <- mat[,matchNames]

  
  cell.annotations <- cell.annotations[matchNames,]
  cell.annotations <- cbind(cell.annotations, hashTable[,6:ncol(hashTable)])

  plate <- if (Batch == 1) {
              ifelse(grepl("^E", cell.annotations$PCR_well), 1, 2)
           } else {
              ifelse(grepl("^F", cell.annotations$PCR_well), 3, 4)
           }
    
  cell.annotations <- cbind(cell.annotations, Plate = as.factor(plate))
  
  cell.annotations$RT_well <- as.factor(cell.annotations$RT_well)
  cell.annotations$Dose <- as.factor(cell.annotations$Dose)
    
  cds = new_cell_data_set(mat, 
                          cell_metadata = cell.annotations, 
                          gene_metadata = gene.annotations)
    
  cds <- cds[row.names(subset(fData(cds), !grepl("MU",id))), cds$Dim == Experiment]
  cds <- detect_genes(cds, min_expr = min_expr_cutoff)
  expressed_genes <- row.names(subset(cds, num_cells_expressed >= num_cells))
  print(paste0("# of expressed genes: ", length(expressed_genes)))
    
  cds <- cds[expressed_genes, ]
    
  return(cds)
}

combine_cds = function(cds_batch_1, cds_batch_2, min_expr_cutoff = 1, num_cells = 0) {
    combined_mat = cbind(counts(cds_batch_1), counts(cds_batch_2))
    combined_colData = rbind(as.data.frame(colData(cds_batch_1)), as.data.frame(colData(cds_batch_2)))

    cds = new_cell_data_set(combined_mat,
                            cell_metadata = combined_colData,
                            gene_metadata = rowData(cds_batch_1))

    cds <- detect_genes(cds, min_expr = min_expr_cutoff)
    expressed_genes <- row.names(subset(cds, num_cells_expressed >= num_cells))
    print(paste0("# of expressed genes: ", length(expressed_genes)))

    cds <- cds[expressed_genes, ]


    return(cds)
}

setwd("../../data/hdaci_qc/")

# HDACi timecourse experiment
## load cds
cds_batch_1 = load.cds('UMI.count.matrix_rep1', '../refdata/gene.annotations.hg38.mm10', 
               'cell.annotations_rep1', 'hashTable_unique_filtered_rep1.txt', 
                Experiment = 1, Batch = 1, min_expr_cutoff = 0, num_cells = 0)

cds_batch_1$Hash_Size_Factor = NULL

dim(as.data.frame(colData(cds_batch_1)))
head(as.data.frame(colData(cds_batch_1)))

## load cds
cds_batch_2 = load.cds('UMI.count.matrix_rep2', '../refdata/gene.annotations.hg38.mm10', 
               'cell.annotations_rep2', 'hashTable_unique_filtered_rep2.txt', 
                Experiment = 1, Batch = 2, min_expr_cutoff = 0, num_cells = 0)

cds_batch_2$Hash_Size_Factor = NULL

dim(as.data.frame(colData(cds_batch_2)))
head(as.data.frame(colData(cds_batch_2)))

saveRDS(cds_batch_1, "cds_dim_1_batch_1.rds")
saveRDS(cds_batch_2, "cds_dim_1_batch_2.rds")

cds = combine_cds(cds_batch_1, cds_batch_2)
head(as.data.frame(colData(cds)))
metadata = as.data.frame(colData(cds))

saveRDS(cds, "cds_dim_1_combined.rds")



# HDACi/Dex experiment
## load cds
cds_batch_1 = load.cds('UMI.count.matrix_rep1', '../refdata/gene.annotations.hg38.mm10', 
               'cell.annotations_rep1', 'hashTable_unique_filtered_rep1.txt', 
                Experiment = 2, Batch = 1, min_expr_cutoff = 0, num_cells = 0)

cds_batch_1$Hash_Size_Factor = NULL

## load cds
cds_batch_2 = load.cds('UMI.count.matrix_rep2', '../refdata/gene.annotations.hg38.mm10', 
               'cell.annotations_rep2', 'hashTable_unique_filtered_rep2.txt', 
                Experiment = 2, Batch = 2, min_expr_cutoff = 0, num_cells = 0)

cds_batch_2$Hash_Size_Factor = NULL

head(as.data.frame(colData(cds_batch_2)))

head(as.data.frame(colData(cds_batch_2)))

saveRDS(cds_batch_1, "cds_dim_2_batch_1.rds")
saveRDS(cds_batch_2, "cds_dim_2_batch_2.rds")


cds = combine_cds(cds_batch_1, cds_batch_2)

cds
head(as.data.frame(colData(cds)))
metadata = as.data.frame(colData(cds))

saveRDS(cds, "cds_dim_2_combined.rds")


