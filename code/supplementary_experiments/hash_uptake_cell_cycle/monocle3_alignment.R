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
    library(VennDiagram)
    library(pheatmap)
    library(piano)
    library(snowfall)
    library(stringr)
    library(cowplot)
    library(eulerr)
    library(tidyr)
})

### load.cds from sci-RNA-seq output
### normalization = c("standard", "hash_size_factor", "hash_ladder", "multinomial")
load.cds = function(mat.path, gene.annotation.path, cell.annotation.path, 
                    hashTable.path) {
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
  
  hashTable <- read.table(hashTable.path, header=T)
  
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
  cell.annotations <- cbind(cell.annotations, hashTable[,c(5:ncol(hashTable))])
  
  batch <- ifelse(grepl("_H0", cell.annotations$PCR_well), 1, 2)
  cell.annotations <- cbind(cell.annotations, Batch = as.factor(batch))
  
  cell.annotations$RT_well <- as.factor(cell.annotations$RT_well)
    
  cds = new_cell_data_set(mat, 
                          cell_metadata = cell.annotations, 
                          gene_metadata = gene.annotations)
  
  return(cds)
}

setwd("../../../data/supplementary_experiments/hash_uptake_cell_cycle")

cds = load.cds('UMI.count.matrix', 'gene.annotations', 'cell.annotations', 'hashTable_unique_filtered.txt')
cds2 <- cds
scale_factors_2 = log(cds2$Total_hash / cds2$Hash_duplication) * cds2$RNA_duplication * cds2$Slope / -cds2$Intercept
cds2$Size_Factor = scale_factors_2 / exp(mean(log(scale_factors_2))) 

saveRDS(cds, "cds_sf.rds")
saveRDS(cds2, "cds_hash.rds")

## Monocle3 alignment
cds = detect_genes(cds, min_expr = 1)
cds2 = detect_genes(cds2, min_expr = 1)
genes = row.names(subset(rowData(cds2), num_cells_expressed > 30))

## Conventional normalization
cds = preprocess_cds(cds, num_dim = 30, use_genes = genes)
cds = align_cds(cds, alignment_k = 20, 
                    residual_model_formula_str = "~log(Total_RNA)")
cds = reduce_dimension(cds, reduction_method = "UMAP", max_components = 2,
                       umap.n_neighbors = 40, umap.min_dist = 0.01)

cds = cluster_cells(cds, resolution=0.0001, partition_qval = 0.5)

cds$Cluster = clusters(cds)
cds$Cluster = factor(cds$Cluster, levels = c('1', '2'))
cds$Cell_type_cluster = ifelse(cds$Cluster == '1', "HEK293T", "A549")

## Hash ladder normalization
cds2 = preprocess_cds(cds2, num_dim = 30, use_genes = genes)
cds2 = align_cds(cds2, alignment_k = 20, 
                    residual_model_formula_str = "~log(Total_RNA)")
cds2 = reduce_dimension(cds2, reduction_method = "UMAP", max_components = 2,
                       umap.n_neighbors = 40, umap.min_dist = 0.01)
cds2 = cluster_cells(cds2, resolution=0.0001, partition_qval = 0.5)

cds2$Cluster = clusters(cds2)
cds2$Cluster = factor(cds2$Cluster, levels = c('1', '2'))
cds2$Cell_type_cluster = ifelse(cds2$Cluster == '1', "HEK293T", "A549")

saveRDS(cds, "cds_sf_cluster_adjusted.rds")
saveRDS(cds2, "cds_hash_cluster_adjusted.rds")
