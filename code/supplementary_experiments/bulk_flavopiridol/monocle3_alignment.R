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
    library(dlookr)
})

setwd("../../../data/supplementary_experiments/bulk_flavopiridol")

### load.cds from sci-RNA-seq output
### normalization = c("standard", "hash_size_factor", "hash_ladder", "multinomial")
load.cds = function(mat.path, gene.annotation.path, cell.annotation.path) {
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
    col.names = c("Sample", "V2"),
    colClasses = c("character", "factor"))
  
    
  
  rownames(gene.annotations) = gene.annotations$id
  rownames(cell.annotations) = cell.annotations$Sample
    
  
  # add a dummy cell to ensure that all genes are included in the matrix
  # even if a gene isn't expressed in any cell
  df = rbind(df, data.frame(
    gene.idx = c(1, nrow(gene.annotations)),
    cell.idx = rep(nrow(cell.annotations)+1, 2),
    count = c(1, 1)))
  

  mat = sparseMatrix(i = df$gene.idx, j = df$cell.idx, x = df$count)
  mat = mat[, 1:(ncol(mat)-1)]
  
  rownames(mat) = gene.annotations$id
  colnames(mat) = cell.annotations$Sample
  
  cell.annotations = cell.annotations %>% rowwise() %>%
                        mutate(RT = as.double(strsplit(Sample, "_")[[1]][4])) %>% ungroup()
  rownames(cell.annotations) = cell.annotations$Sample
    
  cds = new_cell_data_set(mat, 
                          cell_metadata = cell.annotations %>% select(-2), 
                          gene_metadata = gene.annotations)
  
  return(cds)
}

## load cds
cds = load.cds('UMI.count.matrix', 'gene.annotations', 
               'cell.annotations')
cds = detect_genes(cds, min_expr = 0)
cds = cds[rowData(cds)$num_cells_expressed > 0,]
cds

mat = as.matrix(counts(cds))
colnames(mat) = cds$RT

aggregated_mat = vapply(unique(colnames(mat)), function(x) 
      rowSums(mat[,colnames(mat)== x,drop=FALSE], na.rm=TRUE),
                             numeric(nrow(mat)) )

spike_in_table = read.table("ERCC_table_corrected_nb_new.txt", header=T) %>%
                    distinct(RT, .keep_all=T) %>%
                    select(-c(Sample, Count)) %>%
                    mutate(RT = as.factor(RT))
rownames(spike_in_table) = spike_in_table$RT
spike_in_table = spike_in_table[match(rownames(spike_in_table), colnames(aggregated_mat)), ]
head(spike_in_table)


## generate new cds using the spike_in_table
cds = new_cell_data_set(Matrix(aggregated_mat, sparse = TRUE),
                          cell_metadata = spike_in_table, 
                          gene_metadata = gene.annotations)
cds

## Cell count correction; more cells from 1hr timepoint were used for preparing the library
correction_factor = 5/3
mat = counts(cds)
mat[,cds$FP == 1] = round(counts(cds[,cds$FP == 1])/correction_factor)
counts(cds) = mat

RT_correction = data.frame(FP = c(0,1,3,6,12,24), 
                           RT_vol = c(1,correction_factor,1,1,1,1))
RT_df = as.data.frame(colData(cds)) %>%
    left_join(RT_correction, by = "FP") %>%
    mutate(Total_RNA = Total_RNA/RT_vol,
           Total_ERCC = Total_ERCC/RT_vol)

cds$Total_RNA = RT_df$Total_RNA
cds$Total_ERCC = RT_df$Total_ERCC

## Normalization using the ERCC counts
cds2 = cds
scale_factors_2 = log(cds2$Total_ERCC / cds2$ERCC_duplication) * cds2$RNA_duplication * cds2$Slope / -cds2$Intercept
cds2$Size_Factor = scale_factors_2 / exp(mean(log(scale_factors_2))) 

## Supplementary Figure 5a
options(repr.plot.width=3, repr.plot.height=3)
metadata %>%
    ggplot(., aes(y = log10(Total_RNA / Size_Factor), x = as.factor(FP))) +
      geom_point(color = "black", size = 2.5) + 
      geom_point(aes(color = as.factor(FP)), size = 2) + 
      scale_color_brewer(palette = "YlGnBu") + 
      labs(x = "Time", y = "log10 Total RNA count", fill = "Time (hrs)") +
      monocle3:::monocle_theme_opts() + theme(legend.position = "none")

ggsave("R_FP_total_RNA_time.pdf", device = "pdf", width=3, height=3)


saveRDS(cds, "cds_bulk_conventional.rds")
saveRDS(cds2, "cds_bulk_ERCC.rds")
