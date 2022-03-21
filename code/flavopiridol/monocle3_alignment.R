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
    library(tidyr)
})

monocle_theme_opts <- function()
{
  theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
    theme(panel.border = element_blank()) +
    theme(axis.line.x = element_line(size=0.25, color="black")) +
    theme(axis.line.y = element_line(size=0.25, color="black")) +
    theme(panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank()) +
    theme(panel.background = element_rect(fill='white')) +
    theme(legend.key=element_blank())
}

### load.cds from sci-RNA-seq output
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
  cell.annotations <- cbind(cell.annotations, hashTable[,c(7:ncol(hashTable))])
  
  plate <- ifelse(grepl("_B0", cell.annotations$PCR), 1, 2)
  cell.annotations <- cbind(cell.annotations, Plate = as.factor(plate))
  
  cell.annotations$RT_well <- as.factor(cell.annotations$RT_well)
  #rownames(cell.annotations) = cell.annotations$Cell  

  cds = new_cell_data_set(mat, 
                          cell_metadata = cell.annotations, 
                          gene_metadata = gene.annotations)
  
  cds = cds[row.names(subset(rowData(cds), !grepl("MU",id))),]
  return(cds)
}

setwd("../../data/flavopiridol")

## load cds and remove doublets
cds = load.cds('UMI.count.matrix', '../refdata/gene.annotations.hg38.mm10', 
               'cell.annotations', 'hashTable_unique_filtered.txt')

doublets <- read.table("doublets.txt", colClasses = c("character"))
cds <- cds[,-na.omit(match(cds$Cell, doublets[,1]))]
cds


## filter out mito.genes and subtract mito gene counts
mito_genes <- grep(pattern = "^MT-", x = rowData(cds)$gene_short_name, value = TRUE, ignore.case = TRUE)

mito_gene_counts = colSums(counts(cds)[rowData(cds)$gene_short_name %in% mito_genes,])
mito_gene_id = rowData(cds)[rowData(cds)$gene_short_name %in% mito_genes,]$id

cds$Mito_fraction = mito_gene_counts / cds$Total_RNA

cds$Total_RNA = colSums(counts(cds)) - mito_gene_counts
cds = cds[!rowData(cds)$id %in% mito_gene_id, ]
cds$Total_RNA = Matrix::colSums(counts(cds))
cds = estimate_size_factors(cds)

cds = cds[, cds$Total_RNA >= 400]

cds

cds2 <- cds
scale_factors_2 = log(cds2$Total_hash / cds2$Hash_duplication) * cds2$RNA_duplication * cds2$Slope / -cds2$Intercept
cds2$Size_Factor = scale_factors_2 / exp(mean(log(scale_factors_2))) 
summary(cds2$Size_Factor)

saveRDS(cds, "cds_conv.rds")
saveRDS(cds2, "cds_hash_ladder.rds")


metadata = as.data.frame(colData(cds))
original = c("DMSO", "FP_01hr", "FP_03hr", "FP_06hr", "FP_12hr", "FP_24hr")
convertTo = c('0', '1', '3', '6', '12', '24')
metadata = metadata %>% mutate(Time = convertTo[match(metadata$Condition, original)])
metadata = metadata %>% 
         mutate(Time = factor(Time, levels = convertTo))
head(metadata)

#### Proportion of total hash with respect to total RNA
options(repr.plot.width=5, repr.plot.height=3.5)
ggplot(metadata, aes(y = log10(Total_RNA), x = Time)) + 
  geom_boxplot(aes(fill = Time)) + 
  #ylim(c(6.4, 8.4)) + 
  labs(x = "", y = "") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                     panel.border = element_blank(), axis.line = element_line(colour = "black")) + 
  scale_fill_brewer(palette = "YlGnBu")


## align cells with monocle3
# hash ladder normalization
cds = detect_genes(cds, min_expr = 1)
cds2 = detect_genes(cds2, min_expr = 1)
genes = row.names(subset(rowData(cds2), num_cells_expressed > 30))

cds2 = preprocess_cds(cds2, num_dim = 30, use_genes = genes)
cds2 = align_cds(cds2, alignment_k = 20, 
                    residual_model_formula_str = "~log(Total_RNA)")

options(repr.plot.width=4, repr.plot.height=4.5)
cds2 = reduce_dimension(cds2, reduction_method = "UMAP", max_components = 2,
                       umap.n_neighbors = 15, umap.min_dist = 0.01)
plot_cells(cds2, color_cells_by = "Condition", cell_size = 1.2, group_label_size=4)

cds2 = cluster_cells(cds2, resolution=0.001)
plot_cells(cds2, color_cells_by = "cluster", cell_size = 1, group_label_size = 4)

cds2 <- learn_graph(cds2, use_partition = F, close_loop = F,
                    learn_graph_control = list(ncenter=190, maxiter = 10, minimal_branch_len = 10))
options(repr.plot.width=4, repr.plot.height=4)
plot_cells(cds2, color_cells_by = "Condition", cell_size=1, group_label_size = 4)

filter = cds2@reducedDims$UMAP[,1] < -3.7
root_cells <- rownames(as.data.frame(cds2@reducedDims$UMAP[filter, ]))
print(length(root_cells))

cds2 = order_cells(cds2, root_cells = root_cells)
cds2 = order_cells(cds2, 
                root_pr_nodes = unique(cds2@principal_graph_aux[["UMAP"]]$root_pr_nodes))

options(repr.plot.width=5, repr.plot.height=4)
plot_cells(cds2, cell_size = 1.2,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

cds_plot = cds2
cds_plot$pseudotime = pseudotime(cds_plot)

options(repr.plot.width=4, repr.plot.height=3)
ggplot(as.data.frame(colData(cds_plot)), 
       aes(x = pseudotime, y = Condition)) + 
    geom_density_ridges2(aes(height = ..density.. , fill = Condition)) + 
    labs(x = "Pseudotime", y = "Time (hrs)") +
    theme_ridges() + theme(legend.position = "none", axis.text.x = element_blank())

saveRDS(cds2, "cds_hash_ladder_pseudotime.rds")


# conventional normalization
cds = estimate_size_factors(cds)
cds = preprocess_cds(cds, num_dim = 20, use_genes = genes)
cds = align_cds(cds, alignment_k = 10,
                    residual_model_formula_str = "~log(Total_RNA)")

options(repr.plot.width=4, repr.plot.height=4.5)
cds = reduce_dimension(cds, reduction_method = "UMAP", 
                       umap.n_neighbors = 20, umap.min_dist = 0.01)
plot_cells(cds, color_cells_by = "Condition", cell_size = 1, group_label_size=4)

cds = cluster_cells(cds, resolution=0.001)
plot_cells(cds, color_cells_by = "cluster", cell_size = 1, group_label_size = 4)

options(repr.plot.width=4, repr.plot.height=4.5)
cds <- learn_graph(cds, use_partition = F,
                    learn_graph_control = list(ncenter=130, maxiter = 10, minimal_branch_len = 7))
plot_cells(cds, color_cells_by = "Condition", cell_size=1)

filter = cds@reducedDims$UMAP[,1] < -1.6 & cds@reducedDims$UMAP[,2] < -1.2
root_cells <- rownames(as.data.frame(cds@reducedDims$UMAP[filter, ]))
print(length(root_cells))

cds = order_cells(cds, root_cells = root_cells)
cds = order_cells(cds, 
                root_pr_nodes = unique(cds@principal_graph_aux[["UMAP"]]$root_pr_nodes))

options(repr.plot.width=5, repr.plot.height=4)
plot_cells(cds, cell_size = 1.2,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

cds_plot = cds
cds_plot$pseudotime = pseudotime(cds_plot)

options(repr.plot.width=4, repr.plot.height=3)
ggplot(as.data.frame(colData(cds_plot)), 
       aes(x = pseudotime, y = Condition)) + 
    geom_density_ridges2(aes(height = ..density.. , fill = Condition)) + 
    labs(x = "Pseudotime", y = "Time (hrs)") +
    theme_ridges() + theme(legend.position = "none", axis.text.x = element_blank())


saveRDS(cds, "cds_conv_pseudotime.rds")


## differential expression analysis
de_analysis = function(cds_file, cds_hash_file, min_expr_cutoff = 1, 
                           output_folder = ".", full_model_formula = "~1", num_cells = 50,
                           resid_model_formula = "~1", pattern = "", ncores = 1) {
    
    cds = readRDS(cds_file)
    cds_hash = readRDS(cds_hash_file)
    
    cds$Pseudotime = pseudotime(cds)
    cds_hash$Pseudotime = pseudotime(cds_hash)
    
    cds <- detect_genes(cds, min_expr = min_expr_cutoff)
    expressed_genes <- row.names(subset(cds, num_cells_expressed >= num_cells))
    print(paste0("# of expressed genes: ", length(expressed_genes)))
    

    degs <- fit_models(cds[expressed_genes,], 
                   model_formula_str = full_model_formula, 
                   expression_family = "negbinomial",
                   reduction_method="UMAP",
                   cores = ncores,
                   clean_model = TRUE,
                   verbose = FALSE)

    degs_resid <- fit_models(cds[expressed_genes,], 
               model_formula_str = resid_model_formula, 
               expression_family = "negbinomial",
               reduction_method="UMAP",
               cores = ncores,
               clean_model = TRUE,
               verbose = FALSE)

    degs_compared = compare_models(degs, degs_resid)

    print("sf done!")

    degs = degs %>% filter(status == "OK")
    degs_resid = degs_resid %>% filter(status == "OK")

    saveRDS(coefficient_table(degs) %>% select(-c("model", "model_summary")), paste0(output_folder, "/coeff_table_sf.rds"))
    saveRDS(evaluate_fits(degs) %>% select(-c("model", "model_summary")), paste0(output_folder, "/eval_fits_sf.rds"))

    saveRDS(degs_compared, paste0(output_folder, "/degs_sf_LRT.rds"))

    print("sf save done!")
    
    cds_hash = detect_genes(cds_hash, min_expr = min_expr_cutoff)
    
    degs_hash <- fit_models(cds_hash[expressed_genes,], 
                   model_formula_str = full_model_formula, 
                   expression_family = "negbinomial",
                   reduction_method="UMAP",
                   cores = ncores,
                   clean_model = TRUE,
                   verbose = FALSE)
    
    degs_hash_resid <- fit_models(cds_hash[expressed_genes,], 
               model_formula_str = resid_model_formula, 
               expression_family = "negbinomial",
               reduction_method="UMAP",
               cores = ncores,
               clean_model = TRUE,
               verbose = FALSE)
    
    degs_compared_hash = compare_models(degs_hash, degs_hash_resid)
    
    print("Hash done!")
    
    degs_hash = degs_hash %>% filter(status == "OK")
    degs_hash_resid = degs_hash_resid %>% filter(status == "OK")
    
    saveRDS(coefficient_table(degs_hash) %>% select(-c("model", "model_summary")), paste0(output_folder, "/coeff_table_hash.rds"))
    saveRDS(evaluate_fits(degs_hash) %>% select(-c("model", "model_summary")), paste0(output_folder, "/eval_fits_hash.rds"))
    
    saveRDS(degs_compared_hash, paste0(output_folder, "/degs_hash_LRT.rds"))

    print("Hash save done!")
}

de_analysis(cds_file="cds_conv_pseudotime.rds", cds_hash_file="cds_hash_ladder_pseudotime.rds", 
    min_expr_cutoff = 1, num_cells = 30, output_folder = 'de_analysis', 
    full_model_formula = "~ splines::ns(Pseudotime, df=3) + Plate", 
    resid_model_formula = "~1", pattern = "", ncores = 1)


