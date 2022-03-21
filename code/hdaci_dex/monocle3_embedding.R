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
    library(piano)
    library(snowfall)
    library(VennDiagram)
    library(pheatmap)
    library(stringr)
    library(eulerr)
    library(scales)
})

setwd("../../data/hdaci_dex")

## Make cds (filter by total RNA and remove mito genes)
cds_dex = readRDS("cds_dim_2_combined.rds")
metadata = as.data.frame(colData(cds_dex))

cutoff = 1500

mito_genes <- grep(pattern = "^MT-", x = rowData(cds_dex)$gene_short_name, value = TRUE, ignore.case = TRUE)
ribo_genes <- grep(pattern = "^Rp[sl][[:digit:]]", x = rowData(cds_dex)$gene_short_name, value = TRUE)

mito_gene_counts = colSums(counts(cds_dex)[rowData(cds_dex)$gene_short_name %in% mito_genes,])
mito_gene_id = rowData(cds_dex)[rowData(cds_dex)$gene_short_name %in% mito_genes,]$id

cds_dex$Mito_fraction = mito_gene_counts / cds_dex$Total_RNA

cds_dex = cds_dex[!rowData(cds_dex)$id %in% mito_gene_id, ]
cds_dex$Total_RNA = Matrix::colSums(counts(cds_dex))
cds_dex = cds_dex[, cds_dex$Total_RNA >= cutoff]

cds_dex = estimate_size_factors(cds_dex)
cds_dex = cds_dex[, cds_dex$Intercept < 0]
cds_dex

## Hash ladder normalized data
cds_dex_hash = cds_dex

scale_factors = log(cds_dex$Total_hash / cds_dex$Hash_duplication) * cds_dex$Slope / -cds_dex$Intercept *
                cds_dex$RNA_duplication
summary(scale_factors)
cds_dex_hash$Size_Factor = scale_factors / exp(mean(log(scale_factors))) 

summary(cds_dex_hash$Size_Factor)

saveRDS(cds_dex, "cds_dex_sf.rds")
saveRDS(cds_dex_hash, "cds_dex_hash.rds")


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


## Align cells (Conventional size factor normalization)
cds_dex = preprocess_cds(cds_dex, num_dim = 20)
cds_dex = align_cds(cds_dex, alignment_k = 20,
                    residual_model_formula_str = "~ log(Total_RNA)",
                    alignment_group = "Plate")

cds_dex = reduce_dimension(cds_dex, reduction_method = "UMAP", 
                       umap.n_neighbors = 50, umap.min_dist = 0.001)

options(repr.plot.width=4, repr.plot.height=3)
plot_cells(cds_dex, color_cells_by = "Dose", label_cell_groups = F,
           cell_size = 1.2, group_label_size=4)
    
## Identify mislabeled cells
## Here, I define "Acetylation" = 4hrs of HDACi and
## "Metabolic" = 24hrs of HDACi
options(repr.plot.width=3, repr.plot.height=3)
cds_dex = cluster_cells(cds_dex, resolution=0.0001, partition_qval = 0.5)
plot_cells(cds_dex, color_cells_by = "cluster", cell_size = 1)

cds_dex$Cluster = clusters(cds_dex)
cds_dex$Cluster = factor(cds_dex$Cluster, levels = c('3', '1', '2'))
cds_dex$Condition = ifelse(cds_dex$Cluster == '3', "DMSO", 
                          ifelse(cds_dex$Cluster == '1', "Acetylation", "Metabolic"))
cds_dex$Condition = factor(cds_dex$Condition, levels = c("DMSO", "Acetylation", "Metabolic"))

saveRDS(cds_dex, "cds_dex_sf.rds")


## Align cells (Hash ladder normalization)
cds_dex_hash = preprocess_cds(cds_dex_hash, num_dim = 20)
cds_dex_hash = align_cds(cds_dex_hash, alignment_k = 20,
                    residual_model_formula_str = "~ log(Total_RNA)",
                    alignment_group = "Plate")

cds_dex_hash = reduce_dimension(cds_dex_hash, reduction_method = "UMAP", 
                       umap.n_neighbors = 50, umap.min_dist = 0.001)

options(repr.plot.width=3, repr.plot.height=3)
plot_cells(cds_dex_hash, color_cells_by = "Dose", label_cell_groups = F,
           cell_size = 1.2, group_label_size=4) 
    
## Identify mislabeled cells
## Here, I define "Acetylation" = 4hrs of HDACi and
## "Metabolic" = 24hrs of HDACi
options(repr.plot.width=3, repr.plot.height=3)
cds_dex_hash = cluster_cells(cds_dex_hash, resolution=0.0001, partition_qval = 0.5)
plot_cells(cds_dex_hash, color_cells_by = "cluster", cell_size = 1)

cds_dex_hash$Cluster = clusters(cds_dex_hash)
cds_dex_hash$Cluster = factor(cds_dex_hash$Cluster, levels = c('3', '1', '2'))
cds_dex_hash$Condition = ifelse(cds_dex_hash$Cluster == '3', "DMSO", 
                          ifelse(cds_dex_hash$Cluster == '1', "Acetylation", "Metabolic"))
cds_dex_hash$Condition = factor(cds_dex_hash$Condition, levels = c("DMSO", "Acetylation", "Metabolic"))

saveRDS(cds_dex_hash, "cds_dex_hash.rds")


## UMAP plots
cds_dex = readRDS("cds_dex_sf.rds")
cds_dex_hash = readRDS("cds_dex_hash.rds")

cds_dex_hash = detect_genes(cds_dex_hash, 1)
genes = subset(rowData(cds_dex_hash), num_cells_expressed > 50)$id

cds_dex_hash = preprocess_cds(cds_dex_hash, num_dim = 50, use_genes = genes)
cds_dex_hash = align_cds(cds_dex_hash, alignment_k = 20,
                    residual_model_formula_str = "~ log(Total_RNA)")

cds_dex_hash = reduce_dimension(cds_dex_hash, reduction_method = "UMAP", 
                       umap.n_neighbors = 50, umap.min_dist = 0.01)

## Figure 4b
df_hash = data.frame(UMAP1 = cds_dex_hash@reducedDims$UMAP[,1],
                     UMAP2 = cds_dex_hash@reducedDims$UMAP[,2], 
                     Dose = cds_dex_hash$Dose,
                     Dex = cds_dex_hash$Dex,
                     Time = cds_dex_hash$Time)

df_hash[df_mult$Dose == 0, ]$Time = 0
df_hash$Time = as.factor(df_hash$Time)


cbp2 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

g = list()
options(repr.plot.width=9, repr.plot.height=3)
g[[3]] = ggplot(df_hash, aes(UMAP1, UMAP2)) +
    geom_point(aes(color = Dex)) + 
    labs(x = "", y = "")  +
    scale_color_manual(values = c("grey90", "royalblue4")) + 
    monocle_theme_opts() + theme(legend.position = "none", axis.text = element_blank())

g[[2]] = ggplot(df_hash, aes(UMAP1, UMAP2)) +
    geom_point(aes(color = Dose)) + 
    labs(x = "", y = "")  +
    scale_color_manual(values = cbp2) + 
    monocle_theme_opts() + theme(legend.position = "none", axis.text = element_blank())

g[[1]] = ggplot(df_hash, aes(UMAP1, UMAP2)) +
    geom_point(aes(color = Time)) + 
    labs(x = "", y = "")  +
    scale_color_brewer(palette = "Set1") + 
    monocle_theme_opts() + theme(legend.position = "none", axis.text = element_blank())

gg = do.call("grid.arrange", c(g, ncol=3))

ggsave("R_dex_umap_all.png", plot = gg, device = "png", width = 9, height = 3)


# Supplementary Figure 11
get_hash_factors = function(cds) {
    scale_factors = log(cds$Total_hash / cds$Hash_duplication) * 
                    cds$Slope / -cds$Intercept * cds$RNA_duplication
    cds$Size_Factor = scale_factors / exp(mean(log(scale_factors))) 
    return(cds)
}

## 0uM HDACi
cds_dex_hash_0 = cds_dex_hash[, cds_dex_hash$Dose == '0' & cds_dex_hash$Condition == "DMSO"]
cds_dex_hash_0 = get_hash_factors(cds_dex_hash_0)

cds_dex_hash_0 = detect_genes(cds_dex_hash_0, 1)
genes = subset(rowData(cds_dex_hash_0), num_cells_expressed > 20)$id

cds_dex_hash_0 = preprocess_cds(cds_dex_hash_0, num_dim = 50, use_genes = genes)
cds_dex_hash_0 = align_cds(cds_dex_hash_0, alignment_k = 10,
                    residual_model_formula_str = "~ log(Total_RNA)")

cds_dex_hash_0 = reduce_dimension(cds_dex_hash_0, reduction_method = "UMAP", 
                       umap.n_neighbors = 50, umap.min_dist = 0.01)

df_hash = data.frame(UMAP1 = cds_dex_hash_0@reducedDims$UMAP[,1],
                     UMAP2 = cds_dex_hash_0@reducedDims$UMAP[,2], 
                     Dose = cds_dex_hash_0$Dose,
                     Dex = cds_dex_hash_0$Dex)

options(repr.plot.width=3, repr.plot.height=3)
ggplot(df_hash, aes(UMAP1, UMAP2)) +
    geom_point(aes(color = Dex)) + 
    labs(x = "", y = "")  +
    #ggtitle("No HDACi") + 
    scale_color_manual(values = c("grey90", "royalblue4")) + 
    monocle3:::monocle_theme_opts()


## 4hrs HDACi
cds_dex_hash_4 = cds_dex_hash[, cds_dex_hash$Time == 4 & cds_dex_hash$Dose != '0']
cds_dex_hash_4 = get_hash_factors(cds_dex_hash_4)

cds_dex_hash_4 = detect_genes(cds_dex_hash_4, 1)
cds_dex_hash_4 = preprocess_cds(cds_dex_hash_4, num_dim = 50, use_genes = genes)
cds_dex_hash_4 = align_cds(cds_dex_hash_4, alignment_k = 20,
                    residual_model_formula_str = "~ log(Total_RNA)")

cds_dex_hash_4 = reduce_dimension(cds_dex_hash_4, reduction_method = "UMAP", 
                       umap.n_neighbors = 20, umap.min_dist = 0.01)

df_hash = data.frame(UMAP1 = cds_dex_hash_4@reducedDims$UMAP[,1],
                     UMAP2 = cds_dex_hash_4@reducedDims$UMAP[,2], 
                     Dose = cds_dex_hash_4$Dose,
                     Dex = cds_dex_hash_4$Dex)

options(repr.plot.width=3, repr.plot.height=3)
ggplot(df_hash, aes(UMAP1, UMAP2)) +
    geom_point(aes(color = Dex)) + 
    labs(x = "", y = "")  +
    scale_color_manual(values = c("grey90", "royalblue4")) + 
    monocle3:::monocle_theme_opts()

## 1 uM, 4hrs HDACi
cds_dex_hash_4 = cds_dex_hash[, cds_dex_hash$Time == 4 & cds_dex_hash$Dose == '1']
cds_dex_hash_4 = get_hash_factors(cds_dex_hash_4)

cds_dex_hash_4 = detect_genes(cds_dex_hash_4, 1)
cds_dex_hash_4 = preprocess_cds(cds_dex_hash_4, num_dim = 50, use_genes = genes)
cds_dex_hash_4 = align_cds(cds_dex_hash_4, alignment_k = 20,
                    residual_model_formula_str = "~ log(Total_RNA)")

cds_dex_hash_4 = reduce_dimension(cds_dex_hash_4, reduction_method = "UMAP", 
                       umap.n_neighbors = 20, umap.min_dist = 0.01)

df_hash = data.frame(UMAP1 = cds_dex_hash_4@reducedDims$UMAP[,1],
                     UMAP2 = cds_dex_hash_4@reducedDims$UMAP[,2], 
                     Dose = cds_dex_hash_4$Dose,
                     Dex = cds_dex_hash_4$Dex)

options(repr.plot.width=3, repr.plot.height=3)
ggplot(df_hash, aes(UMAP1, UMAP2)) +
    geom_point(aes(color = Dex)) + 
    labs(x = "", y = "")  +
    scale_color_manual(values = c("grey90", "royalblue4")) + 
    monocle3:::monocle_theme_opts() #+ theme(legend.position = "none", axis.text = element_blank())

## 10 uM, 4hrs HDACi
cds_dex_hash_4 = cds_dex_hash[, cds_dex_hash$Time == 4 & cds_dex_hash$Dose == '10']
cds_dex_hash_4 = get_hash_factors(cds_dex_hash_4)

cds_dex_hash_4 = detect_genes(cds_dex_hash_4, 1)
cds_dex_hash_4 = preprocess_cds(cds_dex_hash_4, num_dim = 50, use_genes = genes)
cds_dex_hash_4 = align_cds(cds_dex_hash_4, alignment_k = 20,
                    residual_model_formula_str = "~ log(Total_RNA)")

cds_dex_hash_4 = reduce_dimension(cds_dex_hash_4, reduction_method = "UMAP", 
                       umap.n_neighbors = 20, umap.min_dist = 0.01)

df_hash = data.frame(UMAP1 = cds_dex_hash_4@reducedDims$UMAP[,1],
                     UMAP2 = cds_dex_hash_4@reducedDims$UMAP[,2], 
                     Dose = cds_dex_hash_4$Dose,
                     Dex = cds_dex_hash_4$Dex)

options(repr.plot.width=3, repr.plot.height=3)
ggplot(df_hash, aes(UMAP1, UMAP2)) +
    geom_point(aes(color = Dex)) + 
    labs(x = "", y = "")  +
    scale_color_manual(values = c("grey90", "royalblue4")) + 
    monocle3:::monocle_theme_opts() #+ theme(legend.position = "none", axis.text = element_blank())

## 24hrs HDACi
cds_dex_hash_24 = cds_dex_hash[, cds_dex_hash$Time == 24 & cds_dex_hash$Dose != '0']
cds_dex_hash_24 = get_hash_factors(cds_dex_hash_24)

cds_dex_hash_24 = detect_genes(cds_dex_hash_24, 1)
cds_dex_hash_24 = preprocess_cds(cds_dex_hash_24, num_dim = 50, use_genes = genes)
cds_dex_hash_24 = align_cds(cds_dex_hash_24, alignment_k = 20,
                    residual_model_formula_str = "~ log(Total_RNA)")

cds_dex_hash_24 = reduce_dimension(cds_dex_hash_24, reduction_method = "UMAP", 
                       umap.n_neighbors = 30, umap.min_dist = 0.01)
df_hash = data.frame(UMAP1 = cds_dex_hash_24@reducedDims$UMAP[,1],
                     UMAP2 = cds_dex_hash_24@reducedDims$UMAP[,2], 
                     Dose = cds_dex_hash_24$Dose,
                     Dex = cds_dex_hash_24$Dex)

options(repr.plot.width=3, repr.plot.height=3)
ggplot(df_hash, aes(UMAP1, UMAP2)) +
    geom_point(aes(color = Dex)) + 
    labs(x = "", y = "")  +
    scale_color_manual(values = c("grey90", "royalblue4")) + 
    monocle3:::monocle_theme_opts()

## 1uM, 24hrs HDACi
cds_dex_hash_24 = cds_dex_hash[, cds_dex_hash$Time == 24 & cds_dex_hash$Dose == '1']
cds_dex_hash_24 = get_hash_factors(cds_dex_hash_24)

cds_dex_hash_24 = detect_genes(cds_dex_hash_24, 1)
cds_dex_hash_24 = preprocess_cds(cds_dex_hash_24, num_dim = 50, use_genes = genes)
cds_dex_hash_24 = align_cds(cds_dex_hash_24, alignment_k = 20,
                    residual_model_formula_str = "~ log(Total_RNA)")

cds_dex_hash_24 = reduce_dimension(cds_dex_hash_24, reduction_method = "UMAP", 
                       umap.n_neighbors = 20, umap.min_dist = 0.01)

df_hash = data.frame(UMAP1 = cds_dex_hash_24@reducedDims$UMAP[,1],
                     UMAP2 = cds_dex_hash_24@reducedDims$UMAP[,2], 
                     Dose = cds_dex_hash_24$Dose,
                     Dex = cds_dex_hash_24$Dex,
                     Drug = cds_dex_hash_24$Drug)

options(repr.plot.width=3, repr.plot.height=3)
ggplot(df_hash, aes(UMAP1, UMAP2)) +
    geom_point(aes(color = Dex)) + 
    labs(x = "", y = "")  +
    scale_color_manual(values = c("grey90", "royalblue4")) + 
    monocle3:::monocle_theme_opts() #+ theme(legend.position = "none", axis.text = element_blank())

## 10uM, 24hrs HDACi
cds_dex_hash_24 = cds_dex_hash[, cds_dex_hash$Time == 24 & cds_dex_hash$Dose == '10']
cds_dex_hash_24 = get_hash_factors(cds_dex_hash_24)

cds_dex_hash_24 = detect_genes(cds_dex_hash_24, 1)
cds_dex_hash_24 = preprocess_cds(cds_dex_hash_24, num_dim = 50, use_genes = genes)
cds_dex_hash_24 = align_cds(cds_dex_hash_24, alignment_k = 20,
                    residual_model_formula_str = "~ log(Total_RNA)")

cds_dex_hash_24 = reduce_dimension(cds_dex_hash_24, reduction_method = "UMAP", 
                       umap.n_neighbors = 20, umap.min_dist = 0.01)

df_hash = data.frame(UMAP1 = cds_dex_hash_24@reducedDims$UMAP[,1],
                     UMAP2 = cds_dex_hash_24@reducedDims$UMAP[,2], 
                     Dose = cds_dex_hash_24$Dose,
                     Dex = cds_dex_hash_24$Dex,
                     Drug = cds_dex_hash_24$Drug)

options(repr.plot.width=3, repr.plot.height=3)
ggplot(df_hash, aes(UMAP1, UMAP2)) +
    geom_point(aes(color = Dex)) + 
    labs(x = "", y = "")  +
    scale_color_manual(values = c("grey90", "royalblue4")) + 
    monocle3:::monocle_theme_opts() #+ theme(legend.position = "none", axis.text = element_blank())

