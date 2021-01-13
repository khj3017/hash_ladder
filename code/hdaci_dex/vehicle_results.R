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

convert_gene_to_id <- function(cds, genes) {
    rowData_mat = as.data.frame(rowData(cds))
    return(rowData_mat %>% filter(gene_short_name %in% genes) %>% pull(id))
}

convert_id_to_gene <- function(cds, ids) {
    rowData_mat = as.data.frame(rowData(cds))
    return(rowData_mat %>% filter(id %in% ids) %>% pull(gene_short_name))
}

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

## Load data
setwd("../../data/hdaci_dex")
cds_dex = readRDS("cds_dex_sf.rds")
cds_dex = cds_dex[, cds_dex$Dose == 0 & cds_dex$Condition == "DMSO"]
cds_dex_hash = readRDS("cds_dex_hash.rds")
cds_dex_hash = cds_dex_hash[, cds_dex_hash$Dose == 0 & cds_dex_hash$Condition == "DMSO"]

reddy_degs = read.table("reddy_dex_degs_new.txt", header = T)

coeffs_table_dex <- readRDS("de_analysis/dex_vehicle/coeff_table_sf.rds")
degs_compared = readRDS("de_analysis/dex_vehicle/degs_sf_LRT.rds")

qvalue = 1e-2
num_cells = 20

coeffs_table_dex = coeffs_table_dex %>% 
                    filter(id %in% degs_compared$id) %>%
                    filter(grepl("Dex", term) & q_value < qvalue & num_cells_expressed > num_cells)
dim(coeffs_table_dex)

coeffs_table_dex_hash <- readRDS("de_analysis/dex_vehicle/coeff_table_hash.rds")
degs_compared = readRDS("de_analysis/dex_vehicle/degs_hash_LRT.rds")
coeffs_table_dex_hash = coeffs_table_dex_hash %>% 
                                filter(id %in% degs_compared$id) %>%
                                filter(grepl("Dex", term) & q_value < qvalue & num_cells_expressed > num_cells) 
dim(coeffs_table_dex_hash)


## figure S9A
cds_dex_0 = cds_dex[, cds_dex$Dose == '0' & cds_dex$Condition == "DMSO"]
cds_dex_0 = estimate_size_factors(cds_dex_0)

cds_dex_0 = detect_genes(cds_dex_0, 1)
genes = subset(rowData(cds_dex_0), num_cells_expressed > 20)$id

cds_dex_0 = preprocess_cds(cds_dex_0, num_dim = 50, use_genes = genes)
cds_dex_0 = align_cds(cds_dex_0, alignment_k = 20, 
                    residual_model_formula_str = "~ log(Total_RNA)",
                     alignment_group = "Plate")

cds_dex_0 = reduce_dimension(cds_dex_0, reduction_method = "UMAP", 
                       umap.n_neighbors = 50, umap.min_dist = 0.01)

df_sf = data.frame(UMAP1 = cds_dex_0@reducedDims$UMAP[,1],
                     UMAP2 = cds_dex_0@reducedDims$UMAP[,2], 
                     Dose = cds_dex_0$Dose,
                     Dex = cds_dex_0$Dex)

options(repr.plot.width=3, repr.plot.height=3)
ggplot(df_sf, aes(UMAP1, UMAP2)) +
    geom_point(aes(color = Dex)) + 
    labs(x = "", y = "")  +
    scale_color_manual(values = c("grey90", "royalblue4")) + 
    monocle_theme_opts() + theme(legend.position = "none", axis.text = element_blank())

ggsave("S_conventional_dex_umap_vehicle.png", device = "png", width = 3, height = 3)


## figure S9B
n_common_genes = coeffs_table_dex %>%
    filter(id %in% coeffs_table_dex_hash$id) %>%
    nrow()

n_sf_genes = coeffs_table_dex %>%
    filter(!id %in% coeffs_table_dex_hash$id) %>%
    nrow()

n_hash_genes = coeffs_table_dex_hash %>%
    filter(!id %in% coeffs_table_dex$id) %>%
    nrow()

n_common_genes
n_sf_genes
n_hash_genes

options(repr.plot.width=6, repr.plot.height=3)

g = euler(c("Conventional" = n_sf_genes, "Hash ladder" = n_hash_genes, "Conventional&Hash ladder" = n_common_genes))

pdf("S_venn_dex_vehicle.pdf", width = 6, height = 3)
plot(g, 
     family = "Helveltica",
     fills = list(fill = c("yellow", "darkgreen"), alpha = 0.7),
     quantities = F, 
     labels = NULL, 
     legend = list(labels = c("Conventional", "Hash ladder"), fontsize = 16))
dev.off()


## figure S9C
degs_table_dex = coeffs_table_dex %>% 
            filter(grepl("Dex", term) & q_value < 0.01) %>% 
            count(Dex = term, Type = estimate > 0, name = "Count")

degs_table_dex$Type = ifelse(degs_table_dex$Type, "Upregulated", "Downregulated")
degs_table_dex$Normalization = "Conventional"

df = degs_table_dex


degs_table_dex_hash = coeffs_table_dex_hash %>% 
            filter(grepl("Dex", term) & q_value < 0.01) %>% 
            count(Dex = term, Type = estimate > 0, name = "Count")

degs_table_dex_hash$Type = ifelse(degs_table_dex_hash$Type, "Upregulated", "Downregulated")
degs_table_dex_hash$Normalization = "Hash ladder"

df = rbind(df, degs_table_dex_hash)

options(repr.plot.width=3.5, repr.plot.height=3.5)
bar_order = c("Upregulated", "Downregulated")
ggplot(transform(df, Type = factor(Type, levels = bar_order)),
    aes(x = Dex, y = Count, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +  
  labs(x = "", y = "# of DEGs") + 
  scale_fill_manual(values = c("#e41a1c", "#377eb8")) + 
  facet_wrap(vars(Normalization)) +  
  monocle_theme_opts() + 
  theme(strip.text.x = element_text(size=10), 
        axis.text.x = element_blank(), 
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")

ggsave("S_dex_vehicle_degs.pdf", device = "pdf", width=3.5, height=3.5)


## figure S9D
genes = c("FGD4", "FKBP5", "CYP24A1", "ANGPTL4", "TSC22D3", "AHNAK")
ids = convert_gene_to_id(cds_dex, genes)

options(repr.plot.width=8, repr.plot.height=2)
plot_percent_cells_positive(cds_dex_hash[ids,], min_expr = 1, normalize = T,
                  group_cells_by="Dex", ncol=6) + 
    scale_fill_manual(values = c("grey90", "royalblue4"))

ggsave("S_dex_vehicle_genes.pdf", device = "pdf", width=8, height=2)

