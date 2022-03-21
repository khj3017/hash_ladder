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
    library(readxl)
    library(scales)
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

setwd("../../data/hdaci_dex/")

overall_estimate_df = readRDS("overall_estimate_df_by_time.rds")
no_response_genes = overall_estimate_df %>%
    filter(normalization == "Hash") %>%
    filter(respond == "No") %>%
    mutate(type = ifelse(hdaci_estimate == 0, "Attenuated", 
                    ifelse(sign(hdaci_estimate) == sign(dex_estimate), "Saturated", "Dominated"))) %>%
    transform(type = factor(type, levels = c("Attenuated", "Saturated", "Dominated")))


acetyl_genes = no_response_genes %>%
    filter(time == 4)
meta_genes = no_response_genes %>%
    filter(time != 4)


meta_genes_ordered = meta_genes %>%
    arrange(type, hdaci_estimate) %>%
    select(id, gene_short_name, type, hdaci_estimate, dex_estimate)

meta_gene_matrix = meta_genes_ordered %>%
                        select(hdaci_estimate, dex_estimate)

acetyl_genes_ordered = acetyl_genes %>%
    arrange(type, hdaci_estimate) %>%
    select(id, gene_short_name, type, hdaci_estimate, dex_estimate)

acetyl_gene_matrix = acetyl_genes_ordered %>%
                        select(hdaci_estimate, dex_estimate)

acetyl_genes_ordered %>% 
                count(type) %>%
                arrange(type)
type_counts$nn = c(type_counts$n[1], sum(type_counts$n[1:2]), sum(type_counts$n[1:3]))
bks = seq(-1.5, 1.5, length.out = 100)

## 4hours
options(repr.plot.width=3, repr.plot.height=15)
rownames(acetyl_gene_matrix) = acetyl_genes_ordered$gene_short_name
annotation_row = data.frame(acetyl_genes_ordered$type)
rownames(annotation_row) = rownames(acetyl_gene_matrix)
colnames(annotation_row) = " "

anno_colors <- list(type = RColorBrewer::brewer.pal(3, "Set1"))
names(anno_colors$type) <- unique(acetyl_genes_ordered$type)

ph_res = pheatmap(acetyl_gene_matrix, useRaster = T, cluster_cols = F,
            cluster_rows = F, show_rownames = T, show_colnames = F, silent = TRUE,
            #filename = "R_unresponsive_dex_acetylation_heatmap.pdf", width = 3, height = 15,
            gaps_row = type_counts$nn, breaks = bks)
grid::grid.rect(gp = grid::gpar("fill", col = NA))
grid::grid.draw(ph_res$gtable)


type_counts = meta_genes_ordered %>% 
                count(type) %>%
                arrange(type)
type_counts$nn = c(type_counts$n[1], sum(type_counts$n[1:2]), sum(type_counts$n[1:3]))

## 24hours
options(repr.plot.width=3, repr.plot.height=15)
rownames(meta_gene_matrix) = meta_genes_ordered$gene_short_name
annotation_row = data.frame(type = meta_genes_ordered$type)
rownames(annotation_row) = rownames(meta_gene_matrix)

anno_colors <- list(type = RColorBrewer::brewer.pal(3, "Set1"))
names(anno_colors$type) <- unique(meta_genes_ordered$type)

ph_res = pheatmap(meta_gene_matrix, useRaster = T, cluster_cols = F,
            cluster_rows = F, show_rownames = T, show_colnames = F, silent = TRUE,
            #filename = "R_unresponsive_dex_metabolic_heatmap.pdf", width = 3, height = 15,
            gaps_row = type_counts$nn, breaks = bks)
grid::grid.rect(gp = grid::gpar("fill", col = NA))
grid::grid.draw(ph_res$gtable)
