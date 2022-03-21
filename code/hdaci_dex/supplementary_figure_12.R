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

## Load cds
setwd("../../data/hdaci_dex")
cds_dex = readRDS("cds_dex_sf.rds")
cds_dex = cds_dex[, cds_dex$Dose == 0 & cds_dex$Condition == "DMSO"]
cds_dex_hash = readRDS("cds_dex_hash.rds")
cds_dex_hash = cds_dex_hash[, cds_dex_hash$Dose == 0 & cds_dex_hash$Condition == "DMSO"]

## Load DEGs from the Reddy et al. 2009 Genome Res paper
reddy_degs = read.table("reddy_dex_degs_new.txt", header = T)

## Load de genes
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


# Supplementary Figure 12
## Supplementary figure 12a
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


## Supplementary figure 12b
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


## Supplementary figure 12c
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


## Supplementary figure 12d
## Modified plot_percent_cells_positive function from monocle3
plot_percent_cells_positive_dist <- function(cds_subset,
                                        group_cells_by = NULL,
                                        min_expr = 0,
                                        normalize = TRUE,
                                        bootstrap_samples=100,
                                        conf_int_alpha = .95){

  marker_exprs <- SingleCellExperiment::counts(cds_subset)

  if (normalize) {
    marker_exprs <- Matrix::t(Matrix::t(marker_exprs)/size_factors(cds_subset))
    marker_exprs_melted <- reshape2::melt(round(10000*as.matrix(marker_exprs))/10000)
  } else {
    marker_exprs_melted <- reshape2::melt(as.matrix(marker_exprs))
  }

  colnames(marker_exprs_melted) <- c("f_id", "Cell", "expression")

  marker_exprs_melted <- base::merge(marker_exprs_melted, colData(cds_subset),
                               by.x="Cell", by.y="row.names")
  marker_exprs_melted <- base::merge(marker_exprs_melted, rowData(cds_subset),
                               by.x="f_id", by.y="row.names")

  if (!is.null(marker_exprs_melted$gene_short_name)){
    marker_exprs_melted$feature_label <- marker_exprs_melted$gene_short_name
    marker_exprs_melted$feature_label[
      is.na(marker_exprs_melted$feature_label)] <- marker_exprs_melted$f_id
  } else {
    marker_exprs_melted$feature_label <- marker_exprs_melted$f_id
}


  marker_counts_bootstrap = rsample::bootstraps(marker_exprs_melted, times = bootstrap_samples)

  group_mean_bootstrap <- function(split) {
    rsample::analysis(split) %>%
      dplyr::group_by(!!as.name("feature_label"), !!as.name(group_cells_by)) %>%
      dplyr::summarize(target = sum(expression > min_expr),
                       target_fraction = sum(expression > min_expr)/dplyr::n())
  }
  
    
  marker_counts <-
    marker_counts_bootstrap %>%
    dplyr::mutate(summary_stats = purrr::map(splits, group_mean_bootstrap)) %>%
    tidyr::unnest(summary_stats)
    
  return(marker_counts %>% dplyr::ungroup() %>%
            dplyr::group_by(!!as.name("feature_label"), !!as.name(group_cells_by)) %>%
            select(-splits) %>%
            mutate(target = target * 100,
                   target_fraction = target_fraction * 100))
}


genes = c("FGD4", "FKBP5", "CYP24A1", "ANGPTL4", "TSC22D3", "AHNAK")
ids_to_plot = convert_gene_to_id(cds_dex_hash, genes)

conf_int_alpha = 0.95

hash_plot_data = plot_percent_cells_positive_dist(cds_dex_hash[ids_to_plot,], min_expr = 1,
                                     group_cells_by="Dex")

hash_plot_data_sum = hash_plot_data %>%
                    dplyr::summarize(target_fraction_mean = mean(target_fraction),
                                     target_fraction_low = stats::quantile(target_fraction, (1 - conf_int_alpha) / 2),
                                     target_fraction_high = stats::quantile(target_fraction, 1 - (1 - conf_int_alpha) / 2))

options(repr.plot.width=8, repr.plot.height=2)
ggplot(hash_plot_data_sum, aes(x = as.factor(Dex), y = target_fraction_mean)) + 
        geom_bar(stat="identity", aes(fill = as.factor(Dex))) +
        geom_jitter(data = hash_plot_data, 
                    aes(x = as.factor(Dex), y = target_fraction),
                    position = position_jitter(0.3), size = 0.2,
                    color = "grey50") + 
        geom_linerange(aes(ymin = target_fraction_low, ymax = target_fraction_high)) +
        facet_wrap(~ feature_label, scales = "free", ncol=6) + 
        labs(x = "", y = "") + 
        scale_fill_manual(values = c("grey90", "royalblue4")) + 
        monocle3:::monocle_theme_opts() + theme(legend.position = "none")

ggsave("S_dex_vehicle_genes_jitter.pdf", device = "pdf", width=8, height=2)
