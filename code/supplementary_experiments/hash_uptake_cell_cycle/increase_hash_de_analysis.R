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

hash_normalization = function(cds) {
    scale_factors_2 = log(cds$Total_hash / cds$Hash_duplication) * cds$RNA_duplication * cds$Slope / -cds$Intercept
    cds$Size_Factor = scale_factors_2 / exp(mean(log(scale_factors_2)))
    return(cds)
}

de_analysis = function(cds_hash_file, min_expr_cutoff = 1, cell_type = "A549", scale_factor = 1,
                           output_folder = ".", full_model_formula = "~1", num_cells = 50,
                           resid_model_formula = "~1", pattern = "", ncores = 1) {
    
    cds_hash = readRDS(cds_hash_file)
    cds_hash = cds_hash[, cds_hash$Cell_type_cluster == cell_type & cds_hash$Celltype == cell_type]
    
    metadata = as.data.frame(colData(cds_hash))
    total_hash = metadata[metadata$cycle == "G2M",]$Total_hash
    metadata[metadata$cycle == "G2M",]$Total_hash = total_hash * scale_factor
    cds_hash$Total_hash = metadata$Total_hash
    cds_hash = hash_normalization(cds_hash)
    
    cds_hash <- detect_genes(cds_hash, min_expr = min_expr_cutoff)
    expressed_genes <- row.names(subset(cds_hash, num_cells_expressed >= num_cells))
    print(paste0("# of expressed genes: ", length(expressed_genes)))
 
    
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
    
    degs_hash = degs_hash %>% filter(status == "OK")
    degs_hash_resid = degs_hash_resid %>% filter(status == "OK")
    
    saveRDS(coefficient_table(degs_hash) %>% select(-c("model", "model_summary")), 
            paste0(output_folder, "/coeff_table_hash_", pattern, "_", cell_type, ".rds"))
    saveRDS(evaluate_fits(degs_hash) %>% select(-c("model", "model_summary")), 
            paste0(output_folder, "/eval_fits_hash_", pattern, "_", cell_type, ".rds"))
    
    saveRDS(degs_compared_hash, paste0(output_folder, "/degs_hash_LRT_", pattern, "_", cell_type, ".rds"))

}

de_analysis(cds_mult_file="cds_hash_cluster_adjusted.rds", scale_factor = 5,
            min_expr_cutoff = 1, num_cells = 100, output_folder = 'de_analysis', 
            full_model_formula = "~ cycle + Batch", cell_type = "A549",
            resid_model_formula = "~1", pattern = "5x", ncores = 5)



get_de_genes = function(degs_df, degs_compared, q_value = qvalue, num_cells=30) {
    degs_compared = degs_compared %>% dplyr::filter(q_value < qvalue)
    degs_df = degs_df %>% 
        dplyr::filter(id %in% degs_compared$id & grepl("cycleG2M", term)) %>%
        dplyr::filter(q_value < qvalue & num_cells_expressed > num_cells) %>%
        arrange(id, estimate) %>%
        distinct(id, .keep_all = TRUE) %>% 
        arrange(q_value) %>%
        select(id, gene_short_name, estimate)
    
    return(degs_df)
}

## Retrieve DE genes
degs_control_compared = readRDS("de_analysis/degs_hash_LRT_control_A549.rds")
degs_1.2x_compared = readRDS("de_analysis/degs_hash_LRT_1.2x_A549.rds")
degs_1.5x_compared = readRDS("de_analysis/degs_hash_LRT_1.5x_A549.rds")
degs_2x_compared = readRDS("de_analysis/degs_hash_LRT_2x_A549.rds")
degs_3x_compared = readRDS("de_analysis/degs_hash_LRT_3x_A549.rds")
degs_5x_compared = readRDS("de_analysis/degs_hash_LRT_5x_A549.rds")

degs_control = readRDS("de_analysis/coeff_table_hash_control_A549.rds")
degs_1.2x = readRDS("de_analysis/coeff_table_hash_1.2x_A549.rds")
degs_1.5x = readRDS("de_analysis/coeff_table_hash_1.5x_A549.rds")
degs_2x = readRDS("de_analysis/coeff_table_hash_2x_A549.rds")
degs_3x = readRDS("de_analysis/coeff_table_hash_3x_A549.rds")
degs_5x = readRDS("de_analysis/coeff_table_hash_5x_A549.rds")

qvalue = 1e-2
num_cells = 100

degs_control = get_de_genes(degs_control, degs_control_compared, qvalue, num_cells)
degs_1.2x = get_de_genes(degs_1.2x, degs_1.2x_compared, qvalue, num_cells)
degs_1.5x = get_de_genes(degs_1.5x, degs_1.5x_compared, qvalue, num_cells)
degs_2x = get_de_genes(degs_2x, degs_2x_compared, qvalue, num_cells)
degs_3x = get_de_genes(degs_3x, degs_3x_compared, qvalue, num_cells)
degs_5x = get_de_genes(degs_5x, degs_5x_compared, qvalue, num_cells)

degs_list = list(degs_control, degs_1.2x, degs_1.5x, degs_2x, degs_3x, degs_5x)

## Look at genes in all dfs
joined_dfs = plyr::join_all(degs_list, by = c("id", "gene_short_name"))
joined_dfs = joined_dfs %>% na.omit()

n_overlap = colSums(!is.na(joined_dfs))
n_control = n_overlap[3]
n_overlap = n_overlap[3:length(n_overlap)]

rep_label = c("Control", "1.2x", "1.5x", "2x", "3x", "5x")
n_df = data.frame(n = n_overlap, 
                 increase = rep_label)
n_df$increase = factor(n_df$increase, levels = rep_label)

## Supplementary figure 17a
options(repr.plot.width=3, repr.plot.height=3)
ggplot(n_df, aes(x = increase, y = n / n_control * 100)) +  
    geom_bar(color = "black", stat = "identity") +
    monocle3:::monocle_theme_opts() + theme(legend.position = "none")

ggsave("R_hash_increase_deg_overlap_bar.pdf", device = "pdf", width = 3, height = 3)


fc_1.2x = joined_dfs[,4] / joined_dfs[,3]
fc_1.5x = joined_dfs[,5] / joined_dfs[,3]
fc_2x = joined_dfs[,6] / joined_dfs[,3]
fc_3x = joined_dfs[,7] / joined_dfs[,3]
fc_5x = joined_dfs[,8] / joined_dfs[,3]

rep_label = c("1.2x", "1.5x", "2x", "3x", "5x")
rep_n = c(length(fc_1.2x), length(fc_1.5x), length(fc_2x), length(fc_3x), length(fc_5x))
fc_list = c(fc_1.2x, fc_1.5x, fc_2x, fc_3x, fc_5x)
fc_df = data.frame(fc = fc_list, 
                   increase = rep(rep_label, rep_n))

## Supplementary figure 17b
options(repr.plot.width=3, repr.plot.height=3)
ggplot(fc_df, aes(x = increase, y = fc)) +  
    geom_boxplot(aes(fill = increase)) +
    geom_hline(yintercept = 1) + 
    monocle3:::monocle_theme_opts() + theme(legend.position = "none")

ggsave("R_hash_increase_deg_estimate_fc.pdf", device = "pdf", width = 3, height = 3)
