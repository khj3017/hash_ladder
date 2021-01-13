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

setwd("../../data/hdaci_dex/")

## DEX only vehicle
de_analysis_dex = function(cds_dex_file, cds_dex_hash_file, min_expr_cutoff = 1, 
                           num_cells_expr = 20, output_folder = ".", full_model_formula = "~1", 
                           resid_model_formula = "~1", pattern = "", ncores = 1) {
    
    
    cds_dex = readRDS(cds_dex_file)
    cds_dex = cds_dex[, cds_dex$Dose == 0]
    cds_dex_hash = readRDS(cds_dex_hash_file)
    cds_dex_hash = cds_dex_hash[, cds_dex_hash$Dose == 0]

    
    cds_dex <- detect_genes(cds_dex, min_expr = min_expr_cutoff)
    expressed_genes <- row.names(subset(cds_dex, num_cells_expressed >= num_cells_expr))
    print(paste0("# of expressed genes: ", length(expressed_genes)))
    

    degs_dex <- fit_models(cds_dex[expressed_genes,], 
                   model_formula_str = full_model_formula, 
                   expression_family = "negbinomial",
                   reduction_method="UMAP",
                   cores = ncores,
                   clean_model = TRUE,
                   verbose = FALSE)
    
    degs_dex_resid <- fit_models(cds_dex[expressed_genes,], 
               model_formula_str = resid_model_formula, 
               expression_family = "negbinomial",
               reduction_method="UMAP",
               cores = ncores,
               clean_model = TRUE,
               verbose = FALSE)
    
    degs_compared = compare_models(degs_dex, degs_dex_resid)
    
    print("sf done!")
    
    degs_dex = degs_dex %>% filter(status == "OK")
    degs_dex_resid = degs_dex_resid %>% filter(status == "OK")
    
    saveRDS(coefficient_table(degs_dex) %>% select(-c("model", "model_summary")), paste0(output_folder, "/coeff_table_sf.rds"))
    saveRDS(evaluate_fits(degs_dex) %>% select(-c("model", "model_summary")), paste0(output_folder, "/eval_fits_sf.rds"))
    
    saveRDS(degs_compared, paste0(output_folder, "/degs_sf_LRT.rds"))
    
    print("sf save done!")
    
    cds_dex_hash = detect_genes(cds_dex_hash, min_expr = min_expr_cutoff)
    
    
    degs_dex_hash <- fit_models(cds_dex_hash[expressed_genes,], 
                   model_formula_str = full_model_formula, 
                   expression_family = "negbinomial",
                   reduction_method="UMAP",
                   cores = ncores,
                   clean_model = TRUE,
                   verbose = FALSE)
    
    degs_dex_hash_resid <- fit_models(cds_dex_hash[expressed_genes,], 
               model_formula_str = resid_model_formula, 
               expression_family = "negbinomial",
               reduction_method="UMAP",
               cores = ncores,
               clean_model = TRUE,
               verbose = FALSE)
    
    degs_compared_hash = compare_models(degs_dex_hash, degs_dex_hash_resid)
    
    print("hash done!")
    
    degs_dex_hash = degs_dex_hash %>% filter(status == "OK")
    degs_dex_hash_resid = degs_dex_hash_resid %>% filter(status == "OK")
    
    saveRDS(coefficient_table(degs_dex_hash) %>% select(-c("model", "model_summary")), paste0(output_folder, "/coeff_table_hash.rds"))
    saveRDS(evaluate_fits(degs_dex_hash) %>% select(-c("model", "model_summary")), paste0(output_folder, "/eval_fits_hash.rds"))
    
    saveRDS(degs_compared_hash, paste0(output_folder, "/degs_hash_LRT.rds"))

    print("hash save done!")
}

de_analysis_dex(cds_dex_file="cds_dex_sf.rds", cds_dex_hash_file="cds_dex_hash.rds", 
                min_expr_cutoff = 1, output_folder = "de_analysis/dex_vehicle", 
                num_cells_expr = 20, full_model_formula = "~ Plate + Dex", 
                resid_model_formula = "~1", pattern = "", ncores = 10)



## DEX + HDACi
de_analysis_dex = function(cds_dex_file, cds_dex_hash_file, min_expr_cutoff = 1, 
                           output_folder = ".", full_model_formula = "~1", 
                           resid_model_formula = "~1", pattern = "", ncores = 1) {
    
    cond = "Acetylation"
    time = 4
    
    cds_dex = readRDS(cds_dex_file)
    cds_dex = cds_dex[, cds_dex$Condition == cond & cds_dex$Time == time]
    cds_dex_hash = readRDS(cds_dex_hash_file)
    cds_dex_hash = cds_dex_hash[, cds_dex_hash$Condition == cond & cds_dex_hash$Time == time]
    
    
    cds_dex <- detect_genes(cds_dex, min_expr = min_expr_cutoff)
    dex_genes = readRDS("dex_gene_list_sf.rds")

    degs_dex <- fit_models(cds_dex[dex_genes,], 
                   model_formula_str = full_model_formula, 
                   expression_family = "negbinomial",
                   reduction_method="UMAP",
                   cores = ncores,
                   clean_model = TRUE,
                   verbose = FALSE)
    
    degs_dex_resid <- fit_models(cds_dex[dex_genes,], 
               model_formula_str = resid_model_formula, 
               expression_family = "negbinomial",
               reduction_method="UMAP",
               cores = ncores,
               clean_model = TRUE,
               verbose = FALSE)
    
    degs_compared = compare_models(degs_dex, degs_dex_resid)
    
    print("sf done!")
    
    degs_dex = degs_dex %>% filter(status == "OK")
    degs_dex_resid = degs_dex_resid %>% filter(status == "OK")
    
    saveRDS(coefficient_table(degs_dex) %>% select(-c("model", "model_summary")), paste0(output_folder, "/coeff_table_sf_", pattern, ".rds"))
    saveRDS(evaluate_fits(degs_dex) %>% select(-c("model", "model_summary")), paste0(output_folder, "/eval_fits_sf_", pattern, ".rds"))
    
    saveRDS(degs_compared, paste0(output_folder, "/degs_sf_LRT_", pattern, ".rds"))
    
    print("sf save done!")
    
    cds_dex_hash = detect_genes(cds_dex_hash, min_expr = min_expr_cutoff)
    
    dex_genes = readRDS("dex_gene_list_hash.rds")
    
    degs_dex_hash <- fit_models(cds_dex_hash[dex_genes,], 
                   model_formula_str = full_model_formula, 
                   expression_family = "negbinomial",
                   reduction_method="UMAP",
                   cores = ncores,
                   clean_model = TRUE,
                   verbose = FALSE)
    
    degs_dex_hash_resid <- fit_models(cds_dex_hash[dex_genes,], 
               model_formula_str = resid_model_formula, 
               expression_family = "negbinomial",
               reduction_method="UMAP",
               cores = ncores,
               clean_model = TRUE,
               verbose = FALSE)
    
    degs_compared_hash = compare_models(degs_dex_hash, degs_dex_hash_resid)
    
    print("hash done!")
    
    degs_dex_hash = degs_dex_hash %>% filter(status == "OK")
    degs_dex_hash_resid = degs_dex_hash_resid %>% filter(status == "OK")
    
    saveRDS(coefficient_table(degs_dex_hash) %>% select(-c("model", "model_summary")), paste0(output_folder, "/coeff_table_hash_", pattern, ".rds"))
    saveRDS(evaluate_fits(degs_dex_hash) %>% select(-c("model", "model_summary")), paste0(output_folder, "/eval_fits_hash_", pattern, ".rds"))
    
    saveRDS(degs_dex_hash, paste0(output_folder, "/degs_hash_LRT_", pattern, ".rds"))

    print("hash save done!")
}

de_analysis_dex(cds_dex_file="cds_dex_sf.rds", cds_dex_hash_file="cds_dex_hash.rds", 
                min_expr_cutoff = 1, output_folder = "de_analysis/dex_time", 
                full_model_formula = "~ Plate + Drug + Dose + Dex", 
                resid_model_formula = "~1", pattern = "4hr_dex", ncores = 10)


## Collect DEX vehicle genes
coeffs_table_dex <- readRDS("de_analysis/dex_vehicle/coeff_table_sf.rds")
degs_compared = readRDS("de_analysis/dex_vehicle/degs_sf_LRT.rds")

qvalue = 5e-2
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

saveRDS(coeffs_table_dex$id, "dex_gene_list_sf.rds")
saveRDS(coeffs_table_dex_hash$id, "dex_gene_list_hash.rds")

