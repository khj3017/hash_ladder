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

## load data and filter cells based on RNA UMI/cell
setwd("../../data/hdaci_timecourse/")
cds_timecourse = readRDS("cds_dim_1_24hr_combined.rds")
metadata = as.data.frame(colData(cds_timecourse))

dim(metadata)
head(metadata)

cutoff = 1500

mito_genes <- grep(pattern = "^MT-", x = rowData(cds_timecourse)$gene_short_name, value = TRUE, ignore.case = TRUE)

mito_gene_counts = colSums(counts(cds_timecourse)[rowData(cds_timecourse)$gene_short_name %in% mito_genes,])
mito_gene_id = rowData(cds_timecourse)[rowData(cds_timecourse)$gene_short_name %in% mito_genes,]$id

cds_timecourse$Mito_fraction = mito_gene_counts / cds_timecourse$Total_RNA

cds_timecourse = cds_timecourse[!rowData(cds_timecourse)$id %in% mito_gene_id, ]
cds_timecourse$Total_RNA = Matrix::colSums(counts(cds_timecourse))
cds_timecourse = cds_timecourse[, cds_timecourse$Total_RNA >= cutoff]

cds_timecourse = estimate_size_factors(cds_timecourse)
cds_timecourse

cds_timecourse_hash = cds_timecourse

scale_factors = log(cds_timecourse$Total_hash / cds_timecourse$Hash_duplication) * 
                    cds_timecourse$Slope / -cds_timecourse$Intercept *
                    cds_timecourse$RNA_duplication
cds_timecourse_hash$Size_Factor = scale_factors / exp(mean(log(scale_factors))) 

summary(cds_timecourse_hash$Size_Factor)
metadata = as.data.frame(colData(cds_timecourse_hash))

## Align cells
cds_timecourse = preprocess_cds(cds_timecourse, num_dim = 50)
cds_timecourse = align_cds(cds_timecourse, alignment_k = 20,
                           residual_model_formula_str = "~ log(Total_RNA)",
                           alignment_group = "Plate")

options(repr.plot.width=4, repr.plot.height=4)
cds_timecourse = reduce_dimension(cds_timecourse, reduction_method = "UMAP",  
                       umap.n_neighbors = 50, umap.min_dist = 0.1)
plot_cells(cds_timecourse, color_cells_by = "Time", cell_size = 1.2) + ggtitle("Conventional normalization")

cds_timecourse = cluster_cells(cds_timecourse, resolution=0.0001)
plot_cells(cds_timecourse, color_cells_by = "cluster", cell_size = 1)

cds_timecourse <- learn_graph(cds_timecourse, use_partition = F, 
                              learn_graph_control = list(ncenter=360, maxiter = 20, minimal_branch_len = 20))
plot_cells(cds_timecourse, color_cells_by = "Time", cell_size=1)

filter = cds_timecourse@reducedDims$UMAP[,1] > 1.85 #&
        #cds_timecourse@reducedDims$UMAP[,2] > 2.5
root_cells <- rownames(as.data.frame(cds_timecourse@reducedDims$UMAP[filter, ]))
print(length(root_cells))

cds_timecourse = order_cells(cds_timecourse, root_cells = root_cells)
cds_timecourse = order_cells(cds_timecourse, 
                        root_pr_nodes = unique(cds_timecourse@principal_graph_aux[["UMAP"]]$root_pr_nodes))

options(repr.plot.width=4, repr.plot.height=3)
plot_cells(cds_timecourse, cell_size = 1.2,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

cds_plot = cds_timecourse
cds_plot$pseudotime = pseudotime(cds_plot)

options(repr.plot.width=4, repr.plot.height=3)
ggplot(as.data.frame(colData(cds_plot)), 
       aes(x = pseudotime, y = as.factor(Time))) + 
    geom_density_ridges2(aes(height = ..density.. , fill = Time)) + 
    labs(x = "Pseudotime", y = "Time (hrs)") +
    viridis::scale_fill_viridis(option = "C") +
    theme_ridges() + theme(legend.position = "none", axis.text.x = element_blank())


saveRDS(cds_timecourse, "cds_timecourse_pseudotime.rds")


## Align cells
cds_timecourse_hash = preprocess_cds(cds_timecourse_hash, num_dim = 70)
cds_timecourse_hash = align_cds(cds_timecourse_hash, alignment_k = 20,
                           residual_model_formula_str = "~ log(Total_RNA)",
                           alignment_group = "Plate")


options(repr.plot.width=5, repr.plot.height=4)
cds_timecourse_hash = reduce_dimension(cds_timecourse_hash, reduction_method = "UMAP", 
                       umap.n_neighbors = 100, umap.min_dist = 0.001) # 50: 50, 0.01
plot_cells(cds_timecourse_hash, color_cells_by = "Time", cell_size = 1.2, group_label_size = 4) + 
                ggtitle("Hash ladder normalization")

options(repr.plot.width=5, repr.plot.height=4)
cds_timecourse_hash = reduce_dimension(cds_timecourse_hash, reduction_method = "UMAP", 
                       umap.n_neighbors = 200, umap.min_dist = 0.01) # 50: 200, 0.01
plot_cells(cds_timecourse_hash, color_cells_by = "Time", cell_size = 1.2, group_label_size = 4) + 
                ggtitle("Hash ladder normalization")

options(repr.plot.width=4, repr.plot.height=4)
cds_timecourse_hash = cluster_cells(cds_timecourse_hash, resolution=0.0001)
plot_cells(cds_timecourse_hash, color_cells_by = "cluster", cell_size = 1, group_label_size = 4)

options(repr.plot.width=4, repr.plot.height=3)
cds_timecourse_hash <- learn_graph(cds_timecourse_hash, use_partition = F, close_loop = F,
                                    learn_graph_control = list(ncenter=380, maxiter = 10, minimal_branch_len = 10))
plot_cells(cds_timecourse_hash, color_cells_by = "Time", cell_size=1, alpha = 0.3)

## get root cell names
filter = cds_timecourse_hash@reducedDims$UMAP[,1] > 3 #& # -2.25, -1
            #cds_timecourse_hash@reducedDims$UMAP[,2] < -1
root_cells <- rownames(as.data.frame(cds_timecourse_hash@reducedDims$UMAP[filter, ]))
print(length(root_cells))

## assign pseudotime
cds_timecourse_hash = order_cells(cds_timecourse_hash, root_cells = root_cells)
cds_timecourse_hash = order_cells(cds_timecourse_hash,
                        root_pr_nodes = unique(cds_timecourse_hash@principal_graph_aux[["UMAP"]]$root_pr_nodes))


options(repr.plot.width=4, repr.plot.height=3)
plot_cells(cds_timecourse_hash, cell_size = 1.2,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

cds_plot = cds_timecourse_hash
cds_plot$pseudotime = pseudotime(cds_plot)

options(repr.plot.width=3, repr.plot.height=3)
ggplot(as.data.frame(colData(cds_plot)), 
       aes(x = pseudotime, y = as.factor(Time))) + 
    geom_density_ridges2(aes(height = ..density.. , fill = Time)) + 
    labs(x = "Pseudotime", y = "Time (hrs)") +
    viridis::scale_fill_viridis(option = "C") +
    theme_ridges() + theme(legend.position = "none", axis.text.x = element_blank())


saveRDS(cds_timecourse_hash, "cds_timecourse_hash_pseudotime.rds")


## DE analysis
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

de_analysis(cds_file="cds_timecourse_pseudotime.rds", cds_hash_file="cds_timecourse_hash_pseudotime.rds", 
    min_expr_cutoff = 1, num_cells = 50, output_folder = 'de_analysis', 
    full_model_formula = "~ splines::ns(Pseudotime, df=3) + Plate + RNA_duplication", 
    resid_model_formula = "~1", pattern = "", ncores = 10)
