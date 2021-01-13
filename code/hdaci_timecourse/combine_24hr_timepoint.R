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
})

combine_cds = function(cds_batch_1, cds_batch_2, min_expr_cutoff = 1, num_cells = 0) {
    combined_mat = cbind(counts(cds_batch_1), counts(cds_batch_2))
    combined_colData = rbind(as.data.frame(colData(cds_batch_1)), as.data.frame(colData(cds_batch_2)))
    
    cds = new_cell_data_set(combined_mat, 
                            cell_metadata = combined_colData, 
                            gene_metadata = rowData(cds_batch_1))

    cds <- detect_genes(cds, min_expr = min_expr_cutoff)
    expressed_genes <- row.names(subset(cds, num_cells_expressed >= num_cells))
    print(paste0("# of expressed genes: ", length(expressed_genes)))

    cds <- cds[expressed_genes, ]
    
    return(cds)
}

setwd("../../data/hdaci_qc/")
cds_timecourse = readRDS("cds_dim_1_combined.rds")
cds_dex = readRDS("cds_dim_2_combined.rds")
cds_dex = cds_dex[, cds_dex$Dose == "10" & cds_dex$Time == 24 & cds_dex$Dex == F]

cds_dex

cds = combine_cds(cds_timecourse, cds_dex)
cds

saveRDS(cds, "cds_dim_1_24hr_combined.rds")

