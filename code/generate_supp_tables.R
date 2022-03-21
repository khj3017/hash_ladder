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
    library(writexl)
})

## Supplementary Table 1
fp_data_dir = "../data/flavopiridol/de_analysis"

degs_compared = readRDS(file.path(fp_data_dir, "degs_sf_LRT.rds"))
degs_sf = readRDS(file.path(fp_data_dir, "coeff_table_sf.rds"))
qvalue = 1e-2
num_cells = 30

degs_compared = degs_compared %>% filter(q_value < qvalue)
degs_sf = degs_sf %>% 
    filter(id %in% degs_compared$id & grepl("Pseudo", term)) %>%
    filter(q_value < qvalue & num_cells_expressed > num_cells) %>%
    arrange(id, estimate) %>%
    distinct(id, .keep_all = TRUE) %>% 
    arrange(q_value)
dim(degs_sf)

degs_compared = readRDS(file.path(fp_data_dir, "degs_hash_LRT.rds"))
degs_hash = readRDS(file.path(fp_data_dir, "coeff_table_hash.rds"))

degs_compared = degs_compared %>% filter(q_value < qvalue)
degs_hash = degs_hash %>% 
    filter(id %in% degs_compared$id & grepl("Pseudo", term)) %>%
    filter(q_value < qvalue & num_cells_expressed > num_cells) %>%
    arrange(id, estimate) %>%
    distinct(id, .keep_all = TRUE) %>% 
    arrange(q_value)
dim(degs_hash)

fc_m_sf = readRDS(file.path(fp_data_dir, "fp_sf_fc_log2.rds"))
fc_m_hash = readRDS(file.path(fp_data_dir, "fp_hash_fc_log2.rds"))

cutoff = 1

genes_sf = fc_m_sf[abs(fc_m_sf) > cutoff]
genes_sf = genes_sf[names(genes_sf) %in% degs_sf$id]
degs_sf = degs_sf %>% filter(id %in% names(genes_sf)) %>% mutate(fc = genes_sf)

genes_hash = fc_m_hash[abs(fc_m_hash) > cutoff]
genes_hash = genes_hash[names(genes_hash) %in% degs_hash$id]
degs_hash = degs_hash %>% 
                filter(id %in% names(genes_hash)) %>% 
                mutate(fc = genes_hash)

degs_sf$normalization = "Conventional"
degs_hash$normalization = "Hash ladder"

dim(degs_sf)
dim(degs_hash)

degs_fp = rbind(degs_sf, degs_hash)
dim(degs_fp)

write.table(degs_fp, file = "../supp_tables/supplementary_table_1_fp_degs.txt", col.names = T, row.names = F, quote = F)
write_xlsx(degs_fp, "../supp_tables/supplementary_table_1_fp_degs.xlsx")


## Supplementary Table 2
hdaci_time_data_dir = "../data/hdaci_timecourse/de_analysis"

qvalue = 1e-2
num_cells = 100

degs_compared = readRDS(file.path(hdaci_time_data_dir, "degs_sf_LRT.rds"))
degs_sf = readRDS(file.path(hdaci_time_data_dir, "coeff_table_sf.rds"))

degs_compared = degs_compared %>% filter(q_value < qvalue)
degs_sf = degs_sf %>%
    filter(id %in% degs_compared$id & grepl("Pseudo", term)) %>%
    filter(num_cells_expressed > num_cells & q_value < qvalue) %>%
    distinct(id, .keep_all = TRUE)
dim(degs_sf)

degs_compared = readRDS(file.path(hdaci_time_data_dir, "degs_hash_LRT.rds"))
degs_hash = readRDS(file.path(hdaci_time_data_dir, "coeff_table_hash.rds"))

degs_compared = degs_compared %>% filter(q_value < qvalue)
degs_hash = degs_hash %>%  
    filter(id %in% degs_compared$gene_id & grepl("Pseudo", term)) %>%
    filter(num_cells_expressed > num_cells & q_value < qvalue) %>%
    distinct(id, .keep_all = TRUE) %>%
    select(-gene_id)


degs_sf$normalization = "Conventional"
degs_hash$normalization = "Hash ladder"

degs_hdaci = rbind(degs_sf, degs_hash)
dim(degs_hdaci)

write.table(degs_hdaci, file = "../supp_tables/supplementary_table_2_hdaci_timecourse_degs.txt", col.names = T, row.names = F, quote = F)
write_xlsx(degs_hdaci, "../supp_tables/supplementary_table_2_hdaci_timecourse_degs.xlsx")


## Supplementary Table 3
hdaci_data_dir = "../data/hdaci_timecourse/heatmaps"
go_list_hash = readRDS(file.path(hdaci_data_dir, "gsa_list_go_hash_num_cells_100_qvalue-1e-10_4.rds"))
hallmarks_list_hash = readRDS(file.path(hdaci_data_dir,"gsa_list_hallmarks_hash_num_cells_100_qvalue-1e-10_4.rds"))

num_cells = 100
qvalue = "1e-10"
cluster_n = 4
pattern = paste0(num_cells, "_qvalue_", qvalue, "_", cluster_n)

expectation_sf = readRDS(file.path(hdaci_data_dir, paste0("sf_heatmap_num_cells_", pattern, "_heatmap.rds")))
expectation_hash = readRDS(file.path(hdaci_data_dir, paste0("hash_heatmap_num_cells_", pattern, "_heatmap.rds")))

sf_cluster = readRDS(file.path(hdaci_data_dir, paste0("sf_heatmap_num_cells_", pattern, "_clusters.rds")))
hash_cluster = readRDS(file.path(hdaci_data_dir, paste0("hash_heatmap_num_cells_", pattern, "_clusters.rds")))

show_gsa = function(gsa_list, n) {
    for (i in 1:n) {
        go_df = data.frame(x = gsa_list[[i]]$p.adj) %>% 
                        rownames_to_column() %>% 
                        arrange(x) %>% 
                        filter(x < 0.05)
        display(head(go_df, 30))
    }
}

show_gsa(go_list_mult, 4)
show_gsa(hallmarks_list_hash, 4)

combine_table = function(gsa_list_1, gsa_list_2, n = 4) {
    out_df = list()
    gsa_list = list(gsa_list_1, gsa_list_2)
    genesets = c("GO", "HALLMARKS")
    names(gsa_list) = genesets
    
    for (set in genesets) {
        out_df[[set]] = list()
        for (i in 1:n) {
            out_df[[set]][[i]] = data.frame(x = gsa_list[[set]][[i]]$p.adj) %>% 
                            rownames_to_column() %>% 
                            arrange(x) %>% 
                            filter(x < 0.05)
        }
    }
    return(out_df)
}

out_df = combine_table(go_list_mult, hallmarks_list_hash, 4)

out_df$heatmap = expectation_hash
out_df$cluster = hash_cluster

combine_table = function(gsa_list_1, gsa_list_2, n = 4) {
    out_df = NULL
    gsa_list = list(gsa_list_1, gsa_list_2)
    genesets = c("GO", "HALLMARKS")
    names(gsa_list) = genesets
    
    for (set in genesets) {
        for (i in 1:n) {
            temp = data.frame(x = gsa_list[[set]][[i]]$p.adj) %>% 
                            rownames_to_column(var = "term") %>% 
                            arrange(x) %>% 
                            filter(x < 0.05)
            colnames(temp) = c("term", "q_value")
            temp$geneset = set
            temp$cluster = i
            out_df = rbind(out_df, temp)
        }
    }
    return(out_df)
}

out_df = combine_table(go_list_hash, hallmarks_list_hash, 4)

write.table(out_df, file = "../supp_tables/supplementary_table_3_hdaci_timecourse_gsea.txt", col.names = T, row.names = F, quote = F)
write_xlsx(out_df, "../supp_tables/supplementary_table_3_hdaci_timecourse_gsea.xlsx")


## Supplementary Table 4
dex_data_dir = "../data/hdaci_dex/de_analysis/dex_vehicle"

coeffs_table_dex <- readRDS(file.path(dex_data_dir, "coeff_table_sf.rds"))
degs_compared = readRDS(file.path(dex_data_dir, "degs_sf_LRT.rds"))

qvalue = 5e-2
num_cells = 20

coeffs_table_dex = coeffs_table_dex %>% 
                    filter(id %in% degs_compared$id) %>%
                    filter(grepl("Dex", term) & q_value < qvalue & num_cells_expressed > num_cells)
dim(coeffs_table_dex)

coeffs_table_dex_hash <- readRDS(file.path(dex_data_dir, "coeff_table_hash.rds"))
degs_compared = readRDS(file.path(dex_data_dir, "degs_hash_LRT.rds"))
coeffs_table_dex_hash = coeffs_table_dex_hash %>% 
                                filter(id %in% degs_compared$id) %>%
                                filter(grepl("Dex", term) & q_value < qvalue & num_cells_expressed > num_cells) 
dim(coeffs_table_dex_hash)

coeffs_table_dex$normalization = "Conventional"
coeffs_table_dex_hash$normalization = "Hash ladder"

degs_dex_vehicle = rbind(coeffs_table_dex, coeffs_table_dex_hash)
dim(degs_dex_vehicle)

write.table(degs_dex_vehicle, file = "../supp_tables/supplementary_table_4_dex_vehicle_degs.txt", col.names = T, row.names = F, quote = F)
write_xlsx(degs_dex_vehicle, "../supp_tables/supplementary_table_4_dex_vehicle_degs.xlsx")


## Supplementary Table 5
dex_data_dir = "../data/hdaci_dex"
dex_hdaci_degs = readRDS(file.path(dex_data_dir, "overall_estimate_df_by_time.rds"))

write.table(dex_hdaci_degs, file = "../supp_tables/supplementary_table_5_dex_hdaci_degs.txt", col.names = T, row.names = F, quote = F)
write_xlsx(dex_hdaci_degs, "../supp_tables/supplementary_table_5_dex_hdaci_degs.xlsx")

