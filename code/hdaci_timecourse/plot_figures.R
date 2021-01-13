suppressPackageStartupMessages({
    library(monocle3)
    library(MASS)
    library(reshape2)
    library(VGAM)
    library(tidyverse)
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
setwd("../../data/hdaci_timecourse")
cds_timecourse = readRDS("cds_timecourse_pseudotime.rds")
cds_timecourse_hash = readRDS("cds_timecourse_hash_pseudotime.rds")

cds_timecourse$Pseudotime = pseudotime(cds_timecourse)
cds_timecourse_hash$Pseudotime = pseudotime(cds_timecourse_hash)


## figure 3A
options(repr.plot.width=3, repr.plot.height=3)
g = list()
plot_cells(cds_timecourse, cell_size = 1, alpha = 0.5, 
           label_roots = F,
           trajectory_graph_color = "black", 
           trajectory_graph_segment_size = 1,
           show_trajectory_graph = T,
           color_cells_by = "Time",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5) + 
            theme(legend.position = "none")

ggsave("F3_umap_sf_graph.pdf", device = "pdf", width=3, height=3)

plot_cells(cds_timecourse_hash, cell_size = 1, alpha = 0.5, 
           label_roots = F,
           trajectory_graph_color = "black", 
           trajectory_graph_segment_size = 1,
           show_trajectory_graph = T,
           color_cells_by = "Time",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5) + 
            theme(legend.position = "none")

ggsave("F3_umap_mult_graph.pdf", device = "pdf", width=3, height=3)


## figure 3B
options(repr.plot.width=6, repr.plot.height=3)

g = list()

g[[1]] = ggplot(as.data.frame(colData(cds_timecourse)), aes(x = Pseudotime, y = log10(Total_RNA))) +
    geom_point(aes(colour = Time), alpha = 0.5) +  
    geom_smooth(method = "loess", span = 0.7, colour = "blue") +
    ggtitle("Conventional") + 
    ylim(c(3.2,4.1)) + 
    viridis::scale_color_viridis(option = "C") +
    monocle_theme_opts() + theme(legend.position = "none")

g[[2]] = ggplot(as.data.frame(colData(cds_timecourse_hash)), aes(x = Pseudotime, y = log10(Total_RNA))) +
    geom_point(aes(colour = Time), alpha = 0.5) + 
    geom_smooth(method = "loess", span = 0.7, colour = "blue") + 
    ggtitle("Hash ladder") + 
    ylim(c(3.2,4.1)) + 
    viridis::scale_color_viridis(option = "C") +
    monocle_theme_opts() + theme(legend.position = "none")


gg = do.call("grid.arrange", c(g, ncol=2))
ggsave(plot = gg, "F3_total_RNA_pseudotime.pdf", width = 6, height = 3)

options(repr.plot.width=4, repr.plot.height=3)
cars.lo <- loess(log(Total_RNA) ~ Pseudotime, span = 0.75, 
                 data = as.data.frame(colData(cds_timecourse_hash)))
x = seq(0, max(cds_timecourse_hash$Pseudotime), length.out = 1000)
y = predict(cars.lo, data.frame(Pseudotime = x), se = TRUE)

plot(x,y$fit)

## percentage of total RNA at minimum
exp(y$fit[which(y$fit == min(y$fit))]) / exp(y$fit[1])

cds_plot = cds_timecourse_hash
cds_plot$pseudotime = pseudotime(cds_plot)

options(repr.plot.width=4, repr.plot.height=3)
ggplot(as.data.frame(colData(cds_plot)), 
       aes(x = pseudotime, y = as.factor(Time))) + 
    geom_density_ridges2(aes(height = ..density.. , fill = Time)) + 
    geom_vline(xintercept = x[249], colour = "red", linetype = "longdash") +
    labs(x = "Pseudotime", y = "Time (hrs)") +
    viridis::scale_fill_viridis(option = "C") +
    theme_ridges() + theme(legend.position = "none", axis.text.x = element_blank())


## figure 3C
plot_heatmap <- function(cds_subset, m, cluster_rows = TRUE, bin_size = 100,  
                hclust_method = "ward.D2", min_expr = 0.2, ann_colors = "black",
    num_clusters = 5, hmcols = NULL, add_annotation_row = NULL, annotation_row = NA,
    add_annotation_col = NULL, show_rownames = FALSE, use_gene_short_name = TRUE, 
    norm_method = c("Size_Factor"), scale_max = 3, scale_min = -3, 
    main = NULL, return_heatmap = TRUE, file_name = NA) {
    
    if (use_gene_short_name) {
        rownames(m) = rowData(cds_subset)$gene_short_name
    }
    
    m = m[!apply(m, 1, sd) == 0, ]
    m = Matrix::t(scale(Matrix::t(m), center = TRUE))
    m = m[is.na(row.names(m)) == FALSE, ]
    m[is.nan(m)] = 0
    m[m > scale_max] = scale_max
    m[m < scale_min] = scale_min
    heatmap_matrix <- m
    row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix)))/2)
    row_dist[is.na(row_dist)] <- 1
    if (is.null(hmcols)) {
        bks <- seq(-3.1, 3.1, by = 0.1)
        hmcols <- viridis::inferno(length(bks) - 1)
    } else {
        bks <- seq(-3.1, 3.1, length.out = length(hmcols))
    }
    
    ph <- pheatmap(heatmap_matrix, useRaster = T, cluster_cols = F, 
        cluster_rows = cluster_rows, show_rownames = F, show_colnames = F, 
        clustering_distance_rows = row_dist, clustering_method = hclust_method, 
        cutree_rows = num_clusters, silent = TRUE, filename = NA, 
        breaks = bks, color = hmcols)

        annotation_row <- data.frame(Cluster = factor(cutree(ph$tree_row, 
        num_clusters)))
    
    if (!is.null(add_annotation_row)) {
        old_colnames_length <- ncol(annotation_row)
        annotation_row <- cbind(annotation_row, add_annotation_row[row.names(annotation_row), 
            ])
        colnames(annotation_row)[(old_colnames_length + 1):ncol(annotation_row)] <- colnames(add_annotation_row)
    }
    
    else {
        feature_label <- row.names(heatmap_matrix)
        row_ann_labels <- row.names(annotation_row)
    }
    
    row.names(heatmap_matrix) <- feature_label
    row.names(annotation_row) = feature_label
    colnames(heatmap_matrix) <- c(1:ncol(heatmap_matrix))
    
    ph_res <- pheatmap(heatmap_matrix[, ], useRaster = T, cluster_cols = FALSE, 
        cluster_rows = cluster_rows, show_rownames = show_rownames, 
        show_colnames = F, clustering_distance_rows = row_dist, 
        clustering_method = hclust_method, cutree_rows = num_clusters, 
        annotation_row = annotation_row, treeheight_row = 50, 
        breaks = bks, fontsize = 6, color = hmcols, silent = TRUE, 
        fontsize_row = 3, filename = file_name)
    grid::grid.rect(gp = grid::gpar("fill", col = NA))
    grid::grid.draw(ph_res$gtable)
    
    if (return_heatmap) {
        ph_res$annotation_row = annotation_row
        return(ph_res)
    }
    
}

num_cells = 100
qvalue = "1e-10"
cluster_n = 4
out_dir = "heatmaps"
pattern = paste0(num_cells, "_qvalue_", qvalue, "_", cluster_n)

expectation_sf = readRDS(file.path(out_dir, paste0("sf_heatmap_num_cells_", pattern, "_heatmap.rds")))
expectation_hash = readRDS(file.path(out_dir, paste0("hash_heatmap_num_cells_", pattern, "_heatmap.rds")))

sf_cluster = readRDS(file.path(out_dir, paste0("sf_heatmap_num_cells_", pattern, "_clusters.rds")))
hash_cluster = readRDS(file.path(out_dir, paste0("hash_heatmap_num_cells_", pattern, "_clusters.rds")))

options(repr.plot.width=5.5, repr.plot.height=5)
b = plot_heatmap(cds_timecourse_hash[rownames(expectation_hash),], expectation_hash,
                 num_clusters = 4, min_expr = 0.1, bin_size = 1000)


## figure 3D
gsa_list_hash = readRDS("heatmaps/gsa_list_go_hash_num_cells_100_qvalue-1e-10_4.rds")


convert_gene_to_id <- function(cds, genes) {
    rowData_mat = as.data.frame(rowData(cds))
    return(rowData_mat %>% filter(gene_short_name %in% genes) %>% pull(id))
}

convert_id_to_gene <- function(cds, ids) {
    rowData_mat = as.data.frame(rowData(cds))
    return(rowData_mat %>% filter(id %in% ids) %>% pull(gene_short_name))
}

go_term_names = names(gsa_list_hash[[1]]$gsc)
term_names = "KETONE|CHOLESTEROL|STEROL|STEROID|LIPID|ALCOHOL|CELLULAR|MRNA|FATTY|ACETYL-COA|GLYCO|GLUCOSE"
metabolic_terms = go_term_names[grepl("METABOL|CATABOL", go_term_names) & 
              grepl(term_names, go_term_names)]

metabolic_genes = unique(as.vector(unlist(sapply(metabolic_terms, function(x) {
    gsa_list_hash[[1]]$gsc[[x]]
}))))

metabolic_ids = convert_gene_to_id(cds_timecourse, metabolic_genes)


metabolic_cluster = hash_cluster %>% 
    rownames_to_column() %>%
    filter(rowname %in% metabolic_genes)


hash_cluster_genes = metabolic_cluster %>%
                        filter(Cluster %in% c(1,3)) %>%
                        pull(rowname)

hash_cluster_ids = convert_gene_to_id(cds_timecourse_hash, hash_cluster_genes)
hash_cluster_ids[hash_cluster_ids == "ENSG00000251562.7"] = "ENSG00000278217.1"

scale_min = -3
scale_max = 3

expectation_hash = expectation_hash[!apply(expectation_hash, 1, sd) == 0, ]
expectation_hash = Matrix::t(scale(Matrix::t(expectation_hash), center = TRUE))
expectation_hash = expectation_hash[is.na(row.names(expectation_hash)) == FALSE, ]
expectation_hash[is.nan(expectation_hash)] = 0
expectation_hash[expectation_hash > scale_max] = scale_max
expectation_hash[expectation_hash < scale_min] = scale_min

hash_cluster_mat = expectation_hash[hash_cluster_ids,]

find_closest_point = function(values, target) {
    diff_values = values - target
    x = which(abs(diff_values) == min(abs(diff_values)))
    return(data.frame(x = x, y = values[x]))
}

## point at which normalized expression = 0 
find_closest_index = function(values, target) {
    diff_values = values - target
    index_min = min(which(abs(diff_values) == min(abs(diff_values))))
    return(index_min)
}

find_crisis_point = function(values, target = 0, min_value = 0, min_index = 0) {
    if (min_index == 0) {
        local_min_index = which(diff(sign(diff(values)))==2)+1
        
        local_index = local_min_index
        min_value = min(values)
        
        if (length(local_index) == 0) {
            min_index = find_closest_index(values, min_value)
        } else {
            min_index = min(local_index)
            min_value = values[min_index]
        }
    }     
    
    x = find_closest_index(values, target = target)
   
    if (x > min_index & min_index < 0) {
       return(find_crisis_point(values[-x], target = target, min_value = min_value, min_index = min_index))
    } else if (x < min_index & abs(min_value) > 0.1) {
       return(find_crisis_point(values[-x], target = target, min_value = min_value, min_index = min_index)) 
    } else {
        return(x)
    }
}

omit_ids = c(38,146)
crisis_points = apply(hash_cluster_mat[-omit_ids,], 1, find_crisis_point, target = 0)

cds_plot = cds_timecourse_hash
cds_plot$pseudotime = pseudotime(cds_timecourse_hash)

df_crisis = stack(crisis_points / 1000 * max(cds_plot$pseudotime)) %>% 
    select(2,1) %>%
    rename(id = ind)

options(repr.plot.width=4, repr.plot.height=2)
ggplot(df_crisis, aes(x = values)) + 
    geom_histogram(colour = "black", fill = "white", bins = 40) + 
    scale_x_continuous(breaks=c(0,2,4,6,8,10,12)) + 
    labs(x = "", y = "") + 
    monocle_theme_opts()

ggsave("F3_crisis_hist_norm_expr_0.pdf", device = "pdf", width=4, height=2)


## figure 3E
cds_plot = cds_timecourse_hash
cds_plot$pseudotime = pseudotime(cds_plot)

options(repr.plot.width=5, repr.plot.height=3)
ggplot(as.data.frame(colData(cds_plot)), 
       aes(x = pseudotime, y = as.factor(Time))) + 
    geom_density_ridges2(aes(height = ..density.. , fill = Time)) + 
    geom_vline(xintercept = median(crisis_points)/1000 * max(cds_plot$pseudotime), colour = "red", linetype = "longdash") +
    facet_wrap("Drug") + 
    labs(x = "Pseudotime", y = "Time (hrs)") +
    viridis::scale_fill_viridis(option = "C") +
    theme_ridges() + theme(legend.position = "none", axis.text.x = element_blank())

ggsave("F3_crisis.pdf", device = "pdf", width=5, height=3)


## figure 3F
metabolic_genes = convert_gene_to_id(cds_plot, 
                                     c("ACSL3", "IDH1", "ACLY", "SLC2A3", "CDKN1A", "MKI67"))

options(repr.plot.width=4, repr.plot.height=5)
gg = plot_genes_in_pseudotime(cds_plot[metabolic_genes,],
                         color_cells_by="Time",
                         min_expr=0.1, ncol = 2) + 
        labs(x = "", y = "") + 
        #theme(axis.text = element_blank(), legend.position = "none", strip.text.x = element_blank()) + 
        geom_vline(xintercept = median(crisis_points)/1000 * max(cds_plot$pseudotime), colour = "red", linetype = "longdash")
        
gg
ggsave("F3_metabolic.png", gg, device = "png", width=4, height=5)




## figure S5B
hash_id_df_1 = read.table("../hdaci_qc/hash_id_df.txt_1", 
                          stringsAsFactors = F, header=T)
hash_id_df_2 = read.table("../hdaci_qc/hash_id_df.txt_2",
                          stringsAsFactors = F, header=T)

hash_id_df = rbind(hash_id_df_1, hash_id_df_2) %>% filter(!is.na(top_oligo))

## convert ID to condition
IDinfo <- read.table("../hdaci_qc/hashIDSampleSheet.txt", header=T)
indices <- sapply(hash_id_df$top_oligo, function(x) {
  which(x == IDinfo$ID)
})

temp = IDinfo[indices,-1]
hash_id_df <- cbind(hash_id_df, temp) %>%
                    filter(Dim == 1)
rownames(hash_id_df) <- NULL

g = list()
options(repr.plot.width=7, repr.plot.height=3)

g[[1]] = ggplot(hash_id_df, aes(x = log10(top_hash_umi))) + 
    geom_histogram() + 
    geom_vline(xintercept = log10(30), color = "red") + 
    monocle_theme_opts()

g[[2]] = ggplot(hash_id_df, aes(x = log10(top_to_second_best_ratio))) + 
    geom_histogram() + 
    geom_vline(xintercept = log10(10), color = "red") + 
    monocle_theme_opts()

gg = do.call("grid.arrange", c(g, ncol=2))
ggsave(gg, filename = "S_hashID_filter.pdf", width = 7, height = 3)


## figure S5C
hashTable1  = read.table("../hdaci_qc/hashTable_1.txt", header=T)
hashTable2  = read.table("../hdaci_qc/hashTable_2.txt", header=T)

hashTable = rbind(hashTable1, hashTable2) %>%
                filter(Dim == 1) %>%
                distinct(Cell, .keep_all = T) 

total_ladder1 = readRDS("../hdaci_qc/total_ladder_1.rds")
total_ladder2 = readRDS("../hdaci_qc/total_ladder_2.rds")
total_ladder = rbind(total_ladder1, total_ladder2) %>%
                filter(Dim == 1) %>%
                distinct(Cell, .keep_all = T) 

g = list()
options(repr.plot.width=7, repr.plot.height=3)
g[[1]] = ggplot(total_ladder, aes(x = log10(total_ladder))) + 
    geom_histogram(bins = 40) + 
    geom_vline(xintercept = log10(100), color = "red") + 
    monocle_theme_opts()

g[[2]] = ggplot(hashTable, aes(x = Rsq)) + 
    geom_histogram(bins = 100) + 
    geom_vline(xintercept = 0.7, color = "red") + 
    monocle_theme_opts()

gg = do.call("grid.arrange", c(g, ncol=2))
ggsave(gg, filename = "S_hashLadder_filter.pdf", width = 7, height = 3)


## figure S5D
metadata = as.data.frame(colData(cds_timecourse))
melted_metadata = metadata %>% select(Cell, Total_hash, Total_RNA) %>% melt()

options(repr.plot.width=3.5, repr.plot.height=3)
ggplot(melted_metadata, aes(x = variable, y = log10(value))) + 
    geom_boxplot(aes(fill = variable)) +
    monocle_theme_opts()
ggsave(filename = "S_hash_RNA_count.pdf", width = 3.5, height = 3)


## figure S6A
degs_dir = "de_analysis"
out_dir = "de_analysis"

qvalue = 1e-2
num_cells = 100

degs_compared = readRDS(file.path(degs_dir, "degs_sf_LRT.rds"))
degs_sf = readRDS(file.path(degs_dir, "coeff_table_sf.rds"))

degs_compared = degs_compared %>% filter(q_value < qvalue)
degs_sf = degs_sf %>%
    filter(id %in% degs_compared$id & grepl("Pseudo", term)) %>%
    filter(num_cells_expressed > num_cells & q_value < qvalue) %>%
    distinct(id, .keep_all = TRUE)

degs_compared = readRDS(file.path(degs_dir, "degs_hash_LRT.rds"))
degs_hash = readRDS(file.path(degs_dir, "coeff_table_hash.rds"))

degs_compared = degs_compared %>% filter(q_value < qvalue)
degs_hash = degs_hash %>%  
    filter(id %in% degs_compared$gene_id & grepl("Pseudo", term)) %>%
    filter(num_cells_expressed > num_cells & q_value < qvalue) %>%
    distinct(id, .keep_all = TRUE)

n_sf = degs_sf %>% 
    filter(!id %in% degs_hash$id) %>% nrow()
 
n_hash = degs_hash %>%
    filter(!id %in% degs_sf$id) %>% nrow()

n_overlap = degs_sf %>% 
    filter(id %in% degs_hash$id) %>% nrow()

print(n_sf)
print(n_hash)
print(n_overlap)

options(repr.plot.width=3, repr.plot.height=3)

g = euler(c("Conventional" = n_sf, "Hash ladder" = n_hash, "Conventional&Hash ladder" = n_overlap))

gg = plot(g, 
     family = "Helveltica",
     fills = list(fill = c("yellow", "darkgreen"), alpha = 0.7),
     quantities = F, 
     labels = NULL)

pdf("S_timecourse_venn.pdf", width = 3, height = 3)
gg
dev.off()
gg


## figure S6B
sci_plex_degs = read.delim("aax6234-Srivatsan-Table-S8.txt", header = T, sep = "\t", stringsAsFactors = F) %>%
                    arrange(gene_short_name)

a549_genes = sci_plex_degs %>%
    filter(grepl("A549:", term), q_value < 0.01) %>%
    distinct(id, .keep_all=T)


n_degs_hash = nrow(degs_hash)
n_degs_sf = nrow(degs_sf)

dose_degs_hash = degs_hash %>%
                    filter(id %in% a549_genes$id) %>% 
                    nrow()

dose_degs_sf = degs_sf %>%
                    filter(id %in% a549_genes$id) %>% 
                    nrow()


deg_overlap = rbind(
    data.frame(n = c(dose_degs_hash, n_degs_hash - dose_degs_hash), norm = "Hash ladder", type = c("Yes", "No")),
    data.frame(n = c(dose_degs_sf, n_degs_sf - dose_degs_sf), norm = "Conventional", type = c("Yes", "No")))

deg_overlap$type = factor(deg_overlap$type, levels = c("Yes", "No"))

options(repr.plot.width=4.5, repr.plot.height=3)
ggplot(deg_overlap, aes(x = norm, y = n)) + 
    geom_bar(stat = "identity", position = "dodge", aes(fill = type)) + 
    monocle_theme_opts()
ggsave("S_sci-plex-overlap.pdf", width = 4.5, height = 3)


## figure S6C
model_predictions <- function(model_tbl, new_data, type="response") {
  predict_helper <- function(model, cds){
    tryCatch({
      stats::predict(model, newdata=new_data, type=type)
    }, error = function(e){
      retval = rep_len(NA, nrow(new_data))
      names(retval) = row.names(new_data)
      return(retval)
    })
  }
  model_tbl <- model_tbl %>%
    dplyr::mutate(predictions = purrr::map(model, predict_helper, new_data))
  pred_matrix = t(do.call(cbind, model_tbl$predictions))
  return(pred_matrix)
}

plot_pseudotime_heatmap <- function(cds_subset, cluster_rows = TRUE, bin_size = 100,  
                hclust_method = "ward.D2", min_expr = 0.2, ann_colors = "black",
    num_clusters = 3, hmcols = NULL, add_annotation_row = NULL, annotation_row = NA,
    add_annotation_col = NULL, show_rownames = FALSE, use_gene_short_name = TRUE, 
    norm_method = c("Size_Factor"), scale_max = 3, scale_min = -3, 
    trend_formula = "~ splines::ns(Pseudotime, df=3)", main = "", return_heatmap = FALSE, 
    cores = 1, file_name = NA) 
{
    a = data.frame(table(rowData(cds_subset)$gene_short_name))
    if (nrow(a[a$Freq>1,]) > 0) {
        unique_df = rowData(cds_subset)[rowData(cds_subset)$gene_short_name %in% a[a$Freq>1,]$Var1,]
        unique_id = as.data.frame(unique_df) %>% 
                        dplyr::group_by(gene_short_name) %>% 
                        filter(num_cells_expressed == max(num_cells_expressed)) %>%
                        pull(id)
    } else {
        unique_id = ""
    }
    cds_subset = cds_subset[!(rowData(cds_subset)$id %in%unique_id), ]
    cds_subset$Pseudotime = pseudotime(cds_subset)

    newdata <- data.frame(Pseudotime = seq(min(pData(cds_subset)$Pseudotime), 
            max(pData(cds_subset)$Pseudotime), length.out = bin_size))

    cds_exprs <- SingleCellExperiment::counts(cds_subset)
    cds_exprs <- Matrix::t(Matrix::t(cds_exprs)/size_factors(cds_subset))

    
    new_data <- data.frame(Pseudotime = cds_subset$Pseudotime)
    model_tbl = fit_models(cds_subset, model_formula_str = trend_formula,                    
                           expression_family = "negbinomial",
                           reduction_method="UMAP",
                           cores = cores,
                           clean_model = TRUE,
                           verbose = FALSE)

    m <- model_predictions(model_tbl, new_data = newdata)
    
    #return(m)
    
    if (use_gene_short_name) {
        rownames(m) = rowData(cds_subset)$gene_short_name
    }
    
    m = m[!apply(m, 1, sd) == 0, ]
    m = Matrix::t(scale(Matrix::t(m), center = TRUE))
    m = m[is.na(row.names(m)) == FALSE, ]
    m[is.nan(m)] = 0
    m[m > scale_max] = scale_max
    m[m < scale_min] = scale_min
    heatmap_matrix <- m
    row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix)))/2)
    row_dist[is.na(row_dist)] <- 1
    if (is.null(hmcols)) {
        bks <- seq(-3.1, 3.1, by = 0.1)
        hmcols <- viridis::inferno(length(bks) - 1)
    } else {
        bks <- seq(-3.1, 3.1, length.out = length(hmcols))
    }
    
    hmcols <- viridis::inferno(length(bks) - 1)
    
    ph <- pheatmap(heatmap_matrix, useRaster = T, cluster_cols = F, 
        cluster_rows = cluster_rows, show_rownames = show_rownames, show_colnames = F, 
        silent = TRUE, filename = NA, 
        breaks = bks, color = hmcols)
    
    
    if (!is.null(add_annotation_row)) {
        old_colnames_length <- ncol(annotation_row)
        annotation_row <- cbind(annotation_row, add_annotation_row[row.names(annotation_row), 
            ])
        colnames(annotation_row)[(old_colnames_length + 1):ncol(annotation_row)] <- colnames(add_annotation_row)
    }
    
    else {
        feature_label <- rownames(heatmap_matrix)
        row_ann_labels <- rownames(annotation_row)
    }
    
    
    rownames(heatmap_matrix) <- feature_label
    colnames(heatmap_matrix) <- c(1:ncol(heatmap_matrix))
    
    ph_res <- pheatmap(heatmap_matrix[, ], useRaster = T, cluster_cols = FALSE, 
        cluster_rows = F, show_rownames = show_rownames, 
        show_colnames = F, silent=F, border_color = NA,
        treeheight_row = 0, main = main,
        breaks = bks, fontsize = 10, color = hmcols,
        annotation_names_row = F, annotation_names_col = F,
        fontsize_row = 12, filename = file_name, legend=F)
    grid::grid.rect(gp = grid::gpar("fill", col = NA))
    grid::grid.draw(ph_res$gtable)
    
    if (return_heatmap) {
        ph_res$annotation_row = annotation_row
        return(ph_res)
    }
}

acetyl_coa_genes =
  c("SLC13A3", 
    "SLC2A3", 
    "GLS",
    "IDH1", 
    "ACSS2", 
    "ACSL3", 
    "SIRT2", 
    "SLC25A16", 
    "SLC25A1", 
    "ACLY", 
    "HADHB",
    "HADHA",
    "ACO1",
    "CS")

sf_acetyl_coa_genes = degs_sf %>% filter(gene_short_name %in% acetyl_coa_genes) %>% 
                            pull(id)

hash_acetyl_coa_genes = degs_hash %>% filter(gene_short_name %in% acetyl_coa_genes) %>% 
                            pull(id)

g = list()
options(repr.plot.width=2, repr.plot.height=4)
a =  plot_pseudotime_heatmap(cds_subset = cds_timecourse_hash[hash_acetyl_coa_genes,], num_clusters = 1, 
                                      cluster_rows = F,
                                      return_heatmap = T, bin_size = 500,
                                      show_rownames = T, file_name = NA)
g[[1]] = a[[4]]

options(repr.plot.width=4, repr.plot.height=4)
a = plot_pseudotime_heatmap(cds_subset = cds_timecourse[sf_acetyl_coa_genes,], num_clusters = 1, 
                                      cluster_rows = F,
                                      return_heatmap = T, bin_size = 500,
                                      show_rownames = T, file_name = NA)
g[[2]] = a[[4]]


options(repr.plot.width=6, repr.plot.height=3)
gg = grid.arrange(arrangeGrob(grobs= g,ncol=2))
ggsave(plot = gg, "S_timecourse_acetylcoa_genes.png", width = 6, height = 3)


## figure S6D
slc2a3_id = convert_gene_to_id(cds_timecourse_hash, 
                                     c("SLC2A3"))

options(repr.plot.width=2, repr.plot.height=2)
gg = plot_genes_in_pseudotime_2(cds_timecourse_hash[slc2a3_id,], 
                         color_cells_by="Time",
                         min_expr=0.1, ncol = 1) + 
        labs(x = "", y = "") + 
        #scale_y_log10(limits = c(1,30)) + 
        theme(axis.text = element_blank(), legend.position = "none", strip.text.x = element_blank())
        
gg
ggsave("S_SLC2A3.png", gg, device = "png", width=2, height=2)

slc2a3_expr = data.frame(x = seq(1,1000,by=1), y = hash_cluster_mat[slc2a3_id,])
slc2a3_50_point = find_crisis_point(slc2a3_expr$y, fc = 0.5)
options(repr.plot.width=2.5, repr.plot.height=2)

ggplot(slc2a3_expr, aes(x/1000 * max(pseudotime(cds_timecourse_hash)), y)) + 
    geom_point() +
    geom_point(data = slc2a3_expr[slc2a3_50_point,], aes(x/1000 * max(pseudotime(cds_timecourse_hash)),y), color = "red") + 
    scale_x_continuous() + 
    labs(x = "Pseudotime", y = "Expression") + 
    monocle_theme_opts()
ggsave("S_SLC2A3_pseudotime_fit.pdf", device = "pdf", width=2.5, height=2)


## figure S7A
cds_timecourse_hash = cluster_cells(cds_timecourse_hash, resolution = 0.01, random_seed = 123)

options(repr.plot.width=3, repr.plot.height=3)
plot_cells(cds_timecourse_hash, cell_size = 1, alpha = 0.6, 
           label_groups_by_cluster=F,
           show_trajectory_graph = F,
           color_cells_by = "cluster",  
           label_leaves=FALSE,
           label_branch_points=FALSE,
           group_label_size=4)

ggsave("S_cell_cycle_cluster_umap.pdf", device = "pdf", width=3, height=3)

## figure S7B
cds_timecourse_hash = detect_genes(cds_timecourse_hash, 1)
cds_timecourse_hash = cds_timecourse_hash[rowData(cds_timecourse_hash)$num_cells_expressed > 50,]

cds_hash_graph <- graph_test(cds_timecourse_hash, neighbor_graph="principal_graph", cores=1, verbose = T)

graph_deg_ids <- row.names(subset(cds_hash_graph, q_value < 1e-10))
gene_module_df <- find_gene_modules(cds_timecourse_hash[graph_deg_ids,], resolution=1e-2)

options(repr.plot.width=5, repr.plot.height=4)
cell_group_df <- tibble::tibble(cell=row.names(colData(cds_timecourse_hash)), 
                                cell_group=as.factor(clusters(cds_timecourse_hash)))
agg_mat <- aggregate_gene_expression(cds_timecourse_hash, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
pheatmap::pheatmap(agg_mat, width = 5, height = 4,
                   scale="column", clustering_method="ward.D2")
pheatmap::pheatmap(agg_mat, filename = "S_gene_module.pdf", width = 5, height = 4,
                   scale="column", clustering_method="ward.D2")


## figure S7C
options(repr.plot.width=4.2, repr.plot.height=3)
plot_cells(cds_timecourse_hash, cell_size = 1,
           genes = gene_module_df %>% filter(module == 6),
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE)

ggsave("S_module_6_umap.pdf", device = "pdf", width=4.2, height=3)


## figure S7D
loadGSCSafe <- function (file, type = "auto", addInfo, sep="\t", encoding="latin1") 
{
  if (missing(addInfo)) {
    addUserInfo <- "skip"
    addInfo <- "none"
  }
  else {
    addUserInfo <- "yes"
  }
  tmp <- try(type <- match.arg(type, c("auto", "gmt", "sbml", 
                                       "sif", "data.frame"), several.ok = FALSE), silent = TRUE)
  if (class(tmp) == "try-error") {
    stop("argument type set to unknown value")
  }
  if (type == "auto") {
    if (class(file) == "character") {
      tmp <- unlist(strsplit(file, "\\."))
      type <- tolower(tmp[length(tmp)])
      if (!type %in% c("gmt", "sif", "sbml", "xml")) 
        stop(paste("can not handle .", type, " file extension, read manually using e.g. read.delim() and load as data.frame", 
                   sep = ""))
    }
    else {
      type <- "data.frame"
    }
  }
  if (type == "gmt") {
    con <- file(file, encoding=encoding)
    tmp <- try(suppressWarnings(open(con)), silent = TRUE)
    if (class(tmp) == "try-error") 
      stop("file could not be read")
    if (addUserInfo == "skip") 
      addInfo <- vector()
    gscList <- list()
    i <- 1
    tmp <- try(suppressWarnings(while (length(l <- scan(con, 
                                                        nlines = 1, what = "character", quiet = T, sep=sep)) > 0) {
      if (addUserInfo == "skip") 
        addInfo <- rbind(addInfo, l[1:2])
      tmp <- l[3:length(l)]
      gscList[[l[1]]] <- unique(tmp[tmp != "" & tmp != 
                                      " " & !is.na(tmp)])
      i <- i + 1
    }), silent = TRUE)
    if (class(tmp) == "try-error") 
      stop("file could not be read")
    close(con)
    gsc <- gscList[!duplicated(names(gscList))]
    if (addUserInfo == "skip") 
      addInfo <- unique(addInfo)
  }
  else if (type %in% c("sbml", "xml")) {
    require(rsbml)
    tmp <- try(sbml <- rsbml_read(file))
    if (class(tmp) == "try-error") {
      stop("file could not be read by rsbml_read()")
    }
    gsc <- list()
    for (iReaction in 1:length(reactions(model(sbml)))) {
      metIDs <- names(c(reactants(reactions(model(sbml))[[iReaction]]), 
                        products(reactions(model(sbml))[[iReaction]])))
      geneIDs <- names(modifiers(reactions(model(sbml))[[iReaction]]))
      if (length(geneIDs) > 0) {
        geneNames <- rep(NA, length(geneIDs))
        for (iGene in 1:length(geneIDs)) {
          geneNames[iGene] <- name(species(model(sbml))[[geneIDs[iGene]]])
        }
        for (iMet in 1:length(metIDs)) {
          gsc[[metIDs[iMet]]] <- c(gsc[[metIDs[iMet]]], 
                                   geneNames)
        }
      }
    }
    if (length(gsc) == 0) {
      stop("no gene association found")
    }
    else {
      for (iMet in 1:length(gsc)) {
        tmp1 <- name(species(model(sbml))[[names(gsc)[iMet]]])
        tmp2 <- compartment(species(model(sbml))[[names(gsc)[iMet]]])
        names(gsc)[iMet] <- paste(tmp1, " (", tmp2, ")", 
                                  sep = "")
      }
    }
  }
  else if (type == "sif") {
    tmp <- try(gsc <- as.data.frame(read.delim(file, header = FALSE, 
                                               quote = "", as.is = TRUE), stringsAsFactors = FALSE), 
               silent = TRUE)
    if (class(tmp) == "try-error") {
      stop("argument file could not be read and converted into a data.frame")
    }
    if (ncol(gsc) != 3) {
      stop("sif file should contain three columns")
    }
    if (addUserInfo == "skip") 
      addInfo <- gsc[, c(1, 2)]
    gsc <- gsc[, c(3, 1)]
    tmp <- nrow(gsc)
    gsc <- unique(gsc)
    geneSets <- unique(gsc[, 2])
    gscList <- list()
    for (iGeneSet in 1:length(geneSets)) {
      gscList[[iGeneSet]] <- gsc[gsc[, 2] == geneSets[iGeneSet], 
                                 1]
    }
    names(gscList) <- geneSets
    gsc <- gscList
  }
  else if (type == "data.frame") {
    tmp <- try(gsc <- as.data.frame(file, stringsAsFactors = FALSE), 
               silent = TRUE)
    if (class(tmp) == "try-error") {
      stop("argument file could not be converted into a data.frame")
    }
    for (i in 1:ncol(gsc)) {
      gsc[, i] <- as.character(gsc[, i])
    }
    if (ncol(gsc) != 2) {
      stop("argument file has to contain exactly two columns")
    }
    tmp <- nrow(gsc)
    gsc <- unique(gsc)
    geneSets <- unique(gsc[, 2])
    gscList <- list()
    for (iGeneSet in 1:length(geneSets)) {
      gscList[[iGeneSet]] <- gsc[gsc[, 2] == geneSets[iGeneSet], 
                                 1]
    }
    names(gscList) <- geneSets
    gsc <- gscList
  }
  if (addUserInfo == "yes") {
    tmp <- try(addInfo <- as.data.frame(addInfo, stringsAsFactors = FALSE), 
               silent = TRUE)
    if (class(tmp) == "try-error") {
      stop("failed to convert additional info in argument 'addInfo' into a data.frame")
    }
  }
  if (class(addInfo) == "data.frame") {
    if (ncol(addInfo) != 2) 
      stop("additional info in argument 'file' or 'addInfo' has to contain 2 columns")
    tmp <- nrow(addInfo)
    addInfo <- unique(addInfo[addInfo[, 1] %in% names(gsc), 
                              ])
  }
  else {
  }
  res <- list(gsc, addInfo)
  names(res) <- c("gsc", "addInfo")
  class(res) <- "GSC"
  return(res)
}

gmt_dir = "../gmts"
## Load Gene Set Collections
reactomeGSC<-loadGSCSafe(file=file.path(gmt_dir, "Human_Reactome_August_01_2019_symbol.gmt"))
pantherGSC<-loadGSCSafe(file=file.path(gmt_dir, "Human_Panther_August_01_2019_symbol.gmt"))
KEGGGSC<-loadGSCSafe(file=file.path(gmt_dir, "Human_KEGG_August_01_2019_symbol.gmt"))
NCIGSC<-loadGSCSafe(file=file.path(gmt_dir, "Human_NCI_Nature_August_01_2019_symbol.gmt"))
GOGSC<-loadGSCSafe(file=file.path(gmt_dir, "Human_GO_bp_no_GO_iea_symbol.gmt"))
TFGSC<-loadGSCSafe(file=file.path(gmt_dir, "Human_MSigdb_August_01_2019_symbol.gmt"))
HALLMARKS<-loadGSCSafe(file=file.path(gmt_dir, "h.all.v6.2.symbols.gmt"))


gsa_hyper_clusters <- function(universe, genes, gsc){
    
    GSA_list = 
        tryCatch({
            gsares <- runGSAhyper(genes = genes, universe = universe, gsc=gsc)
            return (gsares)
        }, error = function(e) { print(e); return (NA) } )
    
    return(GSA_list)
}

gsa_list_hash = suppressWarnings({gsa_hyper_clusters(unique(rowData(cds_timecourse_hash)$gene_short_name), 
                                                     gene_module_df %>% 
                                                        filter(module == 6) %>%
                                                        left_join(., as.data.frame(rowData(cds_timecourse_hash))) %>%
                                                        pull(gene_short_name), GOGSC)})

g = gsa_list_hash
go_df = data.frame(x = g$p.adj) %>% 
                    rownames_to_column() %>% 
                    arrange(x) %>% 
                    filter(x < 0.05) %>%
                    mutate(x = -log10(x)) %>%
                    slice(1:6)
options(repr.plot.width=6, repr.plot.height=3)
ggplot(go_df, aes(x = x, y = reorder(rowname, x))) + 
    geom_bar(stat = "identity") + 
    labs(x = "", y = "") + 
    monocle_theme_opts()

ggsave("S_cluster_GSEA.pdf", device = "pdf", width=6, height=3)


## figure S7E
options(repr.plot.width=9, repr.plot.height=3)
plot_cells(cds_timecourse_hash, genes=c("KIF20B", "TUBB4B", "MKI67"), 
           cell_size = 1.2, alpha = 0.8,
           show_trajectory_graph=FALSE, 
           label_cell_groups=FALSE)

ggsave("S_cluster_cell_cycle_genes.pdf", device = "pdf", width=9, height=3)


