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

# Downstream analysis of HDACi timecourse data

## Load cds and de analysis data
setwd("../../data/hdaci_timecourse")
cds_timecourse = readRDS("cds_timecourse_pseudotime.rds")
cds_timecourse_hash = readRDS("cds_timecourse_hash_pseudotime.rds")

cds_timecourse$Pseudotime = pseudotime(cds_timecourse)
cds_timecourse_hash$Pseudotime = pseudotime(cds_timecourse_hash)

degs_dir = "de_analysis"
out_dir = "de_analysis"

qvalue = 1e-2
num_cells = 50

degs_compared = readRDS(file.path(degs_dir, "degs_sf_LRT.rds"))
degs_sf = readRDS(file.path(degs_dir, "coeff_table_sf.rds"))

degs_compared = degs_compared %>% filter(q_value < qvalue)
degs_sf = degs_sf %>%
    filter(id %in% degs_compared$id & grepl("Pseudo", term)) %>%
    filter(num_cells_expressed > num_cells & q_value < qvalue) %>%
    distinct(id, .keep_all = TRUE)
dim(degs_sf)

degs_compared = readRDS(file.path(degs_dir, "degs_hash_LRT.rds"))
degs_hash = readRDS(file.path(degs_dir, "coeff_table_hash.rds"))

degs_compared = degs_compared %>% filter(q_value < qvalue)
degs_hash = degs_hash %>%  
    filter(id %in% degs_compared$gene_id & grepl("Pseudo", term)) %>%
    filter(num_cells_expressed > num_cells & q_value < qvalue) %>%
    distinct(id, .keep_all = TRUE)
dim(degs_hash)


## Hierarchical clustering of DE genes
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
                hclust_method = "ward.D2", min_expr = 0.1, ann_colors = "black",
    num_clusters = 3, hmcols = NULL, add_annotation_row = NULL, annotation_row = NA,
    add_annotation_col = NULL, show_rownames = FALSE, use_gene_short_name = TRUE, 
    norm_method = c("Size_Factor"), scale_max = 3, scale_min = -3, mat_name = NA,
    trend_formula = "~ splines::ns(Pseudotime, df=3)", main = NULL, return_heatmap = FALSE, 
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

    newdata <- data.frame(Plate = "1", Pseudotime = seq(min(pData(cds_subset)$Pseudotime), 
            max(pData(cds_subset)$Pseudotime), length.out = bin_size))
    #print(head(as.data.frame(colData(cds_subset))))

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
    
    if (!is.na(mat_name)) {
        saveRDS(m, mat_name)
    }
    
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


options(repr.plot.width=6, repr.plot.height=8)

num_cells = 100
qval = 1e-10
cluster_n = 5
out_dir = "heatmaps"

expectation = plot_pseudotime_heatmap(cds_timecourse[degs_sf$id,], num_clusters = cluster_n, 
                  return_heatmap = T, bin_size = 1000, cores=1,
                  show_rownames = F, 
                  mat_name = file.path(out_dir, paste0("sf_heatmap_num_cells_", num_cells, "_qvalue_", qval, "_", cluster_n, "_heatmap.rds")),
                  file_name = file.path(out_dir, paste0("sf_heatmap_num_cells_", num_cells, "_qvalue_", qval, "_", cluster_n, ".png")))

expectation_hash = plot_pseudotime_heatmap(cds_timecourse_hash[degs_hash$id,], num_clusters = cluster_n, 
                  return_heatmap = T, bin_size = 1000, cores=1,
                  show_rownames = F, 
                  mat_name = file.path(out_dir, paste0("hash_heatmap_num_cells_", num_cells, "_qvalue_", qval, "_", cluster_n, "_heatmap.rds")),
                  file_name = file.path(out_dir, paste0("hash_heatmap_num_cells_", num_cells, "_qvalue_", qval, "_", cluster_n, ".png")))
                                           
saveRDS(expectation$annotation_row, file.path(out_dir, paste0("sf_heatmap_num_cells_", num_cells, "_qvalue_", qval, "_", cluster_n, "_clusters.rds")))
saveRDS(expectation_hash$annotation_row, file.path(out_dir, paste0("hash_heatmap_num_cells_", num_cells, "_qvalue_", qval, "_", cluster_n, "_clusters.rds")))



## GSEA analysis
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

gsa_hyper_clusters <- function(gene_test_res, gene_clusters, gsc){
    universe = unique(gene_test_res$gene_short_name)
    clusters = as.numeric(unique(unlist(gene_clusters)))
    
    GSA_list = lapply(clusters, function(cluster) {
        print(cluster)
        genes = row.names(gene_clusters[gene_clusters == cluster, , drop=F])
        tryCatch({
            gsares <- runGSAhyper(genes = genes, universe = universe, gsc=gsc)
            return (gsares)
        }, error = function(e) { print(e); return (NA) } )
    })
    
    names(GSA_list) = clusters
    
    return(GSA_list)
}

## Load Gene Set Collections
gmt_dir = "../gmts"
reactomeGSC<-loadGSCSafe(file=file.path(gmt_dir, "Human_Reactome_August_01_2019_symbol.gmt"))
pantherGSC<-loadGSCSafe(file=file.path(gmt_dir, "Human_Panther_August_01_2019_symbol.gmt"))
KEGGGSC<-loadGSCSafe(file=file.path(gmt_dir, "Human_KEGG_August_01_2019_symbol.gmt"))
NCIGSC<-loadGSCSafe(file=file.path(gmt_dir, "Human_NCI_Nature_August_01_2019_symbol.gmt"))
GOGSC<-loadGSCSafe(file=file.path(gmt_dir, "Human_GO_bp_no_GO_iea_symbol.gmt"))
TFGSC<-loadGSCSafe(file=file.path(gmt_dir, "Human_MSigdb_August_01_2019_symbol.gmt"))
HALLMARKS<-loadGSCSafe(file=file.path(gmt_dir, "h.all.v6.2.symbols.gmt"))

## Load clustered DE genes
num_cells = 100
qvalue = "1e-10"
cluster_n = 4
out_dir = "heatmaps"
pattern = paste0(num_cells, "_qvalue_", qvalue, "_", cluster_n)

expectation_sf = readRDS(file.path(out_dir, paste0("sf_heatmap_num_cells_", pattern, "_heatmap.rds")))
expectation_hash = readRDS(file.path(out_dir, paste0("hash_heatmap_num_cells_", pattern, "_heatmap.rds")))

sf_cluster = readRDS(file.path(out_dir, paste0("sf_heatmap_num_cells_", pattern, "_clusters.rds")))
hash_cluster = readRDS(file.path(out_dir, paste0("hash_heatmap_num_cells_", pattern, "_clusters.rds")))

gene_univ = as.data.frame(rowData(cds_timecourse_hash)) %>% 
                filter(num_cells_expressed > 50)


## Perform GSEA analysis
GSC = HALLMARKS

gsa_list_sf = suppressWarnings({gsa_hyper_clusters(gene_univ, sf_cluster, GSC)})
gsa_list_hash = suppressWarnings({gsa_hyper_clusters(gene_univ, hash_cluster, GSC)})

saveRDS(gsa_list_hash, "heatmaps/gsa_list_hallmarks_hash_num_cells_100_qvalue-1e-10_4.rds")
saveRDS(gsa_list_sf, "heatmaps/gsa_list_hallmarks_sf_num_cells_100_qvalue-1e-10_4.rds")

gsa_list_hash = readRDS("heatmaps/gsa_list_go_hash_num_cells_100_qvalue-1e-10_4.rds")

show_gsa = function(gsa_list, num_clusters) {
    for (i in 1:num_clusters) {
        go_df = data.frame(x = gsa_list[[i]]$p.adj) %>% 
                            rownames_to_column() %>% 
                            arrange(x) %>% 
                            filter(x < 0.05)
        display(i)
        display(head(go_df, 30))
    }
}

show_gsa(gsa_list_hash, 4)
show_gsa(gsa_list_sf, 4)

show_gsa = function(gsa_list, num_clusters) {
    for (i in 1:num_clusters) {
        go_df = data.frame(x = gsa_list[[i]]$p.adj) %>% 
                            rownames_to_column() %>% 
                            arrange(x) %>% 
                            filter(x < 0.05)
        display(i)
        display(head(go_df, 30))
    }
}

show_gsa(gsa_list_hash, 4)
show_gsa(gsa_list_sf, 4)

