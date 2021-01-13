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
setwd("../../data/hdaci_dex/")
output_dir = "de_analysis/dex_time" 
sf_file_names = c(file.path(output_dir, "coeff_table_sf_4hr_dex.rds"),
                  file.path(output_dir, "coeff_table_sf_24hr_dex.rds"))

hash_file_names = c(file.path(output_dir, "coeff_table_hash_4hr_dex.rds"),
                  file.path(output_dir, "coeff_table_hash_24hr_dex.rds"))

sf_compared_names = c(file.path(output_dir, "degs_sf_LRT_4hr_dex.rds"),
                      file.path(output_dir, "degs_sf_LRT_24hr_dex.rds"))

hash_compared_names = c(file.path(output_dir, "degs_hash_LRT_4hr_dex.rds"),
                      file.path(output_dir, "degs_hash_LRT_24hr_dex.rds"))


i = 1
qvalue = 5e-2
num_cells = 20
term_name = "Dex"

coeffs_table_dex <- readRDS(sf_file_names[i])
degs_compared = readRDS(sf_compared_names[i])
coeffs_table_dex_1 = coeffs_table_dex %>% 
                    filter(id %in% degs_compared$id) %>%
                    filter(grepl(term_name, term) & q_value < qvalue & num_cells_expressed > num_cells)
dim(coeffs_table_dex_1)

coeffs_table_dex_hash <- readRDS(hash_file_names[i])
degs_compared = readRDS(hash_compared_names[i])
coeffs_table_dex_hash_1 = coeffs_table_dex_hash %>% 
                                filter(id %in% degs_compared$id) %>%
                                filter(grepl(term_name, term) & q_value < qvalue & num_cells_expressed > num_cells) 
dim(coeffs_table_dex_hash_1)

i = 2
coeffs_table_dex <- readRDS(sf_file_names[i])
degs_compared = readRDS(sf_compared_names[i])
coeffs_table_dex_10 = coeffs_table_dex %>% 
                    filter(id %in% degs_compared$id) %>%
                    filter(grepl(term_name, term) & q_value < qvalue & num_cells_expressed > num_cells)
dim(coeffs_table_dex_10)

coeffs_table_dex_hash <- readRDS(hash_file_names[i])
degs_compared = readRDS(hash_compared_names[i])
coeffs_table_dex_hash_10 = coeffs_table_dex_hash %>% 
                                filter(id %in% degs_compared$id) %>%
                                filter(grepl(term_name, term) & q_value < qvalue & num_cells_expressed > num_cells) 
dim(coeffs_table_dex_hash_10)


## make annotations of DEX genes
annotate_helper_1 = function(sub_table) {
    single_annotation_names = c("HDACi")
    grep_types = c("DexTRUE")
    
    if (nrow(sub_table) == 1) {
        out_df = sub_table %>%
                    mutate(annotation = single_annotation_names[which(sapply(grep_types, function(x) grepl(x, term)))])
    } else {
        out_df = sub_table %>%
                    mutate(annotation = single_annotation_names[apply(sapply(grep_types, 
                                                                         function(x) grepl(x, term)), 1, 
                                                                         function(x) which(x))])
    }
                                                                             
    return(out_df)
}


convert_to_df = function(coeffs_table, time, normalization, unique = F) {
    
    grep_types = c("DexTRUE")
    types = c("Dex")
    
    coeffs_table = coeffs_table %>%
                        select(id, gene_short_name, num_cells_expressed, term, estimate, q_value) %>%
                        mutate(time = time, normalization = normalization)  
                                   
    genes = unique(coeffs_table$gene_short_name)
                                                         
    annotation = plyr::rbind.fill(lapply(genes, function(gene) {
        sub_table = coeffs_table %>%
                filter(gene_short_name %in% gene)
        out_df = annotate_helper_1(sub_table)

        return(out_df)           
    }))
}

unique = F
df_sf_1 = convert_to_df(coeffs_table_dex_1, time = 4, normalization = "Conventional", unique = unique)
df_sf_10 = convert_to_df(coeffs_table_dex_10, time = 24, normalization = "Conventional", unique = unique)
df_hash_1 = convert_to_df(coeffs_table_dex_hash_1, time = 4, normalization = "Hash", unique = unique)
df_hash_10 = convert_to_df(coeffs_table_dex_hash_10, time = 24, normalization = "Hash", unique = unique)

df = rbind(df_sf_1, df_sf_10, df_hash_1, df_hash_10)

select_by = c("id", "gene_short_name", "num_cells_expressed", "term", "estimate", "q_value", "normalization", "time", "annotation")
term_estimates = df %>%
                    dplyr::select(select_by) 

final_df = term_estimates %>% 
    tidyr::spread(term, estimate, fill = 0) %>%
    dplyr::rename(hdaci_dex_estimate = DexTRUE)

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

out_df = NULL
for (normal in c("Conventional", "Hash")) {
    if (normal == "Conventional") {
        coeffs_table = coeffs_table_dex
    } else {
        coeffs_table = coeffs_table_dex_hash
    }
    
    for (tim in c(4,24)) {
        temp_df = final_df %>% filter(normalization == normal & time == tim)
        temp_df = coeffs_table %>% filter(id %in% temp_df$id) %>%
                        select(id, estimate) %>%
                        dplyr::rename(dex_estimate = estimate) %>%
                        left_join(temp_df, ., by = "id")
        out_df = rbind(out_df, temp_df)
    }
}

final_df = out_df

for (normal in c("Conventional", "Hash")) {
    if (normal == "Conventional") {
        coeffs_table = coeffs_table_dex
    } else {
        coeffs_table = coeffs_table_dex_hash
    }
    
    for (tim in c(4,24)) {
        sub_df = final_df %>% filter(normalization == normal & time == tim)
        temp_df = coeffs_table %>%
            filter(!id %in% sub_df$id) %>%
            mutate(normalization = normal, time = tim, 
                   hdaci_dex_estimate = 0,
                   annotation = "Dex") %>%
            dplyr::rename(dex_estimate = estimate) %>%
            select(id, gene_short_name, num_cells_expressed, normalization, time, 
                   annotation, hdaci_dex_estimate, dex_estimate, q_value)
        final_df = rbind(final_df, temp_df)
    }
    
}

final_df = final_df %>% select(id:hdaci_dex_estimate, dex_estimate, q_value)

head(final_df)

overall_estimate_df = final_df %>%
        rowwise() %>%
        mutate(respond = if (annotation == "HDACi") {
            "Yes"
        } else {
            "No"
        }) %>%
        as.data.frame() %>%
        arrange(gene_short_name)

saveRDS(overall_estimate_df, "overall_estimate_df_by_time.rds")


## figure 4C
count_estimate_df = overall_estimate_df %>%
        filter(normalization == "Hash") %>%
        as.data.frame() %>%
        mutate(type = ifelse(dex_estimate > 0, "Upregulated", "Downregulated")) %>%
        count(time, respond, type) %>%
        mutate(sample_size = n) %>%
        as.data.frame() %>%
        distinct(time, respond, type, .keep_all = T) %>%
        group_by(time) %>%
        add_tally(name = "m") %>%
        mutate(n = n / m)


count_estimate_df

offset = 0.02
options(repr.plot.width=5.5, repr.plot.height=2.5)
ggplot(count_estimate_df %>%
           transform(respond = factor(respond, levels = c("Yes", "No"))), 
       aes(x = respond, y = n, fill = as.factor(time))) + #reorder(aa, -n))) + #as.factor(dose))) + 
    geom_bar(stat = "identity", position = "dodge", colour = "white") + 
    labs(x = "", y = "Proportion of \ndex response genes", fill = "Effects of HDAC\ninhibition") + 
    scale_fill_manual(values = c("#ffdd55", "#0571b0")) + 
    facet_wrap(~ type, nrow = 1) + 
    geom_text(aes(y = n + offset, label = sample_size), 
             position = position_dodge(width=0.9), size = 3.1) + 
    monocle_theme_opts() # + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

ggsave("F4_anno_bar_new_respond.pdf", device = "pdf", width=5, height=2.5)


## Look at DEX response genes in presence of HDACi
overall_estimate_df = readRDS("overall_estimate_df_by_time.rds")
count_estimate_df = overall_estimate_df %>%
        filter(normalization == "Hash") %>%
        as.data.frame() %>%
        count(time, respond) %>%
        group_by(time) %>%
        mutate(sample_size = n / sum(n)) %>%
        as.data.frame() 
count_estimate_df


## figure 4D
response_gene = overall_estimate_df %>%
    filter(normalization == "Hash") %>%
    filter(respond == "Yes") %>%
    mutate(hdaci_effect = ifelse(time == 4, "Acetylation", "Metabolic"))

acetyl_genes = response_gene %>%
    filter(time == 4)
meta_genes = response_gene %>%
    filter(time == 24)

min_expr = 1
num_cells = 20
cds = readRDS("cds_dex_hash.rds")
cds = detect_genes(cds, min_expr = min_expr)
cds = cds[rowData(cds)$num_cells_expressed > num_cells,]
cds$Align = paste0(substr(cds$Condition, 1,1), "_", substr(cds$Dex, 1,1))
cds$Align = factor(cds$Align, levels = c("D_F", "D_T", "A_F", "A_T", "M_F", "M_T"))

norm_name = "Hash" 
dose_final_df = overall_estimate_df %>%
        filter(normalization == norm_name)

cds = cds[rowData(cds)$id %in% unique(dose_final_df$id),]

factor_list = c("DMSO_TRUE", "Acetylation_TRUE", "Metabolic_TRUE", 
                "DMSO_FALSE", "Acetylation_FALSE", "Metabolic_FALSE")

dose_mat = normalized_counts(cds, norm_method = "size_only")
colnames(dose_mat) = paste0(cds$Condition, "_", cds$Dex)
dose_aggregate <- as.data.frame(vapply(factor_list, function(x) 
      rowMeans(dose_mat[,colnames(dose_mat)== x,drop=FALSE], na.rm=TRUE),
                             numeric(nrow(dose_mat))))


dose_response_df = response_gene

dose_expr = t(apply(dose_response_df, 1, function(row) {
    expr_row = dose_aggregate[row[1], , drop=F]
    c(row[2],
      expr_row[1, "DMSO_TRUE"] / expr_row[1, "DMSO_FALSE"], 
      expr_row[1, paste0(row[11], "_TRUE")] / expr_row[1, paste0(row[11], "_FALSE")])
}))


expr_concat = dose_expr

response_df = data.frame(gene_short_name = expr_concat[,1],
                         Dex = as.double(expr_concat[,2]),
                         HDACi_dex = as.double(expr_concat[,3]),
                         hdaci_effect = dose_response_df[,11])

options(repr.plot.width=4, repr.plot.height=4)
ggplot(response_df, aes(log2(Dex), log2(HDACi_dex))) + 
    geom_point(size = 1.5, aes(color = hdaci_effect)) + 
    scale_color_manual(values = c("#ffdd55", "#0571b0")) + 
    geom_point(data = filter(response_df, gene_short_name == "TSC22D3"), aes(log2(Dex), log2(HDACi_dex)), color = "red") +
    geom_abline(slope = 1, colour = "black", linetype = "dashed", lwd = 0.8) + 
    geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + 
    labs(x = "log2(DEX FC)", y = "log2(DEX+HDACi FC)") + 
    monocle_theme_opts() + theme(legend.position = "none")

ggsave("F4_response_gene_scatter.pdf", device = "pdf", width=4, height=4)


response_df %>%
    mutate(Dex = log(Dex), HDACi_dex = log(HDACi_dex)) %>%
    filter(Dex > HDACi_dex & sign(Dex) == sign(HDACi_dex) & Dex > 0 | 
           Dex < HDACi_dex & sign(Dex) == sign(HDACi_dex) & Dex < 0) %>% 
    nrow

response_df %>%
    nrow

response_df %>%
    mutate(Dex = log2(Dex), HDACi_dex = log2(HDACi_dex)) %>%
    filter(gene_short_name == "TSC22D3")


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

gmt_dir = "../gmts"
## Load Gene Set Collections
reactomeGSC<-loadGSCSafe(file=file.path(gmt_dir, "Human_Reactome_August_01_2019_symbol.gmt"))
pantherGSC<-loadGSCSafe(file=file.path(gmt_dir, "Human_Panther_August_01_2019_symbol.gmt"))
KEGGGSC<-loadGSCSafe(file=file.path(gmt_dir, "Human_KEGG_August_01_2019_symbol.gmt"))
NCIGSC<-loadGSCSafe(file=file.path(gmt_dir, "Human_NCI_Nature_August_01_2019_symbol.gmt"))
GOGSC<-loadGSCSafe(file=file.path(gmt_dir, "Human_GO_bp_no_GO_iea_symbol.gmt"))
TFGSC<-loadGSCSafe(file=file.path(gmt_dir, "Human_MSigdb_August_01_2019_symbol.gmt"))
HALLMARKS<-loadGSCSafe(file=file.path(gmt_dir, "h.all.v6.2.symbols.gmt"))

min_expr = 1
num_cells = 20
cds = readRDS("cds_dex_hash.rds")
cds = detect_genes(cds, min_expr = min_expr)
cds = cds[rowData(cds)$num_cells_expressed > num_cells,]
cds$Align = paste0(substr(cds$Condition, 1,1), "_", substr(cds$Dex, 1,1))
cds$Align = factor(cds$Align, levels = c("D_F", "D_T", "A_F", "A_T", "M_F", "M_T"))
dex_univ = unique(rowData(cds)$gene_short_name)

gsa_hyper_clusters <- function(universe, dose_df, gsc){
    
    genes = dose_df %>% pull(gene_short_name)
    GSA_list = 
        tryCatch({
            gsares <- runGSAhyper(genes = genes, universe = universe, gsc=gsc)
            return (gsares)
        }, error = function(e) { print(e); return (NA) } )
    
    return(GSA_list)
}

GSC = HALLMARKS

gsa_list_hash = suppressWarnings({
                    gsa_hyper_clusters(dex_univ, meta_genes %>%
                                            filter(id %in% acetyl_genes$id), GSC)})


show_gsa = function(gsa_list) {
    go_df = data.frame(x = gsa_list$p.adj) %>% 
                        rownames_to_column() %>% 
                        arrange(x) %>% 
                        filter(x < 0.05)
    display(head(go_df, 30))
}

show_gsa(gsa_list_hash)


## figure S10B
go_df = data.frame(x = gsa_list_hash$p.adj) %>% 
                        rownames_to_column() %>% 
                        arrange(x) %>% 
                        filter(x < 0.05) %>%
                        mutate(x = -log10(x))

options(repr.plot.width=6, repr.plot.height=3)
ggplot(go_df, aes(x = x, y = reorder(rowname, x))) + 
    geom_bar(stat = "identity") + 
    labs(x = "", y = "") + 
    monocle_theme_opts()

ggsave("S_response_GSEA.pdf", device = "pdf", width=6, height=3)


## figure S10A
cds_plot = cds[, cds$Align %in% c("A_F", "A_T", "M_F", "M_T")]
genes_to_plot = convert_gene_to_id(cds_plot, c("SQSTM1", "DUSP1", "TIPARP", "ERRFI1", "SERPINE1", "PTGS2"))

options(repr.plot.width=8, repr.plot.height=2)
plot_genes_violin(cds_plot[genes_to_plot,], 
                  group_cells_by="Align", ncol=6) + 
    scale_fill_manual(values = c("#fff6d5", "#ffdd55", "#92c5de", "#0571b0"))

ggsave("S_response_gene_violin.pdf", device = "pdf", width=8, height=2)


## figure S10C,D
response_df_comp = response_df %>%
                    mutate(Dex = log2(Dex), HDACi_dex = log2(HDACi_dex)) %>%
                    add_count(gene_short_name) %>% filter(n == 2) %>%
                    mutate(fc = HDACi_dex / Dex) %>%
                    arrange(gene_short_name, hdaci_effect) 

response_df_comp$fc_comp = rep(sapply(unique(response_df_comp$gene_short_name), function(gene) {
    sub_acetyl = response_df_comp %>%
        filter(gene_short_name == gene & hdaci_effect == "Acetylation")
    sub_meta = response_df_comp %>%
        filter(gene_short_name == gene & hdaci_effect != "Acetylation")  
    return(sub_acetyl$fc / sub_meta$fc)
}), each = 2)

response_df_comp = response_df_comp %>%
        distinct(gene_short_name, .keep_all = T)

options(repr.plot.width=2, repr.plot.height=2)
ggplot(response_df_comp, aes(x = fc_comp)) +
    geom_histogram(bins = 20, fill = "white", color = "black") + 
    monocle_theme_opts()

ggsave("S_response_hist_4hr_attentuation.pdf", device = "pdf", width=2, height=2)

options(repr.plot.width=2, repr.plot.height=2)
ggplot(response_df_comp %>%
           mutate(type = ifelse(fc_comp < 1, "Acetyl", "Meta")) %>%
           dplyr::count(type), aes(x = type, y = n)) +
    geom_bar(stat = "identity", aes(fill = type)) + 
    monocle_theme_opts() + theme(legend.position = "none")

ggsave("S_response_bar_4hr_attentuation.pdf", device = "pdf", width=2, height=2)

