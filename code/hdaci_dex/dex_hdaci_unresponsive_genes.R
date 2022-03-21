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
overall_estimate_df = readRDS("overall_estimate_df_by_time.rds")

no_response_genes = overall_estimate_df %>%
    filter(normalization == "Hash") %>%
    filter(respond == "No")  %>%
    mutate(type = ifelse(hdaci_estimate == 0, "Attenuated", 
                    ifelse(sign(hdaci_estimate) == sign(dex_estimate), "Saturated", "Dominated"))) %>%
    transform(type = factor(type, levels = c("Attenuated", "Saturated", "Dominated")))

acetyl_genes = no_response_genes %>%
    filter(time == 4)
meta_genes = no_response_genes %>%
    filter(time != 4)


## figure 4e
n_acetyl_genes = acetyl_genes %>%
                    filter(!id %in% meta_genes$id) %>%
                    nrow()

n_meta_genes = meta_genes %>%
                    filter(!id %in% acetyl_genes$id) %>%
                    nrow()

n_common_genes = meta_genes %>%
                    filter(id %in% acetyl_genes$id) %>%
                    nrow()

n_acetyl_genes
n_meta_genes
n_common_genes

options(repr.plot.width=6, repr.plot.height=3)

g = euler(c("Acetyl" = n_acetyl_genes, "Meta" = n_meta_genes, "Acetyl&Meta" = n_common_genes))

pdf("F4_unresponsive_venn.pdf", width = 6, height = 3)
plot(g, 
     family = "Helveltica",
     fills = list(fill = c("#ffdd55", "#0571b0"), alpha = 0.7),
     quantities = F, 
     labels = NULL, 
     legend = list(labels = c("Acetyl", "Meta"), fontsize = 16))
dev.off()


## Supplementary figure 15b 
options(repr.plot.width=4.5, repr.plot.height=2.5)
no_response_genes %>%
    add_count(time, name = "m") %>% 
    count(time, type, m) %>%
    mutate(sample_size = n) %>%
    mutate(n = n / m) %>%
    ggplot(., aes(x = type, y = n, fill = as.factor(time))) +
        geom_bar(stat = "identity", position = "dodge", colour = "white") + 
        geom_text(aes(label = sample_size, y = n + 0.001),
            position = position_dodge(0.9),
            size = 3,
            vjust = 0) + 
        labs(x = "", y = "Proportion of \ndex response genes", fill = "Effects of HDAC\ninhibition") + 
        scale_fill_manual(values = c("#ffdd55", "#0571b0")) + 
        #facet_wrap(~ dose, nrow = 1) + 
        ylim(c(0, 0.6)) + 
        monocle_theme_opts()

ggsave("S_anno_bar_unresponsive_dose_comb.pdf", device = "pdf", width=4.5, height=2.5)

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

## Load Gene Set Collections
gmt_dir = "../gmts"
reactomeGSC<-loadGSCSafe(file=file.path(gmt_dir, "Human_Reactome_August_01_2019_symbol.gmt"))
pantherGSC<-loadGSCSafe(file=file.path(gmt_dir, "Human_Panther_August_01_2019_symbol.gmt"))
KEGGGSC<-loadGSCSafe(file=file.path(gmt_dir, "Human_KEGG_August_01_2019_symbol.gmt"))
NCIGSC<-loadGSCSafe(file=file.path(gmt_dir, "Human_NCI_Nature_August_01_2019_symbol.gmt"))
GOGSC<-loadGSCSafe(file=file.path(gmt_dir, "Human_GO_bp_no_GO_iea_symbol.gmt"))
TFGSC<-loadGSCSafe(file=file.path(gmt_dir, "Human_MSigdb_August_01_2019_symbol.gmt"))
HALLMARKS<-loadGSCSafe(file=file.path(gmt_dir, "h.all.v6.2.symbols.gmt"))

cds_hash = readRDS("cds_dex_hash.rds")
cds_hash = detect_genes(cds_hash, min_expr = 1)
cds_hash = cds_hash[rowData(cds_hash)$num_cells_expressed > 10,]
cds_hash$Align = paste0(substr(cds_hash$Condition, 1,1), "_", substr(cds_hash$Dex, 1,1))
cds_hash$Align = factor(cds_hash$Align, levels = c("D_F", "D_T", "A_F", "A_T", "M_F", "M_T"))


# shared
gsa_hyper_clusters <- function(cds, gene_df, gsc) {
    universe = unique(rowData(cds)$gene_short_name)
    
    genes = gene_df %>%
                add_count(id, name = "m") %>% 
                filter(m == 2) %>%
                distinct(id, .keep_all=T) %>% 
                pull(gene_short_name)
    GSA_res = runGSAhyper(genes = genes, universe = universe, gsc=gsc)

    return(GSA_res)
}

no_response_genes %>%
        add_count(id, name = "m") %>% 
        filter(m == 2) %>%
        distinct(id, .keep_all=T)

# shared genes
GSC = HALLMARKS

gsa_hash = suppressWarnings({gsa_hyper_clusters(cds_hash, no_response_genes, GSC)})


show_gsa = function(gsa_hash) {
    go_df = data.frame(x = gsa_hash$p.adj) %>% 
                     rownames_to_column() %>% 
                      arrange(x) %>% 
                      filter(x < 0.05)
    display(head(go_df, 30))

}

show_gsa(gsa_hash)

# Supplementary Figure 14
## Supplementary figure 14a
go_df = data.frame(x = gsa_hash$p.adj) %>% 
                    rownames_to_column() %>% 
                    arrange(x) %>% 
                    filter(x < 0.05) %>%
                    mutate(x = -log10(x))
options(repr.plot.width=6, repr.plot.height=3)
ggplot(go_df, aes(x = x, y = reorder(rowname, x))) + 
    geom_bar(stat = "identity") + 
    labs(x = "", y = "") + 
    monocle_theme_opts()

ggsave("S_no_response_GSEA.pdf", device = "pdf", width=6, height=3)

no_response_genes %>%
    filter(time != 4) %>%
    filter(gene_short_name %in% sort(GSC$gsc[['HALLMARK_REACTIVE_OXIGEN_SPECIES_PATHWAY']]))

no_response_genes %>%
    filter(time != 4) %>%
    filter(gene_short_name %in% sort(GSC$gsc[['HALLMARK_XENOBIOTIC_METABOLISM']]))


## Supplementary figure 14b
cds_plot = cds_hash
genes_to_plot = convert_gene_to_id(cds_plot, c("ABCC1", "GCLC", "GCLM", "MGST1"))

options(repr.plot.width=4, repr.plot.height=4)
plot_genes_violin(cds_plot[genes_to_plot,], normalize = T,
                  group_cells_by="Align", ncol=2) + 
    scale_fill_manual(values = c("grey90", "grey60", "#fff6d5", "#ffdd55", "#92c5de", "#0571b0"))

ggsave("S_no_response_gene_violin.pdf", device = "pdf", width=4, height=4)


## 4hrs vs 24hrs
gsa_hyper_clusters <- function(cds, gene_df, gsc){
    universe = unique(rowData(cds)$gene_short_name)
    GSA_list = list()
    i = 1
    for (t in c(4,24)) {
         genes = gene_df %>%
                    #add_count(id, name = "m") %>% filter(m == 1) %>%
                    filter(time == t) %>%
                    distinct(id, .keep_all=T) %>% 
                    pull(gene_short_name)
        GSA_list[[i]] = runGSAhyper(genes = genes, universe = universe, gsc=gsc)
        i = i + 1
    }
    return(GSA_list)
}

GSC = HALLMARKS

gsa_hash = suppressWarnings({gsa_hyper_clusters(cds_hash, no_response_genes, GSC)})


show_gsa = function(gsa_list) {
    for (g in gsa_list) {
        go_df = data.frame(x = g$p.adj) %>% 
                    rownames_to_column() %>% 
                    arrange(x) %>% 
                    filter(x < 0.05)
        display(head(go_df, 30))
    }

}

show_gsa(gsa_hash)


no_response_genes %>%
    add_count(id) %>% filter(n == 1 & time == 4) %>%
    arrange(gene_short_name)

no_response_genes %>%
    add_count(id) %>% filter(n == 1 & time == 24) %>%
    arrange(gene_short_name)


## figure 4f
cds_plot = cds_hash[, cds_hash$Align %in% c("D_F", "D_T", "A_F", "A_T")]
genes_to_plot = convert_gene_to_id(cds_plot, c("DGKH", "DOCK4", "CEBPB"))

options(repr.plot.width=5, repr.plot.height=2)
plot_genes_violin(cds_plot[genes_to_plot,], 
                  group_cells_by="Align", ncol=3) + 
    scale_fill_manual(values = c("grey90", "grey60", "#fff6d5", "#ffdd55"))

ggsave("F4_violin_immune_acetyl.pdf", device = "pdf", width=5, height=2)

cds_plot = cds_hash[, cds_hash$Align %in% c("D_F", "D_T", "M_F", "M_T")]
genes_to_plot = convert_gene_to_id(cds_plot, c("ALDH1A1", "BMPR1B", "UGDH"))

options(repr.plot.width=5, repr.plot.height=2)
plot_genes_violin(cds_plot[genes_to_plot,], 
                  group_cells_by="Align", ncol=3) + 
    scale_fill_manual(values = c("grey90", "grey60", "#92c5de", "#0571b0"))

ggsave("F4_violin_fa_meta.pdf", device = "pdf", width=5, height=2)


