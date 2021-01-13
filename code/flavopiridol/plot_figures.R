suppressPackageStartupMessages({
    library(monocle3)
    library(MASS)
    library(reshape2)
    library(VGAM)
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
    library(dplyr)
})

wd = getwd()
wd

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

# load data
setwd("../../data/flavopiridol/")
hashTable = read.table("hashTable.txt", header=T)
hashTable_unique = hashTable %>% distinct(Cell, .keep_all=T)
hashTable_filtered = read.table("hashTable_filtered.txt", header=T)
hashTable_filtered_unique = read.table("hashTable_unique_filtered.txt", header=T)

cds = readRDS("cds_conv_pseudotime.rds")
cds_hash = readRDS("cds_hash_ladder_pseudotime.rds")

original = c("DMSO", "FP_01hr", "FP_03hr", "FP_06hr", "FP_12hr", "FP_24hr")
convertTo = c(0, 1, 3, 6, 12, 24)
cds$Time = convertTo[match(cds$Condition, original)]
cds_hash$Time = convertTo[match(cds_hash$Condition, original)]
metadata = as.data.frame(colData(cds_hash))

metadata %>%
    group_by(Time) %>%
    summarise(a = median(Total_RNA), b = mean(Total_RNA))


## figure 2B
options(repr.plot.width=4, repr.plot.height=3)
ggplot(metadata, aes(y = log10(Total_RNA), x = as.factor(Time))) + 
  geom_boxplot(aes(fill = as.factor(Time))) + 
  labs(x = "", y = "") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                     panel.border = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "none") + 
  scale_fill_brewer(palette = "YlGnBu")
ggsave("F2_fp_total_rna.pdf", device = "pdf", width = 4, height = 3)


## figure 2C
df_sf = data.frame(UMAP1 = cds@reducedDims$UMAP[,1],
                   UMAP2 = cds@reducedDims$UMAP[,2], 
                   Time = cds$Time)

df_hash = data.frame(UMAP1 = cds_hash@reducedDims$UMAP[,1],
                     UMAP2 = cds_hash@reducedDims$UMAP[,2], 
                     Time = cds_hash$Time)

options(repr.plot.width=3, repr.plot.height=3)
ggplot(df_sf, aes(UMAP1, UMAP2)) +
    geom_point(aes(color = as.factor(Time))) + 
    scale_color_brewer(palette = "YlGnBu") + 
    monocle_theme_opts() + theme(legend.position = "none")
ggsave("F2_fp_umap_conv.pdf", device = "pdf", width = 3, height = 3)

options(repr.plot.width=3, repr.plot.height=3)
ggplot(df_hash, aes(UMAP1, UMAP2)) +
    geom_point(aes(color = as.factor(Time))) + 
    scale_color_brewer(palette = "YlGnBu") + 
    monocle_theme_opts() + theme(legend.position = "none")
ggsave("F2_fp_umap_ladder.pdf", device = "pdf", width = 3, height = 3)


## figure 2D
degs_compared = readRDS("de_analysis/degs_sf_LRT.rds")
degs_sf = readRDS("de_analysis/coeff_table_sf.rds")
qvalue = 1e-2
num_cells = 30

degs_compared = degs_compared %>% filter(q_value < qvalue)
degs_sf = degs_sf %>% 
    filter(id %in% degs_compared$id & grepl("Pseudo", term)) %>%
    filter(q_value < qvalue & num_cells_expressed > num_cells) %>%
    arrange(id, estimate) %>%
    distinct(id, .keep_all = TRUE) %>% 
    arrange(q_value)

degs_compared = readRDS("de_analysis/degs_hash_LRT.rds")
degs_hash = readRDS("de_analysis/coeff_table_hash.rds")

degs_compared = degs_compared %>% filter(q_value < qvalue)
degs_hash = degs_hash %>% 
    filter(id %in% degs_compared$id & grepl("Pseudo", term)) %>%
    filter(q_value < qvalue & num_cells_expressed > num_cells) %>%
    arrange(id, estimate) %>%
    distinct(id, .keep_all = TRUE) %>% 
    arrange(q_value)

fc_m_sf = readRDS("de_analysis/fp_sf_fc_log2.rds")
fc_m_hash = readRDS("de_analysis/fp_hash_fc_log2.rds")

cutoff = 1

genes_sf = fc_m_sf[abs(fc_m_sf) > cutoff]
genes_sf = genes_sf[names(genes_sf) %in% degs_sf$id]
degs_sf = degs_sf %>% filter(id %in% names(genes_sf)) %>% mutate(fc = genes_sf)

genes_hash = fc_m_hash[abs(fc_m_hash) > cutoff]
genes_hash = genes_hash[names(genes_hash) %in% degs_hash$id]
degs_hash = degs_hash %>% 
                filter(id %in% names(genes_hash)) %>% 
                mutate(fc = genes_hash)

degs_table_sf = degs_sf

degs_table_sf = degs_table_sf %>%
    dplyr::count(Type = fc > 0, name = "Count") %>%
    group_by(Type)


degs_table_sf$Type = ifelse(degs_table_sf$Type, "Upregulated", "Downregulated")
degs_table_sf$Normalization = "Conventional"

df = degs_table_sf

degs_table_hash = degs_hash

degs_table_hash = degs_table_hash %>%
    dplyr::count(Type = fc > 0, name = "Count") %>%
    group_by(Type)

degs_table_hash$Type = ifelse(degs_table_hash$Type, "Upregulated", "Downregulated")
degs_table_hash$Normalization = "Hash ladder"


df = rbind(df, degs_table_hash)

options(repr.plot.width=3, repr.plot.height=3)
bar_order = c("Upregulated", "Downregulated")

ggplot(base::transform(df, 
                 Type = factor(Type, levels = bar_order)),
    aes(x = Type, y = Count, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +  
  geom_text(aes(label = Count, y = Count),
        position = position_dodge(0.9),
        size = 4, vjust = -0.2) +
  labs(x = "", y = "# of DE genes \n as a function of treatment time", fill = "") + 
  scale_fill_manual(values = c("#be414c", "#40bfaa")) + 
  ylim(c(0,475)) + 
  facet_wrap(vars(Normalization)) +  
  monocle_theme_opts() + 
  theme(strip.text.x = element_text(size=10), 
        axis.text.x = element_blank(), #element_text(angle = 45, hjust = 1) 
        plot.title = element_text(hjust = 0.5), legend.position = "none")

ggsave("F2_degs_bar.pdf", device = "pdf", width = 3, height = 3)


## figure 2E
convert_gene_to_id <- function(cds, genes) {
    rowData_mat = as.data.frame(rowData(cds))
    return(rowData_mat %>% filter(gene_short_name %in% genes) %>% pull(id))
}

convert_id_to_gene <- function(cds, ids) {
    rowData_mat = as.data.frame(rowData(cds))
    return(rowData_mat %>% filter(id %in% ids) %>% pull(gene_short_name))
}

fit_spline = function(data_frame, x, y, n, label, method = "natural") {
    fit = as.data.frame(spline(x, as.vector(unlist(data_frame[y])), n = n, method = method)) %>%
                mutate(feature_label = label)
    return(fit)
}


genes_to_plot = c("LAMC1", "CD151")
ids_to_plot = convert_gene_to_id(cds, genes_to_plot)

options(repr.plot.width=8, repr.plot.height=4)
g = list()
g[[1]] = plot_percent_cells_positive(cds_hash[ids_to_plot,], min_expr = 1,
                                     group_cells_by="Time", ncol=2) +
            theme(axis.text.x=element_text(angle=45, hjust=1))

g[[2]] = plot_percent_cells_positive(cds[ids_to_plot,], min_expr = 1,
                                     group_cells_by="Time", ncol=2) +
                theme(axis.text.x=element_text(angle=45, hjust=1))


a = (g[[1]]$data) %>% filter(feature_label == "CD151")
b = (g[[1]]$data) %>% filter(feature_label != "CD151")

spline_fit_mult = rbind(fit_spline(a, x = 1:6, y = "target_fraction_mean", n = 36, label = "CD151"),
                        fit_spline(b, x = 1:6, y = "target_fraction_mean", n = 36, label = "LAMC1"))

a = (g[[2]]$data) %>% filter(feature_label == "CD151")
b = (g[[2]]$data) %>% filter(feature_label != "CD151")

spline_fit_sf = rbind(fit_spline(a, x = 1:6, y = "target_fraction_mean", n = 36, label = "CD151"),
                      fit_spline(b, x = 1:6, y = "target_fraction_mean", n = 36, label = "LAMC1"))

gg = list()


gg[[1]] = ggplot(g[[1]]$data, aes(x = as.factor(Time), y = target_fraction_mean)) + 
                geom_bar(stat="identity", aes(fill = as.factor(Time))) +
                geom_linerange(aes(ymin = target_fraction_low, ymax = target_fraction_high)) +
                facet_wrap(~ feature_label, scales = "free") + 
                labs(x = "", y = "") + 
                scale_fill_brewer(palette = "YlGnBu") + 
                monocle_theme_opts() + theme(legend.position = "none")
gg[[2]] = ggplot(g[[2]]$data, aes(x = as.factor(Time), y = target_fraction_mean)) + 
                geom_bar(stat="identity", aes(fill = as.factor(Time))) +
                geom_linerange(aes(ymin = target_fraction_low, ymax = target_fraction_high)) +
                facet_wrap(~ feature_label, scales = "free") + 
                labs(y = "", x = "") +
                scale_fill_brewer(palette = "YlGnBu") + 
                monocle_theme_opts() + theme(legend.position = "none")

options(repr.plot.width=4, repr.plot.height=3.5)
cc = do.call("grid.arrange", c(gg, ncol=1))
ggsave("F2_cells_perc_positive_adhesion.pdf", cc, device = "pdf", width = 4, height = 3.5)


## figure 2F
degs_common = inner_join(degs_sf %>% select(id, gene_short_name, estimate, normalized_effect, fc, q_value), 
                         degs_hash %>% select(id, gene_short_name, estimate, normalized_effect, fc, q_value), 
                         by = c("id", "gene_short_name"))

ddf = data.frame(fc = degs_common$fc.y,
                 fcc = degs_common$fc.y / degs_common$fc.x)

ddf$Type = ifelse(ddf$fc < 0, "Downregulated", "Upregulated")

bar_order = c("Upregulated", "Downregulated")

options(repr.plot.width=4, repr.plot.height=2.5)

ggplot(ddf, aes(x = Type, y = fcc)) +  
    geom_hline(yintercept = 1, linetype = 2, colour = "red") +
    geom_violin(adjust = 1.2, draw_quantiles = T, aes(fill = Type)) +
    geom_boxplot(width=.1, outlier.colour=NA, position = "dodge") +
    scale_fill_manual(values = c("#40bfaa", "#be414c")) + 
    coord_flip() +
    labs(x = "", y = "") + 
    monocle_theme_opts() + theme(legend.position = "none")

ggsave("F2_fc_comp.pdf", device = "pdf", width = 4, height = 2.5)

ddf %>% count(Type)



## figure S3A
options(repr.plot.width=4, repr.plot.height=3)
ggplot(metadata, aes(y = Size_Factor, x = as.factor(Time))) + 
    geom_boxplot(aes(fill = as.factor(Time))) + 
    labs(x = "", y = "") +
    scale_fill_brewer(palette = "YlGnBu") +
    monocle_theme_opts() + theme(legend.position = "none")

ggsave("S_fp_hash_size_factor.pdf", device = "pdf", width = 4, height = 3)


## figure S3B
options(repr.plot.width=4, repr.plot.height=3)
ggplot(metadata, aes(y = log10(Total_hash), x = as.factor(Time))) + 
    geom_boxplot(aes(fill = as.factor(Time))) + 
    labs(x = "", y = "") +
    scale_fill_brewer(palette = "YlGnBu") +
    monocle_theme_opts() + theme(legend.position = "none")

ggsave("S_fp_total_hash.pdf", device = "pdf", width = 4, height = 3)


## figure S3C
setwd("/net/trapnell/vol1/khj3017/project/190118_HashLadder_sciRNA/results/sciRNA-FP-timecourse-48-hashesV4-pooled-ID")
hash_id_df = read.table("hash_id_df.txt", header=T)
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


## figure S3D
options(repr.plot.width=3.5, repr.plot.height=3)
ggplot(hashTable_filtered %>%
           dplyr::count(Cell) %>% dplyr::count(n), 
       aes(x = n, y = nn)) + 
      geom_bar(stat = "identity") +
      labs(x = "Number of unique hash oligos", y = "Number of cells") +
      monocle_theme_opts()

ggsave("S_n_hash_plot.pdf", width = 3.5, height = 3)


options(repr.plot.width=6, repr.plot.height=3)
g = list()
g[[1]] = ggplot(hashTable_unique, aes(x = log10(Total_hash))) + 
      geom_histogram(bins = 30) + #, colour = "black", fill = "white") + #ylim(c(0,10)) + 
      geom_vline(xintercept = log10(200), color = "red") + 
      labs(x = "", y = "") +
      monocle_theme_opts()

hashTable_unique %>% dplyr::count(Total_hash > 200)

g[[2]] = ggplot(hashTable_unique, aes(x = Rsq)) + 
      geom_histogram(bins = 30) + #, colour = "black", fill = "white") + #ylim(c(0,10)) + 
      geom_vline(xintercept = 0.7, color = "red") + 
      labs(x = "", y = "") +
      monocle_theme_opts()

hashTable_unique %>% dplyr::count(Rsq > 0.7)

gg = do.call("grid.arrange", c(g, ncol=2))
ggsave(plot = gg, filename = "S_hash_filter.pdf", width = 6, height = 3)


## figure S4A
cds_plot = cds
cds_plot$pseudotime = pseudotime(cds_plot)

g = list()

g[[1]] = ggplot(as.data.frame(colData(cds_plot)), 
       aes(x = pseudotime, y = as.factor(Time))) + 
    geom_density_ridges2(aes(height = ..density.. , fill = as.factor(Time))) + 
    labs(x = "Pseudotime", y = "Time (hrs)") +
    scale_fill_brewer(palette = "YlGnBu") + 
    theme_ridges() + theme(legend.position = "none") #, axis.text.x = element_blank())


cds_plot = cds_hash
cds_plot$pseudotime = pseudotime(cds_plot)

g[[2]] = ggplot(as.data.frame(colData(cds_plot)), 
       aes(x = pseudotime, y = as.factor(Time))) + 
    geom_density_ridges2(aes(height = ..density.. , fill = as.factor(Time))) + 
    labs(x = "Pseudotime", y = "Time (hrs)") +
    scale_fill_brewer(palette = "YlGnBu") + 
    theme_ridges() + theme(legend.position = "none") #, axis.text.x = element_blank())

options(repr.plot.width=5.5, repr.plot.height=3)
gg = do.call("grid.arrange", c(g, ncol=2))
ggsave(plot = gg, "S_FP_pseudo_ordering.pdf", device = "pdf", width = 5.5, height = 3)


## figure S4B
n_hash = degs_hash %>%
    filter(!id %in% degs_sf$id) %>%
    dplyr::count(fc > 0)

n_sf = degs_sf %>%
    filter(!id %in% degs_hash$id) %>%
    dplyr::count(fc > 0)

n_common = degs_sf %>%
    filter(id %in% degs_hash$id) %>%
    dplyr::count(fc > 0)

options(repr.plot.width=6, repr.plot.height=3)
g = list()
gg = list()

g[[1]] = euler(c("Conventional" = n_sf[2,]$n, "Hash ladder" = n_hash[2,]$n, "Conventional&Hash ladder" = n_common[2,]$n))

g[[2]] = euler(c("Conventional" = n_sf[1,]$n, "Hash ladder" = n_hash[1,]$n, "Conventional&Hash ladder" = n_common[1,]$n))


gg[[1]] = plot(g[[1]], 
     family = "Helveltica",
     fills = list(fill = c("yellow", "darkgreen"), alpha = 0.7),
     quantities = F, 
     labels = NULL)

gg[[2]] = plot(g[[2]], 
     family = "Helveltica",
     fills = list(fill = c("yellow", "darkgreen"), alpha = 0.7),
     quantities = F, 
     labels = NULL)

ggg = do.call("grid.arrange", c(gg, ncol=2))
ggsave(plot = ggg, "S_FP_de_overlap.pdf", device = "pdf", width = 6, height = 3)


## figure S4C
cds = detect_genes(cds, min_expr = 0.1)
cds = cds[rowData(cds)$num_cells_expressed > 10, ]
cds_hash = detect_genes(cds_hash, 0.1)
cds_hash = cds_hash[rowData(cds_hash)$num_cells_expressed > 10, ]

goi = rowData(cds)$id 
sf_gene_expr = as.matrix(t(t(counts(cds[goi,])) / cds$Size_Factor))
colnames(sf_gene_expr) = as.factor(cds$Condition)

hash_gene_expr = as.matrix(t(t(counts(cds_hash[goi,])) / cds_hash$Size_Factor))
colnames(hash_gene_expr) = as.factor(cds_hash$Condition)

sf_mean = t(aggregate(t(sf_gene_expr), by = list(Time = colnames(sf_gene_expr)), FUN=mean) %>%
                column_to_rownames(var = "Time"))
hash_mean = t(aggregate(t(hash_gene_expr), by = list(Time = colnames(hash_gene_expr)), FUN=mean) %>%
                column_to_rownames(var = "Time"))
sf_sd = t(aggregate(t(sf_gene_expr), by = list(Time = colnames(sf_gene_expr)), FUN=sd) %>%
                column_to_rownames(var = "Time"))
hash_sd = t(aggregate(t(hash_gene_expr), by = list(Time = colnames(hash_gene_expr)), FUN=sd) %>%
                column_to_rownames(var = "Time"))

sf_fc = data.frame(FP_1hr = sf_mean[,2] / sf_mean[,1],
                   FP_3hr = sf_mean[,3] / sf_mean[,1],
                   FP_6hr = sf_mean[,4] / sf_mean[,1],
                   FP_12hr = sf_mean[,5] / sf_mean[,1],
                   FP_24hr = sf_mean[,6] / sf_mean[,1])

hash_fc = data.frame(FP_1hr = hash_mean[,2] / hash_mean[,1],
                   FP_3hr = hash_mean[,3] / hash_mean[,1],
                   FP_6hr = hash_mean[,4] / hash_mean[,1],
                   FP_12hr = hash_mean[,5] / hash_mean[,1],
                   FP_24hr = hash_mean[,6] / hash_mean[,1])


options(repr.plot.width=8, repr.plot.height=3)
sf_fc_melt = melt(sf_fc)
hash_fc_melt = melt(hash_fc)


sf_fc_melt$type = "Conventional"
hash_fc_melt$type = "hashinomial"
df = rbind(sf_fc_melt, hash_fc_melt)

options(repr.plot.width=7, repr.plot.height=3)
ggplot(df, aes(x = variable, y = log2(value), colour = type)) + 
    geom_hline(yintercept = 0, color = "blue", linetype = "dashed") + 
    geom_violin(aes(fill = type)) + 
    geom_boxplot(width=0.2, position = position_dodge(0.9), outlier.shape = NA) +
    scale_color_manual(values = c("black", "black")) + 
    scale_fill_manual(values = c("#00AFBB", "#FC4E07")) +
    monocle_theme_opts()

ggsave("S_de_time_fc_all_genes.pdf", width = 7, height = 3)


## figure S4D
genes_to_plot = c("ALCAM", "AMOT", "PPP3CA", "EFNA5", "AKT3", "PARD3")
ids_to_plot = convert_gene_to_id(cds, genes_to_plot)

options(repr.plot.width=8, repr.plot.height=4)
g = plot_percent_cells_positive(cds_hash[ids_to_plot,], min_expr = 1,
                                     group_cells_by="Time", ncol=3) +
            theme(axis.text.x=element_text(angle=45, hjust=1))

options(repr.plot.width=5, repr.plot.height=4)
ggplot(g$data, aes(x = as.factor(Time), y = target_fraction_mean)) + 
        geom_bar(stat="identity", aes(fill = as.factor(Time))) +
        geom_linerange(aes(ymin = target_fraction_low, ymax = target_fraction_high)) +
        facet_wrap(~ feature_label, scales = "free") + 
        labs(x = "", y = "") + 
        scale_fill_brewer(palette = "YlGnBu") + 
        monocle_theme_opts() + theme(legend.position = "none")

ggsave("S_gene_neurogenesis_motility.pdf", device = "pdf", width = 5, height = 4)


## figure S4E
goi = degs_sf %>% filter(fc < 0) %>% pull(id)
sf_gene_expr = as.matrix(t(t(counts(cds[goi,])) / cds$Size_Factor))
colnames(sf_gene_expr) = as.factor(cds$Condition)

goi = degs_hash %>% filter(fc < 0) %>% pull(id)
hash_gene_expr = as.matrix(t(t(counts(cds_hash[goi,])) / cds_hash$Size_Factor))
colnames(hash_gene_expr) = as.factor(cds_hash$Condition)

sf_mean = t(aggregate(t(sf_gene_expr), by = list(Time = colnames(sf_gene_expr)), FUN=mean) %>%
                column_to_rownames(var = "Time"))
hash_mean = t(aggregate(t(hash_gene_expr), by = list(Time = colnames(hash_gene_expr)), FUN=mean) %>%
                column_to_rownames(var = "Time"))
sf_sd = t(aggregate(t(sf_gene_expr), by = list(Time = colnames(sf_gene_expr)), FUN=sd) %>%
                column_to_rownames(var = "Time"))
hash_sd = t(aggregate(t(hash_gene_expr), by = list(Time = colnames(hash_gene_expr)), FUN=sd) %>%
                column_to_rownames(var = "Time"))

sf_fc = data.frame(FP_1hr = sf_mean[,2] / sf_mean[,1],
                   FP_3hr = sf_mean[,3] / sf_mean[,1],
                   FP_6hr = sf_mean[,4] / sf_mean[,1],
                   FP_12hr = sf_mean[,5] / sf_mean[,1],
                   FP_24hr = sf_mean[,6] / sf_mean[,1])

hash_fc = data.frame(FP_1hr = hash_mean[,2] / hash_mean[,1],
                   FP_3hr = hash_mean[,3] / hash_mean[,1],
                   FP_6hr = hash_mean[,4] / hash_mean[,1],
                   FP_12hr = hash_mean[,5] / hash_mean[,1],
                   FP_24hr = hash_mean[,6] / hash_mean[,1])

options(repr.plot.width=8, repr.plot.height=3)
sf_fc_melt = melt(sf_fc)
hash_fc_melt = melt(hash_fc)
sf_fc_melt$type = "Conventional"
hash_fc_melt$type = "Hash ladder"
df = rbind(sf_fc_melt, hash_fc_melt)

options(repr.plot.width=5, repr.plot.height=3)
ggplot(df, aes(x = variable, y = log2(value), fill = type)) + 
    geom_hline(yintercept = 0, color = "blue", linetype = "dashed") + 
    geom_boxplot() + 
    monocle_theme_opts()

ggsave("S_de_time_fc_comp_down.pdf", width = 5, height = 3)


goi = degs_sf %>% filter(fc > 0) %>% pull(id)
sf_gene_expr = as.matrix(t(t(counts(cds[goi,])) / cds$Size_Factor))
colnames(sf_gene_expr) = as.factor(cds$Condition)

goi = degs_hash %>% filter(fc > 0) %>% pull(id)
hash_gene_expr = as.matrix(t(t(counts(cds_hash[goi,])) / cds_hash$Size_Factor))
colnames(hash_gene_expr) = as.factor(cds_hash$Condition)

sf_mean = t(aggregate(t(sf_gene_expr), by = list(Time = colnames(sf_gene_expr)), FUN=mean) %>%
                column_to_rownames(var = "Time"))
hash_mean = t(aggregate(t(hash_gene_expr), by = list(Time = colnames(hash_gene_expr)), FUN=mean) %>%
                column_to_rownames(var = "Time"))
sf_sd = t(aggregate(t(sf_gene_expr), by = list(Time = colnames(sf_gene_expr)), FUN=sd) %>%
                column_to_rownames(var = "Time"))
hash_sd = t(aggregate(t(hash_gene_expr), by = list(Time = colnames(hash_gene_expr)), FUN=sd) %>%
                column_to_rownames(var = "Time"))

sf_fc = data.frame(FP_1hr = sf_mean[,2] / sf_mean[,1],
                   FP_3hr = sf_mean[,3] / sf_mean[,1],
                   FP_6hr = sf_mean[,4] / sf_mean[,1],
                   FP_12hr = sf_mean[,5] / sf_mean[,1],
                   FP_24hr = sf_mean[,6] / sf_mean[,1])

hash_fc = data.frame(FP_1hr = hash_mean[,2] / hash_mean[,1],
                   FP_3hr = hash_mean[,3] / hash_mean[,1],
                   FP_6hr = hash_mean[,4] / hash_mean[,1],
                   FP_12hr = hash_mean[,5] / hash_mean[,1],
                   FP_24hr = hash_mean[,6] / hash_mean[,1])

g = list()

options(repr.plot.width=8, repr.plot.height=3)
sf_fc_melt = melt(sf_fc)
hash_fc_melt = melt(hash_fc)
sf_fc_melt$type = "Conventional"
hash_fc_melt$type = "Hash ladder"
df = rbind(sf_fc_melt, hash_fc_melt)

options(repr.plot.width=5, repr.plot.height=3)
ggplot(df, aes(x = variable, y = log2(value), fill = type)) + 
    geom_hline(yintercept = 0, color = "blue", linetype = "dashed") + 
    geom_boxplot() + 
    monocle_theme_opts()

ggsave("S_de_time_fc_comp_up.pdf", width = 5, height = 3)


