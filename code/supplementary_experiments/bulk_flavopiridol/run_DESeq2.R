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
    library(DESeq2)
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

setwd("../../../data/supplementary_experiments/bulk_flavopiridol")


# Run differention expression
cds = readRDS("cds_bulk_conventional.rds")
cds_ERCC = readRDS("cds_bulk_ERCC.rds")
cds$FP = as.factor(cds$FP)
cds_ERCC$FP = as.factor(cds_ERCC$FP)

## Conventional normalization
dds=DESeqDataSetFromMatrix(countData = counts(cds_ERCC), 
                           colData = colData(cds_ERCC), design = ~ FP)
dds = estimateSizeFactors(dds)
log_transformed_counts = log2(counts(dds, normalized=TRUE) + 1)

# shifted log of normalized counts
se <- SummarizedExperiment(log2(counts(dds, normalized=TRUE) + 1),
                           colData=colData(dds))

options(repr.plot.width=4, repr.plot.height=3)
plotPCA( DESeqTransform( se ), intgroup="FP", ntop = 500) +
    theme_minimal()

deseq_res = DESeq(dds)

de_res_all = list()
FP_time = unique(dds$FP)[-1]
for (time in FP_time) {
    res = results(deseq_res, contrast=c("FP",time, "0")) %>%
                as.data.frame() %>%
                na.omit()
    res = res[order(res$padj), , drop=F]
    res = as.data.frame(res) %>%
            rownames_to_column("id") %>%
            filter(padj < 0.01) %>%
            left_join(., as.data.frame(rowData(cds_ERCC)), by = "id") %>%
            select(id, gene_short_name, baseMean:padj) %>%
            na.omit() %>% arrange(padj)
    de_res_all[[time]] = res
}

lapply(de_res_all, dim)
saveRDS(deseq_res, "dds_ERCC_conventional.rds")
saveRDS(de_res_all, "de_res_all_ERCC_conventional.rds")


## ERCC normalization
dds=DESeqDataSetFromMatrix(countData = counts(cds_ERCC), 
                           colData = colData(cds_ERCC), design = ~ FP)
dds = estimateSizeFactors(dds)
log_transformed_counts = log2(counts(dds, normalized=TRUE) + 1)
dds$sizeFactor = dds$Size_Factor

# shifted log of normalized counts
se <- SummarizedExperiment(log2(counts(dds, normalized=TRUE) + 1),
                           colData=colData(dds))

options(repr.plot.width=4, repr.plot.height=3)
plotPCA( DESeqTransform( se ), intgroup="FP", ntop = 500) +
    theme_minimal()

deseq_res = DESeq(dds)

de_res_all = list()
FP_time = unique(dds$FP)[-1]
for (time in FP_time) {
    res = results(deseq_res, contrast=c("FP",time, "0")) %>%
                as.data.frame() %>%
                na.omit()
    res = res[order(res$padj), , drop=F]
    res = as.data.frame(res) %>%
            rownames_to_column("id") %>%
            filter(padj < 0.01) %>%
            left_join(., as.data.frame(rowData(cds_ERCC)), by = "id") %>%
            select(id, gene_short_name, baseMean:padj) %>%
            na.omit() %>% arrange(padj)
    de_res_all[[time]] = res
}

lapply(de_res_all, dim)
saveRDS(deseq_res, "dds_ERCC.rds")
saveRDS(de_res_all, "de_res_all_ERCC.rds")


# DE analysis
de_res_ERCC = readRDS("de_res_all_ERCC.rds")
de_res_ERCC_conv = readRDS("de_res_all_ERCC_conventional.rds")

conv_counts = lapply(seq_along(de_res_ERCC_conv), function(x) {
    de_res_ERCC_conv[[x]] %>% 
        count(type = log2FoldChange > 0) %>%
        mutate(time = names(de_res_ERCC_conv)[[x]],
               normalization = "Conventional") %>%
        rowwise() %>%
        mutate(type = if(type) {
            "Upregulated"
        } else {
            "Downregulated"
        })
})
conv_counts = bind_rows(conv_counts)

ERCC_counts = lapply(seq_along(de_res_ERCC), function(x) {
    de_res_ERCC[[x]] %>% 
        count(type = log2FoldChange > 0) %>%
        mutate(time = names(de_res_ERCC)[[x]],
               normalization = "ERCC") %>%
        rowwise() %>%
        mutate(type = if(type) {
            "Upregulated"
        } else {
            "Downregulated"
        })
})
ERCC_counts = bind_rows(ERCC_counts)


df_comb = rbind(conv_counts, ERCC_counts)
df_comb$time = as.numeric(as.character(df_comb$time))

options(repr.plot.width=3, repr.plot.height=3)
bar_order = c("Upregulated", "Downregulated")

## Supplementary figure 5c
ggplot(base::transform(df_comb, 
                 type = factor(type, levels = bar_order)) %>%
            filter(time == 24),
    aes(x = type, y = n, fill = type)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +  
  geom_text(aes(label = n, y = n),
        position = position_dodge(0.9),
        size = 4, vjust = -0.2) +
  ylim(c(0,15000)) + 
  facet_wrap(~normalization) +  
  monocle_theme_opts() + 
  theme(legend.position = "none")

ggsave("R_DE_24hour_bar.pdf", device = "pdf", width=3, height=3)


# Compare with single cell data
deseq_res = readRDS("dds_ERCC.rds")
de_res_all = list()
FP_time = rev(unique(cds$FP)[-c(1)])
for (time in FP_time) {
    res = results(deseq_res, contrast=c("FP",time, "0")) %>%
                as.data.frame() %>%
                na.omit()
    res = res[order(res$padj), , drop=F]
    res = as.data.frame(res) %>%
            rownames_to_column("id") %>%
            filter(padj < 0.01) %>%
            left_join(., as.data.frame(rowData(cds_ERCC)), by = "id") %>%
            select(id, gene_short_name, baseMean:padj) %>%
            na.omit() %>% arrange(padj)
    de_res_all[[time]] = res
}

degs_compared = readRDS("degs_hash_LRT_discrete_sc.rds")
degs_hash = readRDS("coeff_table_hash_discrete_sc.rds")

qvalue = 5e-2
num_cells = 20

degs_compared = degs_compared %>% filter(q_value < qvalue)
degs_hash = degs_hash %>% 
    filter(id %in% degs_compared$id & grepl("Cond", term)) %>%
    filter(q_value < qvalue & num_cells_expressed > num_cells) %>%
    arrange(id, estimate) %>%
    arrange(q_value)
dim(degs_hash)


## Supplementary figure 5d
times = c("1", "3", "6", "12", "24")
g = list()

for (time in times) {
     g[[time]] = degs_hash %>%
        filter(grepl(paste0(time, "h"), term)) %>%
        filter(id %in% de_res_ERCC[[time]]$id) %>%
        left_join(., de_res_ERCC[[time]] %>% select(id, log2FoldChange)) %>%
        ggplot(aes(x = estimate, y = log2FoldChange)) + 
            geom_point() + 
            geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + 
            geom_smooth(method = "lm") + 
            geom_abline(slope = 1, intercept = 0) + 
            xlim(c(-6,6)) + ylim(c(-6,6)) + 
            ggtitle(time) + 
            monocle_theme_opts()

}

options(repr.plot.width=7, repr.plot.height=5)
gg = do.call("grid.arrange", c(g, ncol=3))
ggsave("R_bulk_sc_correlation.pdf", plot = gg, device = "pdf", width=7, height=5)

## Log2FC vs single cell effect size estimate correlation
for (time in times) {
    b = degs_hash %>%
        filter(grepl(paste0(time, "h"), term)) %>%
        filter(id %in% de_res_ERCC[[time]]$id) %>%
        left_join(., de_res_ERCC[[time]] %>% select(id, log2FoldChange))
    fit_hash = lm("log2FoldChange ~ estimate", b)
    print(time)
    print(summary(fit_hash))
    print(sqrt(summary(fit_hash)$r.squared))
    
}

