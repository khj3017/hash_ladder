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
    library(readxl)
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
setwd("../../data/hdaci_dex")
cds_dex = readRDS("cds_dex_sf.rds")
cds_dex_hash = readRDS("cds_dex_hash.rds")

hash_id_df_1 = read.table("../hdaci_qc/hash_id_df.txt_1", 
                          stringsAsFactors = F, header=T)
hash_id_df_2 = read.table("../hdaci_qc/hash_id_df.txt_2",
                          stringsAsFactors = F, header=T)

hash_id_df = rbind(hash_id_df_1, hash_id_df_2) %>% filter(!is.na(top_oligo))

# convert ID to condition
IDinfo <- read.table("../hdaci_qc/hashIDSampleSheet.txt", header=T)
indices <- sapply(hash_id_df$top_oligo, function(x) {
  which(x == IDinfo$ID)
})

temp = IDinfo[indices,-1]

hash_id_df <- cbind(hash_id_df, temp) %>%
                    filter(Dim == 2)
rownames(hash_id_df) <- NULL
head(hash_id_df)


## figure S8B
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
ggsave(gg, filename = "S_dex_hashID_filter.pdf", width = 7, height = 3)


## figure S8C
hashTable1  = read.table("../hdaci_dex/hashTable_1.txt", header=T)
hashTable2  = read.table("../hdaci_dex/hashTable_2.txt", header=T)

hashTable = rbind(hashTable1, hashTable2) %>%
                filter(Dim == 2) %>%
                distinct(Cell, .keep_all = T) 

total_ladder1 = readRDS("../hdaci_qc/total_ladder_1.rds")
total_ladder2 = readRDS("../hdaci_qc/total_ladder_2.rds")
total_ladder = rbind(total_ladder1, total_ladder2) %>%
                filter(Dim == 2) %>%
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
ggsave(gg, filename = "S_dex_hashLadder_filter.pdf", width = 7, height = 3)

metadata = as.data.frame(colData(cds_dex))
melted_metadata = metadata %>% select(Cell, Total_hash, Total_RNA) %>% melt()

metadata %>%
    pull(Total_hash) %>%
    summary()


## figure S8D
options(repr.plot.width=3.5, repr.plot.height=3)
ggplot(melted_metadata, aes(x = variable, y = log10(value))) + 
    geom_boxplot(aes(fill = variable)) +
    monocle_theme_opts()
ggsave(filename = "S_dex_hash_RNA_count.pdf", width = 3.5, height = 3)

