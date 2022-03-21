suppressPackageStartupMessages({
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
    library(GGally)
})

setwd("../../../data/supplementary_experiments/bulk_flavopiridol")

# Bulk flavopiridol experiment with ERCC spike-ins
ERCC_table = read.table("ERCC.UMI.counts", stringsAsFactors = F)
RT_sample_sheet = read.table("RTSampleSheet.txt", header=T, stringsAsFactors = F)
colnames(ERCC_table) = c("Sample", "ERCC", "Count")

labels = sapply(strsplit(as.character(unlist(ERCC_table$Sample)), "_RT_"), unlist)
ERCC_table = ERCC_table %>% 
                mutate(RT = as.double(labels[2,]), PCR = labels[1,]) %>%
                left_join(., RT_sample_sheet)

## add total RNA UMI counts
cellid_order <- as.data.frame(ERCC_table$Sample)
colnames(cellid_order) <- "Sample"

UMIs_cell <- read.table("UMIs.per.cell.barcode", header=F)
UMIs_cell <- UMIs_cell[,2:3]
colnames(UMIs_cell) <- c("Sample", "Total_RNA")

## filter
UMIs_cell <- UMIs_cell %>% filter(Sample %in% ERCC_table$Sample)

ERCC_table <- ERCC_table %>% left_join(UMIs_cell, by="Sample") %>% na.omit()

## add total hash UMI counts
UMIs_cell <- read.table("UMIs.per.cell.barcode.ERCC", header=F)
UMIs_cell <- UMIs_cell[,2:3]
colnames(UMIs_cell) = c("Sample", "Total_ERCC")

ERCC_table <- ERCC_table %>% left_join(UMIs_cell, by="Sample") %>% na.omit()


RNA_dup = read.table("dup.rate.per.cell", stringsAsFactors = F)
colnames(RNA_dup) = c("Sample", "UMI", "Read", "RNA_duplication")
RNA_dup = RNA_dup %>% filter(Sample %in% ERCC_table$Sample)
RNA_dup = RNA_dup %>% rowwise() %>% mutate(RT = as.double(strsplit(Sample, "_")[[1]][4])) %>%
            group_by(RT) %>% mutate(UMI = sum(UMI), Read = sum(Read)) %>%
            ungroup() %>% distinct(RT, .keep_all = T) %>%
            mutate(RNA_duplication = 1 - UMI / Read)

ladder_dup = read.table("dup.rate.per.cell.ERCC", stringsAsFactors = F)
colnames(ladder_dup) = c("Sample", "UMI", "Read", "ERCC_duplication")
ladder_dup = ladder_dup %>% filter(Sample %in% ERCC_table$Sample)
ladder_dup = ladder_dup %>% rowwise() %>% mutate(RT = as.double(strsplit(Sample, "_")[[1]][4])) %>%
            group_by(RT) %>% mutate(UMI = sum(UMI), Read = sum(Read)) %>%
            ungroup() %>% distinct(RT, .keep_all = T) %>%
            mutate(ERCC_duplication = 1 - UMI / Read)

ERCC_table = ERCC_table %>%
                left_join(., RNA_dup %>% select(RT, RNA_duplication)) %>%
                left_join(., ladder_dup %>% select(RT, ERCC_duplication))


ERCC_table_sub = ERCC_table %>%
                group_by(RT, ERCC) %>%
                mutate(Count = sum(Count)) %>%
                as.data.frame() %>%
                distinct(Sample, RT, .keep_all=T) %>%
                group_by(RT) %>%
                mutate(Total_RNA = sum(Total_RNA),
                       Total_ERCC = sum(Total_ERCC)) %>%
                as.data.frame() %>%
                distinct(Sample, RT, ERCC, .keep_all=T) %>%
                arrange(RT)

ERCC_table = ERCC_table %>% 
                select(-c(Total_RNA, Total_ERCC)) %>%
                left_join(., ERCC_table_sub %>% select(Sample, RT, 
                            Total_RNA, Total_ERCC), by = c("Sample", "RT")) %>%
                group_by(RT, ERCC) %>%
                mutate(Count = sum(Count)) %>%
                as.data.frame()

head(ERCC_table)

## Add ERCC concentration info
ERCC_concentrations = read.table("ERCC.concetrations", stringsAsFactors = F)
colnames(ERCC_concentrations) = c("V1", "ERCC", "V3", "Attomoles_uL", "V5", "V6", "V7")
ERCC_concentrations = ERCC_concentrations %>% select(ERCC, Attomoles_uL)
ERCC_concentrations = ERCC_concentrations %>%
                        mutate(ERCC = gsub("-", "_", ERCC),
                               Mol = Attomoles_uL / 100 * 3 * 6.022e5)
head(ERCC_concentrations)

ERCC_table = ERCC_table %>%
                left_join(., ERCC_concentrations)
head(ERCC_table)


## fit glm line to the calibration curve

hashes_for_fit <- ERCC_concentrations
cellids <- as.vector(unique(ERCC_table$RT))
suppressWarnings(lms <- lapply(cellids, function(x) {
  sub_table <- filter(ERCC_table, RT == x) %>% filter(ERCC %in% hashes_for_fit$ERCC)
  sub_table <- sub_table %>% filter(Count > 1)
  return(
    tryCatch({
      glm.nb(Count ~ log(Mol), data=sub_table, control=glm.control(maxit=100))
    }, error = function(e) {
      print(x)
      glm(Count ~ log(Mol), data=sub_table, control=glm.control(maxit=100),family=poisson)
    })
  )
}))


## extract coefficients
getStats <- sapply(lms, function(x) {
  y <- summary(x)
  #c(y$coefficients[,1][1], y$r.squared)
  c(y$coefficients[,1][2], y$coefficients[,1][1], 1 - (y$deviance/y$null.deviance))
})

colnames(getStats) <- cellids
getStats <- t(getStats)
num_hashes_cell <- data.frame(table(ERCC_table$RT)) %>% filter(Freq > 0)
matched_indices <- match(rownames(getStats),as.character(num_hashes_cell$Var1))
num_hashes_cell <- num_hashes_cell[matched_indices,]

getStatsSlope <- dplyr::slice(as.data.frame(getStats[,1]), 
                              base::rep(1:n(), num_hashes_cell[,2]))
getStatsIntercept <- dplyr::slice(as.data.frame(getStats[,2]), 
                                  base::rep(1:n(), num_hashes_cell[,2]))
getStatsRsq <- dplyr::slice(as.data.frame(getStats[,3]), 
                                  base::rep(1:n(), num_hashes_cell[,2]))

ERCC_table <- ERCC_table %>% filter(RT %in% rownames(getStats))
ERCC_table <- ERCC_table %>% mutate(Slope = getStatsSlope[,1], 
                                  Intercept = getStatsIntercept[,1],
                                  Rsq = getStatsRsq[,1])


## Supplementary figure 5b
## randomly pick a cell and plot calibration curve
options(repr.plot.width=3, repr.plot.height=3)
hashes_for_fit <- ERCC_concentrations
cellids = unique(ERCC_table$RT)

dfplot <- filter(ERCC_table, RT == 0 & Count > 1 & ERCC %in% hashes_for_fit$ERCC)
ggplot(dfplot, aes(log10(Mol), log10(Count))) + 
    geom_point() + #ylim(c(0,10)) + 
    labs(x = "log expected ERCC count", y = "log observed ERCC count") +
    geom_smooth(method = "lm") +
    monocle3:::monocle_theme_opts()

ggsave("R_ERCC_calibration_line.pdf", device = "pdf", width=3, height=3)


write.table(ERCC_table, file="ERCC_table_corrected_nb_new.txt", quote = F, col.names = T, row.names = F)

ERCC_table_unique = ERCC_table %>% 
                    distinct(RT, .keep_all = T)

write.table(ERCC_table_unique, file="ERCC_table_unique_corrected_nb_new.txt", quote = F, col.names = T, row.names = F)
