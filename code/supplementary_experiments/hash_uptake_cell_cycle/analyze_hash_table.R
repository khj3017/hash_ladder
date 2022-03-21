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

setwd("../../../data/supplementary_experiments/hash_uptake_cell_cycle")

hashTable_id <- read.table("./hashRDS/hashTable-labels-filtered.out", header=T) %>%
                      dplyr::select(SampleName, Cell, top_oligo, Axis, top_hash_umi) %>%
                      filter(top_hash_umi > 30)

colnames(hashTable_id) = c("Sample", "Cell", "Hash", "Axis", "Count")

hashTable_ladder = read.table("./hashRDS/hashTable.out", header=T) 
colnames(hashTable_ladder) = c("Sample", "Cell", "Hash", "Axis", "Count")
hashTable_ladder = hashTable_ladder %>%
                        filter(Cell %in% unique(hashTable_id$Cell))

hashTable <- rbind(hashTable_ladder, 
                   hashTable_id) %>%
                arrange(Cell)  %>%
                left_join(., hashTable_id %>% select(Cell, ID = Hash))

head(hashTable)

samples_to_exclude <- unlist(read.table("samples.to.exclude")) # from UMI counts/cell

hashes_ref <- read.table("hash-ladder-exp.txt", header=F)
colnames(hashes_ref) <- c("Hash", "Mol")
hashes_ref[,2] <- hashes_ref[,2]
hashes <- hashes_ref[,1]
hash_moles <- hashes_ref[,2]

## by RT wells
valid_RTs <- paste("_RT_", seq(0,95,1), "$",  sep="")
collapse_RTs <- paste(valid_RTs, collapse="|")
dim(hashTable)
hashTable = hashTable %>% 
    filter(!(Cell %in% samples_to_exclude) & grepl(collapse_RTs, Cell) & Count > 1) %>%
    left_join(hashes_ref, by="Hash") %>%
    group_by(Cell) %>% add_tally(name = "num_hashes") %>% filter(num_hashes > 4) %>%
    ungroup()


## convert ID to condition
IDinfo <- read.table("hashID.txt", header=F)
colnames(IDinfo) <- c("ID", "Celltype")
indices <- sapply(hashTable$ID, function(x) {
  which(x == IDinfo$ID)
})

temp = IDinfo[indices,]

hashTable <- hashTable %>% mutate(Celltype = temp$Celltype, 
                                  Condition = temp$Condition)
hashTable$ID = NULL
hashTable = hashTable %>% na.omit()
rm(temp)


## add RT labels
labels = sapply(strsplit(as.character(unlist(hashTable$Cell)), "_RT_"), unlist)
hashTable <- hashTable %>% mutate(RT_well = labels[2,], PCR_well = labels[1,])
hashTable <- hashTable %>% 
                mutate(PCR_col = as.numeric(substr(hashTable$PCR_well,2,3)))
hashTable <- hashTable %>% 
                dplyr::rowwise() %>%
                dplyr::mutate(cycle = 
                    if (PCR_col <= 6) {
                        "G1"
                    } else if (PCR_col <= 9) {
                        "S"
                    } else {
                        "G2M"
                    }) %>% as.data.frame()
tail(hashTable)

## Fit glm line to the calibration curve
hashes_for_fit <- hashes_ref
cellids <- as.vector(unique(hashTable$Cell))
suppressWarnings(lms <- lapply(cellids, function(x) {
  sub_table <- filter(hashTable, Cell == x) %>% filter(Hash %in% hashes_for_fit$Hash)
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
  c(y$coefficients[,1][2], y$coefficients[,1][1], 1 - (y$deviance/y$null.deviance))
})

colnames(getStats) <- cellids
getStats <- t(getStats)
num_hashes_cell <- data.frame(table(hashTable$Cell)) %>% filter(Freq > 0)
matched_indices <- match(rownames(getStats),as.character(num_hashes_cell$Var1))
num_hashes_cell <- num_hashes_cell[matched_indices,]

getStatsSlope <- dplyr::slice(as.data.frame(getStats[,1]), 
                              base::rep(1:n(), num_hashes_cell[,2]))
getStatsIntercept <- dplyr::slice(as.data.frame(getStats[,2]), 
                                  base::rep(1:n(), num_hashes_cell[,2]))
getStatsRsq <- dplyr::slice(as.data.frame(getStats[,3]), 
                                  base::rep(1:n(), num_hashes_cell[,2]))

hashTable <- hashTable %>% filter(Cell %in% rownames(getStats))
hashTable <- hashTable %>% mutate(Slope = getStatsSlope[,1], 
                                  Intercept = getStatsIntercept[,1],
                                  Rsq = getStatsRsq[,1])

## add total RNA UMI counts
cellid_order <- as.data.frame(hashTable$Cell)
colnames(cellid_order) <- "Cell"

UMIs_cell <- read.table("UMIs.per.cell.barcode", header=F)
UMIs_cell <- UMIs_cell[,2:3]
colnames(UMIs_cell) <- c("Cell", "Total_RNA")

UMIs_cell <- UMIs_cell %>% filter(Cell %in% hashTable$Cell)

hashTable <- hashTable %>% left_join(UMIs_cell, by="Cell") %>% na.omit()

## add total hash UMI counts
cellid_order <- as.data.frame(hashTable$Cell)
colnames(cellid_order) <- "Cell"
UMIs_hash <- aggregate(hashTable$Count, 
                       by=list(Hash = hashTable$Cell),
                       FUN=sum)
colnames(UMIs_hash) = c("Cell", "Total_hash")


hashTable <- hashTable %>% left_join(UMIs_hash, by="Cell") %>% na.omit()

## filter hashTable
hashTable = hashTable %>% filter(Total_hash > 100 & Rsq > 0.6)
dim(hashTable %>% distinct(Cell, .keep_all = T))
head(hashTable)

## add duplication rate
hashTable_unique = hashTable %>% distinct(Cell, .keep_all = T)
RNA_dup = read.table("dup.rate.per.cell")
RNA_dup = RNA_dup %>% filter(V1 %in% hashTable_unique$Cell)

ladder_dup = read.table("hashDupRate.txt")
ladder_dup = ladder_dup %>% filter(V2 %in% hashTable_unique$Cell)

hashTable_unique = hashTable_unique %>%
                    mutate(RNA_duplication = RNA_dup$V4, Hash_duplication = ladder_dup$V5)
head(hashTable_unique)

hashTable_unique_filtered = hashTable_unique %>% filter(Rsq > 0.6)

write.table(hashTable, file="hashTable_filtered.txt", quote = F, col.names = T, row.names = F)
write.table(hashTable_unique_filtered, file="hashTable_unique_filtered.txt", quote = F, col.names = T, row.names = F)


## Supplementary figure 16c
options(repr.plot.width=3.5, repr.plot.height=3)
ggplot(hashTable_unique_filtered, aes(y = log10(Total_hash), x = cycle)) + 
    geom_boxplot(aes(fill = cycle)) + 
    facet_wrap("Celltype") + 
    labs(x = "", y = "log Total Hash") +
    monocle3:::monocle_theme_opts() + 
    theme(legend.position = "none")

ggsave("hash_count_by_cell_cycle.pdf", width = 3.5, height = 3)

## Welch Two Sample t-test
A549_df = hashTable_unique_filtered %>% filter(Celltype == "A549")
HEK293T_df = hashTable_unique_filtered %>% filter(Celltype == "HEK293T")

t.test(A549_df[A549_df$cycle == "G1",]$Total_hash, A549_df[A549_df$cycle == "S",]$Total_hash)
t.test(A549_df[A549_df$cycle == "S",]$Total_hash, A549_df[A549_df$cycle == "G2M",]$Total_hash)
t.test(A549_df[A549_df$cycle == "G1",]$Total_hash, A549_df[A549_df$cycle == "G2M",]$Total_hash)

t.test(HEK293T_df[HEK293T_df$cycle == "G1",]$Total_hash, HEK293T_df[HEK293T_df$cycle == "S",]$Total_hash)
t.test(HEK293T_df[HEK293T_df$cycle == "S",]$Total_hash, HEK293T_df[HEK293T_df$cycle == "G2M",]$Total_hash)
t.test(HEK293T_df[HEK293T_df$cycle == "G1",]$Total_hash, HEK293T_df[HEK293T_df$cycle == "G2M",]$Total_hash)

