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
})

wd = getwd()
wd

setwd("../../data/flavopiridol")

# load hash ID label and ladder tables
hashTable_id <- read.table("hashTable-labels-filtered.out", header=T) %>%
                      dplyr::select(SampleName, Cell, top_oligo, Axis, top_hash_umi) %>%
                      filter(top_hash_umi > 30)

colnames(hashTable_id) = c("Sample", "Cell", "Hash", "Axis", "Count")

hashTable_ladder = read.table("hashTable.out", header=T) 
colnames(hashTable_ladder) = c("Sample", "Cell", "Hash", "Axis", "Count")
hashTable_ladder = hashTable_ladder %>%
                        filter(Cell %in% unique(hashTable_id$Cell))

hashTable <- rbind(hashTable_ladder, 
                   hashTable_id) %>%
                arrange(Cell)  %>%
                left_join(., hashTable_id %>% select(Cell, ID = Hash))

head(hashTable)


# filter out cells based on RNA UMI/cell
samples_to_exclude <- unlist(read.table("samples.to.exclude")) # from UMI counts/cell

hashes_ref <- read.table("hash-ladder-exp.txt", header=F)
colnames(hashes_ref) <- c("Hash", "Mol", "GC")
hashes_ref = hashes_ref %>% select(-GC)

valid_RTs <- paste("_RT_", seq(0,95,1), "$",  sep="")
collapse_RTs <- paste(valid_RTs, collapse="|")

hashTable = hashTable %>% 
    filter(!(Cell %in% samples_to_exclude) & grepl(collapse_RTs, Cell) & Count > 1) %>%
    left_join(hashes_ref, by="Hash") %>%
    group_by(Cell) %>% add_tally(name = "num_hashes") %>% filter(num_hashes > 3) %>%
    ungroup()

dim(hashTable)

## convert ID to condition
IDinfo <- read.table("hashID.txt", header=F)
colnames(IDinfo) <- c("ID", "Celltype", "Condition")
indices <- sapply(hashTable$ID, function(x) {
  which(x == IDinfo$ID)
})

temp = IDinfo[indices,]

hashTable <- hashTable %>% mutate(Celltype = temp$Celltype, 
                                  Condition = temp$Condition)
hashTable$ID = NULL
rm(temp)

## add RT labels
labels = sapply(strsplit(as.character(unlist(hashTable$Cell)), "_RT_"), unlist)
hashTable <- hashTable %>% mutate(RT_well = labels[2,], PCR_well = labels[1,])
tail(hashTable)

## fit glm line to the calibration curve
## Keep in mind that hashes at lower abundance might drop out and
## you might need to omit them for fitting

hashes_for_fit <- hashes_ref[32:48,]
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

# add total RNA UMI counts
cellid_order <- as.data.frame(hashTable$Cell)
colnames(cellid_order) <- "Cell"

UMIs_cell <- read.table("UMIs.per.cell.barcode", header=F)
UMIs_cell <- UMIs_cell[,2:3]
colnames(UMIs_cell) <- c("Cell", "Total_RNA")

UMIs_cell <- UMIs_cell %>% filter(Cell %in% hashTable$Cell)

hashTable <- hashTable %>% left_join(UMIs_cell, by="Cell") %>% na.omit()
head(hashTable)

# add total hash UMI counts
cellid_order <- as.data.frame(hashTable$Cell)
colnames(cellid_order) <- "Cell"
UMIs_hash <- aggregate(hashTable$Count, 
                       by=list(Hash = hashTable$Cell),
                       FUN=sum)
colnames(UMIs_hash) = c("Cell", "Total_hash")


hashTable <- hashTable %>% left_join(UMIs_hash, by="Cell") %>% na.omit()
head(hashTable)

# filter out cells and save
hashTable = hashTable %>% filter(Celltype == "HEK293T")
write.table(hashTable, file=file.path(wd, "hashTable.txt"), quote = F, col.names = T, row.names = F)

# filter out cells and save
hashTable = hashTable %>% filter(Total_hash > 100 & Celltype == "HEK293T" & Rsq > 0.7)
dim(hashTable %>% distinct(Cell, .keep_all = T))
head(hashTable)

write.table(hashTable, file="hashTable_filtered.txt", quote = F, col.names = T, row.names = F)

hashTable_unique = hashTable %>% 
                    distinct(Cell, .keep_all = T)
head(hashTable_unique)
dim(hashTable_unique)

# add duplication rates
RNA_dup = read.table("dup.rate.per.cell")
RNA_dup = RNA_dup %>% filter(V1 %in% hashTable_unique$Cell)

ladder_dup = read.table("hashDupRate.txt")
ladder_dup = ladder_dup %>% filter(V2 %in% hashTable_unique$Cell)

hashTable_unique = hashTable_unique %>%
                    mutate(RNA_duplication = RNA_dup$V4, Hash_duplication = ladder_dup$V5)
head(hashTable_unique)

## hash ladder based size factor calculation
## Compute x_intercept and scale total hash and RNA count by their duplication rate
## (to indirectly adjust for sequencing depth)
x_intercept = -hashTable_unique$Intercept / hashTable_unique$Slope

scale_factors = log(hashTable_unique$Total_hash / x_intercept) * 
                hashTable_unique$RNA_duplication

hashTable_unique$Hash_Size_Factor = scale_factors / exp(mean(log(scale_factors))) 

summary(hashTable_unique$Hash_Size_Factor)

write.table(hashTable_unique, file="hashTable_unique_filtered.txt", quote = F, col.names = T, row.names = F)


