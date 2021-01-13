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
})

wd = getwd()
wd

setwd("../../data/hdaci_qc/")

hashTable_id <- read.table("hashTable-labels-filtered.out_2", header=T) %>%
                      dplyr::select(SampleName, Cell, top_oligo, Axis, top_hash_umi)

colnames(hashTable_id) = c("Sample", "Cell", "Hash", "Axis", "Count")

hashTable_ladder = read.table("hashTable-ladder.out_2", header=T) 
colnames(hashTable_ladder) = c("Sample", "Cell", "Hash", "Axis", "Count")
hashTable_ladder = hashTable_ladder %>%
                        filter(Cell %in% unique(hashTable_id$Cell))

hashTable <- rbind(hashTable_ladder, 
                   hashTable_id) %>%
                arrange(Cell)  %>%
                left_join(., hashTable_id %>% select(Cell, ID = Hash))

head(hashTable)


## define and load features
samples_to_exclude <- as.vector(unlist(read.table("samples.to.exclude_2"))) # from UMI counts/cell

hashes_ref <- read.table("hash-ladder-exp.txt", header=F)
colnames(hashes_ref) <- c("Hash", "Mol")

# filter out based on samples_to_exclude
temp <- hashTable %>% filter(Cell %in% samples_to_exclude)
hashTable <- anti_join(hashTable, temp)
rm(temp)
rm(samples_to_exclude)

dim(hashTable)


## convert ID to condition
IDinfo <- read.table("hashIDSampleSheet.txt", header=T)
indices <- sapply(hashTable$ID, function(x) {
  which(x == IDinfo$ID)
})

temp = IDinfo[indices,-1]
head(temp)

hashTable <- cbind(hashTable, temp)
rownames(hashTable) <- NULL
head(hashTable)

hashTable <- hashTable %>% 
  left_join(hashes_ref, by="Hash")

head(hashTable)

hashTable <- hashTable %>% filter(!(Hash %in% hashes_ref[1:96,1]))
hashTableSub <- hashTable %>% filter(!(Hash %in% hashes_ref[1:96,1]))

# filter by ladder hash count
total_ladder <- aggregate(hashTableSub$Count, 
                          by=list(Cell = hashTableSub$Cell),
                          FUN=sum)
ladder_distinct <- aggregate(hashTableSub$Hash, 
                          by=list(Cell = hashTableSub$Cell),
                          FUN=n_distinct)

hash_filter = total_ladder[total_ladder$x > 100 & ladder_distinct$x > 10,1]
hashTable <- hashTable %>% filter(Cell %in% hash_filter)

hashTable = hashTable %>%
                    left_join(., total_ladder, by = "Cell") %>%
                    dplyr::rename(total_ladder = x)

saveRDS(hashTable, "total_ladder_2.rds")

## add RT labels
labels = sapply(strsplit(as.character(unlist(hashTable$Cell)), "_RT_"), unlist)
hashTable <- hashTable %>% mutate(RT_well = labels[2,], PCR_well = labels[1,])
tail(hashTable)


## glm fit
hashes_for_fit <- hashes_ref[129:144,]
wellids <- as.vector(unique(hashTable$Cell))
suppressWarnings(lms <- lapply(wellids, function(x) {
  sub_table <- filter(hashTable, Cell == x) %>% filter(Hash %in% hashes_for_fit[,1])
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

colnames(getStats) <- wellids
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

head(hashTable)


# add total RNA UMI counts
cellid_order <- as.data.frame(hashTable$Cell)
colnames(cellid_order) <- "Cell"

UMIs_cell <- read.table("UMIs.per.cell.barcode_2", header=F)
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

hashTable_unique = hashTable %>% 
                    distinct(Cell, .keep_all = T)
head(hashTable_unique)
dim(hashTable_unique)

summary(hashTable_unique$Rsq)
summary(hashTable_unique$Total_hash)
hashTable_unique %>% filter(Rsq > 0.7 & Total_hash > 100) %>% dplyr::count(Dim)


## save
write.table(hashTable, file="hashTable_2.txt", quote = F, col.names = T, row.names = F)

dim(hashTable)
hashTable = hashTable %>% filter(Rsq > 0.7 & Total_hash > 100 & Intercept < 0)
dim(hashTable)

write.table(hashTable, file="hashTable_filtered_2.txt", quote = F, col.names = T, row.names = F)

hashTable_unique = read.table("hashTable_filtered_2.txt", header=T) %>%
        distinct(Cell, .keep_all=T)
dim(hashTable_unique)
head(hashTable_unique)

RNA_dup = read.table("dup.rate.per.cell_2")
RNA_dup = RNA_dup %>% filter(V1 %in% hashTable_unique$Cell)

label_dup = read.table("hashDupRate-labels.txt_2")
ladder_dup = read.table("hashDupRate-ladder.txt_2")

label_dup = label_dup %>% filter(V2 %in% hashTable_unique$Cell)
ladder_dup = ladder_dup %>% filter(V2 %in% hashTable_unique$Cell)


hash_dup = data.frame(Cell = label_dup$V2, UMI = label_dup$V3 + ladder_dup$V3, 
                      Reads = label_dup$V4 + ladder_dup$V4)
hash_dup$Duplication_Rate = hash_dup$UMI / hash_dup$Reads

head(hash_dup)

hashTable_unique = hashTable_unique %>% mutate(RNA_duplication = RNA_dup$V4, Hash_duplication = hash_dup$Duplication_Rate)
head(hashTable_unique)

## hash ladder based size factor calculation
## Compute x_intercept and scale total hash and RNA count by their duplication rate
## (to indirectly adjust for sequencing depth)
x_intercept = -hashTable_unique$Intercept / hashTable_unique$Slope

scale_factors = log(hashTable_unique$Total_hash / x_intercept) * 
                hashTable_unique$RNA_duplication

hashTable_unique$Hash_Size_Factor = scale_factors / exp(mean(log(scale_factors))) 

summary(hashTable_unique$Hash_Size_Factor)

write.table(hashTable_unique, file="hashTable_unique_filtered_2.txt", quote = F, col.names = T, row.names = F)

