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
    library(monocle3)
})

wd = getwd()
wd

## Put processed data into the data directory
setwd("../../data/proof_of_concept/")

# load hash table
hashTable <- read.table("hashTable.out", header=F,
                       col.names = c("Sample", "Cell", "Hash", "Dim", "Count"),
                       colClasses = c("character", "character", "factor", "integer", "integer"))[,-c(1,4)]

# filter out cells based on RNA UMI/cell
samples_to_exclude <- unlist(read.table("samples.to.exclude")) # from UMI counts/cell

### Concentration for E5 and F5 flipped when making the ladder
hashes <- c("A5", "B5", "C5", "D5", "E5", "F5", "G5", "H5")
hash_moles <- c(0.2, 0.4, 0.8, 1.6, 6.4, 3.2, 12.8, 25.6)*1000 # femtomoles

N = 2e6
avogadro = 6.022e23
theor_hash_count <- hash_moles*1e-15*avogadro/N

hashes_ref <- data.frame(Hash= hashes, Mol = theor_hash_count)

# filter for cells with valid RT ID
valid_RTs <- paste("_RT_", seq(0,95,1), "$",  sep="")
collapse_RTs <- paste(valid_RTs, collapse="|")

hashTable = hashTable %>% 
    filter(!(Cell %in% samples_to_exclude) & grepl(collapse_RTs, Cell) & Count > 1) %>%
    left_join(hashes_ref, by="Hash") %>%
    group_by(Cell) %>% add_tally(name = "num_hashes") %>% filter(num_hashes > 3) %>%
    ungroup()

labels = sapply(strsplit(as.character(unlist(hashTable$Cell)), "_RT_"), unlist)
hashTable <- hashTable %>% mutate(RT_well = labels[2,], PCR_well = labels[1,])

## fit glm line to the calibration curve
## Keep in mind that hashes at lower abundance might drop out and
## you might need to omit them for fitting

hashes_for_fit <- hashes_ref
cellids <- as.vector(unique(hashTable$Cell))
suppressWarnings(lms <- lapply(cellids, function(x) {
  sub_table <- hashTable %>% filter(Cell == x & Hash %in% hashes_for_fit$Hash)
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

#filter
UMIs_cell <- UMIs_cell %>% filter(Cell %in% hashTable$Cell)

hashTable <- hashTable %>% left_join(UMIs_cell, by="Cell") %>% na.omit()

# add total hash UMI counts
cellid_order <- as.data.frame(hashTable$Cell)
colnames(cellid_order) <- "Cell"
UMIs_hash <- aggregate(hashTable$Count, 
                       by=list(Hash = hashTable$Cell),
                       FUN=sum)
colnames(UMIs_hash) = c("Cell", "Total_hash")


hashTable <- hashTable %>% left_join(UMIs_hash, by="Cell") %>% na.omit()

# filter out cells and save
write.table(hashTable, file="hashTable_not_filtered.txt", quote = F, col.names = T, row.names = F)
hashTable = hashTable %>% filter(Rsq > 0.7 & Total_hash > 100)
hashTable_unique = hashTable %>% distinct(Cell, .keep_all = T)

# add duplication rates
RNA_dup = read.table("dup.rate.per.cell")
RNA_dup = RNA_dup %>% filter(V1 %in% hashTable_unique$Cell)

ladder_dup = read.table("hashDupRate.txt")
ladder_dup = ladder_dup %>% filter(V2 %in% hashTable_unique$Cell)

hashTable_unique = hashTable_unique %>%
                    mutate(RNA_duplication = RNA_dup$V4, Hash_duplication = ladder_dup$V5)

x_intercept = -hashTable_unique$Intercept / hashTable_unique$Slope

scale_factors = log(hashTable_unique$Total_hash) / x_intercept *
                hashTable_unique$RNA_duplication

hashTable_unique$Size_Factor = hashTable_unique$Total_RNA / exp(mean(log(hashTable_unique$Total_RNA)))
hashTable_unique$Hash_Size_Factor = scale_factors / exp(mean(log(scale_factors))) 


# save tables
write.table(hashTable, file="hashTable_filtered.txt", quote = F, col.names = T, row.names = F)
write.table(hashTable_unique, file="hashTable_unique_filtered.txt", quote = F, col.names = T, row.names = F)

### load.cds from sci-RNA-seq output
load.cds = function(mat.path, gene.annotation.path, cell.annotation.path, 
                    hashTable.path) {
  df = read.table(
    mat.path,
    col.names = c("gene.idx", "cell.idx", "count"),
    colClasses = c("integer", "integer", "integer"))
  
  gene.annotations = read.table(
    gene.annotation.path,
    col.names = c("id", "gene_short_name"),
    colClasses = c("character", "character"))
  
  cell.annotations = read.table(
    cell.annotation.path,
    col.names = c("cell", "sample"),
    colClasses = c("character", "factor"))
  
  hashTable <- read.table(hashTable.path, header=T, stringsAsFactors = F)

  
  rownames(gene.annotations) = gene.annotations$id
  rownames(cell.annotations) = cell.annotations$cell
    
  
  # add a dummy cell to ensure that all genes are included in the matrix
  # even if a gene isn't expressed in any cell
  df = rbind(df, data.frame(
    gene.idx = c(1, nrow(gene.annotations)),
    cell.idx = rep(nrow(cell.annotations)+1, 2),
    count = c(1, 1)))
  
  mat = sparseMatrix(i = df$gene.idx, j = df$cell.idx, x = df$count)
  
  mat = mat[, 1:(ncol(mat)-1)]
  
  rownames(mat) = gene.annotations$id
  colnames(mat) = cell.annotations$cell
  
  matchNames <- match(hashTable$Cell, rownames(cell.annotations))
  mat <- mat[,matchNames]

  
  cell.annotations <- cell.annotations[matchNames,]
  cell.annotations <- cbind(cell.annotations, hashTable[,c(4:ncol(hashTable))])
  
    
  cell.annotations$RT_well <- as.factor(cell.annotations$RT_well)
    
  cds = new_cell_data_set(mat, 
                          cell_metadata = cell.annotations, 
                          gene_metadata = gene.annotations)
  
    
  return(cds)
}

## load cds
setwd(wd)
cds = load.cds('UMI.count.matrix', '../refdata/gene.annotations', 
               'cell.annotations', 'hashTable_unique_filtered.txt')


## filter out mito.genes and subtract mito gene counts
mito_genes <- grep(pattern = "^MT-", x = rowData(cds)$gene_short_name, value = TRUE, ignore.case = TRUE)
ribo_genes <- grep(pattern = "^Rp[sl][[:digit:]]", x = rowData(cds)$gene_short_name, value = TRUE)

mito_gene_counts = colSums(counts(cds)[rowData(cds)$gene_short_name %in% mito_genes,])
mito_gene_id = rowData(cds)[rowData(cds)$gene_short_name %in% mito_genes,]$id

cds = estimate_size_factors(cds)

saveRDS(cds, "cds_filtered.rds")
cds

## compute correlations between cells' hash counts
## takes a few hours to run
fit_lines = function(cell1, cell2, hashTable) {
    cell_1_table = hashTable %>%
                    filter(Cell == cell1) %>%
                    select(Cell, Hash, Count)
    cell_2_table = hashTable %>%
                    filter(Cell == cell2) %>%
                    select(Cell, Hash, Count)

    df = inner_join(cell_1_table, cell_2_table, by = "Hash")

    return(summary(lm("log(Count.x) ~ log(Count.y)", data = df))$adj.r.squared)

}

valid_cells = unique(hashTable$Cell)
L = length(valid_cells)
mat = matrix(0L, nrow = L, ncol = L)

for (i in 1:(L-1)) {
    cell_1 = valid_cells[i]
    print(as.character(cell_1))
    flush.console()
    for (j in (i+1):L) {
        cell_2 = valid_cells[j]
        mat[i,j] = fit_lines(cell_1, cell_2, hashTable)
        mat[j,i] = mat[i,j]
    }
}
saveRDS(mat, "corr_mat.rds")
