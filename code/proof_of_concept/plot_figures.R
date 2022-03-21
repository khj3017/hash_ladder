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

# ggplot theme for plotting
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

# load hash tables
setwd("../../data/proof_of_concept/")
hashTable = read.table("hashTable_not_filtered.txt", header=T)
hashTable_unique = read.table("hashTable_unique_filtered.txt", header=T)


# Figure 1

## figure 1b
cellids = unique(hashTable$Cell)
num_hashes = 8
hashes = levels(hashTable$Hash)

hashMatrix <- matrix(0L, nrow = num_hashes, ncol = length(cellids),
                     dimnames = list(hashes, unique(hashTable$Cell)))

for (id in cellids) {
  rowIndices = match(filter(hashTable, Cell == id)$Hash, rownames(hashMatrix))
  colIndices = match(id, colnames(hashMatrix))
  hashMatrix[rowIndices, colIndices] = filter(hashTable, Cell == id)$Count
}

hashes_for_plot <- hashes
hashTable1 <- log10(hashMatrix+1)
hashTable1 <- as.data.frame(hashTable1)
hashTable1 <- hashTable1[match(hashes_for_plot, rownames(hashMatrix)),]

rownames(hashTable1) <- hashes
df <- melt(t(hashTable1))
colnames(df) <- c("Cell", "Hash", "Count")
expected_hash = hashTable %>%
                    distinct(Hash, .keep_all=T) %>%
                    select(Hash, Mol)
df = df %>% left_join(expected_hash, by = "Hash")
        

options(repr.plot.width=3, repr.plot.height=3)
ggplot(df, aes(log10(Mol), Count, fill = as.factor(Mol))) + 
  geom_boxplot(aes(group=Mol), outlier.shape = NA) + 
  scale_fill_brewer(palette = "Reds") + 
  labs(x = "log expected hash count", y = "log observed hash count") +
  monocle_theme_opts() + theme(legend.position = "none")

ggsave("F1_hash_agg.pdf", device = "pdf", width=3, height=3)


## figure 1c
bad_cell = "A06_A06_RT_90"
good_cell = "B02_B02_RT_18"

plotLadderSingleCell <- function(cell) {
    dfplot <- filter(hashTable, Cell == cell & Count > 1)
    g <- ggplot(dfplot, aes(log10(Mol), log10(Count))) + 
        geom_point() +
        labs(x = "log expected hash count", y = "log observed hash count") +
        geom_smooth(method = "lm") +
        monocle_theme_opts()
    print(g)
}

plotLadderSingleCell(bad_cell)
ggsave("F1_bad_cell.pdf", device = "pdf", width=3.5, height=3)
plotLadderSingleCell(good_cell)
ggsave("F1_good_cell.pdf", device = "pdf", width=3.5, height=3)


## figure 1d
b = hashTable_unique
pct = 0.25

b$xintercept = b$Slope / -b$Intercept
b$slope_class = ifelse(hashTable_unique$Slope < quantile(hashTable_unique$Slope, pct), "Bottom 25%", 
                           ifelse(hashTable_unique$Slope > quantile(hashTable_unique$Slope, 1-pct), "Top 25%", "Middle 50%"))
b$intercept_class = ifelse(hashTable_unique$Intercept < quantile(hashTable_unique$Intercept, pct), "Bottom 25%", 
                           ifelse(hashTable_unique$Intercept > quantile(hashTable_unique$Intercept, 1-pct), "Top 25%", "Middle 50%"))
b$xintercept_class = ifelse(b$xintercept < quantile(b$xintercept, pct), "Bottom 25%", 
                           ifelse(b$xintercept > quantile(b$xintercept, 1-pct), "Top 25%", "Middle 50%"))

options(repr.plot.width=4.5, repr.plot.height=3)
ggplot(b, aes(y = log10(Total_RNA), x = as.factor(intercept_class))) + 
    geom_boxplot(aes(fill = as.factor(slope_class))) + 
    labs(x = "", y = "log 10 Total RNA UMI", fill = "Slope class") + 
    monocle_theme_opts()

ggsave("F1_slope_intercept_correlation.pdf", device = "pdf", width=4, height=3)


## figure 1e
size_factor_df = rbind(data.frame(sf = hashTable_unique$Size_Factor, type = "Conventional"),
                       data.frame(sf = hashTable_unique$Hash_Size_Factor, type = "Hash ladder"))

options(repr.plot.width=5, repr.plot.height=3)
ggplot(size_factor_df, aes(x = sf, y = type)) + 
    geom_density_ridges2(aes(height = ..density.., fill = type)) + 
    xlim(c(0,4)) + 
    labs(x = "Size Factor", y = "") + 
    theme_ridges()

ggsave("F1_sf_distribution.pdf", device = "pdf", width=4, height=3)


## figure 1f
cds = readRDS("cds_filtered.rds")
cds2 = cds
cds2$Size_Factor = cds2$Hash_Size_Factor
cds = estimate_size_factors(cds)

compute_cv = function(cds) {
    count_mat = counts(cds) / cds$Size_Factor
    mean_counts = rowMeans(count_mat)
    std_counts = apply(count_mat, 1, sd)
    return(std_counts / mean_counts)
}

compute_mean = function(cds) {
    count_mat = counts(cds) / cds$Size_Factor
    mean_counts = rowMeans(count_mat)
    return(mean_counts)
}

cds = detect_genes(cds, min_expr = 1)
cds = cds[rowData(cds)$num_cells_expressed > 20,]

cds2 = detect_genes(cds2, min_expr = 1)
cds2 = cds2[rowData(cds2)$num_cells_expressed > 20,]

cv_sf = as.data.frame(compute_cv(cds)) %>% rownames_to_column() %>% select(id = 1, cv = 2)
cv_hash = as.data.frame(compute_cv(cds2)) %>% rownames_to_column() %>% select(id = 1, cv = 2)
mean_sf = as.data.frame(compute_mean(cds)) %>% rownames_to_column() %>% select(id = 1, mean = 2)
mean_hash = as.data.frame(compute_mean(cds2)) %>% rownames_to_column() %>% select(id = 1, mean = 2)

df = rbind(left_join(mean_sf, cv_sf, by = "id") %>%
                     mutate(type = "Conventional"),
           left_join(mean_hash, cv_hash, by = "id") %>%
                     mutate(type = "Hash ladder")) %>% na.omit()

options(repr.plot.width=3, repr.plot.height=3)
ggplot(df, aes(x = log10(mean), y = log10(cv), colour = type, alpha = 0.2)) + 
    geom_point() +
    scale_colour_manual(values = c("grey70", "#FF3333")) + 
    labs(x = "", y = "") + 
    monocle_theme_opts() + theme(legend.position = "none")

ggsave("F1_cv_plot.png", device = "png", width=4, height=4)


## figure 1g
df2 = left_join(cv_sf, cv_hash, by = "id") %>%
         dplyr::rename(cv_sf = cv.x, cv_hash = cv.y)

options(repr.plot.width=3, repr.plot.height=3)
ggplot(df2, aes(x = log10(cv_sf), y = log10(cv_hash))) + 
    geom_point(colour = "black") + 
    geom_abline(slope = 1, intercept = 0, lwd = 1.1, colour = "red", linetype = "dashed") + 
    xlim(c(0,1.3)) + ylim(c(0,1.3)) + 
    monocle_theme_opts()

ggsave("F1_cv_comp.png", device = "png", width=3, height=3)


# Supplementary Figure 1
## Supplementary figure 1a
options(repr.plot.width=3.5, repr.plot.height=3)
ggplot(hashTable_unique %>%
        select(Cell, Total_RNA, Total_hash) %>%
        melt(), aes(x = variable, y = log10(value), fill = variable)) + 
  geom_boxplot(aes(outlier.shape = NA)) +
  labs(x = "", y = "log10 observed UMI") +
  monocle_theme_opts()

ggsave("S_RNA_hash_plot.pdf", width = 4, height = 3)


## Supplementary figure 1b
ggplot(hashTable %>% dplyr::count(Cell) %>% dplyr::count(n), 
       aes(x = n, y = nn)) + 
      geom_bar(stat = "identity") +
      labs(x = "Number of unique hash oligos", y = "Number of cells") +
      monocle_theme_opts()

ggsave("S_n_hash_plot.pdf", width = 3.5, height = 3)


## Supplementary figure 1c
valid_cells = unique(hashTable$Cell)
filtered_cells = setdiff(unique(hashTable$Cell), unique(hashTable_unique$Cell))
mat = readRDS("../../data/proof_of_concept/corr_mat.rds")
mat = mat[-match(filtered_cells,valid_cells),-match(filtered_cells,valid_cells)]
corrs = mat[upper.tri(mat, diag = F)]

options(repr.plot.width=3.5, repr.plot.height=3)
ggplot(data.frame(corrs), aes(corrs)) +
    geom_histogram(bins = 200) + 
    xlim(c(0.7,1)) + 
    monocle_theme_opts()

ggsave("S_hash_count_corr.pdf", device = "pdf", width=4, height=4)


## Supplementary figure 1d
hashTable_not_filtered = hashTable %>% distinct(Cell, .keep_all = T)
hashTable_not_filtered$rsq_class = ifelse(hashTable_not_filtered$Rsq < 0.7, "<0.7", 
                ifelse(hashTable_not_filtered$Rsq > 0.7 & hashTable_not_filtered$Rsq < 0.9, "0.7-0.9", 
                       "0.9-1"))

options(repr.plot.width=3.5, repr.plot.height=3)
ggplot(hashTable_not_filtered, aes(y = log10(Total_RNA), x = as.factor(rsq_class))) + 
    geom_boxplot(aes(fill = as.factor(rsq_class))) + 
    labs(x = "", y = "log 10 Total RNA UMI", fill = "Bin") + 
    scale_fill_manual(values = c("#00798c", "#d1495b", "#edae49")) +
    monocle_theme_opts()

ggsave(filename = "S_Rsq_RNA_count.pdf", width = 4, height = 3)


## Supplementary figure 1e
g = list()
g[[1]] = ggplot(hashTable_not_filtered, aes(x = log10(Total_hash))) + 
      geom_histogram(bins = 20, colour = "black", fill = "white") + 
      geom_vline(xintercept = log10(200), linetype = "dashed", colour = "red") + 
      labs(x = "", y = "") +
      monocle_theme_opts()

g[[2]] = ggplot(hashTable_not_filtered, aes(x = Rsq)) + 
      geom_histogram(bins = 20, colour = "black", fill = "white") + 
      geom_vline(xintercept = 0.7, linetype = "dashed", colour = "red") + 
      labs(x = "", y = "") +
      monocle_theme_opts()

options(repr.plot.width=6, repr.plot.height=3)
gg = do.call("grid.arrange", c(g, ncol=2))
ggsave(plot = gg, filename = "S_hash_filter.pdf", width = 6, height = 3)


