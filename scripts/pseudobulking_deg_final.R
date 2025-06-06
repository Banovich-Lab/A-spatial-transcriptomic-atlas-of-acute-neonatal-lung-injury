#### LOAD PACKAGES AND SET ENVIRONMENT ----
library(SeuratObject)
library(tidyverse)
library(Seurat, lib.loc = "/usr/local/lib/R/site-library")
library(edgeR)
library(limma)

# Set seed & working directory
set.seed(0712)
work_dir <- setwd("/scratch/smallapragada/run2")
options(scipen = 0)

# Load data in
run2_final <- readRDS("/scratch/smallapragada/run2/final_complete_cniches_tniches_allsamples_2025_05_26.rds")

#### FILTERING THRESHOLD ----

### Keeping genes that had > 20% expression in all cells

## Pulling raw counts out
cell_expr <- run2_final@assays$RNA$counts
dim(cell_expr)

## Pulling metadata out
cell_object_meta <- run2_final@meta.data

## Appending CTs separately 
sample_celltype_levels <-  paste0(as.character(cell_object_meta$ID),"_",
                                  as.character(cell_object_meta$CT_final))

cell_RNA_expr_perCT <- t(sapply(by(t(cell_expr),sample_celltype_levels,colSums),identity))

dim(cell_RNA_expr_perCT)

## Calculating the proportion of gene expression per cell type
celltype_levels <- cell_object_meta$CT_final
cell_expr_count_per1K <- apply(cell_expr,2,function(x){
  res = x/sum(x)
  res * 1000
})

## Filtering out the proportion of cells with < 3 counts
prop_thresh <- t(sapply(by(t(cell_expr_count_per1K),celltype_levels, 
                           function(x){
                             colSums(x >= 3)/nrow(x) #3
                           }),identity))

test_gene_prop_thresh <- data.frame(prop_thresh,
                                    check.names = FALSE) %>% 
  mutate(CT_final = rownames(prop_thresh)) %>%
  tidyr::pivot_longer(1:343,
                      names_to = "gene",
                      values_to = "prop_cell_expressing_min3")

## Keeping genes and CTs that have expression in at least 20% of that CT
test_gene_prop_thresh_filter0.2 <- test_gene_prop_thresh %>% 
  filter(prop_cell_expressing_min3 >= 0.2) 

#### DEG IN IMMATURE AT2 VS AT2 ----

## Removing immature CT for comparison 
run2_at2 <- subset(run2_final, subset = CT_final %in% c("Immature AT2", "AT2"))

## Pull normalized expression data
expr_matrix <- as.data.frame(as.matrix(run2_at2@assays$RNA@data))
expr_matrix_t <- t(expr_matrix)
expr_matrix_t <- as.data.frame(expr_matrix_t)

## Add sample information to the data
expr_matrix_t$Sample_combined <- run2_at2@meta.data$Sample_combined
expr_matrix_t$CT_final <- run2_at2@meta.data$CT_final

## Summarize by Sample_combined to calculate the mean expression for each gene in each sample
mean_expression_per_sample <- expr_matrix_t %>%
  group_by(Sample_combined, CT_final) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  as.data.frame()

## Assigning rownames to dataframe containing CT and Sample_combined
samp_ct <- paste(mean_expression_per_sample$Sample_combined, mean_expression_per_sample$CT_final, sep = "_") 
rownames(mean_expression_per_sample) <- samp_ct

mean_expression_per_sample$Sample_combined <- NULL
mean_expression_per_sample$CT_final <- NULL

avg_expression_df <- as.data.frame(t(mean_expression_per_sample))

### Limma-Voom

## Grouping samples for design matrix
immature <- avg_expression_df %>%
  select(contains("Immature")) %>%
  colnames()

nonimmature <- avg_expression_df %>%
  select(!contains("Immature")) %>%
  colnames()

## Creating design matrix for transformation
condition <- c(rep("Immature", length(immature)), rep("Nonimmature", length(nonimmature)))
design <- model.matrix(~0 + condition, data = avg_expression_df)

## Creating contrast matrix for linear model
contrast_matrix <- makeContrasts(
  ImmaturevsNonimmature = conditionImmature-conditionNonimmature, 
  NonimmaturevsImmature = conditionNonimmature-conditionImmature,
  levels = colnames(design))

## Identify columns with "immature"
immature_cols <- grep("Immature", colnames(avg_expression_df), value = TRUE)

## Identify other columns
other_cols <- setdiff(colnames(avg_expression_df), immature_cols)

## Reorder dataframe with "immature" columns first
avg_expression_df_reordered <- avg_expression_df[, c(immature_cols, other_cols)]

## Filtering low-exp genes
genes_to_keep <- unique(test_gene_prop_thresh_filter0.2$gene)
avg_expression_df_reordered_keep <- avg_expression_df_reordered[genes_to_keep,]

## Loading DGEList object
d0 <- DGEList(avg_expression_df_reordered_keep)

## Voom transformation
transformed <- voom(d0, design, plot = T)

## Fitting to a linear model
vfit <- lmFit(transformed, design)
model <- contrasts.fit(vfit, contrasts=contrast_matrix)
efit <- eBayes(model)

## Adding the CT as a column
CT_sig_genes <- topTable(efit, coef=2, number = 343)

## Convert row names to a column
CT_sig_genes <- cbind(RowNames = rownames(CT_sig_genes), CT_sig_genes)
rownames(CT_sig_genes) <- NULL

## Save results
write_csv(CT_sig_genes, "/scratch/smallapragada/CT_sig_genes.csv")
