#### LOAD PACKAGES AND SET ENVIRONMENT ----
library(SeuratObject)
library(tidyverse)
library(gplots)
library(rlist)
library(Seurat, lib.loc = "/usr/local/lib/R/site-library")

# Set seed & working directory
set.seed(0712)
work_dir <- setwd("/scratch/smallapragada/run2")

# Read in object & get metadata
run2_final <- readRDS("/scratch/smallapragada/run2/final_complete_cniches_tniches_allsamples_2025_03_26.rds")

#### SUBSETTING OBJECT AND MAKING GENE FILTERING THRESHOLD ----

## Subset
sobj_sub <- run2_final[,run2_final$Sample_combined %in% 
                         c("PDL003", "PDL006", "PDL005",
                           "PDL004", "PDL007", "PDL008", "PDL009",
                           "PDL010", "PDL011", "PDL012", "PDL013")]
## Threshold
min3Per1k_prop20 <- list()

for(ct in unique(sobj_sub$CT_final)){
  message(ct)
  ct_gene_cell <- sobj_sub@assays$RNA@counts[,sobj_sub$CT_final==ct]
  ct_gene_cell_per1k <-  apply(ct_gene_cell, 2, function(x){
    res = x/sum(x)
    res*1000
  })
  keep_per_ct <- rowSums(ct_gene_cell_per1k >=3)/ncol(ct_gene_cell_per1k) >= 0.2
  min3Per1k_prop20[[ct]] <- names(keep_per_ct)[keep_per_ct]
}

saveRDS(min3Per1k_prop20,"/scratch/smallapragada/run2/lm_regression_ga_ls_ruqian_2025_04_11/min3Per1k_prop20_2025_04_17.rds")

## Setting objects up

## Aggregate gene count per cell type per sample
sobh_sub_aggr_matrix <-
  AggregateExpression(sobj_sub, 
                      assays = "RNA",
                      group.by = c("Sample_combined",
                                   "CT_final"))$RNA

sobh_sub_aggr_sce <- 
  SingleCellExperiment(assays = list(counts = sobh_sub_aggr_matrix))
sobh_sub_aggr_sce 

## Pulling sample names
sobh_sub_aggr_sce$sample_id <- sapply(strsplit(colnames(sobh_sub_aggr_sce),"_"),
                                      `[[`,1)

length(unique(sobh_sub_aggr_sce$sample_id))

sobh_sub_aggr_sce$sample_name <- sapply(strsplit(sobh_sub_aggr_sce$sample_id,"-"),
                                        `[[`,1)
length(unique(sobh_sub_aggr_sce$sample_name))

## Pulling CT names
sobh_sub_aggr_sce$CT_final <- sapply(strsplit(colnames(sobh_sub_aggr_sce),"_"),
                                     `[[`,2)
length(unique(sobh_sub_aggr_sce$CT_final))

#### GESTATIONAL AGE AND LIFE SPAN ----
## Pulling metadata 

meta_info <- data.frame(
  Gestational_age_weeks = sobj_sub$Gestational_age_weeks,
  Life_span_weeks = sobj_sub$Life_span_weeks,
  Sample_id = sobj_sub$Sample_combined)


colData(sobh_sub_aggr_sce) <- cbind(colData(sobh_sub_aggr_sce),
                                    meta_info[sobh_sub_aggr_sce$sample_id,])

saveRDS(sobh_sub_aggr_sce,"/scratch/smallapragada/run2/lm_regression_ga_ls_ruqian_2025_04_11/sobh_sub_aggr_sce_2025_04_17.rds")

## Per cell type, test genes with the LM 
# g ~ life_span_weeks + g_age_weeks 

for( test_ct in unique(sobh_sub_aggr_sce$CT_final)){
  
  dir.create(file.path("/scratch/smallapragada/run2/lm_regression_ga_ls_ruqian_2025_04_11/2025-04-17/min3Per1k_prop20/",test_ct),
             recursive = TRUE)
  
  sobh_sub_aggr_sce_ct <- sobh_sub_aggr_sce[,sobh_sub_aggr_sce$CT_final==test_ct]
  life_span_weeks <- as.numeric(sobh_sub_aggr_sce_ct$Life_span_weeks)
  g_age_weeks <- as.numeric(sobh_sub_aggr_sce_ct$Gestational_age_weeks)
  sample_name <- sobh_sub_aggr_sce_ct$sample_name
  
  ## linear terms only 
  test_genes <- min3Per1k_prop20[[test_ct]]
  
  par(mfrow=c(1,2))
  message("length ",test_ct," ", length(test_genes))
  design <- model.matrix(~ life_span_weeks + g_age_weeks)
  v <- voom(sobh_sub_aggr_sce_ct@assays@data$counts[test_genes,], 
            design, plot=F)
  vfit <- lmFit(v, design)
  efit <- eBayes(vfit)
  head(efit$coefficients)
  SE <- as.data.frame(sqrt(efit$s2.post) * efit$stdev.unscaled)
  
  # Pull stats associated with GA
  sig_genes_gage <- topTable(efit,coef =c(3), p.value = 1,
                             number = 300,lfc = 0,
                             sort.by = "logFC")
  sig_genes_gage$gene <- rownames(sig_genes_gage)
  
  # 95% Confidence intervals for LFC
  sig_genes_gage$se_lower <- (sig_genes_gage$logFC - (SE$g_age_weeks * 1.96))
  sig_genes_gage$se_upper <- (sig_genes_gage$logFC + (SE$g_age_weeks * 1.96))
  
  write.csv(sig_genes_gage,file.path("/scratch/smallapragada/run2/lm_regression_ga_ls_ruqian_2025_04_11/2025-04-17/min3Per1k_prop20/",test_ct,
                                     paste0("/test_g_age_genes",
                                            test_ct,"_v_plot.csv")))
  saveRDS(efit, file.path("/scratch/smallapragada/run2/lm_regression_ga_ls_ruqian_2025_04_11/2025-04-17/min3Per1k_prop20/",test_ct,
                          paste0("/efit.",
                                 test_ct,".rds")))
  
  # Pull stats associated with LS
  sig_genes_life_span <- topTable(efit,coef =c(2),
                                  p.value = 1,
                                  number = 300,lfc = 0,
                                  sort.by = "logFC")
  sig_genes_life_span$gene <- rownames(sig_genes_life_span)
  
  # 95% Confidence intervals for LFC
  sig_genes_life_span$se_lower <- (sig_genes_life_span$logFC - (SE$life_span_weeks * 1.96))
  sig_genes_life_span$se_upper <- (sig_genes_life_span$logFC + (SE$life_span_weeks * 1.96))
  
  write.csv(sig_genes_life_span,file.path("/scratch/smallapragada/run2/lm_regression_ga_ls_ruqian_2025_04_11/2025-04-17/min3Per1k_prop20/",test_ct,
                                          paste0("test_life_span_genes",
                                                 test_ct,"_v_plot.csv")))
  
  write.csv(efit$coefficients,
            file.path("/scratch/smallapragada/run2/lm_regression_ga_ls_ruqian_2025_04_11/2025-04-17/min3Per1k_prop20/",test_ct,
                      paste0("no_interaction_coefficient_",
                             test_ct,".csv")))
  
}

## Creating a list of the CT names
id_list <- c("Activated fibroblasts", "Adventitial fibroblasts", "AKR1C1+ & AKR1C2+", "Alveolar macrophages", 
             "Arterial", "AT1", "AT2", "B", "Basal", "Capillaries", "cDC", "Fibroblasts", "Immature AT2", 
             "Lymphatic endothelial", "Mast", "Meg-Ery", "Monocytes", "Multiciliated", "MyoFB", "Neutrophils", 
             "NK & NKT", "pDC", "Pericytes", "Plasma", "PNEC", "Proliferating airway", "Proliferating basal",
             "Proliferating endothelial", "Proliferating Immature alveolar", "Proliferating lymphoid", 
             "Proliferating Meg-Ery", "Proliferating monocytes", "Secretory 3A2+ & 1A1+", "Secretory MUC5B+", 
             "SMC", "SPP1+ macrophages", "T", "Transitional AT2", "Treg", "Venous")

## Loading in all .csv files with stats from LM analysis
ga_test_stats <- list.files(path = "/scratch/smallapragada/run2/lm_regression_ga_ls_ruqian_2025_04_11/2025-04-17/min3Per1k_prop20/", 
                            pattern = "^test_g_age_genes",
                            recursive = TRUE, 
                            full.names = TRUE) 

ls_test_stats <- list.files(path = "/scratch/smallapragada/run2/lm_regression_ga_ls_ruqian_2025_04_11/2025-04-17/min3Per1k_prop20/", 
                            pattern = "^test_life_span_genes",
                            recursive = TRUE, 
                            full.names = TRUE) 

non_int_test_stats <- list.files(path = "/scratch/smallapragada/run2/lm_regression_ga_ls_ruqian_2025_04_11/2025-04-17/min3Per1k_prop20/", 
                                 pattern = "^no_interaction",
                                 recursive = TRUE, 
                                 full.names = TRUE) 

## Parsing through all files and fixing the formats for all
ga_test_stats_list <- lapply(ga_test_stats, function(XX)
  read.delim(XX, sep = ","))

ls_test_stats_list <- lapply(ls_test_stats, function(XX)
  read.delim(XX, sep = ","))

non_int_test_stats_list <- lapply(non_int_test_stats, function(XX)
  read.delim(XX, sep = ","))

## Assigning the loaded-in names to the list of files
names(ga_test_stats_list) <- id_list
names(ls_test_stats_list) <- id_list
names(non_int_test_stats_list) <- id_list

## Binding list together and pulling the values into a consolidated DF
ls_test_stats_list_df <- do.call(rbind, Map(function(df, name) {
  df$celltype <- name
  return(df)
}, ls_test_stats_list, names(ls_test_stats_list))) 

ga_test_stats_list_df <- do.call(rbind, Map(function(df, name) {
  df$celltype <- name
  return(df)
}, ga_test_stats_list, names(ga_test_stats_list)))

non_int_test_stats_list_df <- do.call(rbind, Map(function(df, name) {
  df$celltype <- name
  return(df)
}, non_int_test_stats_list, names(non_int_test_stats_list)))

## Adjusting dfs to include necessary information and join them
all_stats <- non_int_test_stats_list_df %>%
  rename_with(~ "gene", .cols = "X") %>%
  merge(., ls_test_stats_list_df, by = c("gene", "celltype")) %>%
  select(!X) %>%
  select(!logFC) %>%
  select(!AveExpr) %>%
  rename_with(~ "t_ls", .cols = "t") %>%
  rename_with(~ "p_val_ls", .cols = "P.Value") %>%
  rename_with(~ "adj_p_val_ls", .cols = "adj.P.Val") %>%
  rename_with(~ "b_ls", .cols = "B") %>%
  rename_with(~ "se_lower_ls", .cols = "se_lower") %>%
  rename_with(~ "se_upper_ls", .cols = "se_upper") %>%
  merge(., ga_test_stats_list_df, by = c("gene", "celltype")) %>%
  select(!X) %>%
  select(!logFC) %>%
  rename_with(~ "t_ga", .cols = "t") %>%
  rename_with(~ "p_val_ga", .cols = "P.Value") %>%
  rename_with(~ "adj_p_val_ga", .cols = "adj.P.Val") %>%
  rename_with(~ "b_ga", .cols = "B") %>%
  rename_with(~ "se_lower_ga", .cols = "se_lower") %>%
  rename_with(~ "se_upper_ga", .cols = "se_upper") %>%
  relocate(AveExpr, .after = celltype) %>%
  rename_with(~ "ls_coef_lfc", .cols = "life_span_weeks") %>%
  rename_with(~ "ga_coef_lfc", .cols = "g_age_weeks") %>%
  relocate(ga_coef_lfc, .before = t_ga)

write_csv(all_stats, "/scratch/smallapragada/all_stats.csv")

## Adding lineage labels
all_stats$lineage <- "NA"
all_stats$lineage[all_stats$celltype %in% c("AKR1C1+ & AKR1C2+", "AT1", "AT2", "Basal", 
                                            "Immature AT2", "Transitional AT2", 
                                            "Multiciliated", "PNEC", "Proliferating airway", 
                                            "Proliferating basal", "Proliferating Immature alveolar", 
                                            "Secretory 3A2+ & 1A1+", "Secretory MUC5B+")] <- "Epithelial"
all_stats$lineage[all_stats$celltype %in% c("Arterial", "Lymphatic endothelial", 
                                            "Proliferating endothelial", "Venous", 
                                            "Capillaries")] <- "Endothelial"
all_stats$lineage[all_stats$celltype %in% c("Activated fibroblasts", "Adventitial fibroblasts", 
                                            "Fibroblasts", "MyoFB", "Pericytes", "SMC")] <- "Mesenchymal"
all_stats$lineage[all_stats$celltype %in% c("B", "cDC", "Mast", "Meg-Ery", "Monocytes", "Alveolar macrophages", 
                                            "Neutrophils", "NK & NKT", "pDC", "Plasma",
                                            "Proliferating lymphoid",  
                                            "Proliferating Meg-Ery", "Proliferating monocytes",
                                            "SPP1+ macrophages", "T", "Treg")] <- "Immune"

## Adding sublineage labels
all_stats$sublineage <- "NA"
all_stats$sublineage[all_stats$celltype %in% c("AKR1C1+ & AKR1C2+", "AT1", "AT2", 
                                               "Immature AT2", "Transitional AT2", 
                                               "Proliferating Immature alveolar")] <- "Alveolar"
all_stats$sublineage[all_stats$celltype %in% c("Basal", "Multiciliated", "PNEC",
                                               "Proliferating airway", "Proliferating basal", 
                                               "Secretory 3A2+ & 1A1+", "Secretory MUC5B+")] <- "Airway"
all_stats$sublineage[all_stats$celltype %in% c("Arterial", "Lymphatic endothelial", 
                                               "Proliferating endothelial", "Venous", 
                                               "Capillaries")] <- "Endothelial"
all_stats$sublineage[all_stats$celltype %in% c("Activated fibroblasts", "Adventitial fibroblasts", 
                                               "Fibroblasts", "MyoFB", 
                                               "Pericytes", "SMC")] <- "Mesenchymal"
all_stats$sublineage[all_stats$celltype %in% c("Alveolar macrophages", 
                                               "Meg-Ery", "Mast", "Proliferating Meg-Ery",
                                               "Monocytes", "Proliferating monocytes", 
                                               "Neutrophils", "SPP1+ macrophages")] <- "Myeloid"
all_stats$sublineage[all_stats$celltype %in% c("pDC", "NK & NKT", "Treg", "Plasma", 
                                               "cDC", "T", "B", "Proliferating lymphoid")] <- "Lymphoid"

all_stats_ls_sig <- all_stats %>%
  relocate(lineage, .after = celltype) %>% 
  relocate(sublineage, .after = lineage) %>% 
  mutate(gene_ct = paste(gene, "-", celltype)) %>%
  mutate(sig_ls = case_when(adj_p_val_ls < 0.1 ~ "Significant",
                            TRUE ~ "Not significant"))

all_stats_ga_sig <- all_stats %>%
  relocate(lineage, .after = celltype) %>% 
  relocate(sublineage, .after = lineage) %>% 
  mutate(gene_ct = paste(gene, "-", celltype)) %>%
  mutate(sig_ga = case_when(adj_p_val_ga < 0.1 ~ "Significant",
                            TRUE ~ "Not significant"))

#### DX SCORE ----

## Pulling metadata
meta_info <- data.frame(
  dx = sobj_sub$dx,
  Sample_id = sobj_sub$Sample_combined)

colData(sobh_sub_aggr_sce) <- cbind(colData(sobh_sub_aggr_sce),
                                    meta_info[sobh_sub_aggr_sce$sample_id,])

saveRDS(sobh_sub_aggr_sce,"/scratch/smallapragada/run2/lm_regression_ga_ls_ruqian_2025_04_11/sobh_sub_aggr_sce_path_score_2025_04_28.rds")

## Per cell type, test genes with the LM
# g ~ dx
for(test_ct in unique(sobh_sub_aggr_sce$CT_final)) {
  dir.create(file.path("/scratch/smallapragada/run2/lm_regression_ga_ls_ruqian_2025_04_11/2025-04-28/min3Per1k_prop20/",test_ct),
             recursive = TRUE)
  
  sobh_sub_aggr_sce_ct <- sobh_sub_aggr_sce[,sobh_sub_aggr_sce$CT_final==test_ct]
  dx <- as.numeric(sobh_sub_aggr_sce_ct$dx)
  sample_name <- sobh_sub_aggr_sce_ct$sample_name
  
  ## linear terms only 
  test_genes <- min3Per1k_prop20[[test_ct]]
  
  par(mfrow=c(1,2))
  message("length ",test_ct," ", length(test_genes))
  design <- model.matrix(~ dx)
  v <- voom(sobh_sub_aggr_sce_ct@assays@data$counts[test_genes,], 
            design, plot=F)
  vfit <- lmFit(v, design)
  efit <- eBayes(vfit)
  head(efit$coefficients)
  SE <- as.data.frame(sqrt(efit$s2.post) * efit$stdev.unscaled)
  
  ## Pulling stats associated with dx
  sig_genes_dx <- topTable(efit,coef =c(2), p.value = 1,
                           number = 300,lfc = 0,
                           sort.by = "logFC")
  sig_genes_dx$gene <- rownames(sig_genes_dx)
  
  # 95% confidence intervals
  sig_genes_dx$se_lower <- (sig_genes_dx$logFC - (SE$dx * 1.96))
  sig_genes_dx$se_upper <- (sig_genes_dx$logFC + (SE$dx * 1.96))
  
  write.csv(sig_genes_dx, file.path("/scratch/smallapragada/run2/lm_regression_ga_ls_ruqian_2025_04_11/2025-04-28/min3Per1k_prop20/",test_ct,
                                    paste0("/test_dx_genes",
                                           test_ct,"_v_plot.csv")))
  
  saveRDS(efit, file.path("/scratch/smallapragada/run2/lm_regression_ga_ls_ruqian_2025_04_11/2025-04-28/min3Per1k_prop20/",test_ct,
                          paste0("/efit.",
                                 test_ct,".rds")))

## Loading in all .csv files with stats from LM analysis
dx_test_stats <- list.files(path = "/scratch/smallapragada/run2/lm_regression_ga_ls_ruqian_2025_04_11/2025-04-28/min3Per1k_prop20/", 
                            pattern = "^test_dx_genes",
                            recursive = TRUE, 
                            full.names = TRUE) 

## Parsing through all files and fixing the formats for all
dx_test_stats_list <- lapply(dx_test_stats, function(XX)
  read.delim(XX, sep = ","))

## Assigning the loaded-in names to the list of files
names(dx_test_stats_list) <- id_list

## Binding list together and pulling the values into a consolidated DF
dx_test_stats_list_df <- do.call(rbind, Map(function(df, name) {
  df$celltype <- name
  return(df)
}, dx_test_stats_list, names(dx_test_stats_list))) 

## Adjusting df to include necessary information 
dx_stats <- dx_test_stats_list_df %>%
  select(!X) %>%
  rename_with(~ "t_dx", .cols = "t") %>%
  rename_with(~ "p_val_ls", .cols = "P.Value") %>%
  rename_with(~ "adj_p_val_dx", .cols = "adj.P.Val") %>%
  rename_with(~ "b_dx", .cols = "B") %>%
  rename_with(~ "se_lower_dx", .cols = "se_lower") %>%
  rename_with(~ "se_upper_dx", .cols = "se_upper") %>%
  relocate(AveExpr, .after = celltype) %>%
  rename_with(~ "dx_coef_lfc", .cols = "logFC") %>%
  mutate(gene_ct = paste(gene, "-", celltype)) %>%
  mutate(sig_dx = case_when(adj_p_val_dx < 0.1 ~ "Significant",
                            TRUE ~ "Not significant"))

write_csv(dx_stats, "/scratch/smallapragada/dx_stats.csv")

## Adding lineage labels
dx_stats$lineage <- "NA"
dx_stats$lineage[dx_stats$celltype %in% c("AKR1C1+ & AKR1C2+", "AT1", "AT2", "Basal", 
                                          "Immature AT2", "Transitional AT2", 
                                          "Multiciliated", "PNEC", "Proliferating airway", 
                                          "Proliferating basal", "Proliferating Immature alveolar", 
                                          "Secretory 3A2+ & 1A1+", "Secretory MUC5B+")] <- "Epithelial"
dx_stats$lineage[dx_stats$celltype %in% c("Arterial", "Lymphatic endothelial", 
                                          "Proliferating endothelial", "Venous", 
                                          "Capillaries")] <- "Endothelial"
dx_stats$lineage[dx_stats$celltype %in% c("Activated fibroblasts", "Adventitial fibroblasts", 
                                          "Fibroblasts", "MyoFB", "Pericytes", "SMC")] <- "Mesenchymal"
dx_stats$lineage[dx_stats$celltype %in% c("B", "cDC", "Mast", "Meg-Ery", "Monocytes", "Alveolar macrophages", 
                                          "Neutrophils", "NK & NKT", "pDC", "Plasma",
                                          "Proliferating lymphoid",
                                          "Proliferating Meg-Ery", "Proliferating monocytes",
                                          "SPP1+ macrophages", "T", "Treg")] <- "Immune"

#### INTERSECTIONS BETWEEN ALL THREE TERMS ----

## Tabulating results for each term for a two-level p-value cutoff
ls_stats_fdr10 <- all_stats_ls_sig %>%
  filter(adj_p_val_ls < 0.1) 

ls_stats_fdr20 <- all_stats_ls_sig %>%
  filter(adj_p_val_ls < 0.2)

ga_stats_fdr10 <- all_stats_ga_sig %>%
  filter(adj_p_val_ga < 0.1)

ga_stats_fdr20 <- all_stats_ga_sig %>%
  filter(adj_p_val_ga < 0.2)

dx_stats_fdr10 <- dx_stats %>%
  filter(adj_p_val_dx < 0.1)

dx_stats_fdr20 <- dx_stats %>%
  filter(adj_p_val_dx < 0.2)

ga_pairs_10 <- unique(ga_stats_fdr10$gene_ct)
ls_pairs_10 <- unique(ls_stats_fdr10$gene_ct)
dx_pairs_10 <- unique(dx_stats_fdr10$gene_ct)

ga_pairs_20 <- unique(ga_stats_fdr20$gene_ct)
ls_pairs_20 <- unique(ls_stats_fdr20$gene_ct)
dx_pairs_20 <- unique(dx_stats_fdr20$gene_ct)

pairs_to_keep <- unique(c(ga_pairs_10, ls_pairs_10, dx_pairs_10))

ga_pairs_with_thresh <- ga_stats_fdr20 %>%
  filter(gene_ct %in% pairs_to_keep) %>%
  select(gene_ct, ga_coef_lfc) %>%
  rename_with(~ "lfc", .cols = "ga_coef_lfc") %>%
  mutate(origin = "ga") 

ls_pairs_with_thresh <- ls_stats_fdr20 %>%
  filter(gene_ct %in% pairs_to_keep) %>%
  select(gene_ct, ls_coef_lfc) %>%
  rename_with(~ "lfc", .cols = "ls_coef_lfc") %>%
  mutate(origin = "ls") 

dx_pairs_with_thresh <- dx_stats_fdr20 %>%
  filter(gene_ct %in% pairs_to_keep) %>%
  select(gene_ct, dx_coef_lfc) %>%
  rename_with(~ "lfc", .cols = "dx_coef_lfc") %>%
  mutate(origin = "dx") 

## Merge dfs together

all_stats_ga_morphed <- all_stats_ga_sig %>%
  mutate(gene_ct = paste(gene, "-", celltype)) %>%
  select(gene, celltype, gene_ct, ga_coef_lfc, adj_p_val_ga, se_lower_ga, se_upper_ga) %>%
  rename_with(~ "adj_p_val", .cols = "adj_p_val_ga") %>%
  rename_with(~ "lfc", .cols = "ga_coef_lfc") %>%
  rename_with(~ "se_lower", .cols = "se_lower_ga") %>%
  rename_with(~ "se_upper", .cols = "se_upper_ga") %>%
  mutate(origin = "ga") %>%
  mutate(sig = case_when(gene_ct %in% ga_pairs_with_thresh$gene_ct ~ "Significant",
                         TRUE ~ "N.S."))

all_stats_ga_ls_morphed <- all_stats_ls_sig %>%
  mutate(gene_ct = paste(gene, "-", celltype)) %>%
  select(gene, celltype, gene_ct, ls_coef_lfc, adj_p_val_ls, se_lower_ls, se_upper_ls) %>%
  rename_with(~ "adj_p_val", .cols = "adj_p_val_ls") %>%
  rename_with(~ "lfc", .cols = "ls_coef_lfc") %>%
  rename_with(~ "se_lower", .cols = "se_lower_ls") %>%
  rename_with(~ "se_upper", .cols = "se_upper_ls") %>%
  mutate(origin = "ls") %>%
  mutate(sig = case_when(gene_ct %in% ls_pairs_with_thresh$gene_ct ~ "Significant",
                         TRUE ~ "N.S.")) %>%
  rbind(., all_stats_ga_morphed)

all_stats_morphed <- dx_stats %>%
  mutate(gene_ct = paste(gene, "-", celltype)) %>%
  select(gene, celltype, gene_ct, dx_coef_lfc, adj_p_val_dx, se_lower_dx, se_upper_dx) %>%
  rename_with(~ "adj_p_val", .cols = "adj_p_val_dx") %>%
  rename_with(~ "lfc", .cols = "dx_coef_lfc") %>%
  rename_with(~ "se_lower", .cols = "se_lower_dx") %>%
  rename_with(~ "se_upper", .cols = "se_upper_dx") %>%
  mutate(origin = "dx") %>%
  mutate(sig = case_when(gene_ct %in% dx_pairs_with_thresh$gene_ct | adj_p_val < 0.1 ~ "Significant",
                         TRUE ~ "N.S.")) %>%
  rbind(all_stats_ga_ls_morphed)
