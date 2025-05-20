#### LOAD PACKAGES AND SET ENVIRONMENT ----
library(SeuratObject)
library(tidyverse)
library(gplots)
library(rlist)
library(rstatix)
library(Seurat, lib.loc = "/usr/local/lib/R/site-library")

# Set seed & working directory
set.seed(0712)
work_dir <- setwd("/scratch/smallapragada/run2")

#### ALL SAMPLES ----

## Creating a list of the CTs
id_list <- c("Activated fibroblasts", "Adventitial fibroblasts", "AKR1C1+ & AKR1C2+", "Alveolar macrophages", "Arterial",
             "AT1", "AT2", "B", "Basal", "Capillaries", "cDC", "Fibroblasts", "Immature AT2",
             "Lymphatic endothelial", "Mast", "Meg-Ery", "Monocytes", "Multiciliated",
             "MyoFB", "Neutrophils", "NK & NKT", "pDC", "Pericytes", "Plasma", "PNEC", "Proliferating airway", "Proliferating basal",
             "Proliferating endothelial", "Proliferating Immature alveolar", "Proliferating lymphoid", "Proliferating Meg-Ery", 
             "Proliferating monocytes", "Secretory 3A2+ & 1A1+", "Secretory MUC5B+", "SMC", "SPP1+ macrophages", "T", "Transitional AT2",
             "Treg", "Venous")

## Load data and define code paramters
outDir <- '/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/output_stats'
options(contrasts=c("contr.sum","contr.poly"))

## Read in object & get metadata
run2_final <- readRDS("/scratch/smallapragada/run2/final_complete_cniches_tniches_allsamples_2025_03_26.rds")
meta <- run2_final@meta.data

## Add new column in metadata for groups
meta <- meta %>%
  mutate(sample_similarity = case_when(meta$Sample_combined %in% c("PDL001", "PDL002") ~ "Early canalicular",
                                       meta$Sample_combined %in% c("PDL003", "PDL004","PDL005", "PDL006") ~ "Late canalicular",
                                       meta$Sample_combined %in% c("PDL007", "PDL008", "PDL009", "PDL010", "PDL011") ~ "Saccular",
                                       meta$Sample_combined %in% c("PDL012", "PDL013") ~ "Alveolar",
                                       meta$Sample_combined %in% c("PDL014") ~ "Rare disease (CHAOS)",
                                       meta$Sample_combined %in% c("PDL015") ~ "Rare disease (PH + Kidneys)",
                                       meta$Sample_combined %in% c("PDL016") ~ "Rare disease (Infant BPD)",
                                       meta$Sample_combined %in% c("PDL017") ~ "Rare disease (Adult BPD)",
                                       TRUE ~ "help"))

### All samples

## Loop to iterate through each CT and calculate stats
enrich_score <- c()

for(i in id_list) {
  
  celltype_id <- i
  print(celltype_id)
  cell_type_id_out <- gsub("\\/", "\\_", celltype_id)
  
  inpath <-  '/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/output'
  f <- paste0(cell_type_id_out, '_cell_dist_adjdegree.rds')
  df <- readRDS(file.path(inpath, f))
  df <- df[order(df$V1, decreasing = F ),]
  
  # set cell type threshold
  d_threshold <- df %>% 
    group_by(cell.a) %>%
    group_by(idx) %>%
    summarise(d_threshold = mean(V1) + 2*sd(V1)) %>%
    filter(idx == 1) %>%
    pull(d_threshold)
  d_threshold
  
  tmp0 <- df %>% 
    filter(V1 < d_threshold & idx == 1)
  
  tmp0 <- tmp0[!duplicated(tmp0),]
  tmp0 <- tmp0[,c('cell.b', 'celltypeB', 'sid')] ## rm V1, cell.a, degree, angle, idx
  tmp0 <- tmp0[!duplicated(tmp0),]
  
  # proximal cells by cell type
  proximal_cell_counts <- table(tmp0$celltypeB) %>% 
    as.data.frame() %>%
    rename_with(~ "Freq_neighbor_CT", .cols = "Freq") 
  proximal_cell_counts$neighbor_CT_count = sum(proximal_cell_counts$Freq_neighbor_CT)
  
  all_cell_counts <- table(meta$CT_final) %>% 
    as.data.frame()
  all_cell_counts$all_cell_count = sum(all_cell_counts$Freq)
  
  tmp1 <- merge(proximal_cell_counts, all_cell_counts, by = 'Var1', all = T)
  tmp1[is.na(tmp1)] = 0
  
  enrich_score_ct <- c()
  for(j in 1:nrow(tmp1)) {
    ct <- tmp1[j, 1]
    total <- tmp1[j, 5]
    list1 <- tmp1[j, 4]
    list2 <- tmp1[j, 3]
    overlap <- tmp1[j, 2]
    res <- fisher.test(matrix(c((total - (list1 + list2) + overlap),  (list1-overlap),  (list2 - overlap), overlap), nrow=2))
    tmp2 <- data.frame(source_ct = celltype_id, neighboring_ct = ct,
                       p_val = res$p.value, odds_ratio = res$estimate, or05 = res$conf.int[1], or95 = res$conf.int[2])
    enrich_score_ct <- rbind(enrich_score_ct, tmp2)
  }
  
  enrich_score <- rbind(enrich_score_ct, enrich_score)
}

enrich_score <- enrich_score %>% 
  mutate(adj_p_val = p.adjust(p_val, method = 'fdr')) %>%
  mutate(logOR = log(odds_ratio)) %>%
  mutate(logOR_lower = log(or05)) %>%
  mutate(logOR_upper = log(or95)) %>%
  mutate(sample_similarity = "All samples")

# Save for plotting
saveRDS(enrich_score, "/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/proximity_enrichment_30r_2025_04_21.rds")

write_csv(enrich_score, "/scratch/smallapragada/run2/prox_logOR_csvs/enrichment_scores_all_samples_2025_05_08.csv")

#### EARLY CANALICULAR ----

meta_ec <- meta %>%
  filter(sample_similarity == "Early canalicular")

## Loop to iterate through each CT and calculate stats
enrich_score_ec <- c()

for(i in id_list) {
  
  celltype_id <- i
  print(celltype_id)
  cell_type_id_out <- gsub("\\/", "\\_", celltype_id)
  
  inpath <-  '/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/output'
  f <- paste0(cell_type_id_out, '_cell_dist_adjdegree.rds')
  df <- readRDS(file.path(inpath, f))
  df <- df[order(df$V1, decreasing = F ),]
  
  # set cell type threshold
  d_threshold <- df %>%
    filter(idx == 1) %>%
    filter(sid == "s1_PDL001" | sid == "s2_PDL001" | 
             sid == "s1_PDL002" | sid == "s2_PDL002") %>%
    summarise(d_threshold = mean(V1) + 2 * sd(V1)) %>%
    pull(d_threshold)
  
  tmp0 <- df %>%
    filter(V1 < d_threshold & idx == 1) %>%
    filter(sid == "s1_PDL001" | sid == "s2_PDL001" | 
             sid == "s1_PDL002" | sid == "s2_PDL002") %>%
    distinct(cell.b, celltypeB, sid)
  
  tmp0 <- tmp0[!duplicated(tmp0),]
  tmp0 <- tmp0[,c('cell.b', 'celltypeB', 'sid')] ## rm V1, cell.a, degree, angle, idx
  tmp0 <- tmp0[!duplicated(tmp0),]
  
  # proximal cells by cell type
  proximal_cell_counts <- table(tmp0$celltypeB) %>% 
    as.data.frame() %>%
    rename_with(~ "Freq_neighbor_CT", .cols = "Freq") 
  
  if (nrow(proximal_cell_counts) > 0) {
    proximal_cell_counts$neighbor_CT_count = sum(proximal_cell_counts$Freq_neighbor_CT)
  } else {
    message("No proximal cells found for ", celltype_id)
    proximal_cell_counts <- data.frame(Var1 = character(),
                                       Freq_neighbor_CT = integer(),
                                       neighbor_CT_count = integer())
  }
  
  all_cell_counts <- table(meta_ec$CT_final) %>% 
    as.data.frame()
  all_cell_counts$all_cell_count = sum(all_cell_counts$Freq)
  
  tmp1 <- merge(proximal_cell_counts, all_cell_counts, by = 'Var1', all = T)
  tmp1[is.na(tmp1)] = 0
  
  enrich_score_ct <- c()
  
  for (j in 1:nrow(tmp1)) {
    ct <- tmp1[j, 1]
    total <- tmp1[j, 5]
    list1 <- tmp1[j, 4]
    list2 <- tmp1[j, 3]
    overlap <- tmp1[j, 2]
    
    # Safe Fisher matrix calculation
    contingency <- matrix(c(
      total - (list1 + list2) + overlap,
      list1 - overlap,
      list2 - overlap,
      overlap
    ), nrow = 2)
    
    if (any(contingency < 0 | !is.finite(contingency))) {
      res <- list(p.value = NA, estimate = NA)
    } else {
      res <- fisher.test(contingency)
    }
    
    tmp2 <- data.frame(
      source_ct = celltype_id,
      neighboring_ct = ct,
      p_val = res$p.value,
      odds_ratio = as.numeric(res$estimate),
      or05 = res$conf.int[1], 
      or95 = res$conf.int[2]
    )
    
    enrich_score_ct <- rbind(enrich_score_ct, tmp2)
  }
  
  enrich_score_ec <- rbind(enrich_score_ct, enrich_score_ec)
}

enrich_score_ec <- enrich_score_ec %>% 
  mutate(adj_p_val = p.adjust(p_val, method = 'fdr')) %>%
  mutate(logOR = log(odds_ratio)) %>%
  mutate(logOR_lower = log(or05)) %>%
  mutate(logOR_upper = log(or95)) %>%
  mutate(sample_similarity = "Early canalicular")

# Save for plotting
saveRDS(enrich_score_ec, "/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/proximity_enrichment_30r_ec_2025_04_21.rds")

write_csv(enrich_score_ec, "/scratch/smallapragada/run2/prox_logOR_csvs/enrichment_scores_ec_2025_05_08.csv")

#### LATE CANALICULAR ----

meta_lc <- meta %>%
  filter(sample_similarity == "Late canalicular")

## Loop to iterate through each CT and calculate stats
enrich_score_lc <- c()

for(i in id_list) {
  
  celltype_id <- i
  print(celltype_id)
  cell_type_id_out <- gsub("\\/", "\\_", celltype_id)
  
  inpath <-  '/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/output'
  f <- paste0(cell_type_id_out, '_cell_dist_adjdegree.rds')
  df <- readRDS(file.path(inpath, f))
  df <- df[order(df$V1, decreasing = F ),]
  
  # set cell type threshold
  d_threshold <- df %>%
    filter(idx == 1) %>%
    filter(sid == "s1_PDL003" | sid == "s2_PDL003" | 
             sid == "s1_PDL004" | sid == "s2_PDL004" | 
             sid == "s1_PDL005" | sid == "s2_PDL005" |
             sid == "s1_PDL006" | sid == "s2_PDL006") %>%
    summarise(d_threshold = mean(V1) + 2 * sd(V1)) %>%
    pull(d_threshold)
  
  tmp0 <- df %>%
    filter(V1 < d_threshold & idx == 1) %>%
    filter(sid == "s1_PDL003" | sid == "s2_PDL003" | 
             sid == "s1_PDL004" | sid == "s2_PDL004" | 
             sid == "s1_PDL005" | sid == "s2_PDL005" |
             sid == "s1_PDL006" | sid == "s2_PDL006") %>%
    distinct(cell.b, celltypeB, sid)
  
  tmp0 <- tmp0[!duplicated(tmp0),]
  tmp0 <- tmp0[,c('cell.b', 'celltypeB', 'sid')] ## rm V1, cell.a, degree, angle, idx
  tmp0 <- tmp0[!duplicated(tmp0),]
  
  # proximal cells by cell type
  proximal_cell_counts <- table(tmp0$celltypeB) %>% 
    as.data.frame() %>%
    rename_with(~ "Freq_neighbor_CT", .cols = "Freq") 
  
  if (nrow(proximal_cell_counts) > 0) {
    proximal_cell_counts$neighbor_CT_count = sum(proximal_cell_counts$Freq_neighbor_CT)
  } else {
    message("No proximal cells found for ", celltype_id)
    proximal_cell_counts <- data.frame(Var1 = character(),
                                       Freq_neighbor_CT = integer(),
                                       neighbor_CT_count = integer())
  }
  
  all_cell_counts <- table(meta_lc$CT_final) %>% 
    as.data.frame()
  all_cell_counts$all_cell_count = sum(all_cell_counts$Freq)
  
  tmp1 <- merge(proximal_cell_counts, all_cell_counts, by = 'Var1', all = T)
  tmp1[is.na(tmp1)] = 0
  
  enrich_score_ct <- c()
  
  for (j in 1:nrow(tmp1)) {
    ct <- tmp1[j, 1]
    total <- tmp1[j, 5]
    list1 <- tmp1[j, 4]
    list2 <- tmp1[j, 3]
    overlap <- tmp1[j, 2]
    
    # Safe Fisher matrix calculation
    contingency <- matrix(c(
      total - (list1 + list2) + overlap,
      list1 - overlap,
      list2 - overlap,
      overlap
    ), nrow = 2)
    
    if (any(contingency < 0 | !is.finite(contingency))) {
      res <- list(p.value = NA, estimate = NA)
    } else {
      res <- fisher.test(contingency)
    }
    
    tmp2 <- data.frame(
      source_ct = celltype_id,
      neighboring_ct = ct,
      p_val = res$p.value,
      odds_ratio = as.numeric(res$estimate),
      or05 = res$conf.int[1], 
      or95 = res$conf.int[2]
    )
    
    enrich_score_ct <- rbind(enrich_score_ct, tmp2)
  }
  
  enrich_score_lc <- rbind(enrich_score_ct, enrich_score_lc)
}

enrich_score_lc <- enrich_score_lc %>% 
  mutate(adj_p_val = p.adjust(p_val, method = 'fdr')) %>%
  mutate(logOR = log(odds_ratio)) %>%
  mutate(logOR_lower = log(or05)) %>%
  mutate(logOR_upper = log(or95)) %>%
  mutate(sample_similarity = "Late canalicular")

# Save for plotting
saveRDS(enrich_score_lc, "/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/proximity_enrichment_30r_lc_2025_04_21.rds")

write_csv(enrich_score_lc, "/scratch/smallapragada/run2/prox_logOR_csvs/enrichment_scores_lc_2025_05_08.csv")

#### SACCULAR ----

meta_s <- meta %>%
  filter(sample_similarity == "Saccular")

## Loop to iterate through each CT and calculate stats
enrich_score_s <- c()

for(i in id_list) {
  
  celltype_id <- i
  print(celltype_id)
  cell_type_id_out <- gsub("\\/", "\\_", celltype_id)
  
  inpath <-  '/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/output'
  f <- paste0(cell_type_id_out, '_cell_dist_adjdegree.rds')
  df <- readRDS(file.path(inpath, f))
  df <- df[order(df$V1, decreasing = F ),]
  
  # set cell type threshold
  d_threshold <- df %>%
    filter(idx == 1) %>%
    filter(sid == "s1_PDL007" | sid == "s2_PDL007" | 
             sid == "s1_PDL008" | sid == "s2_PDL008" | 
             sid == "s1_PDL009" | sid == "s2_PDL009" |
             sid == "s1_PDL010" | sid == "s2_PDL010" |
             sid == "s1_PDL011" | sid == "s2_PDL011") %>%
    summarise(d_threshold = mean(V1) + 2 * sd(V1)) %>%
    pull(d_threshold)
  
  tmp0 <- df %>%
    filter(V1 < d_threshold & idx == 1) %>%
    filter(sid == "s1_PDL007" | sid == "s2_PDL007" | 
             sid == "s1_PDL008" | sid == "s2_PDL008" | 
             sid == "s1_PDL009" | sid == "s2_PDL009" |
             sid == "s1_PDL010" | sid == "s2_PDL010" |
             sid == "s1_PDL011" | sid == "s2_PDL011") %>%
    distinct(cell.b, celltypeB, sid)
  
  tmp0 <- tmp0[!duplicated(tmp0),]
  tmp0 <- tmp0[,c('cell.b', 'celltypeB', 'sid')] ## rm V1, cell.a, degree, angle, idx
  tmp0 <- tmp0[!duplicated(tmp0),]
  
  # proximal cells by cell type
  proximal_cell_counts <- table(tmp0$celltypeB) %>% 
    as.data.frame() %>%
    rename_with(~ "Freq_neighbor_CT", .cols = "Freq") 
  
  if (nrow(proximal_cell_counts) > 0) {
    proximal_cell_counts$neighbor_CT_count = sum(proximal_cell_counts$Freq_neighbor_CT)
  } else {
    message("No proximal cells found for ", celltype_id)
    proximal_cell_counts <- data.frame(Var1 = character(),
                                       Freq_neighbor_CT = integer(),
                                       neighbor_CT_count = integer())
  }
  
  all_cell_counts <- table(meta_s$CT_final) %>% 
    as.data.frame()
  all_cell_counts$all_cell_count = sum(all_cell_counts$Freq)
  
  tmp1 <- merge(proximal_cell_counts, all_cell_counts, by = 'Var1', all = T)
  tmp1[is.na(tmp1)] = 0
  
  enrich_score_ct <- c()
  
  for (j in 1:nrow(tmp1)) {
    ct <- tmp1[j, 1]
    total <- tmp1[j, 5]
    list1 <- tmp1[j, 4]
    list2 <- tmp1[j, 3]
    overlap <- tmp1[j, 2]
    
    # Safe Fisher matrix calculation
    contingency <- matrix(c(
      total - (list1 + list2) + overlap,
      list1 - overlap,
      list2 - overlap,
      overlap
    ), nrow = 2)
    
    if (any(contingency < 0 | !is.finite(contingency))) {
      res <- list(p.value = NA, estimate = NA)
    } else {
      res <- fisher.test(contingency)
    }
    
    tmp2 <- data.frame(
      source_ct = celltype_id,
      neighboring_ct = ct,
      p_val = res$p.value,
      odds_ratio = as.numeric(res$estimate),
      or05 = res$conf.int[1], 
      or95 = res$conf.int[2]
    )
    
    enrich_score_ct <- rbind(enrich_score_ct, tmp2)
  }
  
  enrich_score_s <- rbind(enrich_score_ct, enrich_score_s)
}

enrich_score_s <- enrich_score_s %>% 
  mutate(adj_p_val = p.adjust(p_val, method = 'fdr')) %>%
  mutate(logOR = log(odds_ratio)) %>%
  mutate(logOR_lower = log(or05)) %>%
  mutate(logOR_upper = log(or95)) %>%
  mutate(sample_similarity = "Saccular")

# Save for plotting
saveRDS(enrich_score_s, "/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/proximity_enrichment_30r_s_2025_04_21.rds")

write_csv(enrich_score_s, "/scratch/smallapragada/run2/prox_logOR_csvs/enrichment_scores_s_2025_05_08.csv")

#### ALVEOLAR ----

meta_a <- meta %>%
  filter(sample_similarity == "Alveolar")

## Loop to iterate through each CT and calculate stats
enrich_score_a <- c()

for(i in id_list) {
  
  celltype_id <- i
  print(celltype_id)
  cell_type_id_out <- gsub("\\/", "\\_", celltype_id)
  
  inpath <-  '/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/output'
  f <- paste0(cell_type_id_out, '_cell_dist_adjdegree.rds')
  df <- readRDS(file.path(inpath, f))
  df <- df[order(df$V1, decreasing = F ),]
  
  # set cell type threshold
  d_threshold <- df %>%
    filter(idx == 1) %>%
    filter(sid == "s1_PDL013" | sid == "s2_PDL013" |
             sid == "s1_PDL012" | sid == "s2_PDL012") %>%
    summarise(d_threshold = mean(V1) + 2 * sd(V1)) %>%
    pull(d_threshold)
  
  tmp0 <- df %>%
    filter(V1 < d_threshold & idx == 1) %>%
    filter(sid == "s1_PDL013" | sid == "s2_PDL013" |
             sid == "s1_PDL012" | sid == "s2_PDL012") %>%
    distinct(cell.b, celltypeB, sid)
  
  tmp0 <- tmp0[!duplicated(tmp0),]
  tmp0 <- tmp0[,c('cell.b', 'celltypeB', 'sid')] ## rm V1, cell.a, degree, angle, idx
  tmp0 <- tmp0[!duplicated(tmp0),]
  
  # proximal cells by cell type
  proximal_cell_counts <- table(tmp0$celltypeB) %>% 
    as.data.frame() %>%
    rename_with(~ "Freq_neighbor_CT", .cols = "Freq") 
  
  if (nrow(proximal_cell_counts) > 0) {
    proximal_cell_counts$neighbor_CT_count = sum(proximal_cell_counts$Freq_neighbor_CT)
  } else {
    message("No proximal cells found for ", celltype_id)
    proximal_cell_counts <- data.frame(Var1 = character(),
                                       Freq_neighbor_CT = integer(),
                                       neighbor_CT_count = integer())
  }
  
  all_cell_counts <- table(meta_a$CT_final) %>% 
    as.data.frame()
  all_cell_counts$all_cell_count = sum(all_cell_counts$Freq)
  
  tmp1 <- merge(proximal_cell_counts, all_cell_counts, by = 'Var1', all = T)
  tmp1[is.na(tmp1)] = 0
  
  enrich_score_ct <- c()
  
  for (j in 1:nrow(tmp1)) {
    ct <- tmp1[j, 1]
    total <- tmp1[j, 5]
    list1 <- tmp1[j, 4]
    list2 <- tmp1[j, 3]
    overlap <- tmp1[j, 2]
    
    # Safe Fisher matrix calculation
    contingency <- matrix(c(
      total - (list1 + list2) + overlap,
      list1 - overlap,
      list2 - overlap,
      overlap
    ), nrow = 2)
    
    if (any(contingency < 0 | !is.finite(contingency))) {
      res <- list(p.value = NA, estimate = NA)
    } else {
      res <- fisher.test(contingency)
    }
    
    tmp2 <- data.frame(
      source_ct = celltype_id,
      neighboring_ct = ct,
      p_val = res$p.value,
      odds_ratio = as.numeric(res$estimate),
      or05 = res$conf.int[1], 
      or95 = res$conf.int[2]
    )
    
    enrich_score_ct <- rbind(enrich_score_ct, tmp2)
  }
  
  enrich_score_a <- rbind(enrich_score_ct, enrich_score_a)
}

enrich_score_a <- enrich_score_a %>% 
  mutate(adj_p_val = p.adjust(p_val, method = 'fdr')) %>%
  mutate(logOR = log(odds_ratio)) %>%
  mutate(logOR_lower = log(or05)) %>%
  mutate(logOR_upper = log(or95)) %>%
  mutate(sample_similarity = "Alveolar")

# Save for plotting
saveRDS(enrich_score_a, "/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/proximity_enrichment_30r_a_2025_04_21.rds")

write_csv(enrich_score_a, "/scratch/smallapragada/run2/prox_logOR_csvs/enrichment_scores_a_2025_05_08.csv")

#### RARE DISEASE (CHAOS) ----

meta_rd_c <- meta %>%
  filter(sample_similarity == "Rare disease (CHAOS)")

## Loop to iterate through each CT and calculate stats
enrich_score_rdc <- c()

for(i in id_list) {
  
  celltype_id <- i
  print(celltype_id)
  cell_type_id_out <- gsub("\\/", "\\_", celltype_id)
  
  inpath <-  '/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/output'
  f <- paste0(cell_type_id_out, '_cell_dist_adjdegree.rds')
  df <- readRDS(file.path(inpath, f))
  df <- df[order(df$V1, decreasing = F ),]
  
  # set cell type threshold
  d_threshold <- df %>%
    filter(idx == 1) %>%
    filter(sid == "s1_PDL014" | sid == "s2_PDL014") %>%
    summarise(d_threshold = mean(V1) + 2 * sd(V1)) %>%
    pull(d_threshold)
  
  tmp0 <- df %>%
    filter(V1 < d_threshold & idx == 1) %>%
    filter(sid == "s1_PDL014" | sid == "s2_PDL014") %>%
    distinct(cell.b, celltypeB, sid)
  
  tmp0 <- tmp0[!duplicated(tmp0),]
  tmp0 <- tmp0[,c('cell.b', 'celltypeB', 'sid')] ## rm V1, cell.a, degree, angle, idx
  tmp0 <- tmp0[!duplicated(tmp0),]
  
  # proximal cells by cell type
  proximal_cell_counts <- table(tmp0$celltypeB) %>% 
    as.data.frame() %>%
    rename_with(~ "Freq_neighbor_CT", .cols = "Freq") 
  
  if (nrow(proximal_cell_counts) > 0) {
    proximal_cell_counts$neighbor_CT_count = sum(proximal_cell_counts$Freq_neighbor_CT)
  } else {
    message("No proximal cells found for ", celltype_id)
    proximal_cell_counts <- data.frame(Var1 = character(),
                                       Freq_neighbor_CT = integer(),
                                       neighbor_CT_count = integer())
  }
  
  all_cell_counts <- table(meta_rd_c$CT_final) %>% 
    as.data.frame()
  all_cell_counts$all_cell_count = sum(all_cell_counts$Freq)
  
  tmp1 <- merge(proximal_cell_counts, all_cell_counts, by = 'Var1', all = T)
  tmp1[is.na(tmp1)] = 0
  
  enrich_score_ct <- c()
  
  for (j in 1:nrow(tmp1)) {
    ct <- tmp1[j, 1]
    total <- tmp1[j, 5]
    list1 <- tmp1[j, 4]
    list2 <- tmp1[j, 3]
    overlap <- tmp1[j, 2]
    
    # Safe Fisher matrix calculation
    contingency <- matrix(c(
      total - (list1 + list2) + overlap,
      list1 - overlap,
      list2 - overlap,
      overlap
    ), nrow = 2)
    
    if (any(contingency < 0 | !is.finite(contingency))) {
      res <- list(p.value = NA, estimate = NA)
    } else {
      res <- fisher.test(contingency)
    }
    
    tmp2 <- data.frame(
      source_ct = celltype_id,
      neighboring_ct = ct,
      p_val = res$p.value,
      odds_ratio = as.numeric(res$estimate),
      or05 = res$conf.int[1], 
      or95 = res$conf.int[2]
    )
    
    enrich_score_ct <- rbind(enrich_score_ct, tmp2)
  }
  
  enrich_score_rdc <- rbind(enrich_score_ct, enrich_score_rdc)
}

enrich_score_rdc <- enrich_score_rdc %>% 
  mutate(adj_p_val = p.adjust(p_val, method = 'fdr')) %>%
  mutate(logOR = log(odds_ratio)) %>%
  mutate(logOR_lower = log(or05)) %>%
  mutate(logOR_upper = log(or95)) %>%
  mutate(sample_similarity = "Rare disease (CHAOS)")

# Save for plotting
saveRDS(enrich_score_rdc, "/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/proximity_enrichment_30r_rdc_2025_04_21.rds")

write_csv(enrich_score_rdc, "/scratch/smallapragada/run2/prox_logOR_csvs/enrichment_scores_rdc_2025_05_08.csv")

#### RARE DISEASE (PH + KIDNEYS) ----

meta_rd_ph <- meta %>%
  filter(sample_similarity == "Rare disease (PH + Kidneys)")

## Loop to iterate through each CT and calculate stats
enrich_score_rdph <- c()

for(i in id_list) {
  
  celltype_id <- i
  print(celltype_id)
  cell_type_id_out <- gsub("\\/", "\\_", celltype_id)
  
  inpath <-  '/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/output'
  f <- paste0(cell_type_id_out, '_cell_dist_adjdegree.rds')
  df <- readRDS(file.path(inpath, f))
  df <- df[order(df$V1, decreasing = F ),]
  
  # set cell type threshold
  d_threshold <- df %>%
    filter(idx == 1) %>%
    filter(sid == "s1_PDL015" | sid == "s2_PDL015") %>%
    summarise(d_threshold = mean(V1) + 2 * sd(V1)) %>%
    pull(d_threshold)
  
  tmp0 <- df %>%
    filter(V1 < d_threshold & idx == 1) %>%
    filter(sid == "s1_PDL015" | sid == "s2_PDL015") %>%
    distinct(cell.b, celltypeB, sid)
  
  tmp0 <- tmp0[!duplicated(tmp0),]
  tmp0 <- tmp0[,c('cell.b', 'celltypeB', 'sid')] ## rm V1, cell.a, degree, angle, idx
  tmp0 <- tmp0[!duplicated(tmp0),]
  
  # proximal cells by cell type
  proximal_cell_counts <- table(tmp0$celltypeB) %>% 
    as.data.frame() %>%
    rename_with(~ "Freq_neighbor_CT", .cols = "Freq") 
  
  if (nrow(proximal_cell_counts) > 0) {
    proximal_cell_counts$neighbor_CT_count = sum(proximal_cell_counts$Freq_neighbor_CT)
  } else {
    message("No proximal cells found for ", celltype_id)
    proximal_cell_counts <- data.frame(Var1 = character(),
                                       Freq_neighbor_CT = integer(),
                                       neighbor_CT_count = integer())
  }
  
  all_cell_counts <- table(meta_rd_ph$CT_final) %>% 
    as.data.frame()
  all_cell_counts$all_cell_count = sum(all_cell_counts$Freq)
  
  tmp1 <- merge(proximal_cell_counts, all_cell_counts, by = 'Var1', all = T)
  tmp1[is.na(tmp1)] = 0
  
  enrich_score_ct <- c()
  
  for (j in 1:nrow(tmp1)) {
    ct <- tmp1[j, 1]
    total <- tmp1[j, 5]
    list1 <- tmp1[j, 4]
    list2 <- tmp1[j, 3]
    overlap <- tmp1[j, 2]
    
    # Safe Fisher matrix calculation
    contingency <- matrix(c(
      total - (list1 + list2) + overlap,
      list1 - overlap,
      list2 - overlap,
      overlap
    ), nrow = 2)
    
    if (any(contingency < 0 | !is.finite(contingency))) {
      res <- list(p.value = NA, estimate = NA)
    } else {
      res <- fisher.test(contingency)
    }
    
    tmp2 <- data.frame(
      source_ct = celltype_id,
      neighboring_ct = ct,
      p_val = res$p.value,
      odds_ratio = as.numeric(res$estimate),
      or05 = res$conf.int[1], 
      or95 = res$conf.int[2]
    )
    
    enrich_score_ct <- rbind(enrich_score_ct, tmp2)
  }
  
  enrich_score_rdph <- rbind(enrich_score_ct, enrich_score_rdph)
}

enrich_score_rdph <- enrich_score_rdph %>% 
  mutate(adj_p_val = p.adjust(p_val, method = 'fdr')) %>%
  mutate(logOR = log(odds_ratio)) %>%
  mutate(logOR_lower = log(or05)) %>%
  mutate(logOR_upper = log(or95)) %>%
  mutate(sample_similarity = "Rare disease (PH + Kidneys)")

# Save for plotting
saveRDS(enrich_score_rdph, "/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/proximity_enrichment_30r_rdph_2025_04_21.rds")

write_csv(enrich_score_rdph, "/scratch/smallapragada/run2/prox_logOR_csvs/enrichment_scores_rdph_2025_05_08.csv")

#### RARE DISEASE (INFANT BPD) ----

meta_rd_i <- meta %>%
  filter(sample_similarity == "Rare disease (Infant BPD)")

## Loop to iterate through each CT and calculate stats
enrich_score_rdi <- c()

for(i in id_list) {
  
  celltype_id <- i
  print(celltype_id)
  cell_type_id_out <- gsub("\\/", "\\_", celltype_id)
  
  inpath <-  '/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/output'
  f <- paste0(cell_type_id_out, '_cell_dist_adjdegree.rds')
  df <- readRDS(file.path(inpath, f))
  df <- df[order(df$V1, decreasing = F ),]
  
  # set cell type threshold
  d_threshold <- df %>%
    filter(idx == 1) %>%
    filter(sid == "s1_PDL016" | sid == "s2_PDL016") %>%
    summarise(d_threshold = mean(V1) + 2 * sd(V1)) %>%
    pull(d_threshold)
  
  tmp0 <- df %>%
    filter(V1 < d_threshold & idx == 1) %>%
    filter(sid == "s1_PDL016" | sid == "s2_PDL016") %>%
    distinct(cell.b, celltypeB, sid)
  
  tmp0 <- tmp0[!duplicated(tmp0),]
  tmp0 <- tmp0[,c('cell.b', 'celltypeB', 'sid')] ## rm V1, cell.a, degree, angle, idx
  tmp0 <- tmp0[!duplicated(tmp0),]
  
  # proximal cells by cell type
  proximal_cell_counts <- table(tmp0$celltypeB) %>% 
    as.data.frame() %>%
    rename_with(~ "Freq_neighbor_CT", .cols = "Freq") 
  
  if (nrow(proximal_cell_counts) > 0) {
    proximal_cell_counts$neighbor_CT_count = sum(proximal_cell_counts$Freq_neighbor_CT)
  } else {
    message("No proximal cells found for ", celltype_id)
    proximal_cell_counts <- data.frame(Var1 = character(),
                                       Freq_neighbor_CT = integer(),
                                       neighbor_CT_count = integer())
  }
  
  all_cell_counts <- table(meta_rd_i$CT_final) %>% 
    as.data.frame()
  all_cell_counts$all_cell_count = sum(all_cell_counts$Freq)
  
  tmp1 <- merge(proximal_cell_counts, all_cell_counts, by = 'Var1', all = T)
  tmp1[is.na(tmp1)] = 0
  
  enrich_score_ct <- c()
  
  for (j in 1:nrow(tmp1)) {
    ct <- tmp1[j, 1]
    total <- tmp1[j, 5]
    list1 <- tmp1[j, 4]
    list2 <- tmp1[j, 3]
    overlap <- tmp1[j, 2]
    
    # Safe Fisher matrix calculation
    contingency <- matrix(c(
      total - (list1 + list2) + overlap,
      list1 - overlap,
      list2 - overlap,
      overlap
    ), nrow = 2)
    
    if (any(contingency < 0 | !is.finite(contingency))) {
      res <- list(p.value = NA, estimate = NA)
    } else {
      res <- fisher.test(contingency)
    }
    
    tmp2 <- data.frame(
      source_ct = celltype_id,
      neighboring_ct = ct,
      p_val = res$p.value,
      odds_ratio = as.numeric(res$estimate),
      or05 = res$conf.int[1], 
      or95 = res$conf.int[2]
    )
    
    enrich_score_ct <- rbind(enrich_score_ct, tmp2)
  }
  
  enrich_score_rdi <- rbind(enrich_score_ct, enrich_score_rdi)
}

enrich_score_rdi <- enrich_score_rdi %>% 
  mutate(adj_p_val = p.adjust(p_val, method = 'fdr')) %>%
  mutate(logOR = log(odds_ratio)) %>%
  mutate(logOR_lower = log(or05)) %>%
  mutate(logOR_upper = log(or95)) %>%
  mutate(sample_similarity = "Rare disease (Infant BPD)")

# Save for plotting
saveRDS(enrich_score_rdi, "/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/proximity_enrichment_30r_rdi_2025_04_21.rds")

write_csv(enrich_score_rdi, "/scratch/smallapragada/run2/prox_logOR_csvs/enrichment_scores_rdi_2025_05_08.csv")

#### RARE DISEASE (ADULT BPD) ----

meta_rd_a <- meta %>%
  filter(sample_similarity == "Rare disease (Adult BPD)")

## Loop to iterate through each CT and calculate stats
enrich_score_rda <- c()

for(i in id_list) {
  
  celltype_id <- i
  print(celltype_id)
  cell_type_id_out <- gsub("\\/", "\\_", celltype_id)
  
  inpath <-  '/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/output'
  f <- paste0(cell_type_id_out, '_cell_dist_adjdegree.rds')
  df <- readRDS(file.path(inpath, f))
  df <- df[order(df$V1, decreasing = F ),]
  
  # set cell type threshold
  d_threshold <- df %>%
    filter(idx == 1) %>%
    filter(sid == "s1_PDL017" | sid == "s2_PDL017") %>%
    summarise(d_threshold = mean(V1) + 2 * sd(V1)) %>%
    pull(d_threshold)
  
  tmp0 <- df %>%
    filter(V1 < d_threshold & idx == 1) %>%
    filter(sid == "s1_PDL017" | sid == "s2_PDL017") %>%
    distinct(cell.b, celltypeB, sid)
  
  tmp0 <- tmp0[!duplicated(tmp0),]
  tmp0 <- tmp0[,c('cell.b', 'celltypeB', 'sid')] ## rm V1, cell.a, degree, angle, idx
  tmp0 <- tmp0[!duplicated(tmp0),]
  
  # proximal cells by cell type
  proximal_cell_counts <- table(tmp0$celltypeB) %>% 
    as.data.frame() %>%
    rename_with(~ "Freq_neighbor_CT", .cols = "Freq") 
  
  if (nrow(proximal_cell_counts) > 0) {
    proximal_cell_counts$neighbor_CT_count = sum(proximal_cell_counts$Freq_neighbor_CT)
  } else {
    message("No proximal cells found for ", celltype_id)
    proximal_cell_counts <- data.frame(Var1 = character(),
                                       Freq_neighbor_CT = integer(),
                                       neighbor_CT_count = integer())
  }
  
  all_cell_counts <- table(meta_rd_a$CT_final) %>% 
    as.data.frame()
  all_cell_counts$all_cell_count = sum(all_cell_counts$Freq)
  
  tmp1 <- merge(proximal_cell_counts, all_cell_counts, by = 'Var1', all = T)
  tmp1[is.na(tmp1)] = 0
  
  enrich_score_ct <- c()
  
  for (j in 1:nrow(tmp1)) {
    ct <- tmp1[j, 1]
    total <- tmp1[j, 5]
    list1 <- tmp1[j, 4]
    list2 <- tmp1[j, 3]
    overlap <- tmp1[j, 2]
    
    # Safe Fisher matrix calculation
    contingency <- matrix(c(
      total - (list1 + list2) + overlap,
      list1 - overlap,
      list2 - overlap,
      overlap
    ), nrow = 2)
    
    if (any(contingency < 0 | !is.finite(contingency))) {
      res <- list(p.value = NA, estimate = NA)
    } else {
      res <- fisher.test(contingency)
    }
    
    tmp2 <- data.frame(
      source_ct = celltype_id,
      neighboring_ct = ct,
      p_val = res$p.value,
      odds_ratio = as.numeric(res$estimate),
      or05 = res$conf.int[1], 
      or95 = res$conf.int[2]
    )
    
    enrich_score_ct <- rbind(enrich_score_ct, tmp2)
  }
  
  enrich_score_rda <- rbind(enrich_score_ct, enrich_score_rda)
}

enrich_score_rda <- enrich_score_rda %>% 
  mutate(adj_p_val = p.adjust(p_val, method = 'fdr')) %>%
  mutate(logOR = log(odds_ratio)) %>%
  mutate(logOR_lower = log(or05)) %>%
  mutate(logOR_upper = log(or95)) %>%
  mutate(sample_similarity = "Rare disease (Adult BPD)")

# Save for plotting
saveRDS(enrich_score_rda, "/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/proximity_enrichment_30r_rda_2025_04_21.rds")

write_csv(enrich_score_rda, "/scratch/smallapragada/run2/prox_logOR_csvs/enrichment_scores_rda_2025_05_08.csv")

#### EACH SAMPLE INDIVIDUALLY (OUTSIDE OF RARE DISEASE) ----

### PDL001
meta_001 <- meta %>%
  filter(Sample_combined == "PDL001")

## Loop to iterate through each CT and calculate stats
enrich_score_001 <- c()

for(i in id_list) {
  
  celltype_id <- i
  print(celltype_id)
  cell_type_id_out <- gsub("\\/", "\\_", celltype_id)
  
  inpath <-  '/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/output'
  f <- paste0(cell_type_id_out, '_cell_dist_adjdegree.rds')
  df <- readRDS(file.path(inpath, f))
  df <- df[order(df$V1, decreasing = F ),]
  
  # set cell type threshold
  d_threshold <- df %>%
    filter(idx == 1) %>%
    filter(sid == "s1_PDL001" | sid == "s2_PDL001") %>%
    summarise(d_threshold = mean(V1) + 2 * sd(V1)) %>%
    pull(d_threshold)
  
  tmp0 <- df %>%
    filter(V1 < d_threshold & idx == 1) %>%
    filter(sid == "s1_PDL001" | sid == "s2_PDL001") %>%
    distinct(cell.b, celltypeB, sid)
  
  tmp0 <- tmp0[!duplicated(tmp0),]
  tmp0 <- tmp0[,c('cell.b', 'celltypeB', 'sid')] ## rm V1, cell.a, degree, angle, idx
  tmp0 <- tmp0[!duplicated(tmp0),]
  
  # proximal cells by cell type
  proximal_cell_counts <- table(tmp0$celltypeB) %>% 
    as.data.frame() %>%
    rename_with(~ "Freq_neighbor_CT", .cols = "Freq") 
  
  if (nrow(proximal_cell_counts) > 0) {
    proximal_cell_counts$neighbor_CT_count = sum(proximal_cell_counts$Freq_neighbor_CT)
  } else {
    message("No proximal cells found for ", celltype_id)
    proximal_cell_counts <- data.frame(Var1 = character(),
                                       Freq_neighbor_CT = integer(),
                                       neighbor_CT_count = integer())
  }
  
  all_cell_counts <- table(meta_001$CT_final) %>% 
    as.data.frame()
  all_cell_counts$all_cell_count = sum(all_cell_counts$Freq)
  
  tmp1 <- merge(proximal_cell_counts, all_cell_counts, by = 'Var1', all = T)
  tmp1[is.na(tmp1)] = 0
  
  enrich_score_ct <- c()
  
  for (j in 1:nrow(tmp1)) {
    ct <- tmp1[j, 1]
    total <- tmp1[j, 5]
    list1 <- tmp1[j, 4]
    list2 <- tmp1[j, 3]
    overlap <- tmp1[j, 2]
    
    # Safe Fisher matrix calculation
    contingency <- matrix(c(
      total - (list1 + list2) + overlap,
      list1 - overlap,
      list2 - overlap,
      overlap
    ), nrow = 2)
    
    if (any(contingency < 0 | !is.finite(contingency))) {
      res <- list(p.value = NA, estimate = NA)
    } else {
      res <- fisher.test(contingency)
    }
    
    tmp2 <- data.frame(
      source_ct = celltype_id,
      neighboring_ct = ct,
      p_val = res$p.value,
      odds_ratio = as.numeric(res$estimate),
      or05 = res$conf.int[1], 
      or95 = res$conf.int[2]
    )
    
    enrich_score_ct <- rbind(enrich_score_ct, tmp2)
  }
  
  enrich_score_001 <- rbind(enrich_score_ct, enrich_score_001)
}

enrich_score_001 <- enrich_score_001 %>% 
  mutate(adj_p_val = p.adjust(p_val, method = 'fdr')) %>%
  mutate(logOR = log(odds_ratio)) %>%
  mutate(logOR_lower = log(or05)) %>%
  mutate(logOR_upper = log(or95)) %>%
  mutate(Sample_combined = "PDL001")

# Save for plotting
saveRDS(enrich_score_001, "/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/proximity_enrichment_30r_001_2025_04_25.rds")

write_csv(enrich_score_001, "/scratch/smallapragada/run2/prox_logOR_csvs/enrichment_scores_001_2025_05_08.csv")

### PDL002
meta_002 <- meta %>%
  filter(Sample_combined == "PDL002")

## Loop to iterate through each CT and calculate stats
enrich_score_002 <- c()

for(i in id_list) {
  
  celltype_id <- i
  print(celltype_id)
  cell_type_id_out <- gsub("\\/", "\\_", celltype_id)
  
  inpath <-  '/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/output'
  f <- paste0(cell_type_id_out, '_cell_dist_adjdegree.rds')
  df <- readRDS(file.path(inpath, f))
  df <- df[order(df$V1, decreasing = F ),]
  
  # set cell type threshold
  d_threshold <- df %>%
    filter(idx == 1) %>%
    filter(sid == "s1_PDL002" | sid == "s2_PDL002") %>%
    summarise(d_threshold = mean(V1) + 2 * sd(V1)) %>%
    pull(d_threshold)
  
  tmp0 <- df %>%
    filter(V1 < d_threshold & idx == 1) %>%
    filter(sid == "s1_PDL002" | sid == "s2_PDL002") %>%
    distinct(cell.b, celltypeB, sid)
  
  tmp0 <- tmp0[!duplicated(tmp0),]
  tmp0 <- tmp0[,c('cell.b', 'celltypeB', 'sid')] ## rm V1, cell.a, degree, angle, idx
  tmp0 <- tmp0[!duplicated(tmp0),]
  
  # proximal cells by cell type
  proximal_cell_counts <- table(tmp0$celltypeB) %>% 
    as.data.frame() %>%
    rename_with(~ "Freq_neighbor_CT", .cols = "Freq") 
  
  if (nrow(proximal_cell_counts) > 0) {
    proximal_cell_counts$neighbor_CT_count = sum(proximal_cell_counts$Freq_neighbor_CT)
  } else {
    message("No proximal cells found for ", celltype_id)
    proximal_cell_counts <- data.frame(Var1 = character(),
                                       Freq_neighbor_CT = integer(),
                                       neighbor_CT_count = integer())
  }
  
  all_cell_counts <- table(meta_002$CT_final) %>% 
    as.data.frame()
  all_cell_counts$all_cell_count = sum(all_cell_counts$Freq)
  
  tmp1 <- merge(proximal_cell_counts, all_cell_counts, by = 'Var1', all = T)
  tmp1[is.na(tmp1)] = 0
  
  enrich_score_ct <- c()
  
  for (j in 1:nrow(tmp1)) {
    ct <- tmp1[j, 1]
    total <- tmp1[j, 5]
    list1 <- tmp1[j, 4]
    list2 <- tmp1[j, 3]
    overlap <- tmp1[j, 2]
    
    # Safe Fisher matrix calculation
    contingency <- matrix(c(
      total - (list1 + list2) + overlap,
      list1 - overlap,
      list2 - overlap,
      overlap
    ), nrow = 2)
    
    if (any(contingency < 0 | !is.finite(contingency))) {
      res <- list(p.value = NA, estimate = NA)
    } else {
      res <- fisher.test(contingency)
    }
    
    tmp2 <- data.frame(
      source_ct = celltype_id,
      neighboring_ct = ct,
      p_val = res$p.value,
      odds_ratio = as.numeric(res$estimate),
      or05 = res$conf.int[1], 
      or95 = res$conf.int[2]
    )
    
    enrich_score_ct <- rbind(enrich_score_ct, tmp2)
  }
  
  enrich_score_002 <- rbind(enrich_score_ct, enrich_score_002)
}

enrich_score_002 <- enrich_score_002 %>% 
  mutate(adj_p_val = p.adjust(p_val, method = 'fdr')) %>%
  mutate(logOR = log(odds_ratio)) %>%
  mutate(logOR_lower = log(or05)) %>%
  mutate(logOR_upper = log(or95)) %>%
  mutate(Sample_combined = "PDL002")

# Save for plotting
saveRDS(enrich_score_002, "/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/proximity_enrichment_30r_002_2025_04_25.rds")

write_csv(enrich_score_002, "/scratch/smallapragada/run2/prox_logOR_csvs/enrichment_scores_002_2025_05_08.csv")

### PDL003
meta_PDL003 <- meta %>%
  filter(Sample_combined == "PDL003")

## Loop to iterate through each CT and calculate stats
enrich_score_PDL003 <- c()

for(i in id_list) {
  
  celltype_id <- i
  print(celltype_id)
  cell_type_id_out <- gsub("\\/", "\\_", celltype_id)
  
  inpath <-  '/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/output'
  f <- paste0(cell_type_id_out, '_cell_dist_adjdegree.rds')
  df <- readRDS(file.path(inpath, f))
  df <- df[order(df$V1, decreasing = F ),]
  
  # set cell type threshold
  d_threshold <- df %>%
    filter(idx == 1) %>%
    filter(sid == "s1_PDL003" | sid == "s2_PDL003") %>%
    summarise(d_threshold = mean(V1) + 2 * sd(V1)) %>%
    pull(d_threshold)
  
  tmp0 <- df %>%
    filter(V1 < d_threshold & idx == 1) %>%
    filter(sid == "s1_PDL003" | sid == "s2_PDL003") %>%
    distinct(cell.b, celltypeB, sid)
  
  tmp0 <- tmp0[!duplicated(tmp0),]
  tmp0 <- tmp0[,c('cell.b', 'celltypeB', 'sid')] ## rm V1, cell.a, degree, angle, idx
  tmp0 <- tmp0[!duplicated(tmp0),]
  
  # proximal cells by cell type
  proximal_cell_counts <- table(tmp0$celltypeB) %>% 
    as.data.frame() %>%
    rename_with(~ "Freq_neighbor_CT", .cols = "Freq") 
  
  if (nrow(proximal_cell_counts) > 0) {
    proximal_cell_counts$neighbor_CT_count = sum(proximal_cell_counts$Freq_neighbor_CT)
  } else {
    message("No proximal cells found for ", celltype_id)
    proximal_cell_counts <- data.frame(Var1 = character(),
                                       Freq_neighbor_CT = integer(),
                                       neighbor_CT_count = integer())
  }
  
  all_cell_counts <- table(meta_PDL003$CT_final) %>% 
    as.data.frame()
  all_cell_counts$all_cell_count = sum(all_cell_counts$Freq)
  
  tmp1 <- merge(proximal_cell_counts, all_cell_counts, by = 'Var1', all = T)
  tmp1[is.na(tmp1)] = 0
  
  enrich_score_ct <- c()
  
  for (j in 1:nrow(tmp1)) {
    ct <- tmp1[j, 1]
    total <- tmp1[j, 5]
    list1 <- tmp1[j, 4]
    list2 <- tmp1[j, 3]
    overlap <- tmp1[j, 2]
    
    # Safe Fisher matrix calculation
    contingency <- matrix(c(
      total - (list1 + list2) + overlap,
      list1 - overlap,
      list2 - overlap,
      overlap
    ), nrow = 2)
    
    if (any(contingency < 0 | !is.finite(contingency))) {
      res <- list(p.value = NA, estimate = NA)
    } else {
      res <- fisher.test(contingency)
    }
    
    tmp2 <- data.frame(
      source_ct = celltype_id,
      neighboring_ct = ct,
      p_val = res$p.value,
      odds_ratio = as.numeric(res$estimate),
      or05 = res$conf.int[1], 
      or95 = res$conf.int[2]
    )
    
    enrich_score_ct <- rbind(enrich_score_ct, tmp2)
  }
  
  enrich_score_PDL003 <- rbind(enrich_score_ct, enrich_score_PDL003)
}

enrich_score_PDL003 <- enrich_score_PDL003 %>% 
  mutate(adj_p_val = p.adjust(p_val, method = 'fdr')) %>%
  mutate(logOR = log(odds_ratio)) %>%
  mutate(logOR_lower = log(or05)) %>%
  mutate(logOR_upper = log(or95)) %>%
  mutate(Sample_combined = "PDL003")

# Save for plotting
saveRDS(enrich_score_PDL003, "/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/proximity_enrichment_30r_PDL003_2025_04_25.rds")

write_csv(enrich_score_PDL003, "/scratch/smallapragada/run2/prox_logOR_csvs/enrichment_scores_PDL003_2025_05_08.csv")

### PDL004
meta_PDL004 <- meta %>%
  filter(Sample_combined == "PDL004")

## Loop to iterate through each CT and calculate stats
enrich_score_PDL004 <- c()

for(i in id_list) {
  
  celltype_id <- i
  print(celltype_id)
  cell_type_id_out <- gsub("\\/", "\\_", celltype_id)
  
  inpath <-  '/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/output'
  f <- paste0(cell_type_id_out, '_cell_dist_adjdegree.rds')
  df <- readRDS(file.path(inpath, f))
  df <- df[order(df$V1, decreasing = F ),]
  
  # set cell type threshold
  d_threshold <- df %>%
    filter(idx == 1) %>%
    filter(sid == "s1_PDL004" | sid == "s2_PDL004") %>%
    summarise(d_threshold = mean(V1) + 2 * sd(V1)) %>%
    pull(d_threshold)
  
  tmp0 <- df %>%
    filter(V1 < d_threshold & idx == 1) %>%
    filter(sid == "s1_PDL004" | sid == "s2_PDL004") %>%
    distinct(cell.b, celltypeB, sid)
  
  tmp0 <- tmp0[!duplicated(tmp0),]
  tmp0 <- tmp0[,c('cell.b', 'celltypeB', 'sid')] ## rm V1, cell.a, degree, angle, idx
  tmp0 <- tmp0[!duplicated(tmp0),]
  
  # proximal cells by cell type
  proximal_cell_counts <- table(tmp0$celltypeB) %>% 
    as.data.frame() %>%
    rename_with(~ "Freq_neighbor_CT", .cols = "Freq") 
  
  if (nrow(proximal_cell_counts) > 0) {
    proximal_cell_counts$neighbor_CT_count = sum(proximal_cell_counts$Freq_neighbor_CT)
  } else {
    message("No proximal cells found for ", celltype_id)
    proximal_cell_counts <- data.frame(Var1 = character(),
                                       Freq_neighbor_CT = integer(),
                                       neighbor_CT_count = integer())
  }
  
  all_cell_counts <- table(meta_PDL004$CT_final) %>% 
    as.data.frame()
  all_cell_counts$all_cell_count = sum(all_cell_counts$Freq)
  
  tmp1 <- merge(proximal_cell_counts, all_cell_counts, by = 'Var1', all = T)
  tmp1[is.na(tmp1)] = 0
  
  enrich_score_ct <- c()
  
  for (j in 1:nrow(tmp1)) {
    ct <- tmp1[j, 1]
    total <- tmp1[j, 5]
    list1 <- tmp1[j, 4]
    list2 <- tmp1[j, 3]
    overlap <- tmp1[j, 2]
    
    # Safe Fisher matrix calculation
    contingency <- matrix(c(
      total - (list1 + list2) + overlap,
      list1 - overlap,
      list2 - overlap,
      overlap
    ), nrow = 2)
    
    if (any(contingency < 0 | !is.finite(contingency))) {
      res <- list(p.value = NA, estimate = NA)
    } else {
      res <- fisher.test(contingency)
    }
    
    tmp2 <- data.frame(
      source_ct = celltype_id,
      neighboring_ct = ct,
      p_val = res$p.value,
      odds_ratio = as.numeric(res$estimate),
      or05 = res$conf.int[1], 
      or95 = res$conf.int[2]
    )
    
    enrich_score_ct <- rbind(enrich_score_ct, tmp2)
  }
  
  enrich_score_PDL004 <- rbind(enrich_score_ct, enrich_score_PDL004)
}

enrich_score_PDL004 <- enrich_score_PDL004 %>% 
  mutate(adj_p_val = p.adjust(p_val, method = 'fdr')) %>%
  mutate(logOR = log(odds_ratio)) %>%
  mutate(logOR_lower = log(or05)) %>%
  mutate(logOR_upper = log(or95)) %>%
  mutate(Sample_combined = "PDL004")

# Save for plotting
saveRDS(enrich_score_PDL004, "/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/proximity_enrichment_30r_PDL004_2025_04_25.rds")

write_csv(enrich_score_PDL004, "/scratch/smallapragada/run2/prox_logOR_csvs/enrichment_scores_PDL004_2025_05_08.csv")

### PDL006
meta_PDL006 <- meta %>%
  filter(Sample_combined == "PDL006")

## Loop to iterate through each CT and calculate stats
enrich_score_PDL006 <- c()

for(i in id_list) {
  
  celltype_id <- i
  print(celltype_id)
  cell_type_id_out <- gsub("\\/", "\\_", celltype_id)
  
  inpath <-  '/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/output'
  f <- paste0(cell_type_id_out, '_cell_dist_adjdegree.rds')
  df <- readRDS(file.path(inpath, f))
  df <- df[order(df$V1, decreasing = F ),]
  
  # set cell type threshold
  d_threshold <- df %>%
    filter(idx == 1) %>%
    filter(sid == "s1_PDL006" | sid == "s2_PDL006") %>%
    summarise(d_threshold = mean(V1) + 2 * sd(V1)) %>%
    pull(d_threshold)
  
  tmp0 <- df %>%
    filter(V1 < d_threshold & idx == 1) %>%
    filter(sid == "s1_PDL006" | sid == "s2_PDL006") %>%
    distinct(cell.b, celltypeB, sid)
  
  tmp0 <- tmp0[!duplicated(tmp0),]
  tmp0 <- tmp0[,c('cell.b', 'celltypeB', 'sid')] ## rm V1, cell.a, degree, angle, idx
  tmp0 <- tmp0[!duplicated(tmp0),]
  
  # proximal cells by cell type
  proximal_cell_counts <- table(tmp0$celltypeB) %>% 
    as.data.frame() %>%
    rename_with(~ "Freq_neighbor_CT", .cols = "Freq") 
  
  if (nrow(proximal_cell_counts) > 0) {
    proximal_cell_counts$neighbor_CT_count = sum(proximal_cell_counts$Freq_neighbor_CT)
  } else {
    message("No proximal cells found for ", celltype_id)
    proximal_cell_counts <- data.frame(Var1 = character(),
                                       Freq_neighbor_CT = integer(),
                                       neighbor_CT_count = integer())
  }
  
  all_cell_counts <- table(meta_PDL006$CT_final) %>% 
    as.data.frame()
  all_cell_counts$all_cell_count = sum(all_cell_counts$Freq)
  
  tmp1 <- merge(proximal_cell_counts, all_cell_counts, by = 'Var1', all = T)
  tmp1[is.na(tmp1)] = 0
  
  enrich_score_ct <- c()
  
  for (j in 1:nrow(tmp1)) {
    ct <- tmp1[j, 1]
    total <- tmp1[j, 5]
    list1 <- tmp1[j, 4]
    list2 <- tmp1[j, 3]
    overlap <- tmp1[j, 2]
    
    # Safe Fisher matrix calculation
    contingency <- matrix(c(
      total - (list1 + list2) + overlap,
      list1 - overlap,
      list2 - overlap,
      overlap
    ), nrow = 2)
    
    if (any(contingency < 0 | !is.finite(contingency))) {
      res <- list(p.value = NA, estimate = NA)
    } else {
      res <- fisher.test(contingency)
    }
    
    tmp2 <- data.frame(
      source_ct = celltype_id,
      neighboring_ct = ct,
      p_val = res$p.value,
      odds_ratio = as.numeric(res$estimate),
      or05 = res$conf.int[1], 
      or95 = res$conf.int[2]
    )
    
    enrich_score_ct <- rbind(enrich_score_ct, tmp2)
  }
  
  enrich_score_PDL006 <- rbind(enrich_score_ct, enrich_score_PDL006)
}

enrich_score_PDL006 <- enrich_score_PDL006 %>% 
  mutate(adj_p_val = p.adjust(p_val, method = 'fdr')) %>%
  mutate(logOR = log(odds_ratio)) %>%
  mutate(logOR_lower = log(or05)) %>%
  mutate(logOR_upper = log(or95)) %>%
  mutate(Sample_combined = "PDL006")

# Save for plotting
saveRDS(enrich_score_PDL006, "/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/proximity_enrichment_30r_PDL006_2025_04_25.rds")

write_csv(enrich_score_PDL006, "/scratch/smallapragada/run2/prox_logOR_csvs/enrichment_scores_PDL006_2025_05_08.csv")

### PDL005 
meta_PDL005 <- meta %>%
  filter(Sample_combined == "PDL005")

## Loop to iterate through each CT and calculate stats
enrich_score_PDL005 <- c()

for(i in id_list) {
  
  celltype_id <- i
  print(celltype_id)
  cell_type_id_out <- gsub("\\/", "\\_", celltype_id)
  
  inpath <-  '/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/output'
  f <- paste0(cell_type_id_out, '_cell_dist_adjdegree.rds')
  df <- readRDS(file.path(inpath, f))
  df <- df[order(df$V1, decreasing = F ),]
  
  # set cell type threshold
  d_threshold <- df %>%
    filter(idx == 1) %>%
    filter(sid == "s1_PDL005" | sid == "s2_PDL005") %>%
    summarise(d_threshold = mean(V1) + 2 * sd(V1)) %>%
    pull(d_threshold)
  
  tmp0 <- df %>%
    filter(V1 < d_threshold & idx == 1) %>%
    filter(sid == "s1_PDL005" | sid == "s2_PDL005") %>%
    distinct(cell.b, celltypeB, sid)
  
  tmp0 <- tmp0[!duplicated(tmp0),]
  tmp0 <- tmp0[,c('cell.b', 'celltypeB', 'sid')] ## rm V1, cell.a, degree, angle, idx
  tmp0 <- tmp0[!duplicated(tmp0),]
  
  # proximal cells by cell type
  proximal_cell_counts <- table(tmp0$celltypeB) %>% 
    as.data.frame() %>%
    rename_with(~ "Freq_neighbor_CT", .cols = "Freq") 
  
  if (nrow(proximal_cell_counts) > 0) {
    proximal_cell_counts$neighbor_CT_count = sum(proximal_cell_counts$Freq_neighbor_CT)
  } else {
    message("No proximal cells found for ", celltype_id)
    proximal_cell_counts <- data.frame(Var1 = character(),
                                       Freq_neighbor_CT = integer(),
                                       neighbor_CT_count = integer())
  }
  
  all_cell_counts <- table(meta_PDL005$CT_final) %>% 
    as.data.frame()
  all_cell_counts$all_cell_count = sum(all_cell_counts$Freq)
  
  tmp1 <- merge(proximal_cell_counts, all_cell_counts, by = 'Var1', all = T)
  tmp1[is.na(tmp1)] = 0
  
  enrich_score_ct <- c()
  
  for (j in 1:nrow(tmp1)) {
    ct <- tmp1[j, 1]
    total <- tmp1[j, 5]
    list1 <- tmp1[j, 4]
    list2 <- tmp1[j, 3]
    overlap <- tmp1[j, 2]
    
    # Safe Fisher matrix calculation
    contingency <- matrix(c(
      total - (list1 + list2) + overlap,
      list1 - overlap,
      list2 - overlap,
      overlap
    ), nrow = 2)
    
    if (any(contingency < 0 | !is.finite(contingency))) {
      res <- list(p.value = NA, estimate = NA)
    } else {
      res <- fisher.test(contingency)
    }
    
    tmp2 <- data.frame(
      source_ct = celltype_id,
      neighboring_ct = ct,
      p_val = res$p.value,
      odds_ratio = as.numeric(res$estimate),
      or05 = res$conf.int[1], 
      or95 = res$conf.int[2]
    )
    
    enrich_score_ct <- rbind(enrich_score_ct, tmp2)
  }
  
  enrich_score_PDL005 <- rbind(enrich_score_ct, enrich_score_PDL005)
}

enrich_score_PDL005 <- enrich_score_PDL005 %>% 
  mutate(adj_p_val = p.adjust(p_val, method = 'fdr')) %>%
  mutate(logOR = log(odds_ratio)) %>%
  mutate(logOR_lower = log(or05)) %>%
  mutate(logOR_upper = log(or95)) %>%
  mutate(Sample_combined = "PDL005")

# Save for plotting
saveRDS(enrich_score_PDL005, "/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/proximity_enrichment_30r_PDL005_2025_04_25.rds")

write_csv(enrich_score_PDL005, "/scratch/smallapragada/run2/prox_logOR_csvs/enrichment_scores_PDL005_2025_05_08.csv")

### PDL007
meta_PDL007 <- meta %>%
  filter(Sample_combined == "PDL007")

## Loop to iterate through each CT and calculate stats
enrich_score_PDL007 <- c()

for(i in id_list) {
  
  celltype_id <- i
  print(celltype_id)
  cell_type_id_out <- gsub("\\/", "\\_", celltype_id)
  
  inpath <-  '/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/output'
  f <- paste0(cell_type_id_out, '_cell_dist_adjdegree.rds')
  df <- readRDS(file.path(inpath, f))
  df <- df[order(df$V1, decreasing = F ),]
  
  # set cell type threshold
  d_threshold <- df %>%
    filter(idx == 1) %>%
    filter(sid == "s1_PDL007" | sid == "s2_PDL007") %>%
    summarise(d_threshold = mean(V1) + 2 * sd(V1)) %>%
    pull(d_threshold)
  
  tmp0 <- df %>%
    filter(V1 < d_threshold & idx == 1) %>%
    filter(sid == "s1_PDL007" | sid == "s2_PDL007") %>%
    distinct(cell.b, celltypeB, sid)
  
  tmp0 <- tmp0[!duplicated(tmp0),]
  tmp0 <- tmp0[,c('cell.b', 'celltypeB', 'sid')] ## rm V1, cell.a, degree, angle, idx
  tmp0 <- tmp0[!duplicated(tmp0),]
  
  # proximal cells by cell type
  proximal_cell_counts <- table(tmp0$celltypeB) %>% 
    as.data.frame() %>%
    rename_with(~ "Freq_neighbor_CT", .cols = "Freq") 
  
  if (nrow(proximal_cell_counts) > 0) {
    proximal_cell_counts$neighbor_CT_count = sum(proximal_cell_counts$Freq_neighbor_CT)
  } else {
    message("No proximal cells found for ", celltype_id)
    proximal_cell_counts <- data.frame(Var1 = character(),
                                       Freq_neighbor_CT = integer(),
                                       neighbor_CT_count = integer())
  }
  
  all_cell_counts <- table(meta_PDL007$CT_final) %>% 
    as.data.frame()
  all_cell_counts$all_cell_count = sum(all_cell_counts$Freq)
  
  tmp1 <- merge(proximal_cell_counts, all_cell_counts, by = 'Var1', all = T)
  tmp1[is.na(tmp1)] = 0
  
  enrich_score_ct <- c()
  
  for (j in 1:nrow(tmp1)) {
    ct <- tmp1[j, 1]
    total <- tmp1[j, 5]
    list1 <- tmp1[j, 4]
    list2 <- tmp1[j, 3]
    overlap <- tmp1[j, 2]
    
    # Safe Fisher matrix calculation
    contingency <- matrix(c(
      total - (list1 + list2) + overlap,
      list1 - overlap,
      list2 - overlap,
      overlap
    ), nrow = 2)
    
    if (any(contingency < 0 | !is.finite(contingency))) {
      res <- list(p.value = NA, estimate = NA)
    } else {
      res <- fisher.test(contingency)
    }
    
    tmp2 <- data.frame(
      source_ct = celltype_id,
      neighboring_ct = ct,
      p_val = res$p.value,
      odds_ratio = as.numeric(res$estimate),
      or05 = res$conf.int[1], 
      or95 = res$conf.int[2]
    )
    
    enrich_score_ct <- rbind(enrich_score_ct, tmp2)
  }
  
  enrich_score_PDL007 <- rbind(enrich_score_ct, enrich_score_PDL007)
}

enrich_score_PDL007 <- enrich_score_PDL007 %>% 
  mutate(adj_p_val = p.adjust(p_val, method = 'fdr')) %>%
  mutate(logOR = log(odds_ratio)) %>%
  mutate(logOR_lower = log(or05)) %>%
  mutate(logOR_upper = log(or95)) %>%
  mutate(Sample_combined = "PDL007")

# Save for plotting
saveRDS(enrich_score_PDL007, "/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/proximity_enrichment_30r_PDL007_2025_04_25.rds")

write_csv(enrich_score_PDL007, "/scratch/smallapragada/run2/prox_logOR_csvs/enrichment_scores_PDL007_2025_05_08.csv")

### PDL008
meta_PDL008 <- meta %>%
  filter(Sample_combined == "PDL008")

## Loop to iterate through each CT and calculate stats
enrich_score_PDL008 <- c()

for(i in id_list) {
  
  celltype_id <- i
  print(celltype_id)
  cell_type_id_out <- gsub("\\/", "\\_", celltype_id)
  
  inpath <-  '/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/output'
  f <- paste0(cell_type_id_out, '_cell_dist_adjdegree.rds')
  df <- readRDS(file.path(inpath, f))
  df <- df[order(df$V1, decreasing = F ),]
  
  # set cell type threshold
  d_threshold <- df %>%
    filter(idx == 1) %>%
    filter(sid == "s1_PDL008" | sid == "s2_PDL008") %>%
    summarise(d_threshold = mean(V1) + 2 * sd(V1)) %>%
    pull(d_threshold)
  
  tmp0 <- df %>%
    filter(V1 < d_threshold & idx == 1) %>%
    filter(sid == "s1_PDL008" | sid == "s2_PDL008") %>%
    distinct(cell.b, celltypeB, sid)
  
  tmp0 <- tmp0[!duplicated(tmp0),]
  tmp0 <- tmp0[,c('cell.b', 'celltypeB', 'sid')] ## rm V1, cell.a, degree, angle, idx
  tmp0 <- tmp0[!duplicated(tmp0),]
  
  # proximal cells by cell type
  proximal_cell_counts <- table(tmp0$celltypeB) %>% 
    as.data.frame() %>%
    rename_with(~ "Freq_neighbor_CT", .cols = "Freq") 
  
  if (nrow(proximal_cell_counts) > 0) {
    proximal_cell_counts$neighbor_CT_count = sum(proximal_cell_counts$Freq_neighbor_CT)
  } else {
    message("No proximal cells found for ", celltype_id)
    proximal_cell_counts <- data.frame(Var1 = character(),
                                       Freq_neighbor_CT = integer(),
                                       neighbor_CT_count = integer())
  }
  
  all_cell_counts <- table(meta_PDL008$CT_final) %>% 
    as.data.frame()
  all_cell_counts$all_cell_count = sum(all_cell_counts$Freq)
  
  tmp1 <- merge(proximal_cell_counts, all_cell_counts, by = 'Var1', all = T)
  tmp1[is.na(tmp1)] = 0
  
  enrich_score_ct <- c()
  
  for (j in 1:nrow(tmp1)) {
    ct <- tmp1[j, 1]
    total <- tmp1[j, 5]
    list1 <- tmp1[j, 4]
    list2 <- tmp1[j, 3]
    overlap <- tmp1[j, 2]
    
    # Safe Fisher matrix calculation
    contingency <- matrix(c(
      total - (list1 + list2) + overlap,
      list1 - overlap,
      list2 - overlap,
      overlap
    ), nrow = 2)
    
    if (any(contingency < 0 | !is.finite(contingency))) {
      res <- list(p.value = NA, estimate = NA)
    } else {
      res <- fisher.test(contingency)
    }
    
    tmp2 <- data.frame(
      source_ct = celltype_id,
      neighboring_ct = ct,
      p_val = res$p.value,
      odds_ratio = as.numeric(res$estimate),
      or05 = res$conf.int[1], 
      or95 = res$conf.int[2]
    )
    
    enrich_score_ct <- rbind(enrich_score_ct, tmp2)
  }
  
  enrich_score_PDL008 <- rbind(enrich_score_ct, enrich_score_PDL008)
}

enrich_score_PDL008 <- enrich_score_PDL008 %>% 
  mutate(adj_p_val = p.adjust(p_val, method = 'fdr')) %>%
  mutate(logOR = log(odds_ratio)) %>%
  mutate(logOR_lower = log(or05)) %>%
  mutate(logOR_upper = log(or95)) %>%
  mutate(Sample_combined = "PDL008")

# Save for plotting
saveRDS(enrich_score_PDL008, "/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/proximity_enrichment_30r_PDL008_2025_04_25.rds")

write_csv(enrich_score_PDL008, "/scratch/smallapragada/run2/prox_logOR_csvs/enrichment_scores_PDL008_2025_05_08.csv")

### PDL009
meta_PDL009 <- meta %>%
  filter(Sample_combined == "PDL009")

## Loop to iterate through each CT and calculate stats
enrich_score_PDL009 <- c()

for(i in id_list) {
  
  celltype_id <- i
  print(celltype_id)
  cell_type_id_out <- gsub("\\/", "\\_", celltype_id)
  
  inpath <-  '/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/output'
  f <- paste0(cell_type_id_out, '_cell_dist_adjdegree.rds')
  df <- readRDS(file.path(inpath, f))
  df <- df[order(df$V1, decreasing = F ),]
  
  # set cell type threshold
  d_threshold <- df %>%
    filter(idx == 1) %>%
    filter(sid == "s1_PDL009" | sid == "s2_PDL009") %>%
    summarise(d_threshold = mean(V1) + 2 * sd(V1)) %>%
    pull(d_threshold)
  
  tmp0 <- df %>%
    filter(V1 < d_threshold & idx == 1) %>%
    filter(sid == "s1_PDL009" | sid == "s2_PDL009") %>%
    distinct(cell.b, celltypeB, sid)
  
  tmp0 <- tmp0[!duplicated(tmp0),]
  tmp0 <- tmp0[,c('cell.b', 'celltypeB', 'sid')] ## rm V1, cell.a, degree, angle, idx
  tmp0 <- tmp0[!duplicated(tmp0),]
  
  # proximal cells by cell type
  proximal_cell_counts <- table(tmp0$celltypeB) %>% 
    as.data.frame() %>%
    rename_with(~ "Freq_neighbor_CT", .cols = "Freq") 
  
  if (nrow(proximal_cell_counts) > 0) {
    proximal_cell_counts$neighbor_CT_count = sum(proximal_cell_counts$Freq_neighbor_CT)
  } else {
    message("No proximal cells found for ", celltype_id)
    proximal_cell_counts <- data.frame(Var1 = character(),
                                       Freq_neighbor_CT = integer(),
                                       neighbor_CT_count = integer())
  }
  
  all_cell_counts <- table(meta_PDL009$CT_final) %>% 
    as.data.frame()
  all_cell_counts$all_cell_count = sum(all_cell_counts$Freq)
  
  tmp1 <- merge(proximal_cell_counts, all_cell_counts, by = 'Var1', all = T)
  tmp1[is.na(tmp1)] = 0
  
  enrich_score_ct <- c()
  
  for (j in 1:nrow(tmp1)) {
    ct <- tmp1[j, 1]
    total <- tmp1[j, 5]
    list1 <- tmp1[j, 4]
    list2 <- tmp1[j, 3]
    overlap <- tmp1[j, 2]
    
    # Safe Fisher matrix calculation
    contingency <- matrix(c(
      total - (list1 + list2) + overlap,
      list1 - overlap,
      list2 - overlap,
      overlap
    ), nrow = 2)
    
    if (any(contingency < 0 | !is.finite(contingency))) {
      res <- list(p.value = NA, estimate = NA)
    } else {
      res <- fisher.test(contingency)
    }
    
    tmp2 <- data.frame(
      source_ct = celltype_id,
      neighboring_ct = ct,
      p_val = res$p.value,
      odds_ratio = as.numeric(res$estimate),
      or05 = res$conf.int[1], 
      or95 = res$conf.int[2]
    )
    
    enrich_score_ct <- rbind(enrich_score_ct, tmp2)
  }
  
  enrich_score_PDL009 <- rbind(enrich_score_ct, enrich_score_PDL009)
}

enrich_score_PDL009 <- enrich_score_PDL009 %>% 
  mutate(adj_p_val = p.adjust(p_val, method = 'fdr')) %>%
  mutate(logOR = log(odds_ratio)) %>%
  mutate(logOR_lower = log(or05)) %>%
  mutate(logOR_upper = log(or95)) %>%
  mutate(Sample_combined = "PDL009")

# Save for plotting
saveRDS(enrich_score_PDL009, "/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/proximity_enrichment_30r_PDL009_2025_04_25.rds")

write_csv(enrich_score_PDL009, "/scratch/smallapragada/run2/prox_logOR_csvs/enrichment_scores_PDL009_2025_05_08.csv")

### PDL010
meta_bl091922 <- meta %>%
  filter(Sample_combined == "PDL010")

## Loop to iterate through each CT and calculate stats
enrich_score_bl091922 <- c()

for(i in id_list) {
  
  celltype_id <- i
  print(celltype_id)
  cell_type_id_out <- gsub("\\/", "\\_", celltype_id)
  
  inpath <-  '/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/output'
  f <- paste0(cell_type_id_out, '_cell_dist_adjdegree.rds')
  df <- readRDS(file.path(inpath, f))
  df <- df[order(df$V1, decreasing = F ),]
  
  # set cell type threshold
  d_threshold <- df %>%
    filter(idx == 1) %>%
    filter(sid == "s1_PDL010" | sid == "s2_PDL010") %>%
    summarise(d_threshold = mean(V1) + 2 * sd(V1)) %>%
    pull(d_threshold)
  
  tmp0 <- df %>%
    filter(V1 < d_threshold & idx == 1) %>%
    filter(sid == "s1_PDL010" | sid == "s2_PDL010") %>%
    distinct(cell.b, celltypeB, sid)
  
  tmp0 <- tmp0[!duplicated(tmp0),]
  tmp0 <- tmp0[,c('cell.b', 'celltypeB', 'sid')] ## rm V1, cell.a, degree, angle, idx
  tmp0 <- tmp0[!duplicated(tmp0),]
  
  # proximal cells by cell type
  proximal_cell_counts <- table(tmp0$celltypeB) %>% 
    as.data.frame() %>%
    rename_with(~ "Freq_neighbor_CT", .cols = "Freq") 
  
  if (nrow(proximal_cell_counts) > 0) {
    proximal_cell_counts$neighbor_CT_count = sum(proximal_cell_counts$Freq_neighbor_CT)
  } else {
    message("No proximal cells found for ", celltype_id)
    proximal_cell_counts <- data.frame(Var1 = character(),
                                       Freq_neighbor_CT = integer(),
                                       neighbor_CT_count = integer())
  }
  
  all_cell_counts <- table(meta_bl091922$CT_final) %>% 
    as.data.frame()
  all_cell_counts$all_cell_count = sum(all_cell_counts$Freq)
  
  tmp1 <- merge(proximal_cell_counts, all_cell_counts, by = 'Var1', all = T)
  tmp1[is.na(tmp1)] = 0
  
  enrich_score_ct <- c()
  
  for (j in 1:nrow(tmp1)) {
    ct <- tmp1[j, 1]
    total <- tmp1[j, 5]
    list1 <- tmp1[j, 4]
    list2 <- tmp1[j, 3]
    overlap <- tmp1[j, 2]
    
    # Safe Fisher matrix calculation
    contingency <- matrix(c(
      total - (list1 + list2) + overlap,
      list1 - overlap,
      list2 - overlap,
      overlap
    ), nrow = 2)
    
    if (any(contingency < 0 | !is.finite(contingency))) {
      res <- list(p.value = NA, estimate = NA)
    } else {
      res <- fisher.test(contingency)
    }
    
    tmp2 <- data.frame(
      source_ct = celltype_id,
      neighboring_ct = ct,
      p_val = res$p.value,
      odds_ratio = as.numeric(res$estimate),
      or05 = res$conf.int[1], 
      or95 = res$conf.int[2]
    )
    
    enrich_score_ct <- rbind(enrich_score_ct, tmp2)
  }
  
  enrich_score_bl091922 <- rbind(enrich_score_ct, enrich_score_bl091922)
}

enrich_score_bl091922 <- enrich_score_bl091922 %>% 
  mutate(adj_p_val = p.adjust(p_val, method = 'fdr')) %>%
  mutate(logOR = log(odds_ratio)) %>%
  mutate(logOR_lower = log(or05)) %>%
  mutate(logOR_upper = log(or95)) %>%
  mutate(Sample_combined = "PDL010")

# Save for plotting
saveRDS(enrich_score_bl091922, "/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/proximity_enrichment_30r_bl091922_2025_04_25.rds")

write_csv(enrich_score_bl091922, "/scratch/smallapragada/run2/prox_logOR_csvs/enrichment_scores_bl091922_2025_05_08.csv")

### PDL012
meta_PDL012 <- meta %>%
  filter(Sample_combined == "PDL012")

## Loop to iterate through each CT and calculate stats
enrich_score_PDL012 <- c()

for(i in id_list) {
  
  celltype_id <- i
  print(celltype_id)
  cell_type_id_out <- gsub("\\/", "\\_", celltype_id)
  
  inpath <-  '/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/output'
  f <- paste0(cell_type_id_out, '_cell_dist_adjdegree.rds')
  df <- readRDS(file.path(inpath, f))
  df <- df[order(df$V1, decreasing = F ),]
  
  # set cell type threshold
  d_threshold <- df %>%
    filter(idx == 1) %>%
    filter(sid == "s1_PDL012" | sid == "s2_PDL012") %>%
    summarise(d_threshold = mean(V1) + 2 * sd(V1)) %>%
    pull(d_threshold)
  
  tmp0 <- df %>%
    filter(V1 < d_threshold & idx == 1) %>%
    filter(sid == "s1_PDL012" | sid == "s2_BL02482023") %>%
    distinct(cell.b, celltypeB, sid)
  
  tmp0 <- tmp0[!duplicated(tmp0),]
  tmp0 <- tmp0[,c('cell.b', 'celltypeB', 'sid')] ## rm V1, cell.a, degree, angle, idx
  tmp0 <- tmp0[!duplicated(tmp0),]
  
  # proximal cells by cell type
  proximal_cell_counts <- table(tmp0$celltypeB) %>% 
    as.data.frame() %>%
    rename_with(~ "Freq_neighbor_CT", .cols = "Freq") 
  
  if (nrow(proximal_cell_counts) > 0) {
    proximal_cell_counts$neighbor_CT_count = sum(proximal_cell_counts$Freq_neighbor_CT)
  } else {
    message("No proximal cells found for ", celltype_id)
    proximal_cell_counts <- data.frame(Var1 = character(),
                                       Freq_neighbor_CT = integer(),
                                       neighbor_CT_count = integer())
  }
  
  all_cell_counts <- table(meta_PDL012$CT_final) %>% 
    as.data.frame()
  all_cell_counts$all_cell_count = sum(all_cell_counts$Freq)
  
  tmp1 <- merge(proximal_cell_counts, all_cell_counts, by = 'Var1', all = T)
  tmp1[is.na(tmp1)] = 0
  
  enrich_score_ct <- c()
  
  for (j in 1:nrow(tmp1)) {
    ct <- tmp1[j, 1]
    total <- tmp1[j, 5]
    list1 <- tmp1[j, 4]
    list2 <- tmp1[j, 3]
    overlap <- tmp1[j, 2]
    
    # Safe Fisher matrix calculation
    contingency <- matrix(c(
      total - (list1 + list2) + overlap,
      list1 - overlap,
      list2 - overlap,
      overlap
    ), nrow = 2)
    
    if (any(contingency < 0 | !is.finite(contingency))) {
      res <- list(p.value = NA, estimate = NA)
    } else {
      res <- fisher.test(contingency)
    }
    
    tmp2 <- data.frame(
      source_ct = celltype_id,
      neighboring_ct = ct,
      p_val = res$p.value,
      odds_ratio = as.numeric(res$estimate),
      or05 = res$conf.int[1], 
      or95 = res$conf.int[2]
    )
    
    enrich_score_ct <- rbind(enrich_score_ct, tmp2)
  }
  
  enrich_score_PDL012 <- rbind(enrich_score_ct, enrich_score_PDL012)
}

enrich_score_PDL012 <- enrich_score_PDL012 %>% 
  mutate(adj_p_val = p.adjust(p_val, method = 'fdr')) %>%
  mutate(logOR = log(odds_ratio)) %>%
  mutate(logOR_lower = log(or05)) %>%
  mutate(logOR_upper = log(or95)) %>%
  mutate(Sample_combined = "PDL012")

# Save for plotting
saveRDS(enrich_score_PDL012, "/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/proximity_enrichment_30r_PDL012_2025_04_25.rds")

write_csv(enrich_score_PDL012, "/scratch/smallapragada/run2/prox_logOR_csvs/enrichment_scores_PDL012_2025_05_08.csv")

### PDL013
meta_PDL013 <- meta %>%
  filter(Sample_combined == "PDL013")

## Loop to iterate through each CT and calculate stats
enrich_score_PDL013 <- c()

for(i in id_list) {
  
  celltype_id <- i
  print(celltype_id)
  cell_type_id_out <- gsub("\\/", "\\_", celltype_id)
  
  inpath <-  '/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/output'
  f <- paste0(cell_type_id_out, '_cell_dist_adjdegree.rds')
  df <- readRDS(file.path(inpath, f))
  df <- df[order(df$V1, decreasing = F ),]
  
  # set cell type threshold
  d_threshold <- df %>%
    filter(idx == 1) %>%
    filter(sid == "s1_PDL013" | sid == "s2_PDL013") %>%
    summarise(d_threshold = mean(V1) + 2 * sd(V1)) %>%
    pull(d_threshold)
  
  tmp0 <- df %>%
    filter(V1 < d_threshold & idx == 1) %>%
    filter(sid == "s1_PDL013" | sid == "s2_PDL013") %>%
    distinct(cell.b, celltypeB, sid)
  
  tmp0 <- tmp0[!duplicated(tmp0),]
  tmp0 <- tmp0[,c('cell.b', 'celltypeB', 'sid')] ## rm V1, cell.a, degree, angle, idx
  tmp0 <- tmp0[!duplicated(tmp0),]
  
  # proximal cells by cell type
  proximal_cell_counts <- table(tmp0$celltypeB) %>% 
    as.data.frame() %>%
    rename_with(~ "Freq_neighbor_CT", .cols = "Freq") 
  
  if (nrow(proximal_cell_counts) > 0) {
    proximal_cell_counts$neighbor_CT_count = sum(proximal_cell_counts$Freq_neighbor_CT)
  } else {
    message("No proximal cells found for ", celltype_id)
    proximal_cell_counts <- data.frame(Var1 = character(),
                                       Freq_neighbor_CT = integer(),
                                       neighbor_CT_count = integer())
  }
  
  all_cell_counts <- table(meta_PDL013$CT_final) %>% 
    as.data.frame()
  all_cell_counts$all_cell_count = sum(all_cell_counts$Freq)
  
  tmp1 <- merge(proximal_cell_counts, all_cell_counts, by = 'Var1', all = T)
  tmp1[is.na(tmp1)] = 0
  
  enrich_score_ct <- c()
  
  for (j in 1:nrow(tmp1)) {
    ct <- tmp1[j, 1]
    total <- tmp1[j, 5]
    list1 <- tmp1[j, 4]
    list2 <- tmp1[j, 3]
    overlap <- tmp1[j, 2]
    
    # Safe Fisher matrix calculation
    contingency <- matrix(c(
      total - (list1 + list2) + overlap,
      list1 - overlap,
      list2 - overlap,
      overlap
    ), nrow = 2)
    
    if (any(contingency < 0 | !is.finite(contingency))) {
      res <- list(p.value = NA, estimate = NA)
    } else {
      res <- fisher.test(contingency)
    }
    
    tmp2 <- data.frame(
      source_ct = celltype_id,
      neighboring_ct = ct,
      p_val = res$p.value,
      odds_ratio = as.numeric(res$estimate),
      or05 = res$conf.int[1], 
      or95 = res$conf.int[2]
    )
    
    enrich_score_ct <- rbind(enrich_score_ct, tmp2)
  }
  
  enrich_score_PDL013 <- rbind(enrich_score_ct, enrich_score_PDL013)
}

enrich_score_PDL013 <- enrich_score_PDL013 %>% 
  mutate(adj_p_val = p.adjust(p_val, method = 'fdr')) %>%
  mutate(logOR = log(odds_ratio)) %>%
  mutate(logOR_lower = log(or05)) %>%
  mutate(logOR_upper = log(or95)) %>%
  mutate(Sample_combined = "PDL013")

# Save for plotting
saveRDS(enrich_score_PDL013, "/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/proximity_enrichment_30r_PDL013_2025_04_25.rds")

write_csv(enrich_score_PDL013, "/scratch/smallapragada/run2/prox_logOR_csvs/enrichment_scores_PDL013_2025_05_08.csv")

### PDL013
meta_PDL011 <- meta %>%
  filter(Sample_combined == "PDL011")

## Loop to iterate through each CT and calculate stats
enrich_score_PDL011 <- c()

for(i in id_list) {
  
  celltype_id <- i
  print(celltype_id)
  cell_type_id_out <- gsub("\\/", "\\_", celltype_id)
  
  inpath <-  '/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/output'
  f <- paste0(cell_type_id_out, '_cell_dist_adjdegree.rds')
  df <- readRDS(file.path(inpath, f))
  df <- df[order(df$V1, decreasing = F ),]
  
  # set cell type threshold
  d_threshold <- df %>%
    filter(idx == 1) %>%
    filter(sid == "s1_PDL011" | sid == "s2_PDL011") %>%
    summarise(d_threshold = mean(V1) + 2 * sd(V1)) %>%
    pull(d_threshold)
  
  tmp0 <- df %>%
    filter(V1 < d_threshold & idx == 1) %>%
    filter(sid == "s1_PDL011" | sid == "s2_PDL011") %>%
    distinct(cell.b, celltypeB, sid)
  
  tmp0 <- tmp0[!duplicated(tmp0),]
  tmp0 <- tmp0[,c('cell.b', 'celltypeB', 'sid')] ## rm V1, cell.a, degree, angle, idx
  tmp0 <- tmp0[!duplicated(tmp0),]
  
  # proximal cells by cell type
  proximal_cell_counts <- table(tmp0$celltypeB) %>% 
    as.data.frame() %>%
    rename_with(~ "Freq_neighbor_CT", .cols = "Freq") 
  
  if (nrow(proximal_cell_counts) > 0) {
    proximal_cell_counts$neighbor_CT_count = sum(proximal_cell_counts$Freq_neighbor_CT)
  } else {
    message("No proximal cells found for ", celltype_id)
    proximal_cell_counts <- data.frame(Var1 = character(),
                                       Freq_neighbor_CT = integer(),
                                       neighbor_CT_count = integer())
  }
  
  all_cell_counts <- table(meta_PDL011$CT_final) %>% 
    as.data.frame()
  all_cell_counts$all_cell_count = sum(all_cell_counts$Freq)
  
  tmp1 <- merge(proximal_cell_counts, all_cell_counts, by = 'Var1', all = T)
  tmp1[is.na(tmp1)] = 0
  
  enrich_score_ct <- c()
  
  for (j in 1:nrow(tmp1)) {
    ct <- tmp1[j, 1]
    total <- tmp1[j, 5]
    list1 <- tmp1[j, 4]
    list2 <- tmp1[j, 3]
    overlap <- tmp1[j, 2]
    
    # Safe Fisher matrix calculation
    contingency <- matrix(c(
      total - (list1 + list2) + overlap,
      list1 - overlap,
      list2 - overlap,
      overlap
    ), nrow = 2)
    
    if (any(contingency < 0 | !is.finite(contingency))) {
      res <- list(p.value = NA, estimate = NA)
    } else {
      res <- fisher.test(contingency)
    }
    
    tmp2 <- data.frame(
      source_ct = celltype_id,
      neighboring_ct = ct,
      p_val = res$p.value,
      odds_ratio = as.numeric(res$estimate),
      or05 = res$conf.int[1], 
      or95 = res$conf.int[2]
    )
    
    enrich_score_ct <- rbind(enrich_score_ct, tmp2)
  }
  
  enrich_score_PDL011 <- rbind(enrich_score_ct, enrich_score_PDL011)
}

enrich_score_PDL011 <- enrich_score_PDL011 %>% 
  mutate(adj_p_val = p.adjust(p_val, method = 'fdr')) %>%
  mutate(logOR = log(odds_ratio)) %>%
  mutate(logOR_lower = log(or05)) %>%
  mutate(logOR_upper = log(or95)) %>%
  mutate(Sample_combined = "PDL011")

# Save for plotting
saveRDS(enrich_score_PDL011, "/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/proximity_enrichment_30r_PDL011_2025_04_25.rds")

write_csv(enrich_score_PDL011, "/scratch/smallapragada/run2/prox_logOR_csvs/enrichment_scores_PDL011_2025_05_08.csv")
