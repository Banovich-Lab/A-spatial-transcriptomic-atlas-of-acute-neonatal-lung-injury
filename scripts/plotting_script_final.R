#### LOAD PACKAGES AND SET ENVIRONMENT ----
library(SeuratObject)
library(tidyverse)
library(gplots)
library(rlist)
library(Seurat, lib.loc = "/usr/local/lib/R/site-library")
library(gprofiler2)
library(ggrepel)
library(speckle)
library(ComplexHeatmap)
library(RColorBrewer)
library(ggpubr)
library(edgeR)
library(limma)
library(patchwork)
library(scCustomize)
library(EnhancedVolcano)
library(ggVennDiagram)
library(forcats)
library(ggrepel)
library(scCustomize)

## Set seed & working directory
set.seed(0712)
work_dir <- setwd("/scratch/smallapragada/run2")

## Load in Seurat object
run2_final <- readRDS("/scratch/smallapragada/run2/final_complete_cniches_tniches_allsamples_2025_05_26.rds")

#### FIGURE 1 ----

### Figure 1A - UMAP grouped by lineage (Biorender figure)

lin_colors <- c("Epithelial" = "#339933",
                "Endothelial" = "darkred",
                "Mesenchymal" = "#A058DD",
                "Immune" = "#5580FF") 

plot <- DimPlot(run2_final,
                reduction = "umap",
                group.by = "Lineage",
                raster = T,
                cols = lin_colors) +
  coord_equal() + 
  NoLegend() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_blank(),
        axis.title = element_blank(),
        line = element_blank()) 

pdf("/home/smallapragada/manuscript_2024/umap_lineage.pdf", width = 2.8, height = 2.2)
plot
dev.off()

### Figure 1B - Cell proportion bar chart

ct_colors <- c(# Epithelial
  "AKR1C1+ & AKR1C2+" = "#006D5B", 
  "Multiciliated" = "#2E8B57", 
  "Basal" = "#3CB371", 
  "PNEC" = "green", 
  "Secretory 3A2+ & 1A1+" = "#73C766",
  "Secretory MUC5B+" = "#b4d382", 
  "Proliferating basal" = "#55BF3B", 
  "Proliferating airway" = "#A0D868", 
  "Immature AT2" = "#004900",
  "AT1" = "#00A36C", 
  "Transitional AT2" = "#008000",
  "AT2" = "#228B22", 
  "Proliferating Immature alveolar" = "#189d75",
  # Endothelial
  "Capillaries" = "#D30000", 
  "Arterial" = "#7C0A02", 
  "Venous" = "#B22222",
  "Lymphatic endothelial" = "#ff2800", 
  "Proliferating endothelial" = "#960019",
  # Mesenchymal
  "Fibroblasts" = "#B894B1", 
  "SMC" = "#6f28d8", 
  "Pericytes" = "#9966cb", 
  "MyoFB" = "#9330A2", 
  "Activated fibroblasts" = "#784bb4", 
  "Adventitial fibroblasts" = "#B200AD",
  # Immune
  "T" = "#6596E7", 
  "cDC" = "#0000ff", 
  "NK & NKT" = "#4160e1", 
  "B" = "#0041c2", 
  "Treg" = "#6495ed", 
  "Proliferating lymphoid" = "#1e90ff", 
  "pDC" = "#0020c2", 
  "Plasma" = "#0014a8",
  "Monocytes" = "#5FA7D0", 
  "Mast" = "#1569c7", 
  "Proliferating Meg-Ery" = "#38acec", 
  "Neutrophils" = "#82cafa", 
  "Meg-Ery" = "#488ac7", 
  "Proliferating monocytes" = "#1f45fc",
  "SPP1+ macrophages" = "#151880",
  "Alveolar macrophages" = "#357ec7")

## Summarize the number of cells per cell type and lineage
cell_counts <- run2_final@meta.data %>%
  group_by(Lineage, CT_final) %>%
  summarise(cell_count = n(), .groups = "drop")

## Reorder cell types (CT_final) by cell count within each lineage
cell_counts <- cell_counts %>%
  group_by(Lineage) %>%
  mutate(CT_final = factor(CT_final, levels = CT_final[order(-cell_count)])) %>%
  ungroup()

## Convert lineage to a factor for proper ordering
cell_counts$Lineage <- factor(cell_counts$Lineage)

## Add formatted cell count labels
cell_counts <- cell_counts %>%
  mutate(label = ifelse(cell_count >= 1000, paste0(round(cell_count / 1000, 1), "k"), as.character(cell_count)))

## Order lineage labels
cell_counts$Lineage = factor(cell_counts$Lineage, levels=c("Epithelial", "Endothelial", "Mesenchymal", "Immune"))

## Create the bar plot with proportional facet widths
plot <- ggplot(cell_counts, aes(x = CT_final, y = cell_count, fill = CT_final)) +
  geom_bar(stat = "identity", width = 0.5) +  # Set bar width
  geom_text(aes(label = label), vjust = -0.5, size = 1) +  # Add labels above bars
  facet_grid(~ Lineage, scales = "free_x", space = "free") +  # Proportional facet widths
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +  # Add space above the bars
  scale_fill_manual(values= ct_colors) + 
  theme_classic() + 
  labs(y = "Number of cells", title = " ", fill = "Lineage") +
  theme(axis.text.x = element_text(color = "black", angle = 45, hjust = 1, size = 5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_line(color = "black"),
        axis.text.y = element_text(color = "black", size = 5),
        axis.title.y = element_text(color = "black", size = 6, face = "bold"),
        strip.background = element_rect(fill = "grey90", color = "black"),
        strip.text = element_text(size = 7, face = "bold"),
        axis.ticks.y = element_line(color = "black"),
        legend.title = element_blank(),
        legend.text = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")

## Change facet colors
plot_facet <- ggplot_gtable(ggplot_build(plot))
striprt <- which(grepl("strip-r", plot_facet$layout$name) | 
                   grepl("strip-t", plot_facet$layout$name))
fills <- c("#339933", "darkred", "#A058DD", "#5580FF") # Lineages & Sample Types
colors <- c(rep("black", 4), rep(NA, 1))
font_colors <- c(rep("white", 4), rep("black", 1))
k <- 1
for (i in striprt) {
  j <- which(grepl("rect", plot_facet$grobs[[i]]$grobs[[1]]$childrenOrder))
  plot_facet$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  plot_facet$grobs[[i]]$grobs[[1]]$children[[j]]$gp$col <- colors[k]
  plot_facet$grobs[[i]]$grobs[[1]]$children[[j+1]]$children[[1]]$gp$col <- font_colors[k]
  k <- k+1
}

pdf("/home/smallapragada/manuscript_2024/ct_props.pdf", width = 5.6, height = 2.25)
gridExtra::grid.arrange(plot_facet)
dev.off()

### Figure 1C - Baby marker dotplot

baby_marker_genes <- rev(c("EPCAM", "AGER", "SFTPC", "SOX9", "KRT8", "SCGB3A2", "AKR1C1",
                           "PECAM1", "COL1A2", "SOD2", "PTPRC", "MARCO", "S100A8", "ELANE", "CPA3", 
                           "FCGR3A", "JCHAIN", "SLC25A37", "MS4A1", "TRAC", "NKG7", "MKI67"))

dotplot <- DotPlot(run2_final, group.by = "CT_final", features = baby_marker_genes)

## Extract DotPlot data
dotplot_data <- dotplot$data

## Scaled expression levels
exp_mat <- dotplot_data %>% 
  select(-pct.exp, -avg.exp) %>%  
  pivot_wider(names_from = features.plot, values_from = avg.exp.scaled) %>% 
  column_to_rownames(var = "id") %>%
  as.matrix() 

## The percentage of cells expressing a feature
percent_mat <- dotplot_data %>% 
  select(-avg.exp, -avg.exp.scaled) %>%  
  pivot_wider(names_from = features.plot, values_from = pct.exp) %>% 
  column_to_rownames(var = "id") %>%
  as.matrix() 

## Ordering cell types
exp_mat_reordered <- exp_mat[c("AT1", "AT2", "Immature AT2", "Transitional AT2", 
                               "Secretory 3A2+ & 1A1+", "Secretory MUC5B+", "Basal", 
                               "Multiciliated", "PNEC", "AKR1C1+ & AKR1C2+", 
                               "Arterial", "Capillaries", "Venous",
                               "Lymphatic endothelial",  
                               "Fibroblasts", "Activated fibroblasts", "Adventitial fibroblasts", 
                               "MyoFB", "SMC", "Pericytes", 
                               "Alveolar macrophages", "SPP1+ macrophages", "Monocytes", "Neutrophils",
                               "Mast", "cDC", "pDC", "Plasma", "Meg-Ery", "B", "T", "Treg",
                               "NK & NKT", "Proliferating Immature alveolar", 
                               "Proliferating airway", "Proliferating basal", "Proliferating endothelial",
                               "Proliferating monocytes", "Proliferating lymphoid", 
                               "Proliferating Meg-Ery"), ]

## Ordering cell types
percent_mat_reordered <- percent_mat[c("AT1", "AT2", "Immature AT2", "Transitional AT2", 
                                       "Secretory 3A2+ & 1A1+", "Secretory MUC5B+", "Basal", 
                                       "Multiciliated", "PNEC", "AKR1C1+ & AKR1C2+", 
                                       "Arterial", "Capillaries", "Venous",
                                       "Lymphatic endothelial",  
                                       "Fibroblasts", "Activated fibroblasts", "Adventitial fibroblasts", 
                                       "MyoFB", "SMC", "Pericytes", 
                                       "Alveolar macrophages", "SPP1+ macrophages", "Monocytes", "Neutrophils",
                                       "Mast", "cDC", "pDC", "Plasma", "Meg-Ery", "B", "T", "Treg",
                                       "NK & NKT", "Proliferating Immature alveolar", 
                                       "Proliferating airway", "Proliferating basal", "Proliferating endothelial",
                                       "Proliferating monocytes", "Proliferating lymphoid", 
                                       "Proliferating Meg-Ery"), ]

## Any value that is greater than 2 will be mapped to bright red
col_fun = circlize::colorRamp2(c(-2, 0, 2), colorspace::diverge_hsv(3))

## Updating lineage order
ct_col_order <- factor(c(rep("Epithelial", 10), 
                         rep("Endothelial", 4),
                         rep("Mesenchymal", 6),
                         rep("Immune", 13),
                         rep("Proliferating", 7)),
                       levels = c("Epithelial", "Endothelial", "Mesenchymal", "Immune", "Proliferating"))

## Creating a layer to add to plot
layer_fun = function(j, i, x, y, w, h, fill) {
  grid.points(x = x, y = y, 
              gp = gpar(col = col_fun(pindex(t(exp_mat_reordered), i, j))), 
              size = pindex(t(percent_mat_reordered), i, j)/100 * unit(2, "mm"),
              pch = 19)
}

## Create legend list for proportion and expression
lgd_list = list(
  Legend(labels = c(0,0.25,0.5,0.75,1), title = "Proportion", type = "points", 
         pch = 19, size = c(0,0.25,0.5,0.75,1) * unit(2, "mm"),
         legend_gp = gpar(col = "black"), direction = "vertical", ncol = 1,
         title_position = "topcenter", background = NA, 
         grid_height = unit(2, "mm"),
         grid_width = unit(2, "mm"),
         title_gp = gpar(fontsize = 4), labels_gp = gpar(fontsize = 4),
         legend_height = unit(5, "mm")))

## Heatmap parameters
hp <- Heatmap(t(exp_mat_reordered),
              column_title_gp = gpar(fontsize = 5),
              column_title_side = "top",
              column_split = ct_col_order,
              col = col_fun,
              rect_gp = gpar(type = "none"),
              layer_fun = layer_fun,
              row_names_gp = gpar(fontsize = 5),
              row_names_side = "right",
              column_names_side = "top",
              column_names_gp = gpar(fontsize = 5),
              column_dend_height = unit(2, "mm"),
              column_dend_side = "bottom",
              border = "black",
              column_names_rot = 90,
              row_names_rot = 180,
              heatmap_legend_param = list(title = "Scaled\nexpression",
                                          labels = c("-2", "-1", "0", "1", "2"),
                                          legend_direction = "vertical",
                                          title_position = "topcenter",
                                          title_gp = gpar(fontsize = 4),
                                          labels_gp = gpar(fontsize = 4),
                                          background = NA,
                                          grid_width = unit(2, "mm"),
                                          grid_height = unit(10, "mm"),
                                          legend_height = unit(10, "mm"),
                                          legend_width = unit(2, "mm")),
              cluster_rows = FALSE, 
              cluster_columns = FALSE,
              cluster_row_slices = FALSE,
              show_heatmap_legend = TRUE)

ht <- draw(hp, 
           heatmap_legend_list = lgd_list,
           heatmap_legend_side = "right",
           annotation_legend_side = "right")

pdf("/home/smallapragada/manuscript_2024/dotplot.pdf", width = 5.65, height = 2.85)
ht
dev.off()

### Figure 1D - HE and Xenium images

## CSV of cell_id and CT for explorer

ct_labels <- run2_final@meta.data %>%
  filter(Slide == "5123") %>%
  filter(grepl("GOLIATH", cell_id)) %>%
  select(cell_id, CT_final) %>%
  rename(group = CT_final)

ct_labels$cell_id <- gsub("s1_[A-Za-z0-9]+_", "", ct_labels$cell_id)

write.csv(ct_labels, "/scratch/smallapragada/ct_labels_run2.csv", row.names = F)

#### FIGURE 2 ----

### Figure 2A - Sample bar chart filled by Sublineage

sample_props <- as.data.frame(table(run2_final$Sublineage, run2_final$Sample_combined)) %>%
  rename_with(~ "Sublineage", .col = "Var1") %>%
  rename_with(~ "Sample_combined", .col = "Var2") %>%
  group_by(Sample_combined) %>% 
  mutate(
    Total = sum(Freq),                
    Proportion = Freq / Total            
  ) %>%
  ungroup() %>%
  as.data.frame()

# Ordering bars
sample_props$Sample_combined <- factor(sample_props$Sample_combined, level = c('PDL001', 'PDL002', 
                                                                               'PDL003', 'PDL006', 'PDL005',
                                                                               'BLPDL008', 'PDL004', 'PDL007', 'PDL009', 
                                                                               'PDL011', 'PDL012', 'PDL010', 'PDL013',
                                                                               'PDL016', 'PDL015', 'PDL014', 'PDL017'))

# Adding a grouping assignment
sample_props <- sample_props %>%
  mutate(groupings = case_when(sample_props$Sample_combined %in% c("PDL001", "PDL002") ~ "Early canalicular",
                               sample_props$Sample_combined %in% c("PDL003", "PDL004","PDL005", "PDL006") ~ "Late canalicular",
                               sample_props$Sample_combined %in% c("PDL007", "BLPDL008", "PDL009", "PDL010", "PDL011") ~ "Saccular",
                               sample_props$Sample_combined %in% c("PDL012", "PDL013") ~ "Alveolar",
                               sample_props$Sample_combined %in% c("PDL014", "PDL015", "PDL016", "PDL017") ~ "Rare disease",
                               TRUE ~ "help"))

sample_props$groupings <- factor(sample_props$groupings, level = c('Early canalicular', 'Late canalicular', "Saccular", "Alveolar", "Rare disease"))

sample_props$Sublineage <- factor(sample_props$Sublineage, level = c("Airway", "Alveolar", "Endothelial", 
                                                                     "Mesenchymal", "Lymphoid", "Myeloid"))

sub_colors <- c("Airway" = "#66cc66", 
                "Alveolar" = "#006600",
                "Endothelial" = "#A42153",
                "Mesenchymal" = "#8d36d6",
                "Myeloid" = "#6699ff", 
                "Lymphoid" = "#0000FF")

sample_plot <- ggplot(sample_props, aes(x = Sample_combined, y = Proportion,
                                        fill = (Sublineage))) +
  scale_y_continuous(expand = c(0, 0)) +  
  geom_bar(position = "stack", stat = "identity", width = 0.75) +
  scale_fill_manual(values = sub_colors) +
  labs(x = "Sample", y = "Proportion of Cells") +
  facet_grid(. ~ groupings, space = "free", scales = "free_x") +
  guides(fill = guide_legend(override.aes = list(size = 2),
                             title.position = "right", nrow = 1,
                             label.theme = element_text(size = 5),
                             by_row = TRUE)) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(color = "black", size = 4),
        axis.title.y = element_text(color = "black", face = "bold", size = 6),
        strip.background = element_rect(fill = "grey90", color = "black"),
        strip.text = element_text(size = 7, face = "bold"),
        axis.ticks.y = element_line(color = "black"),
        axis.title.x = element_text(color = "black", size = 6, face = "bold"),
        plot.title = element_blank(),
        legend.title = element_blank(),
        legend.key.size = unit(2.5, "mm"),
        legend.text = element_text(color = "black", size = 4),
        legend.position = "top") 


pdf("/home/smallapragada/manuscript_2024/sample_barchart.pdf", width = 5, height = 2)
sample_plot
dev.off()

### Figure 2B - Immature vs mature AT2 across timepoint 

## Bar chart - AT2

## Subsetting to samples and CTs of interest

run2_final_subset <- subset(run2_final, subset = dev_stage %in% c("Early canalicular", "Late canalicular",
                                                                  "Saccular", "Alveolar"))

run2_final_imm_at2 <- subset(run2_final_subset, subset = CT_final %in% c("Immature AT2", "AT2"))

## Separating samples into "< 21 weeks" and "> 21 weeks"
run2_final_imm_at2@meta.data <- run2_final_imm_at2@meta.data %>%
  mutate(sample_groupings = case_when(run2_final_imm_at2@meta.data$Sample_combined %in% c("PDL001", "PDL002", "PDL003") ~ "< 21 WPC",
                                      TRUE ~ "> 21 WPC"))

## Summarize the number of cells per cell type and lineage
cell_counts_at2 <- run2_final_imm_at2@meta.data %>%
  group_by(Lineage, CT_final, sample_groupings) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  group_by(sample_groupings) %>%
  mutate(prop = cell_count/sum(cell_count))

## Order 
cell_counts_at2$CT_final <- factor(cell_counts_at2$CT_final, level = c("Immature AT2", "AT2"))

## Create the bar plot with proportional facet widths
plot_at2 <- ggplot(cell_counts_at2, aes(x = sample_groupings, y = prop, fill = CT_final)) +
  geom_bar(position = "stack", stat = "identity", width = 0.5) +  # Set bar width
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_manual(values = ct_colors, limits = c("Immature AT2", "AT2")) + 
  theme_classic() + 
  labs(y = "Proportion of cells", title = " ", fill = "Cell Type") +
  theme(axis.text.x = element_text(color = "black", size = 5, angle = 90),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.key.size = unit(3, "mm"),
        legend.text = element_text(color = "black", size = 4.5),
        axis.ticks.x = element_line(color = "black"),
        axis.text.y = element_text(color = "black", size = 5),
        axis.title.y = element_text(color = "black", size = 6, face = "bold"),
        strip.text = element_blank(),
        axis.ticks.y = element_line(color = "black"),
        axis.title.x = element_blank())

pdf("/home/smallapragada/manuscript_2024/bar_at2.pdf", width = 2.5, height = 1.8)
plot_at2
dev.off()

### Expression - SOX9

## Extract expression data
data <- FetchData(run2_final_imm_at2, vars = c("SOX9", "SFTPC", "CT_final"))

## Reshape data for ggplot
data_melted <- reshape2::melt(data, id.vars = "CT_final", variable.name = "Gene", value.name = "Expression")

## Ordering cell types
data_melted$CT_final <- factor(data_melted$CT_final, level = c("Immature AT2", "AT2"))

plot_exp <- ggplot(data_melted, aes(x = CT_final, y = Expression, fill = Gene)) +
  geom_violin(position = position_dodge(0.6), alpha = 0.6, width = 0.5, scale = "width") +
  scale_fill_manual(values= c("#D1B12D", "#9C0049")) + 
  theme_classic() + 
  labs(y = "Gene expression", title = " ", fill = "Gene") +
  theme(axis.text.x = element_text(color = "black", size = 5, angle = 90),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.key.size = unit(3, "mm"),
        legend.text = element_text(color = "black", size = 4.5),
        axis.ticks.x = element_line(color = "black"),
        axis.text.y = element_text(color = "black", size = 5),
        axis.title.y = element_text(color = "black", size = 6, face = "bold"),
        strip.text = element_blank(),
        axis.ticks.y = element_line(color = "black"),
        axis.title.x = element_blank())

pdf("/home/smallapragada/manuscript_2024/violin_sox9.pdf", width = 2.3, height = 1.95)
plot_exp
dev.off()

### Figure 2C - Proximity enrichment heatmap of the 4 samples

## Loading in each developmental stage's .RDS file
enrich_score_ec <- readRDS("/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/proximity_enrichment_30r_ec_2025_04_21.rds")
enrich_score_lc <- readRDS("/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/proximity_enrichment_30r_lc_2025_04_21.rds")
enrich_score_s <- readRDS("/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/proximity_enrichment_30r_s_2025_04_21.rds")
enrich_score_a <- readRDS("/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/proximity_enrichment_30r_a_2025_04_21.rds")

enrich_score_4_groups <- do.call(rbind, list(enrich_score_ec, enrich_score_lc, enrich_score_s, enrich_score_a))

## Matrix for hclust
enrich_score_mat <- enrich_score_4_groups %>%
  mutate(`source_neighbor` = paste(source_ct, "-", neighboring_ct)) %>%
  select(source_neighbor, logOR, sample_similarity) %>%
  filter(logOR != "-Inf") %>%
  pivot_wider(names_from = sample_similarity, values_from = logOR) %>%
  column_to_rownames(var = "source_neighbor") %>%
  drop_na() %>%
  as.matrix() 

# Binary
enrich_score_mat[enrich_score_mat < 0] <- -1
enrich_score_mat[enrich_score_mat > 0] <- 1

# Distance matrix
d <- dist(enrich_score_mat)

# Hierarchical clustering   
hc <- hclust(d)

plot(as.dendrogram(hc))

hclust_results <- data.frame(cutree(hc, 16)) 

## Plotting all patterns with logOR
hclust_results <- hclust_results %>%
  rownames_to_column() %>%
  rename_with(~ "source_neighbor", .cols = "rowname") 

# Order of patterns follow gestational age
pattern_order <- c("1", "5", "7", "16", "14", "4", "11")

enrich_score_all_groups_binary_hc <- enrich_score_4_groups %>%
  mutate(`source_neighbor` = paste(source_ct, "-", neighboring_ct)) %>%
  inner_join(., hclust_results, by = "source_neighbor") %>%
  filter(cutree.hc..16. == "1" | cutree.hc..16. == "5" | 
           cutree.hc..16. == "7" | cutree.hc..16. == "16" | 
           cutree.hc..16. == "14" | cutree.hc..16. == "4" |
           cutree.hc..16. == "11") %>%
  mutate(cutree.hc..16. = factor(cutree.hc..16., levels = pattern_order)) %>%
  arrange(cutree.hc..16.)

# All significant interaction pairs
int_pairs <- enrich_score_all_groups_binary_hc %>%
  filter(adj_p_val < 0.1) %>%
  distinct(source_neighbor)

enrich_score_all_groups_binary_hc_mat <- enrich_score_all_groups_binary_hc %>%
  select(source_neighbor, logOR, sample_similarity) %>%
  inner_join(., int_pairs) %>%
  pivot_wider(names_from = sample_similarity, values_from = logOR) %>%
  column_to_rownames(var = "source_neighbor") %>%
  as.matrix()

## Draw patterns across the 4 stages
heatmap_colors <- colorRamp2(c(-3, 0, 3), 
                             c("blue", "white", "red"))

heatmap <- Heatmap(enrich_score_all_groups_binary_hc_mat,
                   col = heatmap_colors,
                   name = "logOR", 
                   show_column_names = TRUE,
                   column_names_side = c("top"),
                   column_names_rot = 0,
                   column_title = "Likelihood of cellular proximity",
                   column_title_gp = gpar(fontsize = 8, 
                                          fontface = "bold", 
                                          fontfamily = "Helvetica",
                                          just = "center"),
                   cluster_rows = FALSE,
                   cluster_columns = FALSE,
                   show_row_names = TRUE,
                   column_names_gp = gpar(fontsize = 6, 
                                          fontfamily = "Helvetica"),
                   heatmap_legend_param = list(
                     title = "logOR",
                     direction = "horizontal",
                     legend_width = unit(1.5, "cm"),
                     legend_height = unit(0.1, "cm"),
                     title_gp = gpar(fontsize = 7),
                     labels_gp = gpar(fontsize = 6)
                   ))

pdf("/home/smallapragada/manuscript_2024/prox_heatmap_patterns_labels.pdf", width = 10, height = 23)
draw(heatmap)
dev.off()

 ### Figure 2D - Lollipop plots

## Loading in rare disease .rds files

## Loading in each developmental stage's .RDS file
enrich_score_rda <- readRDS("/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/proximity_enrichment_30r_rda_2025_04_21.rds")
enrich_score_rdc <- readRDS("/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/proximity_enrichment_30r_rdc_2025_04_21.rds")
enrich_score_rdi <- readRDS("/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/proximity_enrichment_30r_rdi_2025_04_21.rds")
enrich_score_rdph <- readRDS("/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/proximity_enrichment_30r_rdph_2025_04_21.rds")

enrich_score_all_groups <- do.call(rbind, list(enrich_score_ec, enrich_score_lc, enrich_score_s, enrich_score_a, 
                                               enrich_score_rda, enrich_score_rdc, enrich_score_rdi, enrich_score_rdph))

sample_order <- c("Early canalicular", "Late canalicular", "Saccular", 
                  "Alveolar", "Rare disease (CHAOS)", "Rare disease (PH + Kidneys)", 
                  "Rare disease (Infant BPD)", "Rare disease (Adult BPD)")

lolli <- enrich_score_all_groups %>%
  mutate(`source_neighbor` = paste(source_ct, "-", neighboring_ct)) %>%
  inner_join(., int_pairs) %>%
  mutate(Significance = case_when(adj_p_val < 0.1 ~ "Significant (FDR < 0.1)",
                                  TRUE ~ "Not significant")) %>%
  filter(source_ct == 'Proliferating endothelial') %>%
  filter(neighboring_ct == 'Alveolar macrophages') %>%
  mutate(across(c(logOR, logOR_lower, logOR_upper), ~ ifelse(is.infinite(.), NA, .))) %>%
  mutate(sample_similarity = fct_relevel(sample_similarity, sample_order)) %>%
  ggplot(aes(x = sample_similarity, y = logOR, color = Significance)) +
  geom_point(size = 0.5) +
  geom_errorbar(aes(ymin = logOR_lower, ymax = logOR_upper), width = 0.5, linewidth = 0.3) +
  theme_bw(base_size = 6.5) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(face = "bold", size = 6, color = "black", hjust = 0.5),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(color = "black", size = 5),
        axis.ticks.y = element_line(color = "black"),
        plot.margin = margin(0,0,0,0, "mm"),
        axis.title.y = element_text(size = 6, color = "black"),
        legend.position = "none") +
  labs(title = 'Prolif. endothelial cells - alveolar macrophages',
       x = NULL,  y = 'log(OR)') +
  scale_color_manual(values = c('gray50', 'hotpink')) +
  geom_hline(yintercept = 0, linetype = 'dashed', linewidth = 0.2)

pdf("/home/smallapragada/manuscript_2024/lolli_prolifendo_alvmacs.pdf", width = 2.03, height = 0.9)
lolli
dev.off()

lolli <- enrich_score_all_groups %>%
  mutate(`source_neighbor` = paste(source_ct, "-", neighboring_ct)) %>%
  inner_join(., int_pairs) %>%
  mutate(Significance = case_when(adj_p_val < 0.1 ~ "Significant (FDR < 0.1)",
                                  TRUE ~ "Not significant")) %>%
  filter(source_ct == 'Secretory 3A2+ & 1A1+') %>%
  filter(neighboring_ct == 'Transitional AT2') %>%
  mutate(across(c(logOR, logOR_lower, logOR_upper), ~ ifelse(is.infinite(.), NA, .))) %>%
  mutate(sample_similarity = fct_relevel(sample_similarity, sample_order)) %>%
  ggplot(aes(x = sample_similarity, y = logOR, color = Significance)) +
  geom_point(size = 0.5) +
  geom_errorbar(aes(ymin = logOR_lower, ymax = logOR_upper), width = 0.5, linewidth = 0.3) +
  theme_bw(base_size = 6.5) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(face = "bold", size = 6, color = "black", hjust = 0.5),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(color = "black", size = 5),
        axis.ticks.y = element_line(color = "black"),
        plot.margin = margin(0,0,0,0, "mm"),
        axis.title.y = element_text(size = 6, color = "black"),
        legend.position = "none") +
  labs(title = 'Secretory - Transitional AT2',
       x = NULL,  y = 'log(OR)') +
  scale_color_manual(values = c('gray50', 'hotpink')) +
  geom_hline(yintercept = 0, linetype = 'dashed', linewidth = 0.2)

pdf("/home/smallapragada/manuscript_2024/lolli_sec3a2_transat2.pdf", width = 2.03, height = 0.9)
lolli
dev.off()

lolli <- enrich_score_all_groups %>%
  mutate(`source_neighbor` = paste(source_ct, "-", neighboring_ct)) %>%
  inner_join(., int_pairs) %>%
  mutate(Significance = case_when(adj_p_val < 0.1 ~ "Significant (FDR < 0.1)",
                                  TRUE ~ "Not significant")) %>%
  filter(source_ct == 'Capillaries') %>%
  filter(neighboring_ct == 'Proliferating endothelial') %>%
  mutate(sample_similarity = fct_relevel(sample_similarity, sample_order)) %>%
  ggplot(aes(x = sample_similarity, y = logOR, color = Significance)) +
  geom_point(size = 0.5) +
  geom_errorbar(aes(ymin = logOR_lower, ymax = logOR_upper), width = 0.5, linewidth = 0.3) +
  theme_bw(base_size = 6.5) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(face = "bold", size = 6, color = "black", hjust = 0.5),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(color = "black", size = 5),
        axis.ticks.y = element_line(color = "black"),
        plot.margin = margin(0,0,0,0, "mm"),
        axis.title.y = element_text(size = 6, color = "black"),
        legend.position = "none") +
  labs(title = 'Capillaries - Prolif. endothelial',
       x = NULL,  y = 'log(OR)') +
  scale_color_manual(values = c('gray50', 'hotpink')) +
  geom_hline(yintercept = 0, linetype = 'dashed', linewidth = 0.2)

pdf("/home/smallapragada/manuscript_2024/lolli_prolifendo_cap.pdf", width = 2.03, height = 0.9)
lolli
dev.off()

lolli <- enrich_score_all_groups %>%
  mutate(`source_neighbor` = paste(source_ct, "-", neighboring_ct)) %>%
  inner_join(., int_pairs) %>%
  mutate(Significance = case_when(adj_p_val < 0.1 ~ "Significant (FDR < 0.1)",
                                  TRUE ~ "Not significant")) %>%
  filter(source_ct == 'AT1') %>%
  filter(neighboring_ct == 'Capillaries') %>%
  mutate(sample_similarity = fct_relevel(sample_similarity, sample_order)) %>%
  ggplot(aes(x = sample_similarity, y = logOR, color = Significance)) +
  geom_point(size = 0.5) +
  geom_errorbar(aes(ymin = logOR_lower, ymax = logOR_upper), width = 0.5, linewidth = 0.3) +
  theme_bw(base_size = 6.5) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(face = "bold", size = 6, color = "black", hjust = 0.5),
        axis.text.x = element_blank(),
        axis.ticks.x = element_line(color = "black"),
        axis.text.y = element_text(color = "black", size = 5),
        axis.ticks.y = element_line(color = "black"),
        plot.margin = margin(0,0,0,0, "mm"),
        axis.title.y = element_text(size = 6, color = "black"),
        legend.position = "none") +
  labs(title = 'Proximity of AT1 cells to Capillaries',
       x = NULL,  y = 'log(OR)') +
  scale_color_manual(values = c('gray50', 'hotpink')) +
  geom_hline(yintercept = 0, linetype = 'dashed', linewidth = 0.2)

pdf("/home/smallapragada/manuscript_2024/lolli_ticks.pdf", width = 2.03, height = 0.9)
lolli
dev.off()

lolli <- enrich_score_all_groups %>%
  mutate(`source_neighbor` = paste(source_ct, "-", neighboring_ct)) %>%
  inner_join(., int_pairs) %>%
  mutate(Significance = case_when(adj_p_val < 0.1 ~ "Significant (FDR < 0.1)",
                                  TRUE ~ "Not significant")) %>%
  filter(source_ct == 'AT1') %>%
  filter(neighboring_ct == 'Capillaries') %>%
  mutate(sample_similarity = fct_relevel(sample_similarity, sample_order)) %>%
  ggplot(aes(x = sample_similarity, y = logOR, color = Significance)) +
  geom_point(size = 0.5) +
  geom_errorbar(aes(ymin = logOR_lower, ymax = logOR_upper), width = 0.5, linewidth = 0.3) +
  theme_bw(base_size = 6.5) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.box.margin = margin(0,0,0,0, "mm"),
        plot.title = element_text(face = "bold", size = 6, color = "black", hjust = 0.5),
        axis.text.x = element_text(color = "black", size = 5, angle = 90, hjust = 1),
        axis.ticks.x = element_line(color = "black"),
        axis.text.y = element_text(color = "black", size = 5),
        axis.ticks.y = element_line(color = "black"),
        plot.margin = margin(0,0,0,0, "mm"),
        axis.title.y = element_text(size = 6, color = "black"),
        legend.position = "bottom", 
        legend.title = element_blank()) +
  labs(title = 'AT1 cells - Capillaries',
       x = NULL,  y = 'log(OR)') +
  scale_color_manual(values = c('gray50', 'hotpink')) +
  geom_hline(yintercept = 0, linetype = 'dashed', linewidth = 0.2)

pdf("/home/smallapragada/manuscript_2024/lolli_bottom_label.pdf", width = 2.1, height = 2.6)
lolli
dev.off()

### Figure 2E - Niche dotplots

## Subset object to samples used
samples <- subset(run2_final, subset = dev_stages %in% c("Late canalicular", "Saccular", "Alveolar"))

## Creating a new dataframe with just cell types and lineage
lineage_df <- samples@meta.data %>%
  dplyr::select(CT_final, Lineage) %>%
  unique()

names(lineage_df) <- c("CT_final", "Lineage")

## Generating row sums for each cell niche assignment
ct_niche_prop_df_cs <- as.data.frame((table(samples$CT_final, samples$CNiches)/rowSums(table(samples$CT_final, samples$CNiches))))

names(ct_niche_prop_df_cs) <- c("CT_final", "CNiche", "Proportion")

## Joining row sums dataframe and lineage dataframe - cell
ct_niche_prop_df_cs_lineage <- ct_niche_prop_df_cs %>%
  full_join(lineage_df) %>%
  mutate(Proportion_Visible = ifelse(Proportion > 0.05, TRUE, FALSE))

## Reorder the levels of lineage
ct_niche_prop_df_cs_lineage$Lineage <- fct_relevel(ct_niche_prop_df_cs_lineage$Lineage, 
                                                   "Epithelial", 
                                                   "Endothelial", 
                                                   "Mesenchymal", 
                                                   "Immune")

## Niche colors
cniche_colors <- c("C1" = "#E63946",
                   "C2" = "#F4A261",
                   "C3" = "#01796f",
                   "C4" = "#4f6389",
                   "C5" = "#D67AB1",
                   "C6" = "#457B9D",
                   "C7" = "#992222",
                   "C8" = "grey70",
                   "C9" = "#6A4C93",
                   "C10" = "#E9C46A",
                   "C11" = "#A8DADC") 

## Order CTs
ct_niche_prop_df_cs_lineage$CT_final <- factor(ct_niche_prop_df_cs_lineage$CT_final,
                                               level = c("AT1", "AT2", "Immature AT2", "Transitional AT2", 
                                                         "Secretory 3A2+ & 1A1+", "Secretory MUC5B+", "Basal", 
                                                         "Multiciliated", "PNEC", "AKR1C1+ & AKR1C2+", "Proliferating airway", 
                                                         "Proliferating basal", "Proliferating Immature alveolar",
                                                         "Arterial", "Capillaries", "Venous",
                                                         "Lymphatic endothelial", "Proliferating endothelial",
                                                         "Fibroblasts", "Activated fibroblasts", "Adventitial fibroblasts", 
                                                         "MyoFB", "SMC", "Pericytes", 
                                                         "Alveolar macrophages", "SPP1+ macrophages", "Monocytes", "Neutrophils",
                                                         "Mast", "cDC", "pDC", "Plasma", "Meg-Ery", "B", "T", "Treg",
                                                         "NK & NKT", "Proliferating monocytes", "Proliferating lymphoid", 
                                                         "Proliferating Meg-Ery"))

# Cell niche order
ct_niche_prop_df_cs_lineage$CNiche <- factor(ct_niche_prop_df_cs_lineage$CNiche, 
                                             levels = c('C1', 'C2', 'C3', 'C4', 'C5', 'C6',
                                                        'C7', 'C8', 'C9', 'C10', 'C11'))

cell_plot <- ggplot(ct_niche_prop_df_cs_lineage, aes(y = reorder(CNiche, dplyr::desc(CNiche)), x = CT_final, 
                                                     size = Proportion, color = as.factor(CNiche), 
                                                     alpha = Proportion_Visible,
                                                     fill = as.factor(CNiche))) +
  geom_point(shape = 21, color = "black") +
  facet_grid(~Lineage, scales = "free", space = "free") +
  theme_bw() +
  scale_fill_manual(values = cniche_colors, guide = "none") +
  scale_alpha_manual(values = c(0, 1), guide = "none") +
  scale_size_continuous(range = c(0, 3), breaks = seq(0.2, 1, 0.2), limits = c(0, 1)) +
  guides(size = guide_legend(direction = "horizontal", title.position = "top")) +
  labs(title = "Cell Types Across Niches", x = "Cell Type", y = "CNiche") +
  theme(axis.text.x = element_blank(),
        plot.title = element_text(color = "black", size = 8, face = "bold", hjust = 0.5),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(color = "grey90", linewidth = 0.2),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(color = "black", size = 5),
        strip.background = element_rect(lin_colors[as.numeric(as.factor(ct_niche_prop_df_cs_lineage$Lineage))], 
                                        color = "black"),
        strip.text.x = element_text(size = 6, face = "bold", lineheight = 0.5),
        axis.ticks.y = element_line(color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "black", face = "bold", size = 6),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.5, "lines"),
        legend.spacing.y = unit(0.2, "cm")) 

## Change facet colors
plot_facet <- ggplot_gtable(ggplot_build(cell_plot))

striprt <- which(grepl("strip-r", plot_facet$layout$name) | 
                   grepl("strip-t", plot_facet$layout$name))
fills <- c("#339933", "darkred", "#A058DD", "#5580FF") # Lineages & Sample Types
colors <- c(rep("black", 4), rep(NA, 1))
font_colors <- c(rep("white", 4), rep("black", 1))
k <- 1
for (i in striprt) {
  j <- which(grepl("rect", plot_facet$grobs[[i]]$grobs[[1]]$childrenOrder))
  plot_facet$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  plot_facet$grobs[[i]]$grobs[[1]]$children[[j]]$gp$col <- colors[k]
  plot_facet$grobs[[i]]$grobs[[1]]$children[[j+1]]$children[[1]]$gp$col <- font_colors[k]
  k <- k+1
}

pdf("/home/smallapragada/manuscript_2024/cniche_dotplot.pdf", width = 6.6, height = 1.85)
gridExtra::grid.arrange(plot_facet)
dev.off()

### Transcript plot

### Generating row sums for each transcript niche assignment
ct_niche_prop_df_ts <- as.data.frame((table(samples$CT_final, samples$TNiches)/rowSums(table(samples$CT_final, samples$TNiches))))

names(ct_niche_prop_df_ts) <- c("CT_final", "TNiches", "Proportion")

## Joining row sums dataframe and lineage dataframe - transcripts
ct_niche_prop_df_ts_lineage <- ct_niche_prop_df_ts %>%
  full_join(lineage_df) %>%
  mutate(Proportion_Visible = ifelse(Proportion > 0.08, TRUE, FALSE))

## Reorder the levels of lineage
ct_niche_prop_df_ts_lineage$Lineage <- fct_relevel(ct_niche_prop_df_ts_lineage$Lineage, 
                                                   "Epithelial", 
                                                   "Endothelial", 
                                                   "Mesenchymal", 
                                                   "Immune")

## Niche colors
tniche_colors <- c("T1" = "#FF5733",
                   "T2" = "#B5651D",
                   "T3" = "#5E8C61",
                   "T4" = "#8A9A5B",
                   "T5" = "#E3DFFF",
                   "T6" = "#FFCC00",
                   "T7" = "#5C3A75",
                   "T8" = "#B39C82",
                   "T9" = "#00416A",
                   "T10" = "#D4A5A5",
                   "T11" = "#FF85A1") 

## Ordering cell types for both niches
ct_niche_prop_df_ts_lineage$CT_final <- factor(ct_niche_prop_df_ts_lineage$CT_final,
                                               level = c("AT1", "AT2", "Immature AT2", "Transitional AT2", 
                                                         "Secretory 3A2+ & 1A1+", "Secretory MUC5B+", "Basal", 
                                                         "Multiciliated", "PNEC", "AKR1C1+ & AKR1C2+", "Proliferating airway", 
                                                         "Proliferating basal", "Proliferating Immature alveolar",
                                                         "Arterial", "Capillaries", "Venous",
                                                         "Lymphatic endothelial", "Proliferating endothelial",
                                                         "Fibroblasts", "Activated fibroblasts", "Adventitial fibroblasts", 
                                                         "MyoFB", "SMC", "Pericytes", 
                                                         "Alveolar macrophages", "SPP1+ macrophages", "Monocytes", "Neutrophils",
                                                         "Mast", "cDC", "pDC", "Plasma", "Meg-Ery", "B", "T", "Treg",
                                                         "NK & NKT", "Proliferating monocytes", "Proliferating lymphoid", 
                                                         "Proliferating Meg-Ery"))

## Plotting both the cell and transcript niches in the same DotPlot

# Transcript
ct_niche_prop_df_ts_lineage$TNiches <- factor(ct_niche_prop_df_ts_lineage$TNiches, 
                                                 levels = c('T1', 'T2', 'T3', 'T4', 'T5',
                                                            'T6', 'T7', 'T8', 'T9', 'T10', 'T11'))

transcript_plot <- ggplot(ct_niche_prop_df_ts_lineage, aes(y = reorder(TNiches, dplyr::desc(TNiches)), x = CT_final, 
                                                           size = Proportion, color = as.factor(TNiches), 
                                                           alpha = Proportion_Visible,
                                                           fill = as.factor(TNiches))) +
  geom_point(shape = 21, color = "black") +
  facet_grid(~Lineage, scales = "free", space = "free") +
  theme_bw() + 
  scale_fill_manual(values = tniche_colors, guide = "none") +
  scale_alpha_manual(values = c(0, 1), guide = "none") +
  scale_size_continuous(range = c(0, 3), breaks = seq(0.2, 1, 0.2), limits = c(0, 1)) +
  labs(x = "Cell Type", y = "TNiche") +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 5, angle = 45, hjust = 1),
        axis.text.y = element_text(color = "black", size = 5),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.ticks.y = element_line(color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "black", face = "bold", size = 6),
        legend.position = "none") 

pdf("/home/smallapragada/manuscript_2024/tniche_dotplot.pdf", width = 5.3, height = 2.1)
transcript_plot
dev.off()

### Figure 2F: Niche bar plot 

### CNiches

## Subset object to samples used
samples <- subset(run2_final, subset = dev_stages %in% c("Early canalicular", "Late canalicular", "Saccular", "Alveolar"))

## Ordering bars
samples$dev_stages <- factor(samples$dev_stages, level = c("Early canalicular", "Late canalicular",
                                                           "Saccular", "Alveolar"))

## Summarize the number of cells per niche and dev stage
cell_counts_cniche <- samples@meta.data %>%
  group_by(CNiches, dev_stages) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  group_by(dev_stages) %>%
  mutate(prop = cell_count/sum(cell_count))

## Ordering CNiches 

cell_counts_cniche$CNiches <- factor(cell_counts_cniche$CNiches, level = c("C1", "C2", "C3", 
                                                                           "C4", "C5", "C6",
                                                                           "C7", "C8", "C9",
                                                                           "C10", "C11"))

cniche_plot <- ggplot(cell_counts_cniche, aes(x = cell_counts_cniche$groupings, y = cell_counts_cniche$prop,
                                              fill = as.factor(cell_counts_cniche$CNiches))) +
  scale_y_continuous(expand = c(0, 0)) +  
  geom_bar(position = "stack", stat = "identity", width = 0.75) +
  scale_fill_manual(values = cniche_colors) +
  labs(title = "Cell niches", x = "", y = "Proportion of Niches") +
  guides(fill = guide_legend(title.position = "right", ncol = 2,
                             label.theme = element_text(size = 5),
                             by_row = TRUE)) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(color = "black", size = 4),
        axis.title.y = element_text(color = "black", size = 6, face = "bold"),
        strip.background = element_rect(fill = "grey90", color = "black"),
        strip.text = element_text(size = 6, face = "bold"),
        axis.ticks.y = element_line(color = "black"),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 6, face = "bold", hjust = 0.5),
        legend.title = element_blank(),
        legend.key.size = unit(3, "mm"),
        legend.text = element_text(color = "black", size = 4),
        legend.position = "right") 

### TNiches

## Summarize the number of cells per niche and dev stage
cell_counts_tniche <- samples@meta.data %>%
  group_by(TNiches, dev_stages) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  group_by(dev_stages) %>%
  mutate(prop = cell_count/sum(cell_count))

## Ordering TNiches 
cell_counts_tniche$TNiches <- factor(cell_counts_tniche$TNiches, level = c("T1", "T2", "T3", 
                                                                                 "T4", "T5", "T6",
                                                                                 "T7", "T8", "T9",
                                                                                 "T10", "T11"))

tniche_plot <- ggplot(cell_counts_tniche, aes(x = cell_counts_tniche$dev_stages, y = cell_counts_tniche$prop,
                                              fill = as.factor(cell_counts_tniche$TNiches))) +
  scale_y_continuous(expand = c(0, 0)) +  
  geom_bar(position = "stack", stat = "identity", width = 0.75) +
  scale_fill_manual(values = tniche_colors) +
  labs(title = "Transcript niches", x = "", y = "Proportion of Niches") +
  guides(fill = guide_legend(title.position = "right", ncol = 2,
                             label.theme = element_text(size = 5),
                             by_row = TRUE)) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_line(color = "black"),
        axis.text.x = element_text(color = "black", size = 5, angle = 90, hjust = 1),
        axis.text.y = element_text(color = "black", size = 4),
        axis.title.y = element_text(color = "black", size = 6, face = "bold"),
        strip.background = element_rect(fill = "grey90", color = "black"),
        strip.text = element_text(size = 6, face = "bold"),
        axis.ticks.y = element_line(color = "black"),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 6, face = "bold", hjust = 0.5),
        legend.title = element_blank(),
        legend.key.size = unit(3, "mm"),
        legend.text = element_text(color = "black", size = 4),
        legend.position = "right")

## Combining plots
c_t_bar <- cniche_plot / tniche_plot

pdf("/home/smallapragada/manuscript_2024/niche_barchart.pdf", width = 2.8, height = 4)
c_t_bar
dev.off()

#### FIGURE 3 ----

### Figure 3B - GA/LS/DX plot

ga_ls_dx_data <- run2_final@meta.data %>%
  filter(dev_stages %in% c("Late canalicular", "Saccular", "Alveolar")) %>%
  select(Sample_combined, Gestational_age_weeks, Life_span_weeks, dx_score) %>%
  distinct() %>%
  rename_with(~ "ga", .cols = "Gestational_age_weeks") %>%
  rename_with(~ "ls", .cols = "Life_span_weeks") 

ga_ls_dx_data$dx_score <- as.factor(ga_ls_dx_data$dx_score)

dx_colors <- c("0" = "#50C878",
               "1" = "#F1C40F",
               "2" = "#E67E22",
               "3"  = "#FF0000")

plot <- ggplot(ga_ls_dx_data, aes(x = as.numeric(ga), y = as.numeric(ls), fill = dx_score)) +
  geom_point(shape = 21, size = 2) +
  geom_text_repel(aes(label = Sample_combined), size = 1.5, max.overlaps = 100) +
  theme_classic() + 
  labs(x = "Gestational age (weeks)", y = "Life span (weeks)", fill = "Disease severity score") +
  scale_fill_manual(values = dx_colors) +
  theme(axis.text.x = element_text(color = "black", size = 5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(color = "black", size = 6, face = "bold"),
        legend.key.size = unit(2.5, "mm"),
        legend.text = element_text(color = "black", size = 4.5),
        axis.ticks.x = element_line(color = "black"),
        axis.text.y = element_text(color = "black", size = 5),
        axis.title.y = element_text(color = "black", size = 6, face = "bold"),
        strip.text = element_blank(),
        axis.ticks.y = element_line(color = "black"),
        axis.title.x = element_text(color = "black", size = 6, face = "bold"))
  
pdf("/home/smallapragada/manuscript_2024/sample_ga_ls_dx.pdf", width = 3.6, height = 2.3)
plot
dev.off()

### Figure 3C - GA LS DX overlap volcano plot

term_colors <- c("ga" = "#0073C2FF", 
                 "ls" = "#EFC000FF", 
                 "dx" = "#CD534CFF")

## Pick the top three per category
top_sig_genes <- all_stats_morphed %>%
  filter(sig == "Significant") %>%
  group_by(origin) %>%
  arrange(desc(lfc), .by_group = TRUE) %>%
  slice_head(n = 3) %>%
  select(origin, gene_ct) %>%
  mutate(gene_col = gene_ct)

## Only label those
test <- all_stats_morphed %>%
  left_join(top_sig_genes %>% select(origin, gene_ct, gene_col),
            by = c("origin", "gene_ct"))

plot <- ggplot(test, aes(x = lfc, y = -log10(adj_p_val), col = origin, alpha = sig)) +
  geom_point(size = 0.5) + 
  geom_label_repel(aes(label = gene_col), 
                   size = 1.8, 
                   max.overlaps = Inf, 
                   box.padding = 0.25, 
                   point.padding = 0.25,
                   segment.size = 0.2,   
                   min.segment.length = 0,
                   show.legend = FALSE) +
  scale_color_manual(values = term_colors, 
                     labels = c("DX score", "Gestational age", "Life span")) + 
  scale_alpha_manual(values = c("Significant" = 1, "N.S." = 0.1)) +
  xlab("logFC") + ylab("-log10(adj_p_val)") +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black", size = 4),
        axis.text.y = element_text(color = "black", size = 4),
        axis.title.x = element_text(color = "black", face = "bold", size = 6),
        axis.title.y = element_text(color = "black", face = "bold", size = 6),
        panel.grid.major.x = element_line(color = "grey80", linewidth = 0.4),
        panel.grid.major.y = element_line(color = "grey80", linewidth = 0.4),
        axis.ticks.x = element_line(color = "black"),
        axis.ticks.y = element_line(color = "black"),
        legend.title = element_blank(),
        legend.position = "top",
        legend.text = element_text(color = "black", size = 6))

pdf("/home/smallapragada/manuscript_2024/volcano_top3.pdf", width = 3.5, height = 3.5)
plot
dev.off()

### Figure 3D - Count heatmap of pairs per CT

ct_order_heat <- c("AT2", "Secretory 3A2+ & 1A1+", "Multiciliated",
                   "Basal", "PNEC", "Proliferating airway", "Arterial",
                   "Fibroblasts", "Alveolar macrophages", "SPP1+ macrophages",
                   "Plasma", "cDC", "Mast", "T", "Meg-Ery", 
                   "Proliferating Meg-Ery", "Proliferating lymphoid")

test <- all_stats_morphed %>%
  filter(sig == "Significant") %>%
  mutate(origin = case_when(origin == "ga" ~ "Gestational age",
                            origin == "ls" ~ "Life span",
                            origin == "dx" ~ "DX score",
                            TRUE ~ "help")) %>%
  group_by(celltype, origin) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = "origin", values_from = "count") %>%
  column_to_rownames("celltype") %>% 
  relocate("DX score", .after = "Life span") %>%
  .[ct_order_heat, , drop = FALSE] %>%  # Reorder rows
  as.matrix()

## Draw patterns across the 4 stages
heatmap_colors <- colorRamp2(c(1, 5, 10, 20), 
                             c("white", "#FF9999", "#FF3333", "#A10000"))

heatmap <- Heatmap(test,
                   col = heatmap_colors,
                   name = "Count", 
                   show_column_names = TRUE,
                   column_names_side = c("bottom"),
                   column_names_rot = 90,
                   column_title = "Count",
                   column_title_gp = gpar(fontsize = 8, 
                                          fontface = "bold", 
                                          fontfamily = "Helvetica",
                                          just = "center"),
                   cluster_rows = FALSE,
                   cluster_columns = FALSE,
                   show_row_names = TRUE,
                   column_names_gp = gpar(fontsize = 5, 
                                          fontfamily = "Helvetica"),
                   row_names_gp = gpar(fontsize = 5,
                                       fontfamily = "Helvetica"),
                   heatmap_legend_param = list(
                     title_gp = gpar(fontsize = 5,
                                     fontface = "bold"),
                     labels_gp = gpar(fontsize = 3),
                     legend_height = unit(2, "cm"),  
                     grid_width = unit(0.3, "cm")   
                   ))

pdf("/home/smallapragada/manuscript_2024/count_heatmap.pdf", width = 2, height = 3.25)
draw(heatmap)
dev.off()

### Figure 3E - Lollipop plots with genes in cell types across GA, LS, and DX

## Picking 5 gene - CT pairs in each overlap

ga_ls_intersect <- intersect(ga_pairs_with_thresh$gene_ct, ls_pairs_with_thresh$gene_ct)
ga_dx_intersect <- intersect(ga_pairs_with_thresh$gene_ct, dx_pairs_with_thresh$gene_ct)
dx_ls_intersect <- intersect(dx_pairs_with_thresh$gene_ct, ls_pairs_with_thresh$gene_ct)

stats_morphed_pairs_chosen <- all_stats_morphed %>%
  filter(gene_ct == 'IL7R - Arterial' |
           gene_ct == 'HIF1A - Multiciliated' |
           gene_ct == 'S100A2 - Basal' |
           gene_ct == 'HMGA1 - Multiciliated' |
           gene_ct == 'COL3A1 - Plasma' |
           gene_ct == 'SEC11C - PNEC') %>%
  mutate(gene_ct = fct_relevel(gene_ct, rev(c('S100A2 - Basal', 'HMGA1 - Multiciliated',
                                          'HIF1A - Multiciliated', 'IL7R - Arterial', 
                                          'SEC11C - PNEC', 'COL3A1 - Plasma'))))

lolli <- stats_morphed_pairs_chosen %>%
  filter(origin == "ga") %>%
  ggplot(aes(x = gene_ct, y = lfc, color = sig)) +
  geom_point(size = 0.5) +
  geom_errorbar(aes(ymin = se_lower, ymax = se_upper), width = 0.5, linewidth = 0.3) +
  theme_bw(base_size = 6.5) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.box.margin = margin(0,0,0,0, "mm"),
        plot.title = element_text(face = "bold", size = 6, color = "black", hjust = 0.5),
        axis.text.x = element_text(color = "black", size = 4),
        axis.ticks.x = element_line(color = "black"),
        axis.text.y = element_text(color = "black", size = 4, angle = 0, hjust = 1),
        axis.ticks.y = element_line(color = "black"),
        plot.margin = margin(0,0,0,0, "mm"),
        axis.title.y = element_text(size = 5, color = "black", face = "bold"),
        legend.position = "bottom", 
        legend.title = element_blank()) +
  labs(title = 'Gestational age',
       x = NULL, y = "logFC") +
  scale_color_manual(values = c('gray50', 'purple')) +
  geom_hline(yintercept = 0, linetype = 'dashed', linewidth = 0.2) +
  coord_flip()

pdf("/home/smallapragada/manuscript_2024/lolli_gene_cts_ga.pdf", width = 2, height = 1.5)
lolli
dev.off()

lolli <- stats_morphed_pairs_chosen %>%
  filter(origin == "ls") %>%
  ggplot(aes(x = gene_ct, y = lfc, color = sig)) +
  geom_point(size = 0.5) +
  geom_errorbar(aes(ymin = se_lower, ymax = se_upper), width = 0.5, linewidth = 0.3) +
  theme_bw(base_size = 6.5) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.box.margin = margin(0,0,0,0, "mm"),
        plot.title = element_text(face = "bold", size = 6, color = "black", hjust = 0.5),
        axis.text.x = element_text(color = "black", size = 4),
        axis.ticks.x = element_line(color = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = margin(0,0,0,0, "mm"),
        axis.title.y = element_text(size = 5, color = "black", face = "bold"),
        legend.position = "bottom", 
        legend.title = element_blank()) +
  labs(title = 'Life span',
       x = NULL, y = "logFC") +
  scale_color_manual(values = c('gray50', 'purple')) +
  geom_hline(yintercept = 0, linetype = 'dashed', linewidth = 0.2) +
  coord_flip()

pdf("/home/smallapragada/manuscript_2024/lolli_gene_cts_ls.pdf", width = 1.43, height = 1.5)
lolli
dev.off()

lolli <- stats_morphed_pairs_chosen %>%
  filter(origin == "dx") %>%
  ggplot(aes(x = gene_ct, y = lfc, color = sig)) +
  geom_point(size = 0.5) +
  geom_errorbar(aes(ymin = se_lower, ymax = se_upper), width = 0.5, linewidth = 0.3) +
  theme_bw(base_size = 6.5) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.box.margin = margin(0,0,0,0, "mm"),
        plot.title = element_text(face = "bold", size = 6, color = "black", hjust = 0.5),
        axis.text.x = element_text(color = "black", size = 4),
        axis.ticks.x = element_line(color = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = margin(0,0,0,0, "mm"),
        axis.title.y = element_text(size = 5, color = "black", face = "bold"),
        legend.position = "bottom", 
        legend.title = element_blank()) +
  labs(title = 'DX score',
       x = NULL, y = "logFC") +
  scale_color_manual(values = c('gray50', 'purple')) +
  geom_hline(yintercept = 0, linetype = 'dashed', linewidth = 0.2) +
  coord_flip()

pdf("/home/smallapragada/manuscript_2024/lolli_gene_cts_dx.pdf", width = 1.43, height = 1.5)
lolli
dev.off()

#### SUPPLEMENTARY FIGURES ----

ct_colors_random <- randomcoloR::distinctColorPalette(40)

### Supplementary figure 1

sample_colors <- c("Early canalicular" = "#F4E6AA",
                   "Late canalicular" = "#F2B342",
                   "Saccular" = "#F07F12",
                   "Alveolar" = "#C43E96",
                   "Rare disease" = "#C03830")

run2_final@meta.data$dev_stages <- factor(run2_final@meta.data$dev_stages, level = c("Early canalicular",
                                                                                     "Late canalicular",
                                                                                     "Saccular",
                                                                                     "Alveolar",
                                                                                     "Rare disease"))

plot <- DimPlot(run2_final,
                reduction = "sp",
                group.by = "dev_stages", 
                raster = T,
                cols = sample_colors) +
  coord_equal() + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_blank(),
        axis.title = element_blank(),
        line = element_blank(),
        legend.key.size = unit(0.3, "cm"),     
        legend.text = element_text(size = 6),      
        legend.spacing.y = unit(0.1, 'cm')         
  ) +
  guides(color = guide_legend(ncol = 1))

pdf("/home/smallapragada/manuscript_2024/sp_both_slides_supp.pdf", width = 7, height = 5.5)
plot
dev.off()

### Supplementary figure 2

## Sample UMAP

run2_final$Sample_combined <- factor(run2_final$Sample_combined, levels = c("PDL001", "PDL002", "PDL003",
                                                                            "PDL004", "PDL005", "PDL006",
                                                                            "PDL007", "PDL008", "PDL009",
                                                                            "PDL010", "PDL011", "PDL012", 
                                                                            "PDL013", "PDL014", "PDL015", 
                                                                            "PDL016", "PDL017"))

plot <- DimPlot(run2_final,
                reduction = "umap",
                group.by = "Sample_combined",
                cols = randomcoloR::distinctColorPalette(17)) +
  coord_equal() + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_blank(),
        axis.title = element_blank(),
        line = element_blank(),
        legend.key.size = unit(0.3, "cm"),     
        legend.text = element_text(size = 6),      
        legend.spacing.y = unit(0.1, 'cm')         
  ) +
  guides(color = guide_legend(ncol = 1))

pdf("/home/smallapragada/manuscript_2024/umap_sample_supp.pdf", width = 3.2, height = 4)
plot
dev.off()

## CT UMAP
run2_final$CT_final <- factor(run2_final$CT_final, levels = c("AT1", "AT2", "Immature AT2", "Transitional AT2", 
                                                              "Secretory 3A2+ & 1A1+", "Secretory MUC5B+", "Basal", 
                                                              "Multiciliated", "PNEC", "AKR1C1+ & AKR1C2+", "Arterial", 
                                                              "Capillaries", "Venous", "Lymphatic endothelial",  
                                                              "Fibroblasts", "Activated fibroblasts", "Adventitial fibroblasts", 
                                                              "MyoFB", "SMC", "Pericytes", 
                                                              "Alveolar macrophages", "SPP1+ macrophages", "Monocytes", 
                                                              "Neutrophils", "Mast", "cDC", "pDC", "Plasma", "Meg-Ery", 
                                                              "B", "T", "Treg", "NK & NKT", "Proliferating Immature alveolar", 
                                                              "Proliferating airway", "Proliferating basal", "Proliferating endothelial",
                                                              "Proliferating monocytes", "Proliferating lymphoid", "Proliferating Meg-Ery"))

plot <- DimPlot(run2_final,
                reduction = "umap",
                group.by = "CT_final",
                cols = ct_colors_random) +
  coord_equal() + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_blank(),
        axis.title = element_blank(),
        line = element_blank(),
        legend.key.size = unit(0.3, "cm"),     
        legend.text = element_text(size = 6),      
        legend.spacing.y = unit(0.1, 'cm')         
  ) +
  guides(color = guide_legend(ncol = 2))

pdf("/home/smallapragada/manuscript_2024/umap_ct_supp.pdf", width = 6, height = 3)
plot
dev.off()

## CT Dotplot

plot <- DotPlot_scCustom(run2_final, group.by = "CT_final",
                         features = c("AGER", "SFTPC", "SOX9", "KRT8", "SCGB3A2", "MUC5B",
                                      "KRT5", "C20orf85", "CALCA", "AKR1C1", "HEY1", "CA4",
                                      "ACKR1", "CCL21", "COL1A2", "SFRP2", "MFAP5", 
                                      "ACTA2", "WNT5A", "PDGFRB", "MARCO", "SPP1", "S100A8", 
                                      "ELANE", "CPA3", "FCGR3A", "JCHAIN", "PIM2", "SLC25A37",  
                                      "BANK1", "TRAC", "CTLA4", "NKG7", "MKI67")) + 
  theme(
    axis.text.x = element_text(size = 8, color = "black", angle = 90),
    axis.text.y = element_text(size = 8, color = "black"),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7),
    legend.key.size = unit(0.3, "cm"),
    legend.spacing.y = unit(0.2, "cm")
  ) + 
  scale_size(range = c(-1, 2.5))

pdf("/home/smallapragada/manuscript_2024/dotplot_supp.pdf", width = 7, height = 5.5)
plot
dev.off()

### Supplementary figure 3

dev_sample_props <- as.data.frame(table(run2_final$Sample_combined, run2_final$dev_stages)) %>%
  rename_with(~ "Sample_combined", .col = "Var1") %>%
  rename_with(~ "dev_stage", .col = "Var2") %>%
  group_by(dev_stage) %>%
  mutate(
    Total = sum(Freq),
    Proportion = Freq / Total
  ) %>%
  ungroup() %>%
  as.data.frame()

# Ordering bars
dev_sample_props$dev_stage <- factor(dev_sample_props$dev_stage , level = c("Early canalicular", "Late canalicular", 
                                                                            "Saccular", "Alveolar", "Rare disease"))

plot <- ggplot(dev_sample_props, aes(x = dev_stage, y = Proportion,
                                        fill = (Sample_combined))) +
  scale_y_continuous(expand = c(0, 0)) +  
  geom_bar(position = "stack", stat = "identity", width = 0.75) +
  scale_fill_manual(values = randomcoloR::distinctColorPalette(17)) +
  labs(x = " ", y = "Proportion of Cells") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color = "black", size = 8),
        axis.text.y = element_text(color = "black", size = 6),
        axis.title.y = element_text(color = "black", face = "bold", size = 8),
        strip.background = element_rect(fill = "grey90", color = "black"),
        axis.ticks.y = element_line(color = "black"),
        axis.title.x = element_blank(),
        plot.title = element_blank(),
        legend.title = element_blank(),
        legend.key.size = unit(4, "mm"),
        legend.text = element_text(color = "black", size = 6),
        legend.position = "right") 


pdf("/home/smallapragada/manuscript_2024/sample_barchart_dev_stage_supp.pdf", width = 7.15, height = 3)
plot
dev.off()

### Supplementary figure 4

sample_props <- as.data.frame(table(run2_final$CT_final, run2_final$Sample_combined)) %>%
  rename_with(~ "CT_final", .col = "Var1") %>%
  rename_with(~ "Sample_combined", .col = "Var2") %>%
  group_by(Sample_combined) %>%
  mutate(
    Total = sum(Freq),
    Proportion = Freq / Total
  ) %>%
  ungroup() %>%
  as.data.frame()

# Groupings
sample_props <- sample_props %>%
  mutate(groupings = case_when(sample_props$Sample_combined %in% c("PDL001", "PDL002") ~ "Early canalicular",
                               sample_props$Sample_combined %in% c("PDL003", "PDL004","PDL005", "PDL006") ~ "Late canalicular",
                               sample_props$Sample_combined %in% c("PDL007", "PDL008", "PDL009", "PDL010", "PDL011") ~ "Saccular",
                               sample_props$Sample_combined %in% c("PDL012", "PDL013") ~ "Alveolar",
                               sample_props$Sample_combined %in% c("PDL014", "PDL015", "PDL016", "PDL017") ~ "Rare disease",
                               TRUE ~ "help"))

# Ordering bars
sample_props$Sample_combined <- factor(sample_props$Sample_combined, level = c('PDL001', 'PDL002',
                                                                               'PDL003', 'PDL006', 'PDL005',
                                                                               'PDL008', 'PDL004', 'PDL007', 'PDL009',
                                                                               'PDL011', 'PDL012', 'PDL010', 'PDL013',
                                                                               'PDL016', 'PDL015', 'PDL014', 'PDL017'))

sample_props$groupings <- factor(sample_props$groupings, level = c("Early canalicular", "Late canalicular", "Saccular",
                                                                   "Alveolar", "Rare disease"))

sample_plot <- ggplot(sample_props, aes(x = Sample_combined, y = Proportion,
                                        fill = (CT_final))) +
  scale_y_continuous(expand = c(0, 0)) +  
  geom_bar(position = "stack", stat = "identity", width = 0.75) +
  scale_fill_manual(values = ct_colors_random) +
  labs(x = "Sample", y = "Proportion of Cells") +
  facet_grid(. ~ groupings, space = "free", scales = "free_x") +
  guides(fill = guide_legend(override.aes = list(size = 2),
                             title.position = "right", nrow = 5,
                             label.theme = element_text(size = 5),
                             by_row = TRUE)) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(color = "black", size = 4),
        axis.title.y = element_text(color = "black", face = "bold", size = 6),
        strip.background = element_rect(fill = "grey90", color = "black"),
        strip.text = element_text(size = 7, face = "bold"),
        axis.ticks.y = element_line(color = "black"),
        axis.title.x = element_text(color = "black", size = 6, face = "bold"),
        plot.title = element_blank(),
        legend.title = element_blank(),
        legend.key.size = unit(2.5, "mm"),
        legend.text = element_text(color = "black", size = 4),
        legend.position = "top") 


pdf("/home/smallapragada/manuscript_2024/sample_barchart_ct_supp.pdf", width = 7.15, height = 3)
sample_plot
dev.off()

### Supplementary figure 5

## Sample barchart of AT2 vs immature AT2

## Summarize the number of cells per cell type and lineage
cell_counts_at2 <- run2_final@meta.data %>%
  filter(dev_stages != "Rare disease") %>%
  filter(CT_final %in% c("Immature AT2", "AT2")) %>%
  group_by(Lineage, CT_final, Sample_combined) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  group_by(Sample_combined) %>%
  mutate(prop = cell_count/sum(cell_count))

## Order 
cell_counts_at2$CT_final <- factor(cell_counts_at2$CT_final, level = c("Immature AT2", "AT2"))

## Create the bar plot with proportional facet widths
plot_at2 <- ggplot(cell_counts_at2, aes(x = Sample_combined, y = prop, fill = CT_final)) +
  geom_bar(position = "stack", stat = "identity", width = 0.5) +  # Set bar width
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_manual(values = c("blue", "red"), limits = c("Immature AT2", "AT2")) + 
  theme_classic() + 
  labs(y = "Proportion of cells", title = " ", fill = "Cell Type") +
  guides(fill = guide_legend(override.aes = list(size = 2),
                             title.position = "right", nrow = 1,
                             label.theme = element_text(size = 6),
                             by_row = TRUE)) +
  theme(axis.text.x = element_text(color = "black", size = 8, angle = 90),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.key.size = unit(3, "mm"),
        legend.text = element_text(color = "black", size = 6),
        legend.position = "top",
        axis.ticks.x = element_line(color = "black"),
        axis.text.y = element_text(color = "black", size = 6),
        axis.title.y = element_text(color = "black", size = 7, face = "bold"),
        strip.text = element_blank(),
        axis.ticks.y = element_line(color = "black"),
        axis.title.x = element_blank())

pdf("/home/smallapragada/manuscript_2024/bar_at2_allcts_supp.pdf", width = 7.15, height = 4.5)
plot_at2
dev.off()

## Violin plot

## Extract expression data
data <- FetchData(run2_final, vars = c("SOX9",  "CT_final"))

## Reshape data for ggplot
data_melted <- reshape2::melt(data, id.vars = "CT_final", variable.name = "Gene", value.name = "Expression")

## Ordering cell types
data_melted$CT_final <- factor(data_melted$CT_final, level = c( "Immature AT2", "Proliferating Immature alveolar",
                                                                "AT1", "AT2", "Transitional AT2", 
                                                               "Secretory 3A2+ & 1A1+", "Secretory MUC5B+", "Proliferating airway", 
                                                               "Basal", "Proliferating basal", 
                                                               "Multiciliated", "PNEC", "AKR1C1+ & AKR1C2+", "Arterial", 
                                                               "Capillaries", "Venous", "Lymphatic endothelial", "Proliferating endothelial",
                                                               "Fibroblasts", "Activated fibroblasts", "Adventitial fibroblasts", 
                                                               "MyoFB", "SMC", "Pericytes", 
                                                               "Alveolar macrophages", "SPP1+ macrophages", "Monocytes", "Proliferating monocytes",
                                                               "Neutrophils", "Mast", "cDC", "pDC", "Plasma", "Meg-Ery", "Proliferating Meg-Ery", 
                                                               "B", "T", "Treg", "NK & NKT", "Proliferating lymphoid"))

## Upper threshold

data_melted <- data_melted %>%
  filter(Expression < 1.5)

plot_exp <- ggplot(data_melted, aes(x = CT_final, y = Expression, fill = Gene)) +
  geom_violin(position = position_dodge(1), alpha = 0.6, scale = "width") +
  theme_classic() + 
  scale_fill_manual(values= c("#D1B12D")) +
  labs(y = "Gene expression", title = " ", fill = "Gene") +
  theme(axis.text.x = element_text(color = "black", size = 8, angle = 90, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.key.size = unit(3, "mm"),
        legend.text = element_text(color = "black", size = 6),
        axis.ticks.x = element_line(color = "black"),
        axis.text.y = element_text(color = "black", size = 6),
        axis.title.y = element_text(color = "black", size = 8, face = "bold"),
        strip.text = element_blank(),
        axis.ticks.y = element_line(color = "black"),
        axis.title.x = element_blank())

pdf("/home/smallapragada/manuscript_2024/violin_sox9_allcts_supp.pdf", width = 7.8, height = 5)
plot_exp
dev.off()

### Supplementary figure 6 - 14

## Proximity heatmaps

enrich_score <- readRDS("/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/proximity_enrichment_30r_2025_04_21.rds")
enrich_score_ec <- readRDS("/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/proximity_enrichment_30r_ec_2025_04_21.rds")
enrich_score_lc <- readRDS("/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/proximity_enrichment_30r_lc_2025_04_21.rds")
enrich_score_s <- readRDS("/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/proximity_enrichment_30r_s_2025_04_21.rds")
enrich_score_a <- readRDS("/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/proximity_enrichment_30r_a_2025_04_21.rds")
enrich_score_rd_c <- readRDS("/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/proximity_enrichment_30r_rdc_2025_04_21.rds")
enrich_score_rd_ph <- readRDS("/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/proximity_enrichment_30r_rdph_2025_04_21.rds")
enrich_score_rd_i <- readRDS("/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/proximity_enrichment_30r_rdi_2025_04_21.rds")
enrich_score_rd_a <- readRDS("/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/proximity_enrichment_30r_rda_2025_04_21.rds")

## All samples
adj_p_val <- reshape2::dcast(source_ct ~ neighboring_ct, value.var = 'adj_p_val', data = enrich_score) 
odds_ratio <- reshape2::dcast(source_ct ~ neighboring_ct, value.var = 'logOR', data = enrich_score)

rownames(adj_p_val) <- adj_p_val$source_ct
adj_p_val$source_ct <- NULL

rownames(odds_ratio) <- odds_ratio$source_ct
odds_ratio$source_ct <- NULL

ct_adj_pval <- adj_p_val
ct_adj_pval[is.na(ct_adj_pval)] = 1

ct_or <- odds_ratio
ct_or[is.na(ct_or)] = 0
ct_or[is.infinite(as.matrix(ct_or))] <- 0

mean(rownames(ct_or) == rownames(ct_adj_pval))
mean(colnames(ct_or) == colnames(ct_adj_pval))

ct_or[ct_adj_pval > 0.1] = 0 # replace non-sig results with 0

# plot odds ratio
min01 <- quantile(as.matrix(ct_or), 0.01)
max99 <- quantile(as.matrix(ct_or), 0.99)
col_fun = colorRamp2(c(min01, 0, max99), c("#377EB8", "white", "#E41A1C"))

hp <- Heatmap(t(ct_or), 
              name = 'log(OR)', 
              row_title = 'Source Cell', 
              row_title_side = 'left', 
              row_names_side = 'left', 
              row_dend_side = 'right',  
              row_title_gp = gpar(fontsize = 8), 
              column_title = 'Nearest Neighbor Cell', 
              column_title_side = "bottom",
              column_names_side = 'top', 
              column_names_rot = 45,
              column_dend_side = 'bottom',  
              column_title_gp = gpar(fontsize = 8),
              row_names_gp = gpar(fontsize = 6),
              column_names_gp = gpar(fontsize = 6),
              rect_gp = gpar(col = "#b8b7b6", lwd = 0.5),
              na_col = 'white', 
              col = col_fun,
              heatmap_legend_param = list(
                title_gp = gpar(fontsize = 6),
                labels_gp = gpar(fontsize = 4),
                legend_height = unit(2, "cm"),  
                grid_width = unit(0.3, "cm")   
              )
)

pdf("/home/smallapragada/manuscript_2024/prox_all_samples_supp.pdf", width = 7.3, height = 5.17)
draw(hp)
dev.off()

## Early human
adj_p_val <- reshape2::dcast(source_ct ~ neighboring_ct, value.var = 'adj_p_val', data = enrich_score_ec) 
odds_ratio <- reshape2::dcast(source_ct ~ neighboring_ct, value.var = 'logOR', data = enrich_score_ec)

rownames(adj_p_val) <- adj_p_val$source_ct
adj_p_val$source_ct <- NULL

rownames(odds_ratio) <- odds_ratio$source_ct
odds_ratio$source_ct <- NULL

ct_adj_pval <- adj_p_val
ct_adj_pval[is.na(ct_adj_pval)] = 1

ct_or <- odds_ratio %>%
  mutate_if(is.numeric, list(~na_if(., -Inf))) %>%
  as.matrix()
# ct_or[is.na(ct_or)] = 0
# ct_or[is.infinite(as.matrix(ct_or))] <- 0

mean(rownames(ct_or) == rownames(ct_adj_pval))
mean(colnames(ct_or) == colnames(ct_adj_pval))

ct_or[ct_adj_pval > 0.1] = 0 # replace non-sig results with 0

# plot odds ratio
min01 <- quantile(as.matrix(ct_or), 0.01)
max99 <- quantile(as.matrix(ct_or), 0.99)
col_fun = colorRamp2(c(min01, 0, max99), c("#377EB8", "white", "#E41A1C"))

hp <- Heatmap(t(ct_or), 
              name = 'log(OR)', 
              row_title = 'Source Cell', 
              row_title_side = 'left', 
              row_names_side = 'left', 
              row_dend_side = 'right',  
              row_title_gp = gpar(fontsize = 8), 
              column_title = 'Nearest Neighbor Cell', 
              column_title_side = "bottom",
              column_names_side = 'top', 
              column_names_rot = 45,
              column_dend_side = 'bottom',  
              column_title_gp = gpar(fontsize = 8),
              row_names_gp = gpar(fontsize = 6),
              column_names_gp = gpar(fontsize = 6),
              rect_gp = gpar(col = "#b8b7b6", lwd = 0.5),
              na_col = 'white', 
              col = col_fun,
              heatmap_legend_param = list(
                title_gp = gpar(fontsize = 6),
                labels_gp = gpar(fontsize = 4),
                legend_height = unit(2, "cm"),  
                grid_width = unit(0.3, "cm")   
              )
)

pdf("/home/smallapragada/manuscript_2024/prox_ec_samples_supp.pdf", width = 7.3, height = 5.17)
draw(hp)
dev.off()

## Late canalicular
adj_p_val <- reshape2::dcast(source_ct ~ neighboring_ct, value.var = 'adj_p_val', data = enrich_score_lc) 
odds_ratio <- reshape2::dcast(source_ct ~ neighboring_ct, value.var = 'logOR', data = enrich_score_lc)

rownames(adj_p_val) <- adj_p_val$source_ct
adj_p_val$source_ct <- NULL

rownames(odds_ratio) <- odds_ratio$source_ct
odds_ratio$source_ct <- NULL

ct_adj_pval <- adj_p_val
ct_adj_pval[is.na(ct_adj_pval)] = 1

ct_or <- odds_ratio
ct_or[is.na(ct_or)] = 0
ct_or[is.infinite(as.matrix(ct_or))] <- 0

mean(rownames(ct_or) == rownames(ct_adj_pval))
mean(colnames(ct_or) == colnames(ct_adj_pval))

ct_or[ct_adj_pval > 0.1] = 0 # replace non-sig results with 0

# plot odds ratio
min01 <- quantile(as.matrix(ct_or), 0.01)
max99 <- quantile(as.matrix(ct_or), 0.99)
col_fun = colorRamp2(c(min01, 0, max99), c("#377EB8", "white", "#E41A1C"))

hp <- Heatmap(t(ct_or), 
              name = 'log(OR)', 
              row_title = 'Source Cell', 
              row_title_side = 'left', 
              row_names_side = 'left', 
              row_dend_side = 'right',  
              row_title_gp = gpar(fontsize = 8), 
              column_title = 'Nearest Neighbor Cell', 
              column_title_side = "bottom",
              column_names_side = 'top', 
              column_names_rot = 45,
              column_dend_side = 'bottom',  
              column_title_gp = gpar(fontsize = 8),
              row_names_gp = gpar(fontsize = 6),
              column_names_gp = gpar(fontsize = 6),
              rect_gp = gpar(col = "#b8b7b6", lwd = 0.5),
              na_col = 'white', 
              col = col_fun,
              heatmap_legend_param = list(
                title_gp = gpar(fontsize = 6),
                labels_gp = gpar(fontsize = 4),
                legend_height = unit(2, "cm"),  
                grid_width = unit(0.3, "cm")   
              )
)

pdf("/home/smallapragada/manuscript_2024/prox_lc_samples_supp.pdf", width = 7.3, height = 5.17)
draw(hp)
dev.off()

## Saccular
adj_p_val <- reshape2::dcast(source_ct ~ neighboring_ct, value.var = 'adj_p_val', data = enrich_score_s) 
odds_ratio <- reshape2::dcast(source_ct ~ neighboring_ct, value.var = 'logOR', data = enrich_score_s)

rownames(adj_p_val) <- adj_p_val$source_ct
adj_p_val$source_ct <- NULL

rownames(odds_ratio) <- odds_ratio$source_ct
odds_ratio$source_ct <- NULL

ct_adj_pval <- adj_p_val
ct_adj_pval[is.na(ct_adj_pval)] = 1

ct_or <- odds_ratio
ct_or[is.na(ct_or)] = 0
ct_or[is.infinite(as.matrix(ct_or))] <- 0

mean(rownames(ct_or) == rownames(ct_adj_pval))
mean(colnames(ct_or) == colnames(ct_adj_pval))

ct_or[ct_adj_pval > 0.1] = 0 # replace non-sig results with 0

calc_heatmap_w <- function(df){width <- ncol(df) * 0.24; return(width)}
calc_heatmap_h <- function(df){height <- nrow(df) * 0.18; return(height)}

# plot odds ratio
min01 <- quantile(as.matrix(ct_or), 0.01)
max99 <- quantile(as.matrix(ct_or), 0.99)
col_fun = colorRamp2(c(min01, 0, max99), c("#377EB8", "white", "#E41A1C"))

hp <- Heatmap(t(ct_or), 
              name = 'log(OR)', 
              row_title = 'Source Cell', 
              row_title_side = 'left', 
              row_names_side = 'left', 
              row_dend_side = 'right',  
              row_title_gp = gpar(fontsize = 8), 
              column_title = 'Nearest Neighbor Cell', 
              column_title_side = "bottom",
              column_names_side = 'top', 
              column_names_rot = 45,
              column_dend_side = 'bottom',  
              column_title_gp = gpar(fontsize = 8),
              row_names_gp = gpar(fontsize = 6),
              column_names_gp = gpar(fontsize = 6),
              rect_gp = gpar(col = "#b8b7b6", lwd = 0.5),
              na_col = 'white', 
              col = col_fun,
              heatmap_legend_param = list(
                title_gp = gpar(fontsize = 6),
                labels_gp = gpar(fontsize = 4),
                legend_height = unit(2, "cm"),  
                grid_width = unit(0.3, "cm")   
              )
)

pdf("/home/smallapragada/manuscript_2024/prox_s_samples_supp.pdf", width = 7.3, height = 5.17)
draw(hp)
dev.off()

## Alveolar
adj_p_val <- reshape2::dcast(source_ct ~ neighboring_ct, value.var = 'adj_p_val', data = enrich_score_a) 
odds_ratio <- reshape2::dcast(source_ct ~ neighboring_ct, value.var = 'logOR', data = enrich_score_a)

rownames(adj_p_val) <- adj_p_val$source_ct
adj_p_val$source_ct <- NULL

rownames(odds_ratio) <- odds_ratio$source_ct
odds_ratio$source_ct <- NULL

ct_adj_pval <- adj_p_val
ct_adj_pval[is.na(ct_adj_pval)] = 1

ct_or <- odds_ratio
ct_or[is.na(ct_or)] = 0
ct_or[is.infinite(as.matrix(ct_or))] <- 0

mean(rownames(ct_or) == rownames(ct_adj_pval))
mean(colnames(ct_or) == colnames(ct_adj_pval))

ct_or[ct_adj_pval > 0.1] = 0 # replace non-sig results with 0

calc_heatmap_w <- function(df){width <- ncol(df) * 0.24; return(width)}
calc_heatmap_h <- function(df){height <- nrow(df) * 0.18; return(height)}

# plot odds ratio
min01 <- quantile(as.matrix(ct_or), 0.01)
max99 <- quantile(as.matrix(ct_or), 0.99)
col_fun = colorRamp2(c(min01, 0, max99), c("#377EB8", "white", "#E41A1C"))

hp <- Heatmap(t(ct_or), 
              name = 'log(OR)', 
              row_title = 'Source Cell', 
              row_title_side = 'left', 
              row_names_side = 'left', 
              row_dend_side = 'right',  
              row_title_gp = gpar(fontsize = 8), 
              column_title = 'Nearest Neighbor Cell', 
              column_title_side = "bottom",
              column_names_side = 'top', 
              column_names_rot = 45,
              column_dend_side = 'bottom',  
              column_title_gp = gpar(fontsize = 8),
              row_names_gp = gpar(fontsize = 6),
              column_names_gp = gpar(fontsize = 6),
              rect_gp = gpar(col = "#b8b7b6", lwd = 0.5),
              na_col = 'white', 
              col = col_fun,
              heatmap_legend_param = list(
                title_gp = gpar(fontsize = 6),
                labels_gp = gpar(fontsize = 4),
                legend_height = unit(2, "cm"),  
                grid_width = unit(0.3, "cm")   
              )
)

pdf("/home/smallapragada/manuscript_2024/prox_a_samples_supp.pdf", width = 7.3, height = 5.17)
draw(hp)
dev.off()

## CHAOS
adj_p_val <- reshape2::dcast(source_ct ~ neighboring_ct, value.var = 'adj_p_val', data = enrich_score_rd_c) 
odds_ratio <- reshape2::dcast(source_ct ~ neighboring_ct, value.var = 'logOR', data = enrich_score_rd_c)

rownames(adj_p_val) <- adj_p_val$source_ct
adj_p_val$source_ct <- NULL

rownames(odds_ratio) <- odds_ratio$source_ct
odds_ratio$source_ct <- NULL

ct_adj_pval <- adj_p_val
ct_adj_pval[is.na(ct_adj_pval)] = 1

ct_or <- odds_ratio
ct_or[is.na(ct_or)] = 0
ct_or[is.infinite(as.matrix(ct_or))] <- 0

mean(rownames(ct_or) == rownames(ct_adj_pval))
mean(colnames(ct_or) == colnames(ct_adj_pval))

ct_or[ct_adj_pval > 0.1] = 0 # replace non-sig results with 0

calc_heatmap_w <- function(df){width <- ncol(df) * 0.24; return(width)}
calc_heatmap_h <- function(df){height <- nrow(df) * 0.18; return(height)}

# plot odds ratio
min01 <- quantile(as.matrix(ct_or), 0.01)
max99 <- quantile(as.matrix(ct_or), 0.99)
col_fun = colorRamp2(c(min01, 0, max99), c("#377EB8", "white", "#E41A1C"))

hp <- Heatmap(t(ct_or), 
              name = 'log(OR)', 
              row_title = 'Source Cell', 
              row_title_side = 'left', 
              row_names_side = 'left', 
              row_dend_side = 'right',  
              row_title_gp = gpar(fontsize = 8), 
              column_title = 'Nearest Neighbor Cell', 
              column_title_side = "bottom",
              column_names_side = 'top', 
              column_names_rot = 45,
              column_dend_side = 'bottom',  
              column_title_gp = gpar(fontsize = 8),
              row_names_gp = gpar(fontsize = 6),
              column_names_gp = gpar(fontsize = 6),
              rect_gp = gpar(col = "#b8b7b6", lwd = 0.5),
              na_col = 'white', 
              col = col_fun,
              heatmap_legend_param = list(
                title_gp = gpar(fontsize = 6),
                labels_gp = gpar(fontsize = 4),
                legend_height = unit(2, "cm"),  
                grid_width = unit(0.3, "cm")   
              )
)

pdf("/home/smallapragada/manuscript_2024/prox_rd_c_samples_supp.pdf", width = 7.3, height = 5.17)
draw(hp)
dev.off()

## PH + kidneys
adj_p_val <- reshape2::dcast(source_ct ~ neighboring_ct, value.var = 'adj_p_val', data = enrich_score_rd_ph) 
odds_ratio <- reshape2::dcast(source_ct ~ neighboring_ct, value.var = 'logOR', data = enrich_score_rd_ph)

rownames(adj_p_val) <- adj_p_val$source_ct
adj_p_val$source_ct <- NULL

rownames(odds_ratio) <- odds_ratio$source_ct
odds_ratio$source_ct <- NULL

ct_adj_pval <- adj_p_val
ct_adj_pval[is.na(ct_adj_pval)] = 1

ct_or <- odds_ratio
ct_or[is.na(ct_or)] = 0
ct_or[is.infinite(as.matrix(ct_or))] <- 0

mean(rownames(ct_or) == rownames(ct_adj_pval))
mean(colnames(ct_or) == colnames(ct_adj_pval))

ct_or[ct_adj_pval > 0.1] = 0 # replace non-sig results with 0

calc_heatmap_w <- function(df){width <- ncol(df) * 0.24; return(width)}
calc_heatmap_h <- function(df){height <- nrow(df) * 0.18; return(height)}

# plot odds ratio
min01 <- quantile(as.matrix(ct_or), 0.01)
max99 <- quantile(as.matrix(ct_or), 0.99)
col_fun = colorRamp2(c(min01, 0, max99), c("#377EB8", "white", "#E41A1C"))

hp <- Heatmap(t(ct_or), 
              name = 'log(OR)', 
              row_title = 'Source Cell', 
              row_title_side = 'left', 
              row_names_side = 'left', 
              row_dend_side = 'right',  
              row_title_gp = gpar(fontsize = 8), 
              column_title = 'Nearest Neighbor Cell', 
              column_title_side = "bottom",
              column_names_side = 'top', 
              column_names_rot = 45,
              column_dend_side = 'bottom',  
              column_title_gp = gpar(fontsize = 8),
              row_names_gp = gpar(fontsize = 6),
              column_names_gp = gpar(fontsize = 6),
              rect_gp = gpar(col = "#b8b7b6", lwd = 0.5),
              na_col = 'white', 
              col = col_fun,
              heatmap_legend_param = list(
                title_gp = gpar(fontsize = 6),
                labels_gp = gpar(fontsize = 4),
                legend_height = unit(2, "cm"),  
                grid_width = unit(0.3, "cm")   
              )
)

pdf("/home/smallapragada/manuscript_2024/prox_rd_ph_samples_supp.pdf", width = 7.3, height = 5.17)
draw(hp)
dev.off()

## Infant BPD
adj_p_val <- reshape2::dcast(source_ct ~ neighboring_ct, value.var = 'adj_p_val', data = enrich_score_rd_i) 
odds_ratio <- reshape2::dcast(source_ct ~ neighboring_ct, value.var = 'logOR', data = enrich_score_rd_i)

rownames(adj_p_val) <- adj_p_val$source_ct
adj_p_val$source_ct <- NULL

rownames(odds_ratio) <- odds_ratio$source_ct
odds_ratio$source_ct <- NULL

ct_adj_pval <- adj_p_val
ct_adj_pval[is.na(ct_adj_pval)] = 1

ct_or <- odds_ratio
ct_or[is.na(ct_or)] = 0
ct_or[is.infinite(as.matrix(ct_or))] <- 0

mean(rownames(ct_or) == rownames(ct_adj_pval))
mean(colnames(ct_or) == colnames(ct_adj_pval))

ct_or[ct_adj_pval > 0.1] = 0 # replace non-sig results with 0

calc_heatmap_w <- function(df){width <- ncol(df) * 0.24; return(width)}
calc_heatmap_h <- function(df){height <- nrow(df) * 0.18; return(height)}

# plot odds ratio
min01 <- quantile(as.matrix(ct_or), 0.01)
max99 <- quantile(as.matrix(ct_or), 0.99)
col_fun = colorRamp2(c(min01, 0, max99), c("#377EB8", "white", "#E41A1C"))

hp <- Heatmap(t(ct_or), 
              name = 'log(OR)', 
              row_title = 'Source Cell', 
              row_title_side = 'left', 
              row_names_side = 'left', 
              row_dend_side = 'right',  
              row_title_gp = gpar(fontsize = 8), 
              column_title = 'Nearest Neighbor Cell', 
              column_title_side = "bottom",
              column_names_side = 'top', 
              column_names_rot = 45,
              column_dend_side = 'bottom',  
              column_title_gp = gpar(fontsize = 8),
              row_names_gp = gpar(fontsize = 6),
              column_names_gp = gpar(fontsize = 6),
              rect_gp = gpar(col = "#b8b7b6", lwd = 0.5),
              na_col = 'white', 
              col = col_fun,
              heatmap_legend_param = list(
                title_gp = gpar(fontsize = 6),
                labels_gp = gpar(fontsize = 4),
                legend_height = unit(2, "cm"),  
                grid_width = unit(0.3, "cm")   
              )
)

pdf("/home/smallapragada/manuscript_2024/prox_rd_i_samples_supp.pdf", width = 7.3, height = 5.17)
draw(hp)
dev.off()

## Adult BPD
adj_p_val <- reshape2::dcast(source_ct ~ neighboring_ct, value.var = 'adj_p_val', data = enrich_score_rd_a) 
odds_ratio <- reshape2::dcast(source_ct ~ neighboring_ct, value.var = 'logOR', data = enrich_score_rd_a)

rownames(adj_p_val) <- adj_p_val$source_ct
adj_p_val$source_ct <- NULL

rownames(odds_ratio) <- odds_ratio$source_ct
odds_ratio$source_ct <- NULL

ct_adj_pval <- adj_p_val
ct_adj_pval[is.na(ct_adj_pval)] = 1

ct_or <- odds_ratio
ct_or[is.na(ct_or)] = 0
ct_or[is.infinite(as.matrix(ct_or))] <- 0

mean(rownames(ct_or) == rownames(ct_adj_pval))
mean(colnames(ct_or) == colnames(ct_adj_pval))

ct_or[ct_adj_pval > 0.1] = 0 # replace non-sig results with 0

calc_heatmap_w <- function(df){width <- ncol(df) * 0.24; return(width)}
calc_heatmap_h <- function(df){height <- nrow(df) * 0.18; return(height)}

# plot odds ratio
min01 <- quantile(as.matrix(ct_or), 0.01)
max99 <- quantile(as.matrix(ct_or), 0.99)
col_fun = colorRamp2(c(min01, 0, max99), c("#377EB8", "white", "#E41A1C"))

hp <- Heatmap(t(ct_or), 
              name = 'log(OR)', 
              row_title = 'Source Cell', 
              row_title_side = 'left', 
              row_names_side = 'left', 
              row_dend_side = 'right',  
              row_title_gp = gpar(fontsize = 8), 
              column_title = 'Nearest Neighbor Cell', 
              column_title_side = "bottom",
              column_names_side = 'top', 
              column_names_rot = 45,
              column_dend_side = 'bottom',  
              column_title_gp = gpar(fontsize = 8),
              row_names_gp = gpar(fontsize = 6),
              column_names_gp = gpar(fontsize = 6),
              rect_gp = gpar(col = "#b8b7b6", lwd = 0.5),
              na_col = 'white', 
              col = col_fun,
              heatmap_legend_param = list(
                title_gp = gpar(fontsize = 6),
                labels_gp = gpar(fontsize = 4),
                legend_height = unit(2, "cm"),  
                grid_width = unit(0.3, "cm")   
              )
)

pdf("/home/smallapragada/manuscript_2024/prox_rd_a_samples_supp.pdf", width = 7.3, height = 5.17)
draw(hp)
dev.off()

### Supplementary figure 15

plot <- DimPlot(cell_niches,
                reduction = "sp",
                group.by = "CNiches", 
                raster = T,
                cols = cniche_colors) +
  coord_equal() + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_blank(),
        axis.title = element_blank(),
        line = element_blank(),
        legend.key.size = unit(0.3, "cm"),     
        legend.text = element_text(size = 6),      
        legend.spacing.y = unit(0.1, 'cm')         
  ) +
  guides(color = guide_legend(ncol = 1))
plot

pdf("/home/smallapragada/manuscript_2024/sp_both_slides_niches_supp.pdf", width = 7, height = 5.5)
plot
dev.off()

run2_final@meta.data$TNiches <- factor(run2_final@meta.data$TNiches, level = c("T1", "T2",
                                                                               "T3", "T4",
                                                                               "T5", "T6",
                                                                               "T7", "T8",
                                                                               "T9", "T10",
                                                                               "T11"))

plot <- DimPlot(run2_final,
                reduction = "sp",
                group.by = "TNiches", 
                raster = T,
                cols = tniche_colors) +
  coord_equal() + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_blank(),
        axis.title = element_blank(),
        line = element_blank(),
        legend.key.size = unit(0.3, "cm"),     
        legend.text = element_text(size = 6),      
        legend.spacing.y = unit(0.1, 'cm')         
  ) +
  guides(color = guide_legend(ncol = 1))

pdf("/home/smallapragada/manuscript_2024/sp_both_slides_niches_transcript_supp.pdf", width = 7, height = 5.5)
plot
dev.off()

### Supplementary figure 16

## Ordering bars
run2_final$Sample_combined <- factor(run2_final$Sample_combined, level = c("PDL001", "PDL002", "PDL003",
                                                                             "PDL004", "PDL005", "PDL006",
                                                                             "PDL007", "PDL008", "PDL009",
                                                                             "PDL010", "PDL011", "PDL012",
                                                                             "PDL013", "PDL014", "PDL015",
                                                                             "PDL016", "PDL017"))

## Summarize the number of cells per niche and dev stage
cell_counts_cniche <- run2_final@meta.data %>%
  group_by(CNiches, Sample_combined) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  group_by(Sample_combined) %>%
  mutate(prop = cell_count/sum(cell_count))

## Ordering CNiches 

cell_counts_cniche$CNiches <- factor(cell_counts_cniche$CNiches, level = c("C1", "C2", "C3", 
                                                                           "C4", "C5", "C6",
                                                                           "C7", "C8", "C9",
                                                                           "C10", "C11"))

cniche_plot <- ggplot(cell_counts_cniche, aes(x = cell_counts_cniche$Sample_combined, y = cell_counts_cniche$prop,
                                              fill = as.factor(cell_counts_cniche$CNiches))) +
  scale_y_continuous(expand = c(0, 0)) +  
  geom_bar(position = "stack", stat = "identity", width = 0.75) +
  scale_fill_manual(values = cniche_colors) +
  labs(title = "Cell niches", x = "Sample", y = "Proportion of Niches") +
  guides(fill = guide_legend(title.position = "right", ncol = 1,
                             label.theme = element_text(size = 5),
                             by_row = TRUE)) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.t3456ext.y = element_text(color = "black", size = 4),
        axis.title.y = element_text(color = "black", size = 6, face = "bold"),
        strip.background = element_rect(fill = "grey90", color = "black"),
        strip.text = element_text(size = 6, face = "bold"),
        axis.ticks.y = element_line(color = "black"),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 6, face = "bold", hjust = 0.5),
        legend.title = element_blank(),
        legend.key.size = unit(3, "mm"),
        legend.text = element_text(color = "black", size = 4),
        legend.position = "right") 

### TNiches

## Summarize the number of cells per niche and dev stage
cell_counts_tniche <- run2_final@meta.data %>%
  group_by(TNiches, Sample_combined) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  group_by(Sample_combined) %>%
  mutate(prop = cell_count/sum(cell_count))

## Ordering TNiches 
cell_counts_tniche$TNiches <- factor(cell_counts_tniche$TNiches, level = c("T1", "T2", "T3", "T4", "T5", "T6",
                                                                         "T7", "T8", "T9", "T10", "T11"))

tniche_plot <- ggplot(cell_counts_tniche, aes(x = cell_counts_tniche$Sample_combined, y = cell_counts_tniche$prop,
                                              fill = as.factor(cell_counts_tniche$TNiches))) +
  scale_y_continuous(expand = c(0, 0)) +  
  geom_bar(position = "stack", stat = "identity", width = 0.75) +
  scale_fill_manual(values = tniche_colors) +
  labs(title = "Transcript niches", x = "Sample", y = "Proportion of Niches") +
  guides(fill = guide_legend(title.position = "right", ncol = 1,
                             label.theme = element_text(size = 5),
                             by_row = TRUE)) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_line(color = "black"),
        axis.text.x = element_text(color = "black", size = 5, hjust = 0.5),
        axis.text.y = element_text(color = "black", size = 4),
        axis.title.y = element_text(color = "black", size = 6, face = "bold"),
        strip.background = element_rect(fill = "grey90", color = "black"),
        strip.text = element_text(size = 6, face = "bold"),
        axis.ticks.y = element_line(color = "black"),
        axis.title.x = element_text(color = "black", size = 6, face = "bold"),
        plot.title = element_text(size = 6, face = "bold", hjust = 0.5),
        legend.title = element_blank(),
        legend.key.size = unit(3, "mm"),
        legend.text = element_text(color = "black", size = 4),
        legend.position = "right")

## Combining plots
c_t_bar <- cniche_plot / tniche_plot

pdf("/home/smallapragada/manuscript_2024/niche_barchart_sample_supp.pdf", width = 7.5, height = 6)
c_t_bar
dev.off()

### Supplementary figure 17 - Non-immune lineages

### Volcano plots

## Pick the top three per category
top_sig_genes <- all_stats_ga_sig %>%
  filter(sig_ga == "Significant") %>%
  arrange(desc(ga_coef_lfc)) %>%
  slice_head(n = 5) %>%
  select(gene_ct) %>%
  mutate(gene_col = gene_ct)

## Only label those
test <- all_stats_ga_sig %>%
  left_join(top_sig_genes %>% select(gene_ct, gene_col),
            by = c("gene_ct"))

plot <- ggplot(test, aes(x = ga_coef_lfc, y = -log10(adj_p_val_ga), alpha = sig_ga)) +
  geom_point(size = 0.5, color = "#0073C2FF") + 
  geom_label_repel(aes(label = gene_col), 
                   size = 1.8, 
                   max.overlaps = Inf, 
                   box.padding = 0.25, 
                   point.padding = 0.25,
                   segment.size = 0.2,   
                   min.segment.length = 0,
                   show.legend = FALSE) +
  scale_alpha_manual(values = c("Significant" = 1, "Not significant" = 0.1)) +
  xlab("logFC") + ylab("-log10(adj_p_val)") +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black", size = 4),
        axis.text.y = element_text(color = "black", size = 4),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "black", face = "bold", size = 6),
        panel.grid.major.x = element_line(color = "grey80", linewidth = 0.4),
        panel.grid.major.y = element_line(color = "grey80", linewidth = 0.4),
        axis.ticks.x = element_line(color = "black"),
        axis.ticks.y = element_line(color = "black"),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(color = "black", size = 6))

pdf("/home/smallapragada/manuscript_2024/volcano_ga_supp_nonimm.pdf", width = 3, height = 3.3)
plot
dev.off()

## Pick the top three per category
top_sig_genes <- all_stats_ls_sig %>%
  filter(sig_ls == "Significant") %>%
  arrange(desc(ls_coef_lfc)) %>%
  slice_head(n = 5) %>%
  select(gene_ct) %>%
  mutate(gene_col = gene_ct)

## Only label those
test <- all_stats_ls_sig %>%
  left_join(top_sig_genes %>% select(gene_ct, gene_col),
            by = c("gene_ct"))

plot <- ggplot(test, aes(x = ls_coef_lfc, y = -log10(adj_p_val_ls), alpha = sig_ls)) +
  geom_point(size = 0.5, color = "#EFC000FF") + 
  geom_label_repel(aes(label = gene_col), 
                   size = 1.8, 
                   max.overlaps = Inf, 
                   box.padding = 0.25, 
                   point.padding = 0.25,
                   segment.size = 0.2,   
                   min.segment.length = 0,
                   show.legend = FALSE) +
  scale_alpha_manual(values = c("Significant" = 1, "Not significant" = 0.1)) +
  xlab("logFC") + ylab("-log10(adj_p_val)") +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black", size = 4),
        axis.text.y = element_text(color = "black", size = 4),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "black", face = "bold", size = 6),
        panel.grid.major.x = element_line(color = "grey80", linewidth = 0.4),
        panel.grid.major.y = element_line(color = "grey80", linewidth = 0.4),
        axis.ticks.x = element_line(color = "black"),
        axis.ticks.y = element_line(color = "black"),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(color = "black", size = 6))

pdf("/home/smallapragada/manuscript_2024/volcano_ls_supp_nonimm.pdf", width = 3, height = 3.3)
plot
dev.off()

## Pick the top three per category
top_sig_genes <- dx_stats %>%
  filter(sig_dx == "Significant") %>%
  arrange(desc(dx_coef_lfc)) %>%
  slice_head(n = 5) %>%
  select(gene_ct) %>%
  mutate(gene_col = gene_ct)

## Only label those
test <- dx_stats %>%
  left_join(top_sig_genes %>% select(gene_ct, gene_col),
            by = c("gene_ct"))

plot <- ggplot(test, aes(x = dx_coef_lfc, y = -log10(adj_p_val_dx), alpha = sig_dx)) +
  geom_point(size = 0.5, color = "#CD534CFF") + 
  geom_label_repel(aes(label = gene_col), 
                   size = 1.8, 
                   max.overlaps = Inf, 
                   box.padding = 0.25, 
                   point.padding = 0.25,
                   segment.size = 0.2,   
                   min.segment.length = 0,
                   show.legend = FALSE) +
  scale_alpha_manual(values = c("Significant" = 1, "Not significant" = 0.1)) +
  xlab("logFC") + ylab("-log10(adj_p_val)") +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black", size = 4),
        axis.text.y = element_text(color = "black", size = 4),
        axis.title.x = element_text(color = "black", face = "bold", size = 6),
        axis.title.y = element_text(color = "black", face = "bold", size = 6),
        panel.grid.major.x = element_line(color = "grey80", linewidth = 0.4),
        panel.grid.major.y = element_line(color = "grey80", linewidth = 0.4),
        axis.ticks.x = element_line(color = "black"),
        axis.ticks.y = element_line(color = "black"),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(color = "black", size = 6))

pdf("/home/smallapragada/manuscript_2024/volcano_dx_supp_nonimm.pdf", width = 3, height = 3.4)
plot
dev.off()

### Venn diagram

pairs_list_20 <- list(`Gestational age` = ga_pairs_with_thresh$gene_ct, 
                      `Life span` = ls_pairs_with_thresh$gene_ct, 
                      `DX score` = dx_pairs_with_thresh$gene_ct)

venn <- ggvenn(pairs_list_20, 
               fill_color = c("#0073C2FF", "#EFC000FF", "#CD534CFF"),
               stroke_size = 0.5, 
               fill_alpha = 0.5,
               set_name_size = 2,
               text_size = 2,
               show_percentage = FALSE)

pdf("/home/smallapragada/manuscript_2024/venn_supp_nonimm.pdf", width = 1.75, height = 1.75)
venn
dev.off()

### Supplementary figure 18 - all lineages

### Volcano plots

## Pick the top three per category
top_sig_genes <- all_stats_ga_sig %>%
  filter(sig_ga == "Significant") %>%
  arrange(desc(ga_coef_lfc)) %>%
  slice_head(n = 5) %>%
  select(gene_ct) %>%
  mutate(gene_col = gene_ct)

## Only label those
test <- all_stats_ga_sig %>%
  left_join(top_sig_genes %>% select(gene_ct, gene_col),
            by = c("gene_ct"))

plot <- ggplot(test, aes(x = ga_coef_lfc, y = -log10(adj_p_val_ga), alpha = sig_ga)) +
  geom_point(size = 0.5, color = "#0073C2FF") + 
  geom_label_repel(aes(label = gene_col), 
                   size = 1.8, 
                   max.overlaps = Inf, 
                   box.padding = 0.25, 
                   point.padding = 0.25,
                   segment.size = 0.2,   
                   min.segment.length = 0,
                   show.legend = FALSE) +
  scale_alpha_manual(values = c("Significant" = 1, "Not significant" = 0.1)) +
  xlab("logFC") + ylab("-log10(adj_p_val)") +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black", size = 4),
        axis.text.y = element_text(color = "black", size = 4),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "black", face = "bold", size = 6),
        panel.grid.major.x = element_line(color = "grey80", linewidth = 0.4),
        panel.grid.major.y = element_line(color = "grey80", linewidth = 0.4),
        axis.ticks.x = element_line(color = "black"),
        axis.ticks.y = element_line(color = "black"),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(color = "black", size = 6))

pdf("/home/smallapragada/manuscript_2024/volcano_ga_supp_all.pdf", width = 3, height = 3.3)
plot
dev.off()

## Pick the top three per category
top_sig_genes <- all_stats_ls_sig %>%
  filter(sig_ls == "Significant") %>%
  arrange(desc(ls_coef_lfc)) %>%
  slice_head(n = 5) %>%
  select(gene_ct) %>%
  mutate(gene_col = gene_ct)

## Only label those
test <- all_stats_ls_sig %>%
  left_join(top_sig_genes %>% select(gene_ct, gene_col),
            by = c("gene_ct"))

plot <- ggplot(test, aes(x = ls_coef_lfc, y = -log10(adj_p_val_ls), alpha = sig_ls)) +
  geom_point(size = 0.5, color = "#EFC000FF") + 
  geom_label_repel(aes(label = gene_col), 
                   size = 1.8, 
                   max.overlaps = Inf, 
                   box.padding = 0.25, 
                   point.padding = 0.25,
                   segment.size = 0.2,   
                   min.segment.length = 0,
                   show.legend = FALSE) +
  scale_alpha_manual(values = c("Significant" = 1, "Not significant" = 0.1)) +
  xlab("logFC") + ylab("-log10(adj_p_val)") +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black", size = 4),
        axis.text.y = element_text(color = "black", size = 4),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "black", face = "bold", size = 6),
        panel.grid.major.x = element_line(color = "grey80", linewidth = 0.4),
        panel.grid.major.y = element_line(color = "grey80", linewidth = 0.4),
        axis.ticks.x = element_line(color = "black"),
        axis.ticks.y = element_line(color = "black"),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(color = "black", size = 6))

pdf("/home/smallapragada/manuscript_2024/volcano_ls_supp_all.pdf", width = 3, height = 3.3)
plot
dev.off()

## Pick the top three per category
top_sig_genes <- dx_stats %>%
  filter(sig_dx == "Significant") %>%
  arrange(desc(dx_coef_lfc)) %>%
  slice_head(n = 5) %>%
  select(gene_ct) %>%
  mutate(gene_col = gene_ct)

## Only label those
test <- dx_stats %>%
  left_join(top_sig_genes %>% select(gene_ct, gene_col),
            by = c("gene_ct"))

plot <- ggplot(test, aes(x = dx_coef_lfc, y = -log10(adj_p_val_dx), alpha = sig_dx)) +
  geom_point(size = 0.5, color = "#CD534CFF") + 
  geom_label_repel(aes(label = gene_col), 
                   size = 1.8, 
                   max.overlaps = Inf, 
                   box.padding = 0.25, 
                   point.padding = 0.25,
                   segment.size = 0.2,   
                   min.segment.length = 0,
                   show.legend = FALSE) +
  scale_alpha_manual(values = c("Significant" = 1, "Not significant" = 0.1)) +
  xlab("logFC") + ylab("-log10(adj_p_val)") +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black", size = 4),
        axis.text.y = element_text(color = "black", size = 4),
        axis.title.x = element_text(color = "black", face = "bold", size = 6),
        axis.title.y = element_text(color = "black", face = "bold", size = 6),
        panel.grid.major.x = element_line(color = "grey80", linewidth = 0.4),
        panel.grid.major.y = element_line(color = "grey80", linewidth = 0.4),
        axis.ticks.x = element_line(color = "black"),
        axis.ticks.y = element_line(color = "black"),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(color = "black", size = 6))

pdf("/home/smallapragada/manuscript_2024/volcano_dx_supp_all.pdf", width = 3, height = 3.4)
plot
dev.off()

### Venn diagram

pairs_list_20 <- list(`Gestational age` = ga_pairs_with_thresh$gene_ct, 
                      `Life span` = ls_pairs_with_thresh$gene_ct, 
                      `DX score` = dx_pairs_with_thresh$gene_ct)

venn <- ggvenn(pairs_list_20, 
               fill_color = c("#0073C2FF", "#EFC000FF", "#CD534CFF"),
               stroke_size = 0.5, 
               fill_alpha = 0.5,
               set_name_size = 2,
               text_size = 2,
               show_percentage = FALSE)

pdf("/home/smallapragada/manuscript_2024/venn_supp_all.pdf", width = 1.75, height = 1.75)
venn
dev.off()
