#### LOAD PACKAGES AND SET ENVIRONMENT ----
library(SeuratObject)
library(tidyverse)
library(gplots)
library(rlist)
library(Seurat, lib.loc = "/usr/local/lib/R/site-library")
library(gprofiler2)
library(ggrepel)
library(ggpubr)
library(edgeR)
library(limma)
library(patchwork)
library(presto)
library(scCustomize) 
library(EnhancedVolcano)
library(rstatix)
library(readxl)

## Set seed & working directory
set.seed(0712)
work_dir <- setwd("/scratch/smallapragada/run2")

## Function to adjust angle for x-axis 
angle <- theme(axis.text.x = element_text(angle = 45, hjust = 1))

## "not in" operand function
'%!in%' <- function(x,y)!('%in%'(x,y))

## Pulling the count and metadata files from Xenium output

# Transcript files contain the gene counts per gene per file for both slides
count_files <- list.files(pattern = "transcripts.csv.gz", recursive = TRUE)

# Cells files contain the metadata for each cell per file for both slides
meta_files <- list.files(pattern = "cells.csv.gz", recursive = TRUE)

## QC metric - checking that all of the correct files were loaded in
print(count_files)
print(meta_files)

## Parsing through all files and fixing the formats for both
count_list <- lapply(count_files, function(XX)
  read.delim(XX, sep = ","))
metadata_list <- lapply(meta_files, function(XX)
  read.delim(XX, sep = ","))

#### ASSIGNING SAMPLE NAMES ----

## Creating a list of the sample name per both slides --> listed file path to assign the correct sample to the correct output
id_list <- c(
  # Slide 1 - 0005123
  s1_PDL006 = "/scratch/smallapragada/run2/0005123/relabeled_output-XETG00048__0005123__PDL006__20231120__192537",
  s1_PDL005 = "/scratch/smallapragada/run2/0005123/relabeled_output-XETG00048__0005123__PDL005__20231120__192537",
  s1_PDL015 = "/scratch/smallapragada/run2/0005123/relabeled_output-XETG00048__0005123__PDL015__20231120__192537",
  s1_PDL009 = "/scratch/smallapragada/run2/0005123/relabeled_output-XETG00048__0005123__PDL009__20231120__192537",
  s1_PDL017 = "/scratch/smallapragada/run2/0005123/relabeled_output-XETG00048__0005123__PDL017__20231120__192537",
  s1_GOLIATH = "/scratch/smallapragada/run2/0005123/relabeled_output-XETG00048__0005123__M_Goliath_B1__20231120__192537",
  # Slide 2 - 0005297
  s2_PDL006  = "/scratch/smallapragada/run2/0005297/relabeled_output-XETG00048__0005297__PDL006__20231120__192537",
  s2_PDL015 = "/scratch/smallapragada/run2/0005297/relabeled_output-XETG00048__0005297__PDL015__20231120__192537",
  s2_PDL009 = "/scratch/smallapragada/run2/0005297/relabeled_output-XETG00048__0005297__PDL009__20231011__185720",
  s2_GOLIATH = "/scratch/smallapragada/run2/0005297/relabeled_output-XETG00048__0005297__M_Goliath_A1__20231120__192537")

## Assigning the loaded-in names to the list of count data and metadata 
ids <- names(id_list)
names(count_list) <- ids
names(metadata_list) <- ids

### Generate gene x cell matrix and use filtering metrics

## Empty lists used in the loop below
percent_negcodes <- list()
percent_negprobes <- list()

## For loop that manipulates these data as necessary

for (i in seq_along(ids)) {
  # Printing to view progress while running
  print(i)  
  
  ### Adjusting count data in various ways
  
  # Adjusting "transcript_id" to the numeric character type and removing scientific notation so the full ID is printed
  count_list[[i]]$transcript_id <- format(as.numeric(count_list[[i]]$transcript_id), scientific = F)
  
  # Filtering out transcripts with no nuclear overlap
  count_list[[i]] <- count_list[[i]][count_list[[i]]$overlaps_nucleus == "1", ]
  
  # Filtering out low-quality transcripts with qv > 20
  count_list[[i]] <- count_list[[i]][count_list[[i]]$qv > 20, ]
  
  # Removing "BLANK" genes
  count_list[[i]] <- count_list[[i]][!grepl(c("BLANK"), count_list[[i]]$feature_name), ]
  
  # Removing "Unassigned" genes
  count_list[[i]] <- count_list[[i]][!grepl(c("Unassigned"), count_list[[i]]$feature_name), ]
  
  count_list[[i]] = count_list[[i]] %>%
    select(cell_id, feature_name) %>%                             # Pull cell ids and genes from all samples within count_list
    table() %>%                                                   # Table() will count the number of instances among the two
    as.data.frame() %>%                                           # Make df so I can use dplyr lol
    pivot_wider(names_from = cell_id, values_from = Freq) %>%     # Makes cell ids column names and genes row names
    as.data.frame()                                               # Make df AGAIN cuz it's easier for me to work with whoops
  
  rownames(count_list[[i]]) <- count_list[[i]]$feature_name       # Make feature_names the column names
  count_list[[i]] <- count_list[[i]][, -c(1)]      # Delete feature_names column
  
  ### Calculating the percent of the number of "Negative Control" instances per cell_id
  
  # Calculating the sum of all genes per cell_id
  total_colsums <- as.data.frame(colSums(count_list[[i]]))
  
  # Calculating the sum of the "NegControlCodeword" instances per cell_id
  negcodes_colsums <- count_list[[i]] %>%
    slice(grep("^NegControlCodeword", row.names(.))) %>%
    colSums()
  
  # Calculating the percent of "NegControlCodeword" instances per cell_id
  percent_negcodes[[i]] <- negcodes_colsums/total_colsums
  
  # Calculating the sum of the "NegControlProbes" instances per cell_id
  negprobes_colsums <- count_list[[i]] %>%
    slice(grep("^NegControlProbe", row.names(.))) %>%
    colSums()
  
  # Calculating the percent of "NegControlProbes" instances per cell_id
  percent_negprobes[[i]] <- negprobes_colsums/total_colsums
  
  # Adding the cell_id column back into count_list data frames from above to join them with the metadata 
  percent_negcodes[[i]]$cell_id <- row.names(percent_negcodes[[i]])
  percent_negprobes[[i]]$cell_id <- row.names(percent_negprobes[[i]])
  
  # Adding the calculated percent data as a column in each data frame in metadata_list
  metadata_list[[i]] <- inner_join(metadata_list[[i]], percent_negcodes[[i]], by ='cell_id')
  metadata_list[[i]] <- inner_join(metadata_list[[i]], percent_negprobes[[i]], by ='cell_id')
  
  # Adjusting metadata to the correct format
  rownames(metadata_list[[i]]) <- metadata_list[[i]]$cell_id      # Make cell_ids the row names
  metadata_list[[i]] <- metadata_list[[i]][, -c(1)]      # Delete cell_ids column
  
  # Renaming to "percent_negcodes" and "percent_negprobes"
  names(metadata_list[[i]])[names(metadata_list[[i]]) == 'colSums(count_list[[i]]).x'] <- 'percent_negcodes'
  names(metadata_list[[i]])[names(metadata_list[[i]]) == 'colSums(count_list[[i]]).y'] <- 'percent_negprobes'
  
  # Change data type to matrix for all data frames in count_list
  count_list[[i]] <- as.matrix(count_list[[i]])
  
  # Remove "NegControl" genes from dataset
  count_list[[i]] <- count_list[[i]][!grepl(c("^NegControl"), row.names(count_list[[i]])), ]
  
}

#### GENERATE LIST OF SEURAT OBJECTS PER FILE AND MERGE INTO ONE ----

## Empty list used in the loop below
obj_list <- list()

## For loop that creates the Seurat objects for each file and adds metadata 

for (i in seq_along(ids)) {
  # Create Seurat objects and add metadata accordingly
  seurat_obj <- CreateSeuratObject(counts = count_list[[i]], assay = "RNA")
  seurat_obj <- AddMetaData(seurat_obj, metadata = metadata_list[[i]])
  
  # Add file name (Sample) to metadata
  seurat_obj$Sample <- ids[[i]]
  
  # Adding file name (Sample) to cell_ids
  seurat_obj <- RenameCells(seurat_obj, add.cell.id = seurat_obj$Sample)
  
  # Filter out cells with nCount_RNA == 0
  seurat_obj <- subset(seurat_obj, subset = nCount_RNA != 0)
  
  # Add each Seurat object into a list of objects
  obj_list[[i]] <- seurat_obj
  
}

## Merging all Seurat objects in obj_list together into one 
merged_obj_unfiltered <- merge(x = obj_list[[1]], y = obj_list[2:length(obj_list)])

#### SPLITTING GOLIATH AND ADDING SAMPLE NAMES TO GOLIATH ----

## Pulling all .csv files per slide that contain the groupings for how "Goliath" is split into the various samples 
s1_cell_files <- list.files(path = "/scratch/smallapragada/run2/0005123", pattern = "s1_", full.names = TRUE, recursive = FALSE)
s2_cell_files <- list.files(path = "/scratch/smallapragada/run2/0005297", pattern = "s2_", full.names = TRUE, recursive = FALSE)

## Parsing through files and fixing the format
s1_cells <- lapply(s1_cell_files, function(XX)
  read.delim(XX, sep = ",", skip = 2))
s2_cells <- lapply(s2_cell_files, function(XX)
  read.delim(XX, sep = ",", skip = 2))

## Creating a list of the sample name per both slides --> listed file path to assign the correct sample to the correct output
s1_id_list <- c(
  s1_PDL008 = "/scratch/smallapragada/run2/0005123/s1_PDL008_06_27_2024_cells_stats.csv",
  s1_PDL003 = "/scratch/smallapragada/run2/0005123/s1_PDL003_06_27_2024_cells_stats.csv", 
  s1_PDL016 = "/scratch/smallapragada/run2/0005123/s1_PDL016_06_27_2024_cells_stats.csv",
  s1_PDL011 = "/scratch/smallapragada/run2/0005123/s1_PDL011_06_27_2024_cells_stats.csv",
  s1_PDL012 = "/scratch/smallapragada/run2/0005123/s1_PDL012_cells_stats_06_27_2024.csv", 
  s1_PDL004 = "/scratch/smallapragada/run2/0005123/s1_PDL004_06_27_2024_cells_stats.csv", 
  s1_PDL010 = "/scratch/smallapragada/run2/0005123/s1_PDL010_06_27_2024_cells_stats.csv",
  s1_PDL013 = "/scratch/smallapragada/run2/0005123/s1_PDL013_06_27_2024_cells_stats.csv",
  s1_PDL007 = "/scratch/smallapragada/run2/0005123/s1_PDL007_06_27_2024_cells_stats.csv",
  s1_PDL014 = "/scratch/smallapragada/run2/0005123/s1_PDL014_06_27_2024_cells_stats.csv",
  s1_PDL001 = "/scratch/smallapragada/run2/0005123/s1_PDL001_06_27_2024_cells_stats.csv",
  s1_PDL002 = "/scratch/smallapragada/run2/0005123/s1_PDL002_06_27_2024_cells_stats.csv")

s2_id_list <- c(
  s2_PDL008 = "/scratch/smallapragada/run2/0005297/s2_PDL008_06_27_2024_cells_stats.csv",
  s2_PDL003 = "/scratch/smallapragada/run2/0005297/s2_PDL003_06_27_2024_cells_stats.csv",
  s2_PDL016 = "/scratch/smallapragada/run2/0005297/s2_PDL016_06_27_2024_cells_stats.csv",
  s2_PDL011 = "/scratch/smallapragada/run2/0005297/s2_PDL011_06_27_2024_cells_stats.csv",
  s2_PDL012 = "/scratch/smallapragada/run2/0005297/s2_PDL012_06_27_2024_cells_stats.csv",
  s2_PDL004 = "/scratch/smallapragada/run2/0005297/s2_PDL004_06_27_2024_cells_stats.csv",
  s2_PDL010 = "/scratch/smallapragada/run2/0005297/s2_PDL010_06_27_2024_cells_stats.csv",
  s2_PDL013 = "/scratch/smallapragada/run2/0005297/s2_PDL013_06_27_2024_cells_stats.csv",
  s2_PDL005 = "/scratch/smallapragada/run2/0005297/s2_PDL005_06_27_2024_cells_stats.csv",
  s2_PDL007 = "/scratch/smallapragada/run2/0005297/s2_PDL007_06_27_2024_cells_stats.csv",
  s2_PDL014 = "/scratch/smallapragada/run2/0005297/s2_PDL014_06_27_2024_cells_stats.csv",
  s2_PDL017 = "/scratch/smallapragada/run2/0005297/s2_PDL017_06_27_2024_cells_stats.csv",
  s2_PDL001 = "/scratch/smallapragada/run2/0005297/s2_PDL001_06_27_2024_cells_stats.csv",
  s2_PDL002 = "/scratch/smallapragada/run2/0005297/s2_PDL002_06_27_2024_cells_stats.csv")

## Assigning the loaded-in names to "s1_cells" and "s2_cells"
s1_ids <- names(s1_id_list)
names(s1_cells) <- s1_ids

s2_ids <- names(s2_id_list)
names(s2_cells) <- s2_ids

## Pasting sample names to the beginning of cell_id for both slides to ensure uniqueness
s1_cells$s1_PDL008 <- paste0("s1_GOLIATH_", s1_cells$s1_PDL008$Cell.ID)
s1_cells$s1_PDL003 <- paste0("s1_GOLIATH_", s1_cells$s1_PDL003$Cell.ID)
s1_cells$s1_PDL016 <- paste0("s1_GOLIATH_", s1_cells$s1_PDL016$Cell.ID)
s1_cells$s1_PDL011 <- paste0("s1_GOLIATH_", s1_cells$s1_PDL011$Cell.ID)
s1_cells$s1_PDL012 <- paste0("s1_GOLIATH_", s1_cells$s1_PDL012$Cell.ID)
s1_cells$s1_PDL004 <- paste0("s1_GOLIATH_", s1_cells$s1_PDL004$Cell.ID)
s1_cells$s1_PDL010 <- paste0("s1_GOLIATH_", s1_cells$s1_PDL010$Cell.ID)
s1_cells$s1_PDL013 <- paste0("s1_GOLIATH_", s1_cells$s1_PDL013$Cell.ID)
s1_cells$s1_PDL007 <- paste0("s1_GOLIATH_", s1_cells$s1_PDL007$Cell.ID)
s1_cells$s1_PDL014 <- paste0("s1_GOLIATH_", s1_cells$s1_PDL014$Cell.ID)
s1_cells$s1_PDL017 <- paste0("s1_PDL017_", s1_cells$s1_PDL017$Cell.ID)
s1_cells$s1_PDL001 <- paste0("s1_GOLIATH_", s1_cells$s1_PDL001$Cell.ID)
s1_cells$s1_PDL002 <- paste0("s1_GOLIATH_", s1_cells$s1_PDL002$Cell.ID)

s2_cells$s2_PDL008 <- paste0("s2_GOLIATH_", s2_cells$s2_PDL008$Cell.ID)
s2_cells$s2_PDL003 <- paste0("s2_GOLIATH_", s2_cells$s2_PDL003$Cell.ID)
s2_cells$s2_PDL016 <- paste0("s2_GOLIATH_", s2_cells$s2_PDL016$Cell.ID)
s2_cells$s2_PDL011 <- paste0("s2_GOLIATH_", s2_cells$s2_PDL011$Cell.ID)
s2_cells$s2_PDL012 <- paste0("s2_GOLIATH_", s2_cells$s2_PDL012$Cell.ID)
s2_cells$s2_PDL004 <- paste0("s2_GOLIATH_", s2_cells$s2_PDL004$Cell.ID)
s2_cells$s2_PDL010 <- paste0("s2_GOLIATH_", s2_cells$s2_PDL010$Cell.ID)
s2_cells$s2_PDL013 <- paste0("s2_GOLIATH_", s2_cells$s2_PDL013$Cell.ID)
s2_cells$s2_PDL005 <- paste0("s2_GOLIATH_", s2_cells$s2_PDL005$Cell.ID)
s2_cells$s2_PDL007 <- paste0("s2_GOLIATH_", s2_cells$s2_PDL007$Cell.ID)
s2_cells$s2_PDL014 <- paste0("s2_GOLIATH_", s2_cells$s2_PDL014$Cell.ID)
s2_cells$s2_PDL017 <- paste0("s2_GOLIATH_", s2_cells$s2_PDL017$Cell.ID)
s2_cells$s2_PDL001 <- paste0("s2_GOLIATH_", s2_cells$s2_PDL001$Cell.ID)
s2_cells$s2_PDL002 <- paste0("s2_GOLIATH_", s2_cells$s2_PDL002$Cell.ID)

## Making a new column named "Sample" that assigns cells to their individual samples 
merged_obj_unfiltered@meta.data <- merged_obj_unfiltered@meta.data %>%
  # Slide 1 - 0005123
  mutate(Sample = case_when(colnames(merged_obj_unfiltered) %in% s1_cells[["s1_PDL008"]] ~ "s1_PDL008",
                            colnames(merged_obj_unfiltered) %in% s1_cells[["s1_PDL003"]] ~ "s1_PDL003",
                            colnames(merged_obj_unfiltered) %in% s1_cells[["s1_PDL016"]] ~ "s1_PDL016",
                            colnames(merged_obj_unfiltered) %in% s1_cells[["s1_PDL011"]] ~ "s1_PDL011",
                            colnames(merged_obj_unfiltered) %in% s1_cells[["s1_PDL012"]] ~ "s1_PDL012",
                            colnames(merged_obj_unfiltered) %in% s1_cells[["s1_PDL004"]] ~ "s1_PDL004",
                            colnames(merged_obj_unfiltered) %in% s1_cells[["s1_PDL010"]] ~ "s1_PDL010",
                            colnames(merged_obj_unfiltered) %in% s1_cells[["s1_PDL013"]] ~ "s1_PDL013",
                            colnames(merged_obj_unfiltered) %in% s1_cells[["s1_PDL007"]] ~ "s1_PDL007",
                            colnames(merged_obj_unfiltered) %in% s1_cells[["s1_PDL014"]] ~ "s1_PDL014",
                            colnames(merged_obj_unfiltered) %in% s1_cells[["s1_PDL001"]] ~ "s1_PDL001",
                            colnames(merged_obj_unfiltered) %in% s1_cells[["s1_PDL002"]] ~ "s1_PDL002",
                            startsWith(colnames(merged_obj_unfiltered), "s1_PDL006") ~ "s1_PDL006",
                            startsWith(colnames(merged_obj_unfiltered), "s1_PDL015") ~ "s1_PDL015",
                            startsWith(colnames(merged_obj_unfiltered), "s1_PDL009") ~ "s1_PDL009",
                            startsWith(colnames(merged_obj_unfiltered), "s1_PDL017") ~ "s1_PDL017",
                            startsWith(colnames(merged_obj_unfiltered), "s1_PDL005") ~ "s1_PDL005",
                            
                            # Slide 2 - 0005297
                            colnames(merged_obj_unfiltered) %in% s2_cells[["s2_PDL008"]] ~ "s2_PDL008",
                            colnames(merged_obj_unfiltered) %in% s2_cells[["s2_PDL003"]] ~ "s2_PDL003",
                            colnames(merged_obj_unfiltered) %in% s2_cells[["s2_PDL016"]] ~ "s2_PDL016",
                            colnames(merged_obj_unfiltered) %in% s2_cells[["s2_PDL011"]] ~ "s2_PDL011",
                            colnames(merged_obj_unfiltered) %in% s2_cells[["s2_PDL012"]] ~ "s2_PDL012",
                            colnames(merged_obj_unfiltered) %in% s2_cells[["s2_PDL004"]] ~ "s2_PDL004",
                            colnames(merged_obj_unfiltered) %in% s2_cells[["s2_PDL010"]] ~ "s2_PDL010",
                            colnames(merged_obj_unfiltered) %in% s2_cells[["s2_PDL013"]] ~ "s2_PDL013",
                            colnames(merged_obj_unfiltered) %in% s2_cells[["s2_PDL005"]] ~ "s2_PDL005",
                            colnames(merged_obj_unfiltered) %in% s2_cells[["s2_PDL007"]] ~ "s2_PDL007",
                            colnames(merged_obj_unfiltered) %in% s2_cells[["s2_PDL014"]] ~ "s2_PDL014",
                            colnames(merged_obj_unfiltered) %in% s2_cells[["s2_PDL017"]] ~ "s2_PDL017",
                            colnames(merged_obj_unfiltered) %in% s2_cells[["s2_PDL001"]] ~ "s2_PDL001",
                            colnames(merged_obj_unfiltered) %in% s2_cells[["s2_PDL002"]] ~ "s2_PDL002",
                            startsWith(colnames(merged_obj_unfiltered), "s2_PDL009") ~ "s2_PDL009",
                            startsWith(colnames(merged_obj_unfiltered), "s2_PDL006") ~ "s2_PDL006",
                            startsWith(colnames(merged_obj_unfiltered), "s2_PDL015") ~ "s2_PDL015",
                            TRUE ~ "dropped"))

### Adding demographic information into new columns in Seurat metadata

## Slide number
merged_obj_unfiltered$Slide <- "NA"
merged_obj_unfiltered$Slide[grepl("^s1", colnames(merged_obj_unfiltered))] <- "0005123"
merged_obj_unfiltered$Slide[grepl("^s2", colnames(merged_obj_unfiltered))] <- "0005297"

## Sample timepoint - weeks gestation
merged_obj_unfiltered@meta.data$Timepoint <- NA
merged_obj_unfiltered@meta.data$Timepoint[merged_obj_unfiltered@meta.data$Sample == "s1_PDL001" | 
                                            merged_obj_unfiltered@meta.data$Sample == "s1_PDL002" | merged_obj_unfiltered@meta.data$Sample == "s1_PDL003" | 
                                            merged_obj_unfiltered@meta.data$Sample == "s2_PDL001" | merged_obj_unfiltered@meta.data$Sample == "s2_PDL002" | 
                                            merged_obj_unfiltered@meta.data$Sample == "s2_PDL003"] <- "16 - 21 weeks"
merged_obj_unfiltered@meta.data$Timepoint[merged_obj_unfiltered@meta.data$Sample == "s1_PDL004" |
                                            merged_obj_unfiltered@meta.data$Sample == "s1_PDL006" | merged_obj_unfiltered@meta.data$Sample == "s1_PDL005" | 
                                            merged_obj_unfiltered@meta.data$Sample == "s1_PDL007" | merged_obj_unfiltered@meta.data$Sample == "s1_PDL008" | 
                                            merged_obj_unfiltered@meta.data$Sample == "s2_PDL004" | merged_obj_unfiltered@meta.data$Sample == "s2_PDL006" | 
                                            merged_obj_unfiltered@meta.data$Sample == "s2_PDL005" | merged_obj_unfiltered@meta.data$Sample == "s2_PDL007" | 
                                            merged_obj_unfiltered@meta.data$Sample == "s2_PDL008"] <- "23 - 26 weeks"
merged_obj_unfiltered@meta.data$Timepoint[merged_obj_unfiltered@meta.data$Sample == "s1_PDL009" |
                                            merged_obj_unfiltered@meta.data$Sample == "s1_PDL014" | merged_obj_unfiltered@meta.data$Sample == "s2_PDL009" |
                                            merged_obj_unfiltered@meta.data$Sample == "s2_PDL014"] <- "28 - 29 weeks"
merged_obj_unfiltered@meta.data$Timepoint[merged_obj_unfiltered@meta.data$Sample == "s1_PDL010" | 
                                            merged_obj_unfiltered@meta.data$Sample == "s1_PDL015" | merged_obj_unfiltered@meta.data$Sample == "s2_PDL015" |
                                            merged_obj_unfiltered@meta.data$Sample == "s1_PDL011" | merged_obj_unfiltered@meta.data$Sample == "s1_PDL016" |
                                            merged_obj_unfiltered@meta.data$Sample == "s1_PDL013" | merged_obj_unfiltered@meta.data$Sample == "s2_PDL010" |
                                            merged_obj_unfiltered@meta.data$Sample == "s2_PDL011" | merged_obj_unfiltered@meta.data$Sample == "s2_PDL016" |
                                            merged_obj_unfiltered@meta.data$Sample == "s2_PDL013" | merged_obj_unfiltered@meta.data$Sample == "s1_PDL012" | 
                                            merged_obj_unfiltered@meta.data$Sample == "s2_PDL012"] <- "33 weeks - term"
merged_obj_unfiltered@meta.data$Timepoint[merged_obj_unfiltered@meta.data$Sample == "s1_PDL017" |
                                            merged_obj_unfiltered@meta.data$Sample == "s2_PDL017"] <- "30 years old"

## Disease phenotype
merged_obj_unfiltered@meta.data$Phenotype <- NA
merged_obj_unfiltered@meta.data$Phenotype[merged_obj_unfiltered@meta.data$Sample == "s1_PDL001" | 
                                            merged_obj_unfiltered@meta.data$Sample == "s1_PDL002" | 
                                            merged_obj_unfiltered@meta.data$Sample == "s2_PDL001" | 
                                            merged_obj_unfiltered@meta.data$Sample == "s2_PDL002"] <- "Fetal"
merged_obj_unfiltered@meta.data$Phenotype[merged_obj_unfiltered@meta.data$Sample == "s1_PDL003" | 
                                            merged_obj_unfiltered@meta.data$Sample == "s2_PDL003" | 
                                            merged_obj_unfiltered@meta.data$Sample == "s1_PDL006" | 
                                            merged_obj_unfiltered@meta.data$Sample == "s2_PDL006" |
                                            merged_obj_unfiltered@meta.data$Sample == "s1_PDL005" | 
                                            merged_obj_unfiltered@meta.data$Sample == "s2_PDL005" ] <- "Preterm - uninjured"

merged_obj_unfiltered@meta.data$Phenotype[merged_obj_unfiltered@meta.data$Sample == "s1_PDL004" | 
                                            merged_obj_unfiltered@meta.data$Sample == "s2_PDL004" |
                                            merged_obj_unfiltered@meta.data$Sample == "s1_PDL007" | 
                                            merged_obj_unfiltered@meta.data$Sample == "s2_PDL007" |
                                            merged_obj_unfiltered@meta.data$Sample == "s1_PDL008" | 
                                            merged_obj_unfiltered@meta.data$Sample == "s2_PDL008" |
                                            merged_obj_unfiltered@meta.data$Sample == "s1_PDL009" | 
                                            merged_obj_unfiltered@meta.data$Sample == "s2_PDL009"] <- "Preterm - injured"

merged_obj_unfiltered@meta.data$Phenotype[merged_obj_unfiltered@meta.data$Sample == "s1_PDL010" | 
                                            merged_obj_unfiltered@meta.data$Sample == "s2_PDL010" | 
                                            merged_obj_unfiltered@meta.data$Sample == "s1_PDL011" | 
                                            merged_obj_unfiltered@meta.data$Sample == "s2_PDL011" |
                                            merged_obj_unfiltered@meta.data$Sample == "s1_PDL013" | 
                                            merged_obj_unfiltered@meta.data$Sample == "s2_PDL013" | 
                                            merged_obj_unfiltered@meta.data$Sample == "s1_PDL012" | 
                                            merged_obj_unfiltered@meta.data$Sample == "s2_PDL012"] <- "Term - uninjured"
merged_obj_unfiltered@meta.data$Phenotype[merged_obj_unfiltered@meta.data$Sample == "s1_PDL016" | 
                                            merged_obj_unfiltered@meta.data$Sample == "s2_PDL016" |
                                            merged_obj_unfiltered@meta.data$Sample == "s1_PDL017" | 
                                            merged_obj_unfiltered@meta.data$Sample == "s2_PDL017" | 
                                            merged_obj_unfiltered@meta.data$Sample == "s1_PDL014" | 
                                            merged_obj_unfiltered@meta.data$Sample == "s2_PDL014" | 
                                            merged_obj_unfiltered@meta.data$Sample == "s1_PDL015" | 
                                            merged_obj_unfiltered@meta.data$Sample == "s2_PDL015"] <- "Rare disease"

## Saving the unfiltered object
saveRDS(merged_obj_unfiltered, "0005123_0005297_unfiltered_split_2024_08_15.rds")

#### FILTERING AND QC ----

## nCount_RNA vs nFeature_RNA

smoothScatter(merged_obj_unfiltered@meta.data$nCount_RNA,
              merged_obj_unfiltered@meta.data$nFeature_RNA,
              cex = 0.5, pch = 16)

VlnPlot(merged_obj_unfiltered@meta.data, features = c("nFeature_RNA", "nCount_RNA")) +
  NoLegend()

## nCount_RNA vs percent_negcodes
smoothScatter(merged_obj_unfiltered@meta.data$nCount_RNA,
              merged_obj_unfiltered@meta.data$percent_negcodes,
              cex = 0.5, pch = 16)

## nCount_RNA vs percent_negprobes
smoothScatter(merged_obj_unfiltered@meta.data$nCount_RNA,
              merged_obj_unfiltered@meta.data$percent_negprobes,
              cex = 0.5, pch = 16)

## nFeature_RNA vs. percent_negcodes
smoothScatter(merged_obj_unfiltered@meta.data$nFeature_RNA,
              merged_obj_unfiltered@meta.data$percent_negcodes,
              cex = 0.5, pch = 16)

## nFeature_RNA vs. percent_negprobes
smoothScatter(merged_obj_unfiltered@meta.data$nFeature_RNA,
              merged_obj_unfiltered@meta.data$percent_negprobes,
              cex = 0.5, pch = 16)

## nucleus_area vs. nCount_RNA 
smoothScatter(merged_obj_unfiltered@meta.data$nucleus_area,
              merged_obj_unfiltered@meta.data$nCount_RNA,
              cex = 0.5, pch = 16)

## nucleus_area vs. nFeature_RNA 
smoothScatter(merged_obj_unfiltered@meta.data$nucleus_area,
              merged_obj_unfiltered@meta.data$nFeature_RNA,
              cex = 0.5, pch = 16)

### Filtering
merged_obj <- subset(merged_obj_unfiltered, subset = nCount_RNA >= 10
                     & nFeature_RNA >= 5
                     & nucleus_area <= 75
                     & percent_negcodes < 0.15
                     & percent_negprobes == 0)

## Number of cells per sample before vs after filtering
dim(merged_obj_unfiltered)
dim(merged_obj)

## Saving filtered object
saveRDS(merged_obj, "0005123_0005297_filtered_split_2024_08_15.rds")

#### ADJUSTING SAMPLE COORDINATES FOR VISUALIZATION ----

## Subsetting out dropped cells - cells that weren't grouped into any sample
merged_obj_filtered <- subset(merged_obj, subset = Sample != "dropped")

## Adding a new column in the metadata that adjusts the x_centroid and y_centroid 

# X coords
merged_obj_filtered@meta.data <- merged_obj_filtered@meta.data %>%
  mutate(adj_x_centroid = case_when(merged_obj_filtered$Sample == "s1_PDL002" ~ (-merged_obj_filtered$x_centroid + 400),
                                    merged_obj_filtered$Sample == "s1_PDL009" ~ (-merged_obj_filtered$x_centroid - 7500),
                                    merged_obj_filtered$Sample == "s1_PDL006" ~ (-merged_obj_filtered$x_centroid + 500),
                                    merged_obj_filtered$Sample == "s1_PDL015" ~ (-merged_obj_filtered$x_centroid - 3000),
                                    merged_obj_filtered$Sample == "s1_PDL012" ~ (-merged_obj_filtered$x_centroid - 100),
                                    merged_obj_filtered$Sample == "s1_PDL010" ~ (-merged_obj_filtered$x_centroid),
                                    merged_obj_filtered$Sample == "s1_PDL001" ~ (-merged_obj_filtered$x_centroid),
                                    merged_obj_filtered$Sample == "s1_PDL017" ~ (-merged_obj_filtered$x_centroid - 3000),
                                    merged_obj_filtered$Sample == "s1_PDL014" ~ (-merged_obj_filtered$x_centroid + 500),
                                    merged_obj_filtered$Sample == "s1_PDL007" ~ (-merged_obj_filtered$x_centroid - 300),
                                    merged_obj_filtered$Sample == "s1_PDL005" ~ (-merged_obj_filtered$x_centroid),
                                    merged_obj_filtered$Sample == "s1_PDL013" ~ (-merged_obj_filtered$x_centroid),
                                    merged_obj_filtered$Sample == "s1_PDL004" ~ (-merged_obj_filtered$x_centroid),
                                    merged_obj_filtered$Sample == "s1_PDL011" ~ (-merged_obj_filtered$x_centroid),
                                    merged_obj_filtered$Sample == "s1_PDL016" ~ (-merged_obj_filtered$x_centroid),
                                    merged_obj_filtered$Sample == "s1_PDL003" ~ (-merged_obj_filtered$x_centroid),
                                    merged_obj_filtered$Sample == "s1_PDL008" ~ (merged_obj_filtered$x_centroid - 3000),
                                    merged_obj_filtered$Sample == "s2_PDL002" ~ (merged_obj_filtered$x_centroid + 750 + 5000),
                                    merged_obj_filtered$Sample == "s2_PDL009" ~ (merged_obj_filtered$x_centroid + 400 + 5000),
                                    merged_obj_filtered$Sample == "s2_PDL006" ~ (merged_obj_filtered$x_centroid + 11000),
                                    merged_obj_filtered$Sample == "s2_PDL015" ~ (merged_obj_filtered$x_centroid + 4000 + 5000),
                                    merged_obj_filtered$Sample == "s2_PDL012" ~ (merged_obj_filtered$x_centroid + 6000),
                                    merged_obj_filtered$Sample == "s2_PDL010" ~ (merged_obj_filtered$x_centroid + 200 + 5000),
                                    merged_obj_filtered$Sample == "s2_PDL001" ~ (merged_obj_filtered$x_centroid + 5200),
                                    merged_obj_filtered$Sample == "s2_PDL017" ~ (merged_obj_filtered$x_centroid + 5000),
                                    merged_obj_filtered$Sample == "s2_PDL014" ~ (merged_obj_filtered$x_centroid + 5000),
                                    merged_obj_filtered$Sample == "s2_PDL007" ~ (merged_obj_filtered$x_centroid + 5000),
                                    merged_obj_filtered$Sample == "s2_PDL005" ~ (merged_obj_filtered$x_centroid + 5000),
                                    merged_obj_filtered$Sample == "s2_PDL013" ~ (merged_obj_filtered$x_centroid + 5300),
                                    merged_obj_filtered$Sample == "s2_PDL004" ~ (merged_obj_filtered$x_centroid + 5000),
                                    merged_obj_filtered$Sample == "s2_PDL011" ~ (merged_obj_filtered$x_centroid + 5000),
                                    merged_obj_filtered$Sample == "s2_PDL016" ~ (merged_obj_filtered$x_centroid + 5000),
                                    merged_obj_filtered$Sample == "s2_PDL003" ~ (merged_obj_filtered$x_centroid + 5000),
                                    merged_obj_filtered$Sample == "s2_PDL008" ~ (merged_obj_filtered$x_centroid + 5000),
                                    TRUE ~ merged_obj_filtered$x_centroid))

# Y coords
merged_obj_filtered@meta.data <- merged_obj_filtered@meta.data %>%
  mutate(adj_y_centroid = case_when(merged_obj_filtered$Sample == "s1_PDL002" ~ (-merged_obj_filtered$y_centroid + 14000),
                                    merged_obj_filtered$Sample == "s1_PDL006" ~ (-merged_obj_filtered$y_centroid - 500),
                                    merged_obj_filtered$Sample == "s1_PDL015" ~ (-merged_obj_filtered$y_centroid - 1000),
                                    merged_obj_filtered$Sample == "s1_PDL009" ~ (-merged_obj_filtered$y_centroid - 750),
                                    merged_obj_filtered$Sample == "s1_PDL010" ~ (-merged_obj_filtered$y_centroid + 14000),
                                    merged_obj_filtered$Sample == "s1_PDL004" ~ (-merged_obj_filtered$y_centroid + 14000),
                                    merged_obj_filtered$Sample == "s1_PDL012" ~ (-merged_obj_filtered$y_centroid + 14000),
                                    merged_obj_filtered$Sample == "s1_PDL016" ~ (-merged_obj_filtered$y_centroid + 14000),
                                    merged_obj_filtered$Sample == "s1_PDL008" ~ (merged_obj_filtered$y_centroid - 1000),
                                    merged_obj_filtered$Sample == "s1_PDL003" ~ (-merged_obj_filtered$y_centroid + 14000),
                                    merged_obj_filtered$Sample == "s1_PDL011" ~ (-merged_obj_filtered$y_centroid + 14000),
                                    merged_obj_filtered$Sample == "s1_PDL007" ~ (merged_obj_filtered$y_centroid + 10500),
                                    merged_obj_filtered$Sample == "s1_PDL004" ~ (-merged_obj_filtered$y_centroid + 14000),
                                    merged_obj_filtered$Sample == "s1_PDL013" ~ (-merged_obj_filtered$y_centroid + 14000),
                                    merged_obj_filtered$Sample == "s1_PDL005" ~ (-merged_obj_filtered$y_centroid + 2000),
                                    merged_obj_filtered$Sample == "s1_PDL014" ~ (-merged_obj_filtered$y_centroid + 13500),
                                    merged_obj_filtered$Sample == "s1_PDL017" ~ (-merged_obj_filtered$y_centroid + 750),
                                    merged_obj_filtered$Sample == "s1_PDL001" ~ (-merged_obj_filtered$y_centroid + 14000),
                                    merged_obj_filtered$Sample == "s2_PDL002" ~ (merged_obj_filtered$y_centroid - 800),
                                    merged_obj_filtered$Sample == "s2_PDL006" ~ (merged_obj_filtered$y_centroid - 6500),
                                    merged_obj_filtered$Sample == "s2_PDL015" ~ (merged_obj_filtered$y_centroid - 4000),
                                    merged_obj_filtered$Sample == "s2_PDL009" ~ (merged_obj_filtered$y_centroid - 3300),
                                    merged_obj_filtered$Sample == "s2_PDL010" ~ (merged_obj_filtered$y_centroid - 200),
                                    merged_obj_filtered$Sample == "s2_PDL004" ~ (merged_obj_filtered$y_centroid - 350),
                                    merged_obj_filtered$Sample == "s2_PDL012" ~ (merged_obj_filtered$y_centroid - 750),
                                    TRUE ~ merged_obj_filtered$y_centroid))

## Add in spatial information as dimension reduction objects
position_xy <- cbind(merged_obj_filtered$adj_x_centroid, merged_obj_filtered$adj_y_centroid)
row.names(position_xy) <- row.names(merged_obj_filtered@meta.data)
colnames(position_xy) <- c("SP_1", "SP_2")
merged_obj_filtered[["sp"]] <- CreateDimReducObject(embeddings = position_xy, key = "SP_",
                                                    assay = DefaultAssay(merged_obj_filtered))

## Spatial visualization of all samples from both slides
DimPlot(merged_obj_filtered,
        reduction = "sp",
        group.by = "Sample",
        raster = T,
        label = T,
        cols = randomcoloR::distinctColorPalette(34)) +
  NoLegend()

## Saving the filtered and split object
saveRDS(merged_obj_filtered, "run2_filtered_split_adj_2024_08_15.rds")

################## JUMP TO `seurat_to_anndata_final.ipynb` FOR CLUSTERING ##################
                   
########################### RETURN HERE FOR CELL TYPE ANNOTATION ###########################

#### CELL TYPE ANNOTATION - SUBLINEAGE OBJECT SEPARATION ----

## Importing object post-RAPIDS processing
run2 <- readRDS("/scratch/smallapragada/rapids_pipeline/run2_filtered_split_adj_2024_08_15.rds")

## Setting ident for object
Idents(run2) <- 'leiden_brl_0.5'

## Assigning hallmark genes to lineages 
other_features <- c("EPCAM", "PECAM1", "PTPRC", "COL1A2", "MKI67", "TOP2A", "CENPF", "CDK1", "SOX9", "SOX2")

epithelial_features <- c("SFTPC", "KRT18", "COL4A3", "NAPSA", "MMP7", "DMBT1", "AGER", "RTKN2", "DUOX1", "NKX2-1", "KRT8", 
                         "PGC", "EGFR", "MUC5B", "MUC5AC", "TP63", "KRT5", "S100A2", "KRT15", 
                         "KRT14", "KRT17", "BPIFA1", "SAA2", "CXCL9", "CCNA1", "GCLM", "C20orf85", "FOXJ1", "TP73", 
                         "SCGB1A1", "CCL18", "SCGB3A2", "ERN2", "CFTR", "LTF", "AKR1C1", "AKR1C2", "CEACAM5", "CEACAM6",
                         "FCGBP", "CALCA", "CHGB", "SCG2")

endothelial_features <- c("EPAS1", "KDR", "GNG11", "CLDN5", "RAMP2", "CD34", "BMPR2", 
                          "FCN3", "APLN", "APLNR", "SPP1", "CA4", "CCL21", "SLC7A11", 
                          "IL7R", "HEY1", "PLVAP", "ACKR1", "POSTN", "PTGDS", "FGF2", 
                          "GCLM", "COL15A1", "ZEB1", "RNASE1", "HIF1A")

immune_features <- c("PPARG", "CXCR4", "HLA-DRA", "HLA-DQA1", "HAVCR2", "C1QC", "CD86", "KIT", "CDKN2A", 
                     "CD27", "BCL2L11", "PDIA6", "ISG20", "PDIA4", "HERPUD1", "GPR183", "GZMB", "LILRA4", 
                     "IRF7", "HYOU1", "CD79B", "PIM2", "CCR7", "BANK1", "MS4A1", "LTB", "TNFRSF13C", "CD19", 
                     "CD79A", "TCL1A", "BCL2L1", "HIST1H1C", "IL7R", "CD69", "TRAC", "CD3E", "CD3D", "CD2", 
                     "KLRB1", "CD8A", "CD3G", "CD68", "PLIN2", "MARCO", "FCER1G", "MS4A7",  "CD52", 
                     "FABP4", "MCEMP1", "S100A12", "S100A8", "S100A9", "LCK", "FGFBP2", "GNLY", "KLRG1", 
                     "BCL2", "IL2RA", "ITGAM", "TNFRSF9", "CPA3", "FASLG", "SLC25A37", "ITM2C", 
                     "TPSAB1", "SNCA", "CTLA4", "FOXP3")

mesenchymal_features <- c("SPARCL1", "ELN", "CTHRC1", "HAS1", "TGFB3", "MFAP5", "PI16", 
                          "SFRP2", "FAP", "MEG3", "AXL", "LGR6", "HAS2", "LUM","COL1A1",
                          "COL3A1", "DCN", "FN1", "SNAI2", "WNT2", "ITM2C", 
                          "WNT5A", "CCL2", "FAS", "SLC25A4", "CSPG4", "ACTA2", "LYZ", 
                          "PDGFRB", "PDGFRA", "CD4", "SOD2")

## Epithelial lineage dotplot
DotPlot(run2, group.by = "leiden_brl_0.5", features = c(other_features, epithelial_features)) + angle

## Endothelial lineage dotplot
DotPlot(run2, group.by = "leiden_brl_0.5", features = c(other_features, endothelial_features)) + angle

## Immune lineage dotplot
DotPlot(run2, group.by = "leiden_brl_0.5", features = c(other_features, immune_features)) + angle

## Mesenchymal lineage dotplot
DotPlot(run2, group.by = "leiden_brl_0.5", features = c(other_features, mesenchymal_features)) + angle

## UMAP of all cells
DimPlot(run2,
        reduction = "umap",
        group.by = "leiden_brl_0.5", 
        raster = T,
        label = F, 
        cols = randomcoloR::distinctColorPalette(50)) 

## Top markers in all clusters
run2_res1 <- FindAllMarkers(run2)
run2_res1 %>%
  mutate(pct.diff = abs(pct.1 - pct.2))

## Looking at the nCount_RNA and the nFeature_RNA distribution across clusters to determine low-quality clusters
ggplot(run2@meta.data, aes(factor(run2@meta.data$leiden_brl_0.5), (run2@meta.data$nCount_RNA))) + geom_violin()
ggplot(run2@meta.data, aes(factor(run2@meta.data$leiden_brl_0.5), (run2@meta.data$nFeature_RNA))) + geom_violin()

## Separating sublineages out via clusters based on top genes and hallmark genes
airway <- subset(run2, subset = leiden_brl_0.5 %in% c(0, 5, 13, 1))       
alv <- subset(run2, subset = leiden_brl_0.5 %in% c(2, 4, 8, 14))
endo <- subset(run2, subset = leiden_brl_0.5 %in% c(7, 15))            
mes <- subset(run2, subset = leiden_brl_0.5 %in% c(6, 18, 19))  
lymph <- subset(run2, subset = leiden_brl_0.5 %in% c(12, 14, 16))      
mye <- subset(run2, subset = leiden_brl_0.5 %in% c(9, 10, 11, 17))      

## Saving objects for RAPIDS processing
saveRDS(airway, "/scratch/smallapragada/run2/run2_airway_2024_08_15.rds")
saveRDS(alv, "/scratch/smallapragada/run2/run2_alv_2024_08_15.rds")
saveRDS(endo, "/scratch/smallapragada/run2/run2_endo_2024_08_15.rds")
saveRDS(mes, "/scratch/smallapragada/run2/run2_mes_2024_08_15.rds")
saveRDS(lymph, "/scratch/smallapragada/run2/run2_lymph_2024_08_15.rds")
saveRDS(mye, "/scratch/smallapragada/run2/run2_mye_2024_08_15.rds")

#### CELL TYPE ANNOTATION - AIRWAY ---- 

## Reading in airway object
airway <- readRDS("/scratch/smallapragada/rapids_pipeline/run2_airway_2024_08_15.rds")

## Looking at the nCount_RNA and the nFeature_RNA distribution across clusters to determine low-quality clusters
ggplot(airway@meta.data, aes(factor(airway@meta.data$airway_brl_0.4), (airway@meta.data$nCount_RNA))) + geom_violin()
ggplot(airway@meta.data, aes(factor(airway@meta.data$airway_brl_0.4), (airway@meta.data$nFeature_RNA))) + geom_violin()

## UMAP visualization of clusters in airway object
DimPlot(airway,
        reduction = "umap",
        group.by = "airway_brl_0.4", 
        raster = T,
        label = T, 
        cols = randomcoloR::distinctColorPalette(10)) 

## UMAP visualization of phenotype in airway object
DimPlot(airway,
        reduction = "umap",
        group.by = "Phenotype", 
        raster = T,
        label = T, 
        cols = randomcoloR::distinctColorPalette(10)) 

## A summary of how many cells per cluster belong to a specific disease phenotype
table(airway$airway_brl_0.4, airway$Phenotype)

## Top markers for airway
Idents(airway) = "airway_brl_0.4"
airway_markers <- FindAllMarkers(airway)
airway_markers %>%
  mutate(pct.diff = abs(pct.1 - pct.2)) %>%
  filter(avg_log2FC >= 0) %>%
  View()

## Airway dotplot
DotPlot(airway, group.by = "airway_brl_0.4", features = c(other_features, epithelial_features)) + angle

## Cluster 0 - AKR1C1+/1C2+
DotPlot(airway, features = c("SLC7A11", "GCLM", "HMOX1", "HMGA1", "PLIN2", "AKR1C1", "MGST1",
                             "AKR1C2", "GDF15", "ATF4", "MYC", "GSR")) + angle 

FeaturePlot(airway, features = c("SLC7A11", "GCLM", "HMOX1", "HMGA1", "PLIN2", "AKR1C1", "MGST1",
                                 "AKR1C2", "GDF15", "ATF3", "MYC", "GSR"))

## Cluster 1 - Basal/Proliferating Basal (Subcluster)
DotPlot(airway, features = c("KRT5", "KRT15", "KRT6A", "S100A2", "TP63", "KRT17", "MEG3", 
                             "KRT14", "AKR1C1", "CDH26", "KRT18", "MMP10", "MKI67", "TOP2A")) + angle 

FeaturePlot(airway, features = c("KRT5", "KRT15", "KRT6A", "S100A2", "TP63", "KRT17", "MEG3", 
                                 "KRT14", "AKR1C1", "CDH26", "KRT18", "MMP10", "MKI67", "TOP2A"))

## Cluster 2 - Ciliated/secretory (Recluster)
obj <- subset(airway, subset = airway_brl_0.4 == 2)

FeaturePlot(obj, features = c("AGR3", "C20orf85", "FOXJ1", "NUCB2", "MGST1", "CCNA1", "KRT8", "SCGB1A1", 
                              "DUOX1", "PDIA4", "GCLM", "UGDH", "SCGB3A2", "BCL2L1", "LMAN1", "SOX2"))

DotPlot(airway, features = c(epithelial_features, "AGR3", "NUCB2", "MGST1", "PDIA4",
                             "UGDH", "BCL2L1", "LMAN1", "SOX2")) + angle

## Cluster 3 - AT2 (Recluster)
obj <- subset(airway, subset = airway_brl_0.4 == 3)

FeaturePlot(obj, features = c("PGC", "SFTPC", "MSLN", "NAPSA", "RTKN2", "SFTPD", "LAMP3", "AGER", 
                              "RNASE1", "PTPRC", "KRT14", "KRT17", "S100A8"))

DotPlot(airway, features = c("PGC", "SFTPC", "MSLN", "NAPSA", "RTKN2", "SFTPD", "LAMP3", "AGER", 
                             "RNASE1", "PTPRC", "KRT14", "KRT17", "S100A8")) + angle

## Cluster 4 - FB/SMC (Subcluster)
obj <- subset(airway_subclustering, subset = airway_brl_0.4 == 4)

FeaturePlot(obj, features = c("COL1A2", "COL3A1", "MEG3", "FN1", "COL1A1", "LUM", "PDGFRA", "SOD2", 
                              "VIM", "POSTN", "SPARCL1", "PDGFRB", "ACTA2", "SNAI2", "HAS2", "PDGFRA"))

DotPlot(airway, features = c(mesenchymal_features)) + angle

## Cluster 5 - Capillaries (Recluster)
obj <- subset(airway, subset = airway_brl_0.4 == 5)

FeaturePlot(obj, features = c("PECAM1", "EPAS1", "CD34", "KDR", "CA4", "MEG3", "GNG11", "ACKR1", 
                              "VIM", "FCN3", "CCL21", "IL7R", "BMPR2", "CLDN5", "RAMP2", "KIT", "APLN", "APLNR"))

DotPlot(airway, features = c(endothelial_features)) + angle

## Cluster 6 - Secretory (Recluster)
obj <- subset(airway, subset = airway_brl_0.4 == 6)

FeaturePlot(obj, features = c("SAA2", "WFDC2", "HLA-DRA", "MUC5B", "SOD2", "STAT1", "XBP1", "NUCB2", 
                              "SCGB3A2", "HYOU1", "MGST1", "PDIA4", "AGR3", "C20orf85", "KRT8", "IFIT3"))

DotPlot(airway, group.by = "airway_brl_0.4", features = c("SAA2", "WFDC2", "HLA-DRA", "MUC5B", "SOD2", "STAT1", "XBP1", "NUCB2", 
                                                          "SCGB3A2", "HYOU1", "MGST1", "PDIA4", "AGR3", "C20orf85", "KRT8", "IFIT3")) + angle

## Cluster 7 - PNEC (Subcluster)
obj <- subset(airway, subset = airway_brl_0.4 == 7)

FeaturePlot(obj, features = c("CALCA", "CHGB", "SEC11C", "SCG2", "MEG3", "ASCL1", "ATP2A3", "ITM2C", 
                              "SNCA", "DIRAS3", "BCL2", "CGA", "SCGB3A2", "HYOU1", "XBP1", "SCGB1A1"))

DotPlot(airway, group.by = "airway_brl_0.4", features = c(epithelial_features, "CALCA", "CHGB", "SEC11C", "SCG2", "MEG3", 
                                                          "ASCL1", "ATP2A3", "ITM2C", "SNCA", "DIRAS3", "BCL2", "CGA")) + angle

### Subclustering of clusters 1 and 4

airway_subclustering <- readRDS("/scratch/smallapragada/rapids_pipeline/run2_airway_subclustering_2024_08_15.rds")

## Cluster 1

## UMAP visualization of clusters in subclustered airway object
DimPlot(airway_subclustering,
        reduction = "umap",
        group.by = "airway_sc_0.4_c1_0.2", 
        raster = T,
        label = T, 
        cols = randomcoloR::distinctColorPalette(20)) 

## UMAP visualization of phenotypes in subclustered airway object
DimPlot(airway_subclustering,
        reduction = "umap",
        group.by = "Phenotype", 
        raster = T,
        label = F, 
        cols = randomcoloR::distinctColorPalette(5)) 

## A summary of how many cells per cluster belong to a specific disease phenotype
table(airway_subclustering$airway_sc_0.4_c1_0.2, airway_subclustering$Phenotype)

## Looking at the nCount_RNA and the nFeature_RNA distribution across clusters to determine low-quality clusters
ggplot(airway_subclustering@meta.data, aes(factor(airway_subclustering@meta.data$airway_sc_0.4_c1_0.1), 
                                           (airway_subclustering@meta.data$nCount_RNA))) + geom_violin()
ggplot(airway_subclustering@meta.data, aes(factor(airway_subclustering@meta.data$airway_sc_0.4_c1_0.1), 
                                           (airway_subclustering@meta.data$nFeature_RNA))) + geom_violin()

## Temp object to visualize clusters of interest
obj <- subset(airway_subclustering, subset = airway_sc_0.4_c1_0.2 == "1,0" | airway_sc_0.4_c1_0.2 == "1,1" |
                airway_sc_0.4_c1_0.2 == "1,2" | airway_sc_0.4_c1_0.2 == "1,3")

DotPlot(obj, group.by = "airway_sc_0.4_c1_0.2", features = c(other_features, epithelial_features)) + angle

## Top markers of temp object
Idents(obj) = "airway_sc_0.4_c1_0.2"
airway_markers <- FindAllMarkers(obj)
airway_markers %>%
  mutate(pct.diff = abs(pct.1 - pct.2)) %>%
  filter(avg_log2FC >= 0) %>%
  View()

FeaturePlot(obj, features = c("DMBT1", "LTF", "KRT5", "KRT15", "KRT17", "TP63", "KRT6A", "BPIFA1", 
                              "MUC5B", "AGR3", "C20orf85", "SCGB1A1", "FOXJ1", "MUC5AC", "MKI67"))

DotPlot(run2_old, features = c("DMBT1", "LTF", "MUC5B", "SAA2", "LYZ", "WFDC2", "SOX9", "RNASE1", "ACTA2", "MMP7", "SCGB3A2")) + angle

# 1,0 = Proliferating basal
# 1,1 = AT2
# 1,2 = Multiciliated
# 1,3 = Basal

## Cluster 4

## UMAP visualization of clusters in subclustered airway object
DimPlot(airway_subclustering,
        reduction = "umap",
        group.by = "airway_sc_0.4_c4_0.2", 
        raster = T,
        label = T, 
        cols = randomcoloR::distinctColorPalette(10)) 

## UMAP visualization of phenotype in subclustered airway object
DimPlot(airway_subclustering,
        reduction = "umap",
        group.by = "Phenotype", 
        raster = T,
        label = F, 
        cols = randomcoloR::distinctColorPalette(5)) 

## A summary of how many cells per cluster belong to a specific disease phenotype
table(airway_subclustering$airway_sc_0.4_c4_0.2, airway_subclustering$Phenotype)

## Looking at the nCount_RNA and the nFeature_RNA distribution across clusters to determine low-quality clusters
ggplot(airway_subclustering@meta.data, aes(factor(airway_subclustering@meta.data$airway_sc_0.4_c4_0.2), 
                                           (airway_subclustering@meta.data$nCount_RNA))) + geom_violin()
ggplot(airway_subclustering@meta.data, aes(factor(airway_subclustering@meta.data$airway_sc_0.4_c4_0.2), 
                                           (airway_subclustering@meta.data$nFeature_RNA))) + geom_violin()

## Temp object to visualize clusters of interest
obj <- subset(airway_subclustering, subset = airway_sc_0.4_c4_0.2 == "4,0" | 
                airway_sc_0.4_c4_0.2 == "4,1" | airway_sc_0.4_c4_0.2 == "4,2")

DotPlot(obj, group.by = "airway_sc_0.4_c4_0.2", features = c(other_features, epithelial_features)) + angle

## Top markers for temp object
Idents(obj) = "airway_sc_0.4_c4_0.2"
airway_markers <- FindAllMarkers(obj)
airway_markers %>%
  mutate(pct.diff = abs(pct.1 - pct.2)) %>%
  filter(avg_log2FC >= 0) %>%
  View()

FeaturePlot(obj, features = c("SNAI2", "WNT2", "LUM", "HAS2", "SOD2", "CD44", "HIF1A", "CCN2", 
                              "ACTA2", "PDGFRA", "SCG2", "CALCA", "CHGB", "ASCL1"))

DotPlot(obj, features = c("SNAI2", "WNT2", "LUM", "HAS2", "SOD2", "CD44", "HIF1A", "CCN2", 
                          "ACTA2", "PDGFRA", "SCG2", "CALCA", "CHGB", "ASCL1")) + angle

# 4,0 = Fibroblasts
# 4,1 = SMC
# 4,2 = PNEC

### Reclustering of clusters 2, 3, 5, and 6

cluster2 <- subset(airway, subset = airway_brl_0.4 == 2)
cluster3 <- subset(airway, subset = airway_brl_0.4 == 3)
cluster5 <- subset(airway, subset = airway_brl_0.4 == 5)
cluster6 <- subset(airway, subset = airway_brl_0.4 == 6)

saveRDS(cluster2, "/scratch/smallapragada/run2/run2_airway_reclustering_c2_2024_08_15.rds")
saveRDS(cluster3, "/scratch/smallapragada/run2/run2_airway_reclustering_c3_2024_08_15.rds")
saveRDS(cluster5, "/scratch/smallapragada/run2/run2_airway_reclustering_c5_2024_08_15.rds")
saveRDS(cluster6, "/scratch/smallapragada/run2/run2_airway_reclustering_c6_2024_08_15.rds")

## Cluster 2 

cluster2 <- readRDS("/scratch/smallapragada/rapids_pipeline/run2_airway_reclustering_c2_2024_08_15.rds")

## UMAP visualization of clusters in reclustered airway object
DimPlot(cluster2,
        reduction = "umap",
        group.by = "airway_c2_brl_0.3", 
        raster = T,
        label = T, 
        cols = randomcoloR::distinctColorPalette(10)) 

## UMAP visualization of phenotype in reclustered airway object
DimPlot(cluster2,
        reduction = "umap",
        group.by = "Phenotype", 
        raster = T,
        label = F, 
        cols = randomcoloR::distinctColorPalette(5)) 

## A summary of how many cells per cluster belong to a specific disease phenotype
table(cluster2$airway_c2_brl_0.3, cluster2$Phenotype)

## Looking at the nCount_RNA and the nFeature_RNA distribution across clusters to determine low-quality clusters
ggplot(cluster2@meta.data, aes(factor(cluster2@meta.data$airway_c2_brl_0.3), 
                               (cluster2@meta.data$nCount_RNA))) + geom_violin()
ggplot(cluster2@meta.data, aes(factor(cluster2@meta.data$airway_c2_brl_0.3), 
                               (cluster2@meta.data$nFeature_RNA))) + geom_violin()

DotPlot(cluster2, group.by = "airway_c2_brl_0.3", features = c(other_features, epithelial_features)) + angle

## Top markers for reclustered cluster
Idents(cluster2) = "airway_c2_brl_0.3"
airway_markers <- FindAllMarkers(cluster2)
airway_markers %>%
  mutate(pct.diff = abs(pct.1 - pct.2)) %>%
  filter(avg_log2FC >= 0) %>%
  View()

FeaturePlot(cluster2, features = c("MGST1", "KRT8", "CTNNB1", "KRT18", "PDIA4", "NUCB2", "UGDH", "RHOA", 
                                   "C20orf85", "NUTF2", "SAA2", "MUC5B", "WFDC2", "LTF", "BPIFA1", "VEGFA", 
                                   "SCGB1A1", "HES1", "VIM", "CD44", "TOP2A", "SFTPC", "SCGB3A2", "AGER"))

DotPlot(cluster2, features = c("MGST1", "KRT8", "CTNNB1", "KRT18", "PDIA4", "NUCB2", "UGDH", "RHOA", 
                               "C20orf85", "NUTF2", "SAA2", "MUC5B", "WFDC2", "LTF", "BPIFA1", "VEGFA", 
                               "SCGB1A1", "HES1", "VIM", "CD44", "TOP2A", "SFTPC", "SCGB3A2", "AGER")) + angle

# Cluster 0 - Multiciliated
# Cluster 1 - Secretory MUC5B
# Cluster 2 - Secretory SCGB3A2+/1A1+
# Cluster 3 - Secretory SCGB3A2+/1A1+

## Cluster 3 

cluster3 <- readRDS("/scratch/smallapragada/rapids_pipeline/run2_airway_reclustering_c3_2024_08_15.rds")

## UMAP visualization of clusters in reclustered airway object
DimPlot(cluster3,
        reduction = "umap",
        group.by = "airway_c3_brl_0.4", 
        raster = T,
        label = T, 
        cols = randomcoloR::distinctColorPalette(10)) 

## UMAP visualization of phenotype in reclustered airway object
DimPlot(cluster3,
        reduction = "umap",
        group.by = "Phenotype", 
        raster = T,
        label = F, 
        cols = randomcoloR::distinctColorPalette(5)) 

## A summary of how many cells per cluster belong to a specific disease phenotype
table(cluster3$airway_c3_brl_0.4, cluster3$Phenotype)

## Looking at the nCount_RNA and the nFeature_RNA distribution across clusters to determine low-quality clusters
ggplot(cluster3@meta.data, aes(factor(cluster3@meta.data$airway_c3_brl_0.4), 
                               (cluster3@meta.data$nCount_RNA))) + geom_violin()
ggplot(cluster3@meta.data, aes(factor(cluster3@meta.data$airway_c3_brl_0.4), 
                               (cluster3@meta.data$nFeature_RNA))) + geom_violin()

DotPlot(cluster3, group.by = "airway_c3_brl_0.4", features = c(other_features, epithelial_features)) + angle

## Top markers
Idents(cluster3) = "airway_c3_brl_0.4"
airway_markers <- FindAllMarkers(cluster3)
airway_markers %>%
  mutate(pct.diff = abs(pct.1 - pct.2)) %>%
  filter(avg_log2FC >= 0) %>%
  View()

FeaturePlot(cluster3, features = c("MS4A1", "CPA3", "TPSAB1", "S100A2", "MEG3", "COL1A2", "COL3A1", "VEGFA", 
                                   "ICAM1", "MSLN", "AGER", "RTKN2", "SCGB3A2", "S100A8", "S100A9", "S100A12"))

FeaturePlot(cluster3, features = c("ITGAX", "XBP1", "MUC5B", "SLC25A37", "SAA2", "SFTPC", "C20orf85", "AGR3", "TCL1A", "KRT17",
                                   "SPCS3", "UBE2J1", "ACKR1", "SLC7A11", "CSPG4", "TOP2A", "FCER1G", "FCGR3A", "LAMP3",
                                   "SFTPC", "RNASE1", "NAPSA"))

DotPlot(cluster3, features = c("MS4A1", "CPA3", "TPSAB1", "S100A2", "MEG3", "COL1A2", "COL3A1", "VEGFA", 
                               "ICAM1", "MSLN", "AGER", "RTKN2", "SCGB3A2", "S100A8", "S100A9", "S100A12", 
                               "ITGAX", "XBP1", "MUC5B", "SLC25A37", "SAA2", "SFTPC", "C20orf85", "AGR3", "TCL1A", "KRT17",
                               "SPCS3", "UBE2J1", "ACKR1", "SLC7A11", "CSPG4", "TOP2A", "FCER1G", "FCGR3A", "LAMP3",
                               "RNASE1", "NAPSA", "SPP1")) + angle

# Cluster 0 - Secretory MUC5B+
# Cluster 1 - Mast
# Cluster 2 - AT1
# Cluster 3 - Low quality
# Cluster 4 - B
# Cluster 5 - cDC
# Cluster 6 - AT2

## Cluster 5 

cluster5 <- readRDS("/scratch/smallapragada/rapids_pipeline/run2_airway_reclustering_c5_2024_08_15.rds")

## UMAP visualization of clusters in reclustered airway object
DimPlot(cluster5,
        reduction = "umap",
        group.by = "airway_c5_brl_0.3", 
        raster = T,
        label = T, 
        cols = randomcoloR::distinctColorPalette(10)) 

## UMAP visualization of phenotype in reclustered airway object
DimPlot(cluster5,
        reduction = "umap",
        group.by = "Phenotype", 
        raster = T,
        label = F, 
        cols = randomcoloR::distinctColorPalette(5)) 

## A summary of how many cells per cluster belong to a specific disease phenotype
table(cluster5$airway_c5_brl_0.3, cluster5$Phenotype)

## Looking at the nCount_RNA and the nFeature_RNA distribution across clusters to determine low-quality clusters
ggplot(cluster5@meta.data, aes(factor(cluster5@meta.data$airway_c5_brl_0.3), 
                               (cluster5@meta.data$nCount_RNA))) + geom_violin()
ggplot(cluster5@meta.data, aes(factor(cluster5@meta.data$airway_c5_brl_0.3), 
                               (cluster5@meta.data$nFeature_RNA))) + geom_violin()

DotPlot(cluster5, group.by = "airway_c5_brl_0.3", features = c(other_features, epithelial_features)) + angle

## Top markers for reclustered object
Idents(cluster5) = "airway_c5_brl_0.3"
airway_markers <- FindAllMarkers(cluster5)
airway_markers %>%
  mutate(pct.diff = abs(pct.1 - pct.2)) %>%
  filter(avg_log2FC >= 0) %>%
  View()

FeaturePlot(cluster5, features = c("KRT8", "KRT18", "MGST1", "HSPA5", "SFTPC", "LAMP3", "CD44"))      # Cluster 0

FeaturePlot(cluster5, features = c("FCN3", "KRT15", "SLC7A11", "ACKR1", "PGC", "IL4R", "KRT14", "SNAI2",
                                   "BPIFA1", "MUC5AC", "RACGAP1"))      # Cluster 1

FeaturePlot(cluster5, features = c("C20orf85", "FOXJ1", "CCNA1", "SCGB3A2", "AGR3", "SCGB1A1", 
                                   "SOX2", "SFTA2", "CTHRC1", "CEACAM6"))     # Cluster 2

FeaturePlot(cluster5, features = c("CCL21", "MRC1", "TP63", "COL3A1", "GNG11", "PPARG", 
                                   "CTHRC1"))     # Cluster 3

FeaturePlot(cluster5, features = c("CD34", "ACKR1", "APLNR", "POSTN", "RAMP2"))       # Cluster 4


DotPlot(cluster5, features = c("KRT8", "KRT18", "MGST1", "HSPA5", "SFTPC", "LAMP3", "CD44",
                               "FCN3", "KRT15", "SLC7A11", "ACKR1", "PGC", "IL4R", "KRT14", "SNAI2",
                               "BPIFA1", "MUC5AC", "RACGAP1", "C20orf85", "FOXJ1", "CCNA1", "SCGB3A2", 
                               "AGR3", "SCGB1A1", "SOX2", "SFTA2", "CTHRC1", "CEACAM6", 
                               "CCL21", "MRC1", "TP63", "COL3A1", "GNG11", "PPARG", 
                               "CTHRC1", "CD34", "ACKR1", "APLNR", "POSTN", "RAMP2")) + angle

# Cluster 0 - Transitional AT2
# Cluster 1 - Venous
# Cluster 2 - Multiciliated
# Cluster 3 - Lymphatic endothelial
# Cluster 4 - Capillaries

## Cluster 6 

cluster6 <- readRDS("/scratch/smallapragada/rapids_pipeline/run2_airway_reclustering_c6_2024_08_15.rds")

## UMAP visualization of clusters in reclustered airway object
DimPlot(cluster6,
        reduction = "umap",
        group.by = "airway_c6_brl_0.4", 
        raster = T,
        label = T, 
        cols = randomcoloR::distinctColorPalette(10)) 

## UMAP visualization of phenotype in reclustered airway object
DimPlot(cluster6,
        reduction = "umap",
        group.by = "Phenotype", 
        raster = T,
        label = F, 
        cols = randomcoloR::distinctColorPalette(5)) 

## A summary of how many cells per cluster belong to a specific disease phenotype
table(cluster6$airway_c6_brl_0.4, cluster6$Phenotype)

## Looking at the nCount_RNA and the nFeature_RNA distribution across clusters to determine low-quality clusters
ggplot(cluster6@meta.data, aes(factor(cluster6@meta.data$airway_c6_brl_0.4), 
                               (cluster6@meta.data$nCount_RNA))) + geom_violin()
ggplot(cluster6@meta.data, aes(factor(cluster6@meta.data$airway_c6_brl_0.4), 
                               (cluster6@meta.data$nFeature_RNA))) + geom_violin()

DotPlot(cluster6, group.by = "airway_c6_brl_0.4", features = c(other_features, epithelial_features)) + angle

## Top markers for reclustered object
Idents(cluster6) = "airway_c6_brl_0.4"
airway_markers <- FindAllMarkers(cluster6)
airway_markers %>%
  mutate(pct.diff = abs(pct.1 - pct.2)) %>%
  filter(avg_log2FC >= 0) %>%
  View()


FeaturePlot(cluster6, features = c("MUC5B", "LTF", "BPIFA1", "SAA2", "SCGB3A2"))    # Cluster 0

FeaturePlot(cluster6, features = c("MUC5B", "KRT8", "KDR", "LTF", "CFTR", "CXCL9"))     # Cluster 1 & 2

FeaturePlot(cluster6, features = c("TOP2A", "MKI67", "CCNB2", "GDF15", "CFTR", "SCGB3A2", "SCGB1A1"))     # Cluster 3

FeaturePlot(cluster6, features = c("CEACAM6", "RNASE1", "SCGB1A1", "VEGFA", "GDF15", "IDH1"))     # Cluster 4

FeaturePlot(cluster6, features = c("ITGB1", "NKX2-1", "SPINK1", "CD44", "HMGA1", "AXL", 
                                   "SEC11C", "PDGFRB"))     # Cluster 5

FeaturePlot(cluster6, features = c("IFIT3", "CCN2", "C20orf85", "OAS2", "SAA2", "IRF1"))      # Cluster 6

FeaturePlot(cluster6, features = c("DMBT1", "NAPSA", "LAMP3", "SFTPC", "SFTPD", "RTKN2", "AGER"))     # Cluster 7

FeaturePlot(cluster6, features = c("CCL22", "IL7R", "CD68", "IL2RA", "TNFRSF9", "GPR183", "FCER1G"))  # Cluster 8

DotPlot(cluster6, features = c("MUC5B", "LTF", "BPIFA1", "SAA2", "SCGB3A2", "KRT8", 
                               "KDR", "LTF", "CFTR", "CXCL9", "TOP2A", "MKI67", "CCNB2", 
                               "GDF15", "CFTR", "SCGB1A1", "CEACAM6", "RNASE1", "VEGFA", 
                               "IDH1", "ITGB1", "NKX2-1", "SPINK1", "CD44", "HMGA1", "AXL", 
                               "SEC11C", "PDGFRB", "IFIT3", "CCN2", "C20orf85", "OAS2", "IRF1",
                               "DMBT1", "NAPSA", "LAMP3", "SFTPC", "SFTPD", "RTKN2", "AGER",
                               "CCL22", "IL7R", "CD68", "IL2RA", "TNFRSF9", "GPR183", "FCER1G")) + angle

# Cluster 0 - Secretory MUC5B+
# Cluster 1 - Fetal secretory MUC5B+
# Cluster 2 - Low quality
# Cluster 3 - Proliferating secretory SCGB3A2+
# Cluster 4 - Secretory SCGB3A2+
# Cluster 5 - Pericytes
# Cluster 6 - Multiciliated
# Cluster 7 - AT2
# Cluster 8 - pDC

## Adding cell type labels to each of the reclustered objects
cluster2$CT <- ""
cluster2$CT[cluster2$airway_c2_brl_0.3 == 0] <- "Multiciliated"
cluster2$CT[cluster2$airway_c2_brl_0.3 == 1] <- "Secretory MUC5B+"
cluster2$CT[cluster2$airway_c2_brl_0.3 == 2 | cluster2$airway_c2_brl_0.3 == 3] <- "Secretory 3A2+ & 1A1+"

cluster3$CT <- ""
cluster3$CT[cluster3$airway_c3_brl_0.4 == 0] <- "Secretory MUC5B+"
cluster3$CT[cluster3$airway_c3_brl_0.4 == 1] <- "Mast"
cluster3$CT[cluster3$airway_c3_brl_0.4 == 2] <- "AT2"
cluster3$CT[cluster3$airway_c3_brl_0.4 == 3] <- "Low quality"
cluster3$CT[cluster3$airway_c3_brl_0.4 == 4] <- "B"
cluster3$CT[cluster3$airway_c3_brl_0.4 == 5] <- "cDC"
cluster3$CT[cluster3$airway_c3_brl_0.4 == 6] <- "AT2"

cluster5$CT <- ""
cluster5$CT[cluster5$airway_c5_brl_0.3 == 0] <- "Transitional AT2"
cluster5$CT[cluster5$airway_c5_brl_0.3 == 1] <- "Venous"
cluster5$CT[cluster5$airway_c5_brl_0.3 == 2] <- "Multiciliated"
cluster5$CT[cluster5$airway_c5_brl_0.3 == 3] <- "Lymphatic endothelial"
cluster5$CT[cluster5$airway_c5_brl_0.3 == 4] <- "Capillaries"

cluster6$CT <- ""
cluster6$CT[cluster6$airway_c6_brl_0.4 == 0 | cluster6$airway_c6_brl_0.4 == 1] <- "Secretory MUC5B+"
cluster6$CT[cluster6$airway_c6_brl_0.4 == 2] <- "Low quality"
cluster6$CT[cluster6$airway_c6_brl_0.4 == 3] <- "Proliferating airway"
cluster6$CT[cluster6$airway_c6_brl_0.4 == 4] <- "Secretory 3A2+ & 1A1+"
cluster6$CT[cluster6$airway_c6_brl_0.4 == 5] <- "Pericytes"
cluster6$CT[cluster6$airway_c6_brl_0.4 == 6] <- "Multiciliated"
cluster6$CT[cluster6$airway_c6_brl_0.4 == 7] <- "AT2"
cluster6$CT[cluster6$airway_c6_brl_0.4 == 8] <- "pDC"

## Combining objects together and merging into one
obj_list <- list(cluster2, cluster3, cluster5, cluster6)
firstpass <- merge(x = obj_list[[1]], y = obj_list[2:length(obj_list)])

## Adding labels from the merged reclustered object to the whole airway object + all other CT labels 
airway_subclustering$CT <- ""
airway_subclustering$CT <- firstpass$CT
airway_subclustering$CT[airway_subclustering$airway_brl_0.4 == 0] <- "AKR1C1+ & AKR1C2+"
airway_subclustering$CT[airway_subclustering$airway_sc_0.4_c1_0.2 == "1,0"] <- "Proliferating basal"
airway_subclustering$CT[airway_subclustering$airway_sc_0.4_c1_0.2 == "1,1"] <- "AT2"
airway_subclustering$CT[airway_subclustering$airway_sc_0.4_c1_0.2 == "1,2"] <- "Multiciliated"
airway_subclustering$CT[airway_subclustering$airway_sc_0.4_c1_0.2 == "1,3"] <- "Basal"
airway_subclustering$CT[airway_subclustering$airway_sc_0.4_c4_0.2 == "4,0"] <- "Fibroblasts"
airway_subclustering$CT[airway_subclustering$airway_sc_0.4_c4_0.2 == "4,1"] <- "SMC"
airway_subclustering$CT[airway_subclustering$airway_sc_0.4_c4_0.2 == "4,2"] <- "PNEC"
airway_subclustering$CT[airway_subclustering$airway_brl_0.4 == 7] <- "PNEC"

## Save airway object
saveRDS(airway_subclustering, "/scratch/smallapragada/run2/final_airway_CT_2024_08_15.rds")

#### CELL TYPE ANNOTATION - ALVEOLAR ----

## Reading in alveolar object
alv <- readRDS("/scratch/smallapragada/rapids_pipeline/run2_alv_2024_08_15.rds")

## Looking at the nCount_RNA and the nFeature_RNA distribution across clusters to determine low-quality clusters
ggplot(alv@meta.data, aes(factor(alv@meta.data$alv_brl_0.4), (alv@meta.data$nCount_RNA))) + geom_violin()
ggplot(alv@meta.data, aes(factor(alv@meta.data$alv_brl_0.4), (alv@meta.data$nFeature_RNA))) + geom_violin()

## UMAP visualization of clusters in alveolar object
DimPlot(alv,
        reduction = "umap",
        group.by = "alv_brl_0.4", 
        raster = T,
        label = T, 
        cols = randomcoloR::distinctColorPalette(10)) 

## UMAP visualization of phenotype in airway object
DimPlot(alv,
        reduction = "umap",
        group.by = "Phenotype", 
        raster = T,
        label = T, 
        cols = randomcoloR::distinctColorPalette(10)) 

## A summary of how many cells per cluster belong to a specific disease phenotype
table(alv$alv_brl_0.4, alv$Phenotype)

## Top markers for alveolar
Idents(alv) = "alv_brl_0.4"
alv_markers <- FindAllMarkers(alv)
alv_markers %>%
  mutate(pct.diff = abs(pct.1 - pct.2)) %>%
  filter(avg_log2FC >= 0) %>%
  View()

## Alveolar dotplot
DotPlot(alv, group.by = "alv_brl_0.4", features = c(other_features, epithelial_features)) + angle

## Cluster 0 - AT2

# Temp object
obj <- subset(alv, subset = alv_brl_0.4 == 0)

FeaturePlot(obj, features = c("SFTPD", "NAPSA", "DUOX1", "SFTPC", "LAMP3", "PGC", "MGST1", "RNASE1", "NKX2-1",
                              "CD44", "ITGB6"))

DotPlot(alv, features = c(other_features, epithelial_features)) + angle

## Cluster 1 - Alveolar macrophages (Recluster)

# Temp object
obj <- subset(alv, subset = alv_brl_0.4 == 1)

FeaturePlot(obj, features = c("CCL18", "MS4A7", "MRC1", "CD68", "PTGDS", "PTPRC", "TNFRSF17", "MARCO", "FABP4",
                              "JCHAIN", "LYZ", "CD14", "TP73", "SFRP2"))

DotPlot(alv, features = c("CCL18", "MS4A7", "MRC1", "CD68", "PTGDS", "PTPRC", "TNFRSF17", "MARCO", "FABP4",
                          "JCHAIN", "LYZ", "CD14", "TP73", "SFRP2")) + angle

## Cluster 2 - PDGFRA+ FB

# Temp object
obj <- subset(alv, subset = alv_brl_0.4 == 2)

FeaturePlot(obj, features = c("PDGFRA", "MEG3", "COL1A2", "COL3A1", "SOD2", "POSTN", "COL1A1", "DMBT1", "DCN",
                              "FN1", "HYOU1", "LUM", "PDIA4", "VIM", "MMP7"))

DotPlot(alv, features = c(mesenchymal_features)) + angle

## Cluster 3 - Fetal-enriched transitional AT2

# Temp object
obj <- subset(alv, subset = alv_brl_0.4 == 3)

FeaturePlot(obj, features = c("MMP7", "HYOU1", "PDIA4", "SOD2", "DMBT1", "MGST1", "XBP1", "HSPA5",
                              "LAMP3", "SFTPC", "PGC", "KRT18", "KRT8", "SPCS2", "NUCB2", "AGER"))

DotPlot(alv, features = c("MMP7", "HYOU1", "PDIA4", "SOD2", "DMBT1", "MGST1", "XBP1", "HSPA5",
                          "LAMP3", "SFTPC", "PGC", "KRT18", "KRT8", "SPCS2", "NUCB2", "AGER")) + angle

## Cluster 4 - Monocytes (Recluster)

# Temp object
obj <- subset(alv, subset = alv_brl_0.4 == 4)

FeaturePlot(obj, features = c("S100A8", "S100A9", "S100A12", "ICAM1", "CCL2", "SOD2", "IRF1", "LYZ",
                              "EPAS1", "LUM", "SFTPD", "COL3A1"))

DotPlot(alv, features = c(mesenchymal_features)) + angle

## Cluster 5 - Fetal-enriched AT2

# Temp object
obj <- subset(alv, subset = alv_brl_0.4 == 5)

FeaturePlot(obj, features = c("HES1", "SOX9", "CD44", "WNT7B", "LGR5", "SOX4", "TGFB3", "LPAR2",
                              "FCGBP", "VEGFA", "BMP4", "CXCL14", "SPRY2", "CFTR", "MKI67", "SFTPC"))

DotPlot(alv, features = c("HES1", "SOX9", "CD44", "WNT7B", "LGR5", "SOX4", "TGFB3", "LPAR2",
                          "FCGBP", "VEGFA", "BMP4", "CXCL14", "SPRY2", "CFTR", "MKI67", "SFTPC")) + angle

## Cluster 6 - AT1 

# Temp object
obj <- subset(alv, subset = alv_brl_0.4 == 6)

FeaturePlot(obj, features = c("RTKN2", "ICAM1", "COL4A3", "AGER", "CEACAM6", "CEACAM5", "SCGB3A2", "KRT15", "SLC2A1",
                              "KRT17", "MSLN", "ITGA3", "KRT14", "VEGFA", "FGF2", "AKR1C1", "AKR1C2", "BCL2L1",
                              "SNCA", "KRT8"))

DotPlot(alv, group.by = "alv_brl_0.4", features = c("RTKN2", "GDF15", "COL4A3", "AGER", "CEACAM6", "CEACAM5", 
                                                    "KRT15", "SLC2A1", "KRT17", "MSLN", "ICAM1", "KRT14", "AGR3", 
                                                    "FGF2", "AKR1C1", "AKR1C2")) + angle

## Cluster 7 - Proliferating Fetal-enriched AT2/secretory (Recluster)

# Temp object
obj <- subset(alv, subset = alv_brl_0.4 == 7)

FeaturePlot(obj, features = c("MKI67", "TOP2A", "RTKN2", "HMGA1", "KRT8", "SFTPC", 
                              "LAMP3", "KRT5", "FOXJ1", "CALCA", "SCGB3A2", "C20orf85"))

DotPlot(alv, group.by = "alv_brl_0.4", features = c(epithelial_features)) + angle

### Reclustering of clusters 1 and 7 

cluster1 <- subset(alv, subset = alv_brl_0.4 == 1)
cluster7 <- subset(alv, subset = alv_brl_0.4 == 7)

saveRDS(cluster1, "/scratch/smallapragada/run2/run2_alv_reclustering_c1_2024_08_15.rds")
saveRDS(cluster7, "/scratch/smallapragada/run2/run2_alv_reclustering_c7_2024_08_15.rds")

## Cluster 1

cluster1 <- readRDS("/scratch/smallapragada/rapids_pipeline/run2_alv_reclustering_c1_2024_08_15.rds")

## UMAP visualization of clusters in reclustered alveolar object
DimPlot(cluster1,
        reduction = "umap",
        group.by = "alv_c1_brl_0.4", 
        raster = T,
        label = T, 
        cols = randomcoloR::distinctColorPalette(10)) 

## UMAP visualization of phenotype in reclustered alveolar object
DimPlot(cluster1,
        reduction = "umap",
        group.by = "Phenotype", 
        raster = T,
        label = F, 
        cols = randomcoloR::distinctColorPalette(5)) 

## A summary of how many cells per cluster belong to a specific disease phenotype
table(cluster1$alv_c1_brl_0.4, cluster1$Phenotype)

## Looking at the nCount_RNA and the nFeature_RNA distribution across clusters to determine low-quality clusters
ggplot(cluster1@meta.data, aes(factor(cluster1@meta.data$alv_c1_brl_0.4), 
                               (cluster1@meta.data$nCount_RNA))) + geom_violin()
ggplot(cluster1@meta.data, aes(factor(cluster1@meta.data$alv_c1_brl_0.4), 
                               (cluster1@meta.data$nFeature_RNA))) + geom_violin()

DotPlot(cluster1, group.by = "alv_c1_brl_0.4", features = c(other_features, epithelial_features)) + angle

## Top markers for reclustered object
Idents(cluster1) = "alv_c1_brl_0.4"
alv_markers <- FindAllMarkers(cluster1)
alv_markers %>%
  mutate(pct.diff = abs(pct.1 - pct.2)) %>%
  filter(avg_log2FC >= 0) %>%
  View()

FeaturePlot(cluster1, features = c("PGC", "EPAS1", "PECAM1", "SPARCL1", "FCN3", "IL7R", 
                                   "PTPRC", "CD34", "NOX4", "CA4", "BMPR2", "MEG3", "DUOX1"))     ## Cluster 0 

FeaturePlot(cluster1, features = c("CD8A", "CD28", "HLA-DRA", "RETN"))       ## Cluster 1

FeaturePlot(cluster1, features = c("KRT8", "KRT18", "SAA2", "AGR3" ,"LTF", "WFDC2", 
                                   "SCGB3A2", "SFTPC", "AGER"))     ## Cluster 2

FeaturePlot(cluster1, features = c("SCGB3A2", "VEGFA", "RTKN2", "SLC2A1" ,"MSLN", 
                                   "AGER", "HES1", "SFTPD", "NKX2-1", "CD44", "KRT8"))      ## Cluster 3

FeaturePlot(cluster1, features = c("HSPA5", "PDIA4", "MGST1", "LAMP3", "KRT18", 
                                   "HYOU1", "PDIA3", "SOD2", "HIF1A", "KRT8", 
                                   "MMP7", "HMGA1", "DMBT1", "XBP1", "SFTPC", "MKI67"))     ## Cluster 4

FeaturePlot(cluster1, features = c("COL3A1", "COL1A2", "COL1A1", "FN1", "LUM", 
                                   "VIM", "ACTA2", "POSTN", "DCN", "MEG3", "PTGDS", 
                                   "TP73", "CCNA1"))      ## Cluster 5

FeaturePlot(cluster1, features = c("HLA-DRA", "LYZ", "CD68", "NKX2-1", "DUOX1", 
                                   "ICAM1", "MSLN", "CCL18", "VIM", "NAPSA", "CD44", 
                                   "ITGA3", "VEGFA", "SOD2", "FCER1G"))       ## Cluster 6 

FeaturePlot(cluster1, features = c("IL11", "VIM", "NAPSA", "CD44", "LYZ", "CCL18", 
                                   "MSLN", "DUOX1", "LAMP3"))       ## Cluster 7

# Cluster 0 - Capillaries
# Cluster 1 - T
# Cluster 2 - Transitional AT2
# Cluster 3 - Fetal-enriched transitional AT2
# Cluster 4 - Transitional AT2
# Cluster 5 - FB
# Cluster 6 - Low quality
# Cluster 7 - Low quality

## Cluster 7

cluster7 <- readRDS("/scratch/smallapragada/rapids_pipeline/run2_alv_reclustering_c7_2024_08_15.rds")

## UMAP visualization of clusters in reclustered alveolar object
DimPlot(cluster7,
        reduction = "umap",
        group.by = "alv_c7_brl_0.4", 
        raster = T,
        label = T, 
        cols = randomcoloR::distinctColorPalette(10)) 

## UMAP visualization of phenotype in reclustered alveolar object
DimPlot(cluster7,
        reduction = "umap",
        group.by = "Phenotype", 
        raster = T,
        label = F, 
        cols = randomcoloR::distinctColorPalette(5)) 

## A summary of how many cells per cluster belong to a specific disease phenotype
table(cluster7$alv_c7_brl_0.4, cluster7$Phenotype)

## Looking at the nCount_RNA and the nFeature_RNA distribution across clusters to determine low-quality clusters
ggplot(cluster7@meta.data, aes(factor(cluster7@meta.data$alv_c7_brl_0.4), 
                               (cluster7@meta.data$nCount_RNA))) + geom_violin()
ggplot(cluster7@meta.data, aes(factor(cluster7@meta.data$alv_c7_brl_0.4), 
                               (cluster7@meta.data$nFeature_RNA))) + geom_violin()

DotPlot(cluster7, group.by = "alv_c7_brl_0.4", features = c(other_features, epithelial_features)) + angle

## Top markers for reclustered object
Idents(cluster7) = "alv_c7_brl_0.4"
alv_markers <- FindAllMarkers(cluster7)
alv_markers %>%
  mutate(pct.diff = abs(pct.1 - pct.2)) %>%
  filter(avg_log2FC >= 0) %>%
  View()

FeaturePlot(cluster7, features = c("SCGB1A1", "SCGB3A2", "KDR", "HES1", "VEGFA", "XBP1", 
                                   "AGR3", "GDF15", "SPARCL1", "HSPA5", "WWTR1", 
                                   "WNT5A", "CFTR", "FGF2", "CDH26"))        ## Cluster 0 

FeaturePlot(cluster7, features = c("TOP2A", "MKI67", "SOX9", "SOX4", "SFTPC", "KRT8", 
                                   "AGER"))       ## Cluster 1

FeaturePlot(cluster7, features = c("PDCD1", "LY6D", "SMAD4", "ITM2C", "SOX4", "CTNNB1", 
                                   "LMAN1", "BMP4", "SOX9", "AXL", "HMGA1", "CD44", 
                                   "SOX2"))       ## Cluster 2

FeaturePlot(cluster7, features = c("PGC", "NAPSA", "DUOX1", "MGST1", "EPAS1", "SFTPC", 
                                   "SFTPD", "LAMP3"))      ## Cluster 3

FeaturePlot(cluster7, features = c("SCGB1A1", "SCGB3A2", "C20orf85", "KDR"))       ## Cluster 4

FeaturePlot(cluster7, features = c("C20orf85", "COL1A2", "COL3A1", "MEG3", "COL1A1", 
                                   "SPARCL1", "FN1", "LUM", "POSTN", "PDGFRA", "DCN", 
                                   "FGF7", "ZEB1", "COL15A1", "ELN"))       ## Cluster 5

FeaturePlot(cluster7, features = c("SOD2", "DMBT1", "PGC", "XBP1", "LAMP3", "EPAS1", 
                                   "MMP7", "ICAM1", "SFTPC", "KRT8", "KRT5", "AGER", 
                                   "MKI67", "TOP2A"))        ## Cluster 6

FeaturePlot(cluster7, features = c("ASCL1", "CENFP", "ATF3", "SOX2", "KRT15", "SNCA", 
                                   "BANK1", "TP63", "S100A7", "TRAC", "AGER", "MKI67", 
                                   "TOP2A"))                 ## Cluster 7

FeaturePlot(cluster7, features = c("PGC", "NAPSA", "SFTPC", "SFTPD", "DMBT1", "EPAS1"))      ## Cluster 8

## Adding cell type labels to each of the reclustered objects
cluster1$CT <- ""
cluster1$CT[cluster1$alv_c1_brl_0.4 == 0] <- "Capillaries"
cluster1$CT[cluster1$alv_c1_brl_0.4 == 1] <- "T"
cluster1$CT[cluster1$alv_c1_brl_0.4 == 2] <- "Transitional AT2"
cluster1$CT[cluster1$alv_c1_brl_0.4 == 3] <- "Fetal-enriched transitional AT2"
cluster1$CT[cluster1$alv_c1_brl_0.4 == 4] <- "Transitional AT2"
cluster1$CT[cluster1$alv_c1_brl_0.4 == 5] <- "Fibroblasts"
cluster1$CT[cluster1$alv_c1_brl_0.4 == 6] <- "Low quality"
cluster1$CT[cluster1$alv_c1_brl_0.4 == 7] <- "Low quality"

cluster7$CT <- ""
cluster7$CT[cluster7$alv_c7_brl_0.4 == 0 | cluster7$alv_c7_brl_0.4 == 4] <- "Secretory 3A2+ & 1A1+"
cluster7$CT[cluster7$alv_c7_brl_0.4 == 1 | cluster7$alv_c7_brl_0.4 == 7 | cluster7$alv_c7_brl_0.4 == 2] <- "Proliferating Fetal-enriched AT2"
cluster7$CT[cluster7$alv_c7_brl_0.4 == 3 | cluster7$alv_c7_brl_0.4 == 8] <- "AT2"
cluster7$CT[cluster7$alv_c7_brl_0.4 == 5] <- "Multiciliated"
cluster7$CT[cluster7$alv_c7_brl_0.4 == 6] <- "Proliferating AT2"

## Combining objects together and merging into one
firstpass <- merge(cluster1, cluster7)

## Adding labels from the merged reclustered object to the whole alveolar object + all other CT labels 
alv$CT <- ""
alv$CT <- firstpass$CT
alv$CT[alv$alv_brl_0.4 == 0] <- "AT2"
alv$CT[alv$alv_brl_0.4 == 2] <- "Fetal-enriched fibroblasts"
alv$CT[alv$alv_brl_0.4 == 3] <- "Fetal-enriched transitional AT2"
alv$CT[alv$alv_brl_0.4 == 4] <- "Fibroblasts"
alv$CT[alv$alv_brl_0.4 == 5] <- "Fetal-enriched AT2"
alv$CT[alv$alv_brl_0.4 == 6] <- "AT1"

## Save alveolar object
saveRDS(alv, "/scratch/smallapragada/run2/final_alv_CT_2024_08_15.rds")

#### CELL TYPE ANNOTATION - ENDOTHELIAL ----

## Reading in endothelial object
endo <- readRDS("/scratch/smallapragada/rapids_pipeline/run2_endo_2024_08_15.rds")

## Looking at the nCount_RNA and the nFeature_RNA distribution across clusters to determine low-quality clusters
ggplot(endo@meta.data, aes(factor(endo@meta.data$endo_brl_0.3), (endo@meta.data$nCount_RNA))) + geom_violin()
ggplot(endo@meta.data, aes(factor(endo@meta.data$endo_brl_0.3), (endo@meta.data$nFeature_RNA))) + geom_violin()

## UMAP visualization of clusters in endothelial object
DimPlot(endo,
        reduction = "umap",
        group.by = "endo_brl_0.3", 
        raster = T,
        label = T, 
        cols = randomcoloR::distinctColorPalette(10)) 

## UMAP visualization of phenotype in endothelial object
DimPlot(endo,
        reduction = "umap",
        group.by = "Phenotype", 
        raster = T,
        label = F, 
        cols = randomcoloR::distinctColorPalette(5)) 

## A summary of how many cells per cluster belong to a specific disease phenotype
table(endo$endo_brl_0.3, endo$Phenotype)

## Top markers for endothelial
Idents(endo) = "endo_brl_0.3"
endo_markers <- FindAllMarkers(endo)
endo_markers %>%
  mutate(pct.diff = abs(pct.1 - pct.2)) %>%
  filter(avg_log2FC >= 0) %>%
  View()

# Endothelial dotplot
DotPlot(endo, group.by = "endo_brl_0.3", features = c(other_features, endothelial_features)) + angle

## Cluster 0 - low quality

# Temp object
obj <- subset(endo, subset = endo_brl_0.3 == 0)

FeaturePlot(obj, features = c("PECAM1"))

DotPlot(endo, features = c(other_features, endothelial_features)) + angle

## Cluster 1 - Meg-Ery

# Temp object
obj <- subset(endo, subset = endo_brl_0.3 == 1)

FeaturePlot(obj, features = c("SLC25A37", "SNCA", "SLC2A1", "BCL2L1", "PIM2", "TOP1", 
                              "BCL2L11", "GCLM", "PLVAP"))

DotPlot(endo, features = c(other_features, endothelial_features, 
                           "SLC25A37", "SNCA", "SLC2A1", "BCL2L1", "PIM2", "TOP1")) + angle

## Cluster 2 - Unsure (Recluster)

# Temp object
obj <- subset(endo, subset = endo_brl_0.3 == 2)

FeaturePlot(obj, features = c("SLC7A11", "ATF4", "GCLM", "HMGA1", "COL1A1", "PLIN2", "COL1A2", "COL3A1", 
                              "CCN2", "AKR1C1"))

DotPlot(endo, features = c("SLC7A11", "ATF4", "GCLM", "HMGA1", "COL1A1", "PLIN2", "COL1A2", "COL3A1", 
                           "CCN2", "AKR1C1")) + angle

## Cluster 3 - Venous (recluster)

# Temp object
obj <- subset(endo, subset = endo_brl_0.3 == 3)

FeaturePlot(obj, features = c("ACKR1", "ICAM1", "CCL2", "PTGDS", "HIF1A", "CCN2", "PLVAP", "IRF1", 
                              "COL15A1", "ITGAV"))

DotPlot(endo, features = c("ACKR1", "ICAM1", "CCL2", "PTGDS", "HIF1A", "CCN2", "PLVAP", "IRF1", 
                           "COL15A1", "ITGAV")) + angle

## Cluster 4 - Monocytes (recluster)

# Temp object
obj <- subset(endo, subset = endo_brl_0.3 == 4)

FeaturePlot(obj, features = c("S100A8", "S100A9", "S100A12", "SPP1", "LYZ", "FCER1G", "AIF1"))

DotPlot(endo, features = c("S100A8", "S100A9", "S100A12", "SPP1", "LYZ", "FCER1G", "AIF1")) + angle

## Cluster 5 - Lymphatic endothelial

# Temp object
obj <- subset(endo, subset = endo_brl_0.3 == 5)

FeaturePlot(obj, features = c("CCL21", "MRC1", "GNG11", "ITGB1", "CTHRC1", "COL1A1"))

DotPlot(endo, features = c("S100A8", "S100A9", "S100A12", "SPP1", "LYZ", "FCER1G", "AIF1")) + angle

## Cluster 6 - Capillaries

# Temp object
obj <- subset(endo, subset = endo_brl_0.3 == 6)

FeaturePlot(obj, features = c("CD34", "CA4", "APLN", "APLNR", "KDR", "GNG11", "BMPR2", "CLDN5"))

DotPlot(endo, features = c(endothelial_features)) + angle

## Cluster 7 - Fetal myoFB

# Temp object
obj <- subset(endo, subset = endo_brl_0.3 == 7)

FeaturePlot(obj, features = c("SOD2", "HSPA5", "HYOU1", "PDIA4", "PGC", "KRT18", "KRT8", "LAMP3", "SFTPC",
                              "ELANE", "DCN", "LUM", "PDGFRA", "HAS2", "VEGFA"))

DotPlot(endo, features = c("SOD2", "HSPA5", "HYOU1", "PDIA4", "PGC", "KRT18", "KRT8", "LAMP3", "SFTPC",
                           "ELANE")) + angle

### Reclustering of clusters 2, 3, and 4

cluster2 <- subset(endo, subset = endo_brl_0.3 == 2)
cluster3 <- subset(endo, subset = endo_brl_0.3 == 3)
cluster4 <- subset(endo, subset = endo_brl_0.3 == 4)

saveRDS(cluster2, "/scratch/smallapragada/run2/run2_endo_reclustering_c2_2024_08_15.rds")
saveRDS(cluster3, "/scratch/smallapragada/run2/run2_endo_reclustering_c3_2024_08_15.rds")
saveRDS(cluster4, "/scratch/smallapragada/run2/run2_endo_reclustering_c4_2024_08_15.rds")

## Cluster 2

cluster2 <- readRDS("/scratch/smallapragada/rapids_pipeline/run2_endo_reclustering_c2_2024_08_15.rds")

## UMAP visualization of clusters in reclustered endothelial object
DimPlot(cluster2,
        reduction = "umap",
        group.by = "endo_c2_brl_0.1", 
        raster = T,
        label = T, 
        cols = randomcoloR::distinctColorPalette(10)) 

## UMAP visualization of phenotype in reclustered endothelial object
DimPlot(cluster2,
        reduction = "umap",
        group.by = "Phenotype", 
        raster = T,
        label = F, 
        cols = randomcoloR::distinctColorPalette(5)) 

## A summary of how many cells per cluster belong to a specific disease phenotype
table(cluster2$endo_c2_brl_0.3, cluster2$Phenotype)

## Looking at the nCount_RNA and the nFeature_RNA distribution across clusters to determine low-quality clusters
ggplot(cluster2@meta.data, aes(factor(cluster2@meta.data$endo_c2_brl_0.3), 
                               (cluster2@meta.data$nCount_RNA))) + geom_violin()
ggplot(cluster2@meta.data, aes(factor(cluster2@meta.data$endo_c2_brl_0.3), 
                               (cluster2@meta.data$nFeature_RNA))) + geom_violin()

DotPlot(cluster2, group.by = "endo_c2_brl_0.3", features = c(other_features, endothelial_features)) + angle

## Top markers for endothelial reclustered object
Idents(cluster2) = "endo_c2_brl_0.3"
endo_markers <- FindAllMarkers(cluster2)
endo_markers %>%
  mutate(pct.diff = abs(pct.1 - pct.2)) %>%
  filter(avg_log2FC >= 0) %>%
  View()

FeaturePlot(cluster2, features = c("HIF1A", "ACTA2", "MGST1", "PLIN2", "SPARCL1", "KRT18", "FN1", "SOD2", "AKR1C1",
                                   "CA4", "CTNNB1", "KIT", "UBE2J1"))      ## Cluster 0, 2 and 3

DotPlot(cluster2, features = c("HIF1A", "ACTA2", "MGST1", "PLIN2", "SPARCL1", "KRT18", "FN1", "SOD2", "AKR1C1",
                               "CA4", "CTNNB1", "KIT", "UBE2J1")) + angle

# Cluster 0 - Capillaries
# Cluster 1 - Low quality
# Cluster 2 - Capillaries
# Cluster 3 - Capillaries

## Cluster 3

cluster3 <- readRDS("/scratch/smallapragada/rapids_pipeline/run2_endo_reclustering_c3_2024_08_15.rds")

## UMAP visualization of clusters in reclustered endothelial object
DimPlot(cluster3,
        reduction = "umap",
        group.by = "endo_c3_brl_0.4", 
        raster = T,
        label = T, 
        cols = randomcoloR::distinctColorPalette(10)) 

## UMAP visualization of phenotype in reclustered endothelial object
DimPlot(cluster3,
        reduction = "umap",
        group.by = "Phenotype", 
        raster = T,
        label = F, 
        cols = randomcoloR::distinctColorPalette(5)) 

## A summary of how many cells per cluster belong to a specific disease phenotype
table(cluster3$endo_c3_brl_0.4, cluster3$Phenotype)

## Looking at the nCount_RNA and the nFeature_RNA distribution across clusters to determine low-quality clusters
ggplot(cluster3@meta.data, aes(factor(cluster3@meta.data$endo_c3_brl_0.4), 
                               (cluster3@meta.data$nCount_RNA))) + geom_violin()
ggplot(cluster3@meta.data, aes(factor(cluster3@meta.data$endo_c3_brl_0.4), 
                               (cluster3@meta.data$nFeature_RNA))) + geom_violin()

DotPlot(cluster3, group.by = "endo_c3_brl_0.4", features = c(other_features, endothelial_features)) + angle

## Top markers for reclustered endothelial object
Idents(cluster3) = "endo_c3_brl_0.4"
endo_markers <- FindAllMarkers(cluster3)
endo_markers %>%
  mutate(pct.diff = abs(pct.1 - pct.2)) %>%
  filter(avg_log2FC >= 0) %>%
  View()

FeaturePlot(cluster3, features = c("CCL2", "SOD2", "IRF1", "ICAM1", "IL7R", "CA4", 
                                   "SPARCL1", "HIF1A", "CD44", "IL32", "S100A8", "LYZ", "CLDN5"))     ## Cluster 0 

FeaturePlot(cluster3, features = c("MMP10", "MEG3", "NHSL2", "ITGB1", "HAS2", "POSTN", 
                                   "SPCS2", "MYC", "HSPA5", "PDIA6", "HYOU1", "TGFB1", 
                                   "PTGDS", "TGFB2", "VEGFA"))        ## Cluster 1 

FeaturePlot(cluster3, features = c("COL15A1", "PLVAP", "KDR", "CD34", "COL1A2", "COL3A1", 
                                   "COL1A1", "POSTN", "CCN2", "FABP4", "CTHRC1", "SFRP2", 
                                   "ACKR1", "APLNR"))                 ## Cluster 2

FeaturePlot(cluster3, features = c("KRT17", "IL32", "SELENOS", "FCER1G", "TP63", "WNT7B", 
                                   "BMP4", "EPAS1", "TNFRSF13C"))     ## Cluster 3

FeaturePlot(cluster3, features = c("FCN3", "ACKR1", "PGC", "BMPR2", "XBP1", "SLC7A11", 
                                   "ATF4", "MYC", "IL4R"))            ## Cluster 4

FeaturePlot(cluster3, features = c("S100A8", "CCN2", "S100A9", "CCL2", "LYZ", "CD44", 
                                   "HSPA5"))                           ## Cluster 5

FeaturePlot(cluster3, features = c("PTGDS", "ACKR1", "APLNR", "RNASE1", "MEG3", 
                                   "POSTN", "ACTA2", "FN1"))          ## Cluster 6

FeaturePlot(cluster3, features = c("LILRA4", "GZMB", "BMPR2", "RHOA", "KDR", 
                                   "FN1", "GNLY", "AXL"))             ## Cluster 7

DotPlot(cluster3, features = c("CCL2", "SOD2", "IRF1", "ICAM1", "IL7R", "CA4", 
                               "SPARCL1", "HIF1A", "CD44", "IL32", "S100A8", "LYZ", 
                               "CLDN5", "MMP10", "MEG3", "NHSL2", "ITGB1", "HAS2", 
                               "POSTN", "SPCS2", "MYC", "HSPA5", "PDIA6", "HYOU1", 
                               "TGFB1", "TGFB2", "VEGFA", "COL15A1", "PLVAP", "KDR", 
                               "CD34", "COL1A2", "COL3A1", "COL1A1", "FABP4", "CTHRC1", 
                               "SFRP2", "ACKR1", "APLNR", "KRT17", "SELENOS", "FCER1G", 
                               "TP63", "WNT7B", "BMP4", "EPAS1", "TNFRSF13C", "FCN3", 
                               "PGC", "BMPR2", "XBP1", "SLC7A11", "ATF4", "MYC", "IL4R",
                               "CCN2", "S100A9", "PTGDS",  "RNASE1", "ACTA2", "FN1", 
                               "LILRA4", "GZMB", "RHOA",  "GNLY", "AXL")) + angle
# Cluster 0 - Capillaries
# Cluster 1 - Basal
# Cluster 2 - Capillaries
# Cluster 3 - Basal
# Cluster 4 - Venous
# Cluster 5 - Monocytes
# Cluster 6 - Capillaries
# Cluster 7 - NK & NKT

#### Cluster 4 - reclustered

cluster4 <- readRDS("/scratch/smallapragada/rapids_pipeline/run2_endo_reclustering_c4_2024_08_15.rds")

## UMAP visualization of clusters in reclustered endothelial object
DimPlot(cluster4,
        reduction = "umap",
        group.by = "endo_c4_brl_0.4", 
        raster = T,
        label = T, 
        cols = randomcoloR::distinctColorPalette(10)) 

## UMAP visualization of clusters in reclustered endothelial object
DimPlot(cluster4,
        reduction = "umap",
        group.by = "Phenotype", 
        raster = T,
        label = F, 
        cols = randomcoloR::distinctColorPalette(5)) 

## A summary of how many cells per cluster belong to a specific disease phenotype
table(cluster4$endo_c4_brl_0.3, cluster4$Phenotype)

## Looking at the nCount_RNA and the nFeature_RNA distribution across clusters to determine low-quality clusters
ggplot(cluster4@meta.data, aes(factor(cluster4@meta.data$endo_c4_brl_0.3), 
                               (cluster4@meta.data$nCount_RNA))) + geom_violin()
ggplot(cluster4@meta.data, aes(factor(cluster4@meta.data$endo_c4_brl_0.3), 
                               (cluster4@meta.data$nFeature_RNA))) + geom_violin()

DotPlot(cluster4, group.by = "endo_c4_brl_0.3", features = c(other_features, endothelial_features)) + angle

## Top markers
Idents(cluster4) = "endo_c4_brl_0.3"
endo_markers <- FindAllMarkers(cluster4)
endo_markers %>%
  mutate(pct.diff = abs(pct.1 - pct.2)) %>%
  filter(avg_log2FC >= 0) %>%
  View()

FeaturePlot(cluster4, features = c("CALCA", "SCG2", "CHGB", "FOXI1", "STAT6"))     ## Cluster 0 

FeaturePlot(cluster4, features = c("ELANE", "NUCB2", "LYZ", "CEACAM6", "HMGA1", "PDIA3"))      ## Cluster 1

FeaturePlot(cluster4, features = c("FCER1G", "FCGR3A", "PTPRC", "AIF1", "MS4A7", "CD44"))        ## Cluster 2

FeaturePlot(cluster4, features = c("FOXJ1", "ATF4", "PGC", "JCHAIN", "CST3", "SNAI1"))      ## Cluster 3

FeaturePlot(cluster4, features = c("S100A12", "S100A8", "S100A9", "PECAM1", "CD34", "IL7R"))     ## Cluster 4

FeaturePlot(cluster4, features = c("COL1A2", "PDGFRB", "COL3A1", "FN1", "LUM", "AXL"))       ## Cluster 5

DotPlot(cluster4, features = c("CALCA", "SCG2", "CHGB", "FOXI1", "STAT6", "ELANE", 
                               "NUCB2", "LYZ", "CEACAM6", "HMGA1", "PDIA3", "FCER1G", 
                               "FCGR3A", "PTPRC", "AIF1", "MS4A7", "CD44", "FOXJ1", 
                               "ATF4", "PGC", "JCHAIN", "CST3", "SNAI1", "S100A12", 
                               "S100A8", "S100A9", "PECAM1", "CD34", "IL7R", 
                               "COL1A2", "PDGFRB", "COL3A1", "FN1", "LUM", "AXL")) + angle

# Cluster 0 - PNEC
# Cluster 1 - Neutrophils
# Cluster 2 - cDC
# Cluster 3 - Plasma
# Cluster 4 - Monocytes
# Cluster 5 - Pericytes

## Adding cell type labels to each of the reclustered objects
cluster2$CT <- ""
cluster2$CT[cluster2$endo_c2_brl_0.3 == 1] <- "Low quality"
cluster2$CT[cluster2$endo_c2_brl_0.3 == 0 | 
              cluster2$endo_c2_brl_0.3 == 2 | 
              cluster2$endo_c2_brl_0.3 == 3 | 
              cluster2$endo_c2_brl_0.3 == 4] <- "Capillaries"

cluster3$CT <- ""
cluster3$CT[cluster3$endo_c3_brl_0.4 == 0 | 
              cluster3$endo_c3_brl_0.4 == 2 | 
              cluster3$endo_c3_brl_0.4 == 6] <- "Capillaries"
cluster3$CT[cluster3$endo_c3_brl_0.4 == 1 | 
              cluster3$endo_c3_brl_0.4 == 3] <- "Basal"
cluster3$CT[cluster3$endo_c3_brl_0.4 == 4] <- "Venous"
cluster3$CT[cluster3$endo_c3_brl_0.4 == 5] <- "Monocytes"
cluster3$CT[cluster3$endo_c3_brl_0.4 == 7] <- "NK & NKT"

cluster4$CT <- ""
cluster4$CT[cluster4$endo_c4_brl_0.3 == 0] <- "PNEC"
cluster4$CT[cluster4$endo_c4_brl_0.3 == 1] <- "Neutrophils"
cluster4$CT[cluster4$endo_c4_brl_0.3 == 2] <- "cDC"
cluster4$CT[cluster4$endo_c4_brl_0.3 == 3] <- "Plasma"
cluster4$CT[cluster4$endo_c4_brl_0.3 == 4] <- "Monocytes"
cluster4$CT[cluster4$endo_c4_brl_0.3 == 5] <- "Pericytes"

## Combining objects together and merging into one
obj_list <- list(cluster2, cluster3, cluster4)
firstpass <- merge(x = obj_list[[1]], y = obj_list[2:length(obj_list)])

## Adding labels from the merged reclustered object to the whole endothelial object + all other CT labels 
endo$CT <- ""
endo$CT <- firstpass$CT
endo$CT[endo$endo_brl_0.3 == 0] <- "Low quality"
endo$CT[endo$endo_brl_0.3 == 1] <- "Meg-Ery"
endo$CT[endo$endo_brl_0.3 == 5] <- "Lymphatic endothelial"
endo$CT[endo$endo_brl_0.3 == 6] <- "Capillaries"
endo$CT[endo$endo_brl_0.3 == 7] <- "Fetal-enriched SOD2+ myoFB"

## Save endothelial object
saveRDS(endo, "/scratch/smallapragada/run2/final_endo_CT_2024_08_15.rds")

#### CELL TYPE ANNOTATION - MESENCHYMAL ----

## Reading in mesenchymal object
mes <- readRDS("/scratch/smallapragada/rapids_pipeline/run2_mes_2024_08_15.rds")

## Looking at the nCount_RNA and the nFeature_RNA distribution across clusters to determine low-quality clusters
ggplot(mes@meta.data, aes(factor(mes@meta.data$mes_brl_0.3), (mes@meta.data$nCount_RNA))) + geom_violin()
ggplot(mes@meta.data, aes(factor(mes@meta.data$mes_brl_0.3), (mes@meta.data$nFeature_RNA))) + geom_violin()

## UMAP visualization of clusters in mesenchymal object
DimPlot(mes,
        reduction = "umap",
        group.by = "mes_brl_0.3", 
        raster = T,
        label = T, 
        cols = randomcoloR::distinctColorPalette(20)) 

## UMAP visualization of phenotypes in mesenchymal object
DimPlot(mes,
        reduction = "umap",
        group.by = "Phenotype", 
        raster = T,
        label = F, 
        cols = randomcoloR::distinctColorPalette(5)) 

## A summary of how many cells per cluster belong to a specific disease phenotype
table(mes$mes_brl_0.3, mes$Phenotype)

DotPlot(mes, group.by = "mes_brl_0.3", features = c(other_features, mesenchymal_features)) + angle

## Top markers for mesenchymal
Idents(mes) = "mes_brl_0.3"
mes_markers <- FindAllMarkers(mes)
mes_markers %>%
  mutate(pct.diff = abs(pct.1 - pct.2)) %>%
  filter(avg_log2FC >= 0) %>%
  View()

## Cluster 0 - FB (subcluster)

# Temp object
obj <- subset(mes, subset = mes_brl_0.3 == 0)

FeaturePlot(obj, features = c("SFRP2", "PTGDS", "POSTN", "LUM", "CTHRC1", "DCN", 
                              "FGF7", "FAP", "SFRP4", "ANKRD28", "COL3A1"))

DotPlot(mes, features = c("SFRP2", "PTGDS", "POSTN", "LUM", "CTHRC1", "DCN", 
                          "FGF7", "FAP", "SFRP4", "ANKRD28", "COL3A1")) + angle

## Cluster 1 - Adventitial FB

# Temp object
obj <- subset(mes, subset = mes_brl_0.3 == 1)

FeaturePlot(obj, features = c("PI16", "ICAM1", "CCL2", "SFRP4", "SFRP2", "HYOU1", 
                              "DCN", "HSPA5", "MFAP5", "FGF7", "HIF1A"))

DotPlot(mes, features = c("PI16", "ICAM1", "CCL2", "SFRP4", "SFRP2", "HYOU1", "DCN", 
                          "HSPA5", "MFAP5", "FGF7", "HIF1A")) + angle

## Cluster 2 - Pericytes/general endothelial (subcluster with 13)

# Temp object
obj <- subset(mes, subset = mes_brl_0.3 == 2)

FeaturePlot(obj, features = c("PDGFRB", "PECAM1", "CSPG4", "AXL", "ITGB1", "ITM2C", 
                              "CD34", "SPARCL1", "CTNNB1", "KDR", "CD4", "ACTA2"))

DotPlot(mes, features = c("PDGFRB", "PECAM1", "CSPG4", "AXL", "ITGB1", "ITM2C", 
                          "CD34", "SPARCL1", "CTNNB1", "KDR", "CD4", "BCL2L11")) + angle

# Plot to spatially visualize where certain cells are located
samp <- subset(mes, subset = Sample == "s1_PDL003")

DimPlot(samp,
        reduction = "SP",
        group.by = "mes_brl_0.3", 
        raster = T,
        label = F,
        cols = c(`2` = "red"), na.value = "grey90") + coord_equal()

## Cluster 3 - SOD2+ FB

# Temp object
obj <- subset(mes, subset = mes_brl_0.3 == 3)

FeaturePlot(obj, features = c("ICAM1", "CCL2", "SOD2", "CD44", "LUM", "HIF1A", 
                              "FAS", "PTGDS", "FGF7", "DCN", "HYOU1", "HLA-DRA"))

DotPlot(mes, features = c(mesenchymal_features)) + angle

## Cluster 4 - Adventitial FB

# Temp object
obj <- subset(mes, subset = mes_brl_0.3 == 4)

FeaturePlot(obj, features = c("MFAP5", "PI16", "CCN2", "DCN", "CST3", "CXCL14", 
                              "CD34", "CCL2", "COL3A1", "SFRP2", "SFRP4", "HAS1"))

DotPlot(mes, features = c("MFAP5", "PI16", "CCN2", "DCN", "CST3", "CXCL14", "CD34", 
                          "CCL2", "COL3A1", "SFRP2", "SFRP4", "HAS1")) + angle

## Cluster 5 - AKR1C1+/1C2+

# Temp object
obj <- subset(mes, subset = mes_brl_0.3 == 5)

FeaturePlot(obj, features = c("AKR1C1", "AKR1C2", "IL11", "HMGA1", "MGST1", "PKM", "POSTN", "CD44", 
                              "ATF4", "HAS2", "HIF1A", "PLIN2", "MYC", "FKBP11", "XBP1", "PDIA4"))

DotPlot(mes, features = c("AKR1C1", "AKR1C2", "IL11", "HMGA1", "MGST1", "PKM", "POSTN", "CD44", 
                          "ATF4", "HAS2", "HIF1A", "PLIN2", "MYC", "FKBP11", "XBP1", "PDIA4")) + angle

## Cluster 6 - Low quality

# Temp object
obj <- subset(mes, subset = mes_brl_0.3 == 6)

FeaturePlot(obj, features = c("COL1A1", "LUM", "COL3A1"))

DotPlot(obj, features = c(mesenchymal_features)) + angle

## Cluster 7 - FB (Subcluster)

# Temp object
obj <- subset(mes, subset = mes_brl_0.3 == 7)

FeaturePlot(obj, features = c("LUM", "POSTN", "DCN", "FGF7", "HAS2", "PDIA3", 
                              "CD44", "YAP1", "FAP", "FAS", "SPARCL1", "PDIA4", 
                              "PDGFRA", "ACTA2", "PLIN2"))

DotPlot(mes, features = c(mesenchymal_features)) + angle

## Cluster 8 - PLIN2+ FB

# Temp object
obj <- subset(mes, subset = mes_brl_0.3 == 8)

FeaturePlot(obj, features = c("GDF15", "PLIN2", "HIF1A", "VEGFA", "DDIT3", "HSPA5", 
                              "FGF7", "ICAM1", "FAP", "POSTN", "SOX4", "XBP1", 
                              "LUM", "FGF2", "HERPUD1", "HMGA1"))

DotPlot(mes, features = c("GDF15", "PLIN2", "HIF1A", "VEGFA", "DDIT3", "HSPA5", 
                          "FGF7", "ICAM1", "FAP", "POSTN", "SOX4", "XBP1", "LUM", 
                          "FGF2", "HERPUD1", "HMGA1", epithelial_features)) + angle

## Cluster 9 - SMC

# Temp object
obj <- subset(mes, subset = mes_brl_0.3 == 9)

FeaturePlot(obj, features = c("ACTA2", "WNT5A", "LGR6", "SPARCL1", "CXCL14", "COL15A1", 
                              "SOX4", "CCN2", "ITGB1", "TGFB2"))

DotPlot(mes, features = c("ACTA2", "WNT5A", "LGR6", "SPARCL1", "CXCL14", "COL15A1", 
                          "SOX4", "CCN2", "ITGB1", "TGFB2")) + angle

# Plot to spatially visualize where certain cells are located
samp <- subset(mes, subset = Sample == "s1_PDL001")

DimPlot(samp,
        reduction = "SP",
        group.by = "mes_brl_0.3", 
        raster = T,
        label = F,
        cols = c(`9` = "red"), na.value = "grey90") + coord_equal()

# Plot to spatially visualize the presence of specific genes
FeaturePlot(samp, reduction = "SP", features = c("ACTA2", "CXCL14"))

## Cluster 10 - Adventitial FB
obj <- subset(mes, subset = mes_brl_0.3 == 10 | mes_brl_0.3 == 4 | mes_brl_0.3 == 1)

FeaturePlot(obj, features = c("MFAP5", "PI16", "DCN", "CD34", "CCN2", "CST3", "LUM", "SFRP2", "SPRY2", "HERPUD1"))

DotPlot(obj, features = c("MFAP5", "PI16", "DCN", "CD34", "CCN2", "CST3", "LUM", "SFRP2", "SPRY2", "HERPUD1")) + angle

## Cluster 11 - Monocytes (Recluster)
obj <- subset(mes, subset = mes_brl_0.3 == 11)

FeaturePlot(obj, features = c("S100A8", "S100A9", "S100A12", "LYZ", "FCER1G", "MCEMP1", "FCN1", "AIF1", "RETN", "LTF"))

DotPlot(mes, features = c("S100A8", "S100A9", "S100A12", "LYZ", "FCER1G", "MCEMP1", "FCN1", "AIF1", "RETN", "LTF")) + angle

## Cluster 12 - SPP1+ macrophages/pDC (Recluster)
obj <- subset(mes, subset = mes_brl_0.3 == 12)

FeaturePlot(obj, features = c("HLA-DRA", "FCER1G", "SPP1", "MS4A7", "CD68", "MRC1", "AIF1", "CD44", "C1QC", "GPR183", "HMOX1",
                              "IDH1", "RNASE1"))

DotPlot(mes, features = c("HLA-DRA", "FCER1G", "SPP1", "MS4A7", "CD68", "MRC1", "AIF1", "CD44", "C1QC", "GPR183", "HMOX1",
                          "IDH1", "RNASE1", "S100A8", "S100A9", "S100A12", "MARCO")) + angle

## Cluster 13 - Pericytes/general endothelial (Absorb into 2)
obj <- subset(mes, subset = mes_brl_0.3 == 13)

FeaturePlot(obj, features = c("HIF1A", "CCL2", "EPAS1", "AXL", "ICAM1", "ACTA2", "SOD2", "COL15A1", "SLC25A4", "ITGB1", "HSPA5",
                              "CSPG4", "PDGFRB"))

DotPlot(mes, features = c("HIF1A", "CCL2", "EPAS1", "AXL", "ICAM1", "ACTA2", "SOD2", "COL15A1", "SLC25A4", "ITGB1", "HSPA5",
                          "CSPG4", "PDGFRB")) + angle

## Cluster 14 - PDGFRA+ FB (Absorb into 15 - Subcluster)
obj <- subset(mes, subset = mes_brl_0.3 == 14)

FeaturePlot(obj, features = c("SOD2", "PDGFRA", "COL3A1", "LUM", "DCN", "CD44", "STAT1", "HIF1A", "HIST1H1C", "HSPA5", "PDIA4", "YAP1", "PDIA3", "AXIN2",
                              "HYOU1", "PDGFRB"))

DotPlot(mes, features = c("SOD2", "PDGFRA", "COL3A1", "LUM", "DCN", "CD44", "STAT1", "HIF1A", "HIST1H1C", "HSPA5", "PDIA4", "YAP1", "PDIA3", "AXIN2",
                          "HYOU1", "PDGFRB")) + angle

## Cluster 15 - General mesenchymal (Absorb into 14 - Subcluster)
obj <- subset(mes, subset = mes_brl_0.3 == 15)

FeaturePlot(obj, features = c("SOD2", "PDGFRA", "COL3A1", "LUM", "DCN", "CD44", "STAT1", "HIF1A", "HIST1H1C", "HSPA5", "PDIA4", "YAP1", "PDIA3", "AXIN2",
                              "HYOU1", "PDGFRB"))

DotPlot(mes, features = c("SOD2", "PDGFRA", "COL3A1", "LUM", "DCN", "CD44", "STAT1", "HIF1A", "HIST1H1C", "HSPA5", "PDIA4", "YAP1", "PDIA3", "AXIN2",
                          "HYOU1", "PDGFRB")) + angle

### Subclustering of clusters 0, 2, and 7

mes_subclustering <- readRDS("/scratch/smallapragada/rapids_pipeline/run2_mes_subclustering_2024_08_15.rds")

## Cluster 0 

## UMAP visualization of clusters in subclustered mesenchymal object
DimPlot(mes_subclustering,
        reduction = "umap",
        group.by = "mes_sc_0.3_c0_0.5", 
        raster = T,
        label = T, 
        cols = randomcoloR::distinctColorPalette(20)) 

## UMAP visualization of phenotypes in subclustered mesenchymal object
DimPlot(mes_subclustering,
        reduction = "umap",
        group.by = "Phenotype", 
        raster = T,
        label = F, 
        cols = randomcoloR::distinctColorPalette(5)) 

## Temp object for the specific clusters of interest
obj <- subset(mes_subclustering, subset = mes_sc_0.3_c0_0.5 == "0,0" | 
                mes_sc_0.3_c0_0.5 == "0,1" | mes_sc_0.3_c0_0.5 == "0,2" |
                mes_sc_0.3_c0_0.5 == "0,3" | mes_sc_0.3_c0_0.5 == "0,4" | 
                mes_sc_0.3_c0_0.5 == "0,5")

## A summary of how many cells per cluster belong to a specific disease phenotype
table(obj$mes_sc_0.3_c0_0.5, obj$Phenotype)

## Looking at the nCount_RNA and the nFeature_RNA distribution across clusters to determine low-quality clusters
ggplot(obj@meta.data, aes(factor(obj@meta.data$mes_sc_0.3_c0_0.5), (obj@meta.data$nCount_RNA))) + geom_violin()
ggplot(obj@meta.data, aes(factor(obj@meta.data$mes_sc_0.3_c0_0.5), (obj@meta.data$nFeature_RNA))) + geom_violin()

DotPlot(obj, group.by = "mes_sc_0.3_c0_0.5", features = c(other_features, mesenchymal_features)) + angle

## Top markers for temp object
Idents(obj) = "mes_sc_0.3_c0_0.5"
mes_markers <- FindAllMarkers(obj)
mes_markers %>%
  mutate(pct.diff = abs(pct.1 - pct.2)) %>%
  filter(avg_log2FC >= 0) %>%
  View()

FeaturePlot(obj, features = c("SFRP2", "PTGDS", "POSTN", "LUM", "CTHRC1", "DCN", 
                              "FGF7", "ATF4", "SFRP4", "PDGFRA", "ACTA2", "HLA-DRA", 
                              "MS4A7", "ITGAX", "MRC1", "FCER1G", "GPR183", "PLIN2", 
                              "JCHAIN", "COL15A1", "WNT5A", "YAP1", "MEG3", "SOD2", 
                              "KLRB1", "CD2", "CD8A", "GZMA", "TRAC"))

DotPlot(obj, features = c("SFRP2", "PTGDS", "POSTN", "LUM", "CTHRC1", "DCN", 
                          "FGF7", "ATF4", "SFRP4", "PDGFRA", "ACTA2", "HLA-DRA", 
                          "MS4A7", "ITGAX", "MRC1", "FCER1G", "GPR183", "PLIN2", 
                          "JCHAIN", "COL15A1", "WNT5A", "YAP1", "MEG3", "SOD2", 
                          "KLRB1", "CD2", "CD8A", "GZMA", "TRAC")) + angle

# 0,0 - Activated FB
# 0,1 - Low quality
# 0,2 - Activated FB
# 0,3 - Activated FB
# 0,4 - FB
# 0,5 - T

## Clusters 2 and 13 (combined)

## UMAP visualization of clusters in subclustered mesenchymal object
DimPlot(mes_subclustering,
        reduction = "umap",
        group.by = "mes_sc_0.3_c2_13_0.3", 
        raster = T,
        label = T, 
        cols = randomcoloR::distinctColorPalette(20)) 

## UMAP visualization of phenotypes in subclustered mesenchymal object
DimPlot(mes_subclustering,
        reduction = "umap",
        group.by = "Phenotype", 
        raster = T,
        label = F, 
        cols = randomcoloR::distinctColorPalette(5)) 

## Temp object for the specific clusters of interest
obj <- subset(mes_subclustering, subset = mes_sc_0.3_c2_13_0.2 == "2-13,0" | 
                mes_sc_0.3_c2_13_0.2 == "2-13,1" | 
                mes_sc_0.3_c2_13_0.2 == "2-13,2")

## A summary of how many cells per cluster belong to a specific disease phenotype
table(obj$mes_sc_0.3_c2_13_0.2, obj$Phenotype)

## Looking at the nCount_RNA and the nFeature_RNA distribution across clusters to determine low-quality clusters
ggplot(obj@meta.data, aes(factor(obj@meta.data$mes_sc_0.3_c2_13_0.2), (obj@meta.data$nCount_RNA))) + geom_violin()
ggplot(obj@meta.data, aes(factor(obj@meta.data$mes_sc_0.3_c2_13_0.2), (obj@meta.data$nFeatures_RNA))) + geom_violin()

DotPlot(obj, group.by = "mes_sc_0.3_c2_13_0.2", features = c(other_features, mesenchymal_features)) + angle

## Top markers for temp object
Idents(obj) = "mes_sc_0.3_c2_13_0.2"
mes_markers <- FindAllMarkers(obj)
mes_markers %>%
  mutate(pct.diff = abs(pct.1 - pct.2)) %>%
  filter(avg_log2FC >= 0) %>%
  View()

FeaturePlot(obj, features = c("PDGFRB", "CSPG4", "AXL", "ITM2C", "PECAM1", "KDR", 
                              "CD34", "PDGFRA", "CA4", "SOD2", "KIT", "FCN3", "CCL2", 
                              "ICAM1", "COL15A1", "ACTA2", "SLC25A4", "FGF2"))

DotPlot(obj, features = c("PDGFRB", "CSPG4", "AXL", "ITM2C", "PECAM1", "KDR", "CD34", 
                          "PDGFRA", "CA4", "SOD2", "KIT", "FCN3", "CCL2", "ICAM1", 
                          "COL15A1", "ACTA2", "SLC25A4", "FGF2")) + angle

# 2-13,0 - Pericytes
# 2-13,1 - SOD2+ FB
# 2-13,2 - Pericytes

## Cluster 7 

## UMAP visualization of clusters in subclustered mesenchymal object
DimPlot(mes_subclustering,
        reduction = "umap",
        group.by = "mes_sc_0.3_c7_0.2", 
        raster = T,
        label = T, 
        cols = randomcoloR::distinctColorPalette(20)) 

## UMAP visualization of phenotypes in subclustered mesenchymal object
DimPlot(mes_subclustering,
        reduction = "umap",
        group.by = "Phenotype", 
        raster = T,
        label = F, 
        cols = randomcoloR::distinctColorPalette(5)) 

## Temp object for the specific clusters of interest
obj <- subset(mes_subclustering, subset = mes_sc_0.3_c7_0.2 == "7,0" | 
                mes_sc_0.3_c7_0.2 == "7,1" | 
                mes_sc_0.3_c7_0.2 == "7,2")

## A summary of how many cells per cluster belong to a specific disease phenotype
table(obj$mes_sc_0.3_c7_0.2, obj$Phenotype)

## Looking at the nCount_RNA and the nFeature_RNA distribution across clusters to determine low-quality clusters
ggplot(obj@meta.data, aes(factor(obj@meta.data$mes_sc_0.3_c7_0.2), (obj@meta.data$nCount_RNA))) + geom_violin()
ggplot(obj@meta.data, aes(factor(obj@meta.data$mes_sc_0.3_c7_0.2), (obj@meta.data$nFeature_RNA))) + geom_violin()

DotPlot(obj, group.by = "mes_sc_0.3_c7_0.2", features = c(other_features, mesenchymal_features, "SOD2")) + angle

## Top markers for temp object
Idents(obj) = "mes_sc_0.3_c7_0.2"
mes_markers <- FindAllMarkers(obj)
mes_markers %>%
  mutate(pct.diff = abs(pct.1 - pct.2)) %>%
  filter(avg_log2FC >= 0) %>%
  View()

FeaturePlot(obj, features = c("PDGFRA", "PTGDS", "SOD2", "SOX9", "PGC", "FCN3", 
                              "HMGA1", "MGST1", "PLIN2", "AKR1C1", "PKM", "CD44", 
                              "POSTN", "ATF4", "IL11", "AKR1C2", "MYC", "PDIA3", 
                              "PDIA4", "HYOU1"))

DotPlot(obj, features = c("PDGFRB", "CSPG4", "AXL", "ITM2C", "PECAM1", "KDR", "CD34", 
                          "PDGFRA", "CA4", "SOD2", "KIT", "FCN3", "CCL2", "ICAM1", 
                          "COL15A1", "ACTA2", "SLC25A4", "FGF2")) + angle

# 7,0 - SOD2+ FB
# 7,1 - AKR1C1+/1C2+
# 7,2 - SOD2+ FB

### Reclustering of clusters 11, 12, 14, and 15
cluster11 <- subset(mes, subset = mes_brl_0.3 == 11)
cluster12 <- subset(mes, subset = mes_brl_0.3 == 12)
cluster14_15 <- subset(mes, subset = mes_brl_0.3 == 14 | mes_brl_0.3 == 15)

saveRDS(cluster11, "/scratch/smallapragada/run2/run2_mes_reclustering_c11_2024_08_15.rds")
saveRDS(cluster12, "/scratch/smallapragada/run2/run2_mes_reclustering_c12_2024_08_15.rds")
saveRDS(cluster14_15, "/scratch/smallapragada/run2/run2_mes_reclustering_c14_c15_2024_08_15.rds")

## Cluster 11

cluster11 <- readRDS("/scratch/smallapragada/rapids_pipeline/run2_mes_reclustering_c11_2024_08_15.rds")

## UMAP visualization of clusters in reclustered mesenchymal object
DimPlot(cluster11,
        reduction = "umap",
        group.by = "mes_c11_brl_0.4", 
        raster = T,
        label = T, 
        cols = randomcoloR::distinctColorPalette(20)) 

## UMAP visualization of phenotype in reclustered mesenchymal object
DimPlot(cluster11,
        reduction = "umap",
        group.by = "Phenotype", 
        raster = T,
        label = F, 
        cols = randomcoloR::distinctColorPalette(5)) 

## A summary of how many cells per cluster belong to a specific disease phenotype
table(cluster11$mes_c11_brl_0.4, cluster11$Phenotype)

## Looking at the nCount_RNA and the nFeature_RNA distribution across clusters to determine low-quality clusters
ggplot(cluster11@meta.data, aes(factor(cluster11$mes_c11_brl_0.3), (cluster11$nCount_RNA))) + geom_violin()
ggplot(cluster11@meta.data, aes(factor(cluster11$mes_c11_brl_0.3), (cluster11$nFeature_RNA))) + geom_violin()

DotPlot(cluster11, group.by = "mes_c11_brl_0.4", features = c(other_features, mesenchymal_features, "SOD2")) + angle

## Top markers
Idents(cluster11) = "mes_c11_brl_0.4"
markers <- FindAllMarkers(cluster11)
markers %>%
  mutate(pct.diff = abs(pct.1 - pct.2)) %>%
  filter(avg_log2FC >= 0) %>%
  View()

FeaturePlot(cluster11, features = c("LCK", "KRT17", "S100A8", "S100A9", "PDGFRA", "SOD2", 
                                    "LYZ", "S100A12", "EGFR"))          ## Cluster 0

FeaturePlot(cluster11, features = c("PDGFRB", "EPAS1", "AXL", "CSPG4", "ITGB1", "ITM2C", 
                                    "SPARCL1", "BCL2L1", "TOP2A")).     ## Cluster 1

FeaturePlot(cluster11, features = c("ACTA2", "WNT5A", "SPARCL1", "POSTN", "AXIN2", "ELN", 
                                    "CXCL14", "LGR6", "CTNNB1"))        ## Cluster 2

FeaturePlot(cluster11, features = c("CD8B", "MMP10", "VPREB3", "S100A12", "LMAN1", "SLC1A3", 
                                    "MCEMP1", "LTF", "SPP1"))           ## Cluster 3

FeaturePlot(cluster11, features = c("LUM", "SPARCL1", "DCN", "POSTN", "MEG3", "FN1", 
                                    "PDGFRA", "PLIN2"))                 ## Cluster 4

FeaturePlot(cluster11, features = c("LUM", "SPARCL1", "DCN", "POSTN", "MEG3", "FN1", 
                                    "PDGFRA", "PLIN2", "SOD2"))         ## Cluster 5

FeaturePlot(cluster11, features = c("ITGAX", "PTPRC", "POSTN", "SLC25A37", "FCER1G", 
                                    "S100A9", "FCGR3A", "MCEMP1", "ITGAM"))      ## Cluster 6

FeaturePlot(cluster11, features = c("SOD2", "CCL2", "ICAM1", "PTGDS", "FGF7", "DCN", 
                                    "PDGFRA", "PLIN2"))                 ## Cluster 7

DotPlot(cluster11, features = c("LCK", "KRT17", "S100A8", "S100A9", "PDGFRA", "SOD2", 
                                "LYZ", "S100A12", "EGFR", "PDGFRB", "EPAS1", "AXL", 
                                "CSPG4", "ITGB1", "ITM2C", "SPARCL1", "BCL2L1", "TOP2A",
                                "ACTA2", "WNT5A",  "POSTN", "AXIN2", "ELN", 
                                "CXCL14", "LGR6", "CTNNB1", "CD8B", "MMP10", "VPREB3", 
                                "LMAN1", "SLC1A3", "MCEMP1", "LTF", "SPP1",
                                "LUM", "DCN", "MEG3", "FN1", "PLIN2", "ITGAX", "PTPRC", 
                                "SLC25A37", "FCER1G", "FCGR3A", "ITGAM", "CCL2", 
                                "ICAM1", "PTGDS", "FGF7")) + angle

# Cluster 0 - Monocytes
# Cluster 1 - Pericytes
# Cluster 2 - SMC
# Cluster 3 - T
# Cluster 4 - FB
# Cluster 5 - Fetal SOD2+ FB
# Cluster 6 - Monocytes
# Cluster 7 - Fetal SOD2+ FB

## Cluster 12

cluster12 <- readRDS("/scratch/smallapragada/rapids_pipeline/run2_mes_reclustering_c12_2024_08_15.rds")

## UMAP visualization of clusters in reclustered mesenchymal object
DimPlot(cluster12,
        reduction = "umap",
        group.by = "mes_c12_brl_0.3", 
        raster = T,
        label = T, 
        cols = randomcoloR::distinctColorPalette(20)) 

## UMAP visualization of clusters in reclustered mesenchymal object
DimPlot(cluster12,
        reduction = "umap",
        group.by = "Phenotype", 
        raster = T,
        label = F, 
        cols = randomcoloR::distinctColorPalette(5)) 

## A summary of how many cells per cluster belong to a specific disease phenotype
table(cluster12$mes_c12_brl_0.3, cluster12$Phenotype)

## Looking at the nCount_RNA and the nFeature_RNA distribution across clusters to determine low-quality clusters
ggplot(cluster12@meta.data, aes(factor(cluster12$mes_c12_brl_0.3), (cluster12$nCount_RNA))) + geom_violin()
ggplot(cluster12@meta.data, aes(factor(cluster12$mes_c12_brl_0.3), (cluster12$nFeature_RNA))) + geom_violin()

DotPlot(cluster12, group.by = "mes_c11_brl_0.4", features = c(other_features, mesenchymal_features)) + angle

## Top markers
Idents(cluster12) = "mes_c12_brl_0.3"
markers <- FindAllMarkers(cluster12)
markers %>%
  mutate(pct.diff = abs(pct.1 - pct.2)) %>%
  filter(avg_log2FC >= 0) %>%
  View()

FeaturePlot(cluster12, features = c("POSTN", "SPARCL1", "HAS2", "ACTA2", "LYZ", 
                                    "CD68", "ITGAX", "PTGDS", "CCL18", "SPP1", 
                                    "PDGFRA", "FCER1G", "SOD2", "AXIN2", "GPR183",
                                    "STAT1"))           # Clusters 0 and 1

DotPlot(cluster12, features = c("POSTN", "SPARCL1", "HAS2", "ACTA2", "LYZ", 
                                "CD68", "ITGAX", "PTGDS", "CCL18", "SPP1", 
                                "PDGFRA", "FCER1G", "SOD2", "AXIN2", "GPR183",
                                "STAT1")) + angle

# Cluster 0 - FB
# Cluster 1 - Fetal-enriched SOD2+ myoFB

## Cluster 14 and 15

cluster14_15 <- readRDS("/scratch/smallapragada/rapids_pipeline/run2_mes_reclustering_c14_c15_2024_08_15.rds")

## UMAP visualization of clusters in reclustered mesenchymal object
DimPlot(cluster14_15,
        reduction = "umap",
        group.by = "mes_c14_15_brl_0.3", 
        raster = T,
        label = T, 
        cols = randomcoloR::distinctColorPalette(20)) 

## UMAP visualization of phenotype in reclustered mesenchymal object
DimPlot(cluster14_15,
        reduction = "umap",
        group.by = "Phenotype", 
        raster = T,
        label = F, 
        cols = randomcoloR::distinctColorPalette(5)) 

## A summary of how many cells per cluster belong to a specific disease phenotype
table(cluster14_15$mes_c14_15_brl_0.3, cluster14_15$Phenotype)

## Looking at the nCount_RNA and the nFeature_RNA distribution across clusters to determine low-quality clusters
ggplot(cluster14_15@meta.data, aes(factor(cluster14_15$mes_c14_15_brl_0.3), (cluster14_15$nCount_RNA))) + geom_violin()
ggplot(cluster14_15@meta.data, aes(factor(cluster14_15$mes_c14_15_brl_0.3), (cluster14_15$nFeature_RNA))) + geom_violin()

DotPlot(cluster14_15, group.by = "mes_c14_15_brl_0.3", features = c(other_features, mesenchymal_features, "SOD2")) + angle

## Top markers
Idents(cluster14_15) = "mes_c14_15_brl_0.3"
markers <- FindAllMarkers(cluster14_15)
markers %>%
  mutate(pct.diff = abs(pct.1 - pct.2)) %>%
  filter(avg_log2FC >= 0) %>%
  View()

FeaturePlot(cluster14_15, features = c("CCN2", "STAT1", "CXCL9", "CXCL14", "ZEB1", "CTNNB1", 
                                       "IFIT3", "CST3", "SNAI2", "IRF1"))        ## Cluster 0 

FeaturePlot(cluster14_15, features = c("HLA-DRA", "FCER1G", "AIF1", "MS4A7", "SPP1", "CD68", "FCGR3A", "GPR183", "MRC1", 
                                       "S100A9", "CD14", "C1QC"))                ## Cluster 1

FeaturePlot(cluster14_15, features = c("CCL21", "GNG11", "KDR", "PECAM1", "CCL2", "MRC1", 
                                       "MYC", "CTHRC1", "SFRP2"))                ## Cluster 2

FeaturePlot(cluster14_15, features = c("SPARCL1", "SCGB3A2", "PTGDS", "ACTA2", "SFTPC", 
                                       "SFTPD", "SCGB1A1", "KRT8", "COL15A1"))   ## Cluster 3

# Plot to spatially visualize where certain cells are located
samp <- subset(cluster14_15, subset = Sample == "s1_PDL002")

DimPlot(samp,
        reduction = "SP",
        group.by = "mes_c14_15_brl_0.3", 
        raster = T,
        label = F,
        cols = c(`3` = "green"), na.value = "red") + coord_equal()

FeaturePlot(cluster14_15, features = c("PGC", "LAMP3", "SFTPC", "KRT18", "RNASE1", "KRT8", 
                                       "DMBT1", "MMP7", "EPCAM"))                ## Cluster 4

FeaturePlot(cluster14_15, features = c("TCL1A", "CD79A", "S100A9", "FCER1G", "ATG7", "LMAN1",
                                       "PKM", "GDF15", "CD8A"))                  ## Cluster 5

FeaturePlot(cluster14_15, features = c("PDGFRA", "LUM", "DCN", "SOD2", "AXIN2", "YAP1", 
                                       "CD44", "HIF1A", "HAS2"))                 ## Cluster 6

## Adding cell type labels to each of the reclustered objects
cluster11$CT <- ""
cluster11$CT[cluster11$mes_c11_brl_0.4 == 0] <- "Monocytes"
cluster11$CT[cluster11$mes_c11_brl_0.4 == 1] <- "Pericytes"
cluster11$CT[cluster11$mes_c11_brl_0.4 == 2] <- "SMC"
cluster11$CT[cluster11$mes_c11_brl_0.4 == 3] <- "T"
cluster11$CT[cluster11$mes_c11_brl_0.4 == 4] <- "Fibroblasts"
cluster11$CT[cluster11$mes_c11_brl_0.4 == 5] <- "Fetal-enriched SOD2+ myoFB"
cluster11$CT[cluster11$mes_c11_brl_0.4 == 6] <- "Monocytes"
cluster11$CT[cluster11$mes_c11_brl_0.4 == 7] <- "Fetal-enriched SOD2+ myoFB"

cluster12$CT <- ""
cluster12$CT[cluster12$mes_c12_brl_0.3 == 0] <- "Fibroblasts"
cluster12$CT[cluster12$mes_c12_brl_0.3 == 1] <- "Fetal-enriched SOD2+ myoFB"

cluster14_15$CT <- ""
cluster14_15$CT[cluster14_15$mes_c14_c15_brl_0.3 == 0] <- "Fetal-enriched fibroblasts"
cluster14_15$CT[cluster14_15$mes_c14_c15_brl_0.3 == 1] <- "pDC"
cluster14_15$CT[cluster14_15$mes_c14_c15_brl_0.3 == 2] <- "Lymphatic endothelial"
cluster14_15$CT[cluster14_15$mes_c14_c15_brl_0.3 == 3] <- "SMC"
cluster14_15$CT[cluster14_15$mes_c14_c15_brl_0.3 == 4] <- "Fetal-enriched transitional AT2"
cluster14_15$CT[cluster14_15$mes_c14_c15_brl_0.3 == 5] <- "T"
cluster14_15$CT[cluster14_15$mes_c14_c15_brl_0.3 == 6] <- "Fetal-enriched SOD2+ myoFB"

## Combining objects together and merging into one
obj_list <- list(cluster11, cluster12, cluster14_15)
firstpass <- merge(x = obj_list[[1]], y = obj_list[2:length(obj_list)])

## Adding labels from the merged reclustered object to the whole mesenchymal object + all other CT labels 
mes_subclustering$CT <- ""
mes_subclustering$CT <- firstpass$CT
mes_subclustering$CT[mes_subclustering$mes_sc_0.3_c0_0.5 == "0,0"] <- "Activated fibroblasts"
mes_subclustering$CT[mes_subclustering$mes_sc_0.3_c0_0.5 == "0,1"] <- "Low quality"
mes_subclustering$CT[mes_subclustering$mes_sc_0.3_c0_0.5 == "0,2"] <- "Activated fibroblasts"
mes_subclustering$CT[mes_subclustering$mes_sc_0.3_c0_0.5 == "0,3"] <- "Activated fibroblasts"
mes_subclustering$CT[mes_subclustering$mes_sc_0.3_c0_0.5 == "0,4"] <- "Fibroblasts" 
mes_subclustering$CT[mes_subclustering$mes_sc_0.3_c0_0.5 == "0,5"] <- "T" 
mes_subclustering$CT[mes_subclustering$mes_brl_0.3 == 1] <- "Adventitial fibroblasts"
mes_subclustering$CT[mes_subclustering$mes_sc_0.3_c2_c13_0.2 == "2-13,0"] <- "Pericytes"
mes_subclustering$CT[mes_subclustering$mes_sc_0.3_c2_c13_0.2 == "2-13,1"] <- "Fetal-enriched SOD2+ myoFB"
mes_subclustering$CT[mes_subclustering$mes_sc_0.3_c2_c13_0.2 == "2-13,2"] <- "Pericytes"
mes_subclustering$CT[mes_subclustering$mes_brl_0.3 == 3] <- "Fetal-enriched SOD2+ myoFB"
mes_subclustering$CT[mes_subclustering$mes_brl_0.3 == 4] <- "Adventitial fibroblasts"
mes_subclustering$CT[mes_subclustering$mes_brl_0.3 == 5] <- "AKR1C1+ & AKR1C2+"
mes_subclustering$CT[mes_subclustering$mes_brl_0.3 == 6] <- "Low quality"
mes_subclustering$CT[mes_subclustering$mes_sc_0.3_c7_0.2 == "7,0"] <- "Fetal-enriched SOD2+ myoFB"
mes_subclustering$CT[mes_subclustering$mes_sc_0.3_c7_0.2 == "7,1"] <- "AKR1C1+ & AKR1C2+"
mes_subclustering$CT[mes_subclustering$mes_sc_0.3_c7_0.2 == "7,2"] <- "Fetal-enriched SOD2+ myoFB"
mes_subclustering$CT[mes_subclustering$mes_brl_0.3 == 8] <- "PLIN2+ fibroblasts"
mes_subclustering$CT[mes_subclustering$mes_brl_0.3 == 9] <- "SMC"
mes_subclustering$CT[mes_subclustering$mes_brl_0.3 == 10] <- "Adventitial fibroblasts"

## Save mesenchymal object
saveRDS(mes_subclustering, "/scratch/smallapragada/run2/final_mes_CT_2024_08_15.rds")

#### CELL TYPE ANNOTATION - LYMPHOID ----

## Reading in lymphoid object
lymph <- readRDS("/scratch/smallapragada/rapids_pipeline/run2_lymph_2024_08_15.rds")

## Looking at the nCount_RNA and the nFeature_RNA distribution across clusters to determine low-quality clusters
ggplot(lymph@meta.data, aes(factor(lymph@meta.data$lymph_brl_0.3), (lymph@meta.data$nCount_RNA))) + geom_violin()
ggplot(lymph@meta.data, aes(factor(lymph@meta.data$lymph_brl_0.3), (lymph@meta.data$nFeature_RNA))) + geom_violin()

## UMAP visualization of clusters in lymphoid object
DimPlot(lymph,
        reduction = "umap",
        group.by = "lymph_brl_0.2", 
        raster = T,
        label = T, 
        cols = randomcoloR::distinctColorPalette(20)) 

## UMAP visualization of phenotype in lymphoid object
DimPlot(lymph,
        reduction = "umap",
        group.by = "Phenotype", 
        raster = T,
        label = F, 
        cols = randomcoloR::distinctColorPalette(5)) 

## A summary of how many cells per cluster belong to a specific disease phenotype
table(lymph$lymph_brl_0.2, lymph$Phenotype)

DotPlot(lymph, group.by = "lymph_brl_0.2", features = c(other_features, immune_features)) + angle

## Top markers for lymphoid
Idents(lymph) = "lymph_brl_0.2"
lymph_markers <- FindAllMarkers(lymph)
lymph_markers %>%
  mutate(pct.diff = abs(pct.1 - pct.2)) %>%
  filter(avg_log2FC >= 0) %>%
  View()

## Cluster 0 - SMC/PDGFRA+ FB (Subcluster)

# Temp object
obj <- subset(lymph, subset = lymph_brl_0.2 == 0)

FeaturePlot(obj, features = c("COL1A2", "COL3A1", "MEG3", "LUM", "FN1", "COL1A1", 
                              "DCN", "PDGFRA", "POSTN", "PDGFRB", "SNAI2", "HAS2", 
                              "ACTA2", "WNT2", "CSPG4", "FGF10", "CD4", "PTGDS"))

DotPlot(lymph, features = c("COL1A2", "COL3A1", "MEG3", "LUM", "FN1", "COL1A1", 
                            "DCN", "PDGFRA", "POSTN", "PDGFRB", "SNAI2", "HAS2", 
                            "ACTA2", "WNT2", "CSPG4", "FGF10", "CD4", "PTGDS")) + angle

## Cluster 1 - Endothelial/Low quality (Subcluster)

# Temp object
obj <- subset(lymph, subset = lymph_brl_0.2 == 1)

FeaturePlot(obj, features = c("PECAM1", "CD34", "FCN3", "KDR", "SFTPD", "CA4", 
                              "LYZ", "GNG11", "RAMP2", "CLDN5", "RTKN2", "HEY1", 
                              "CEACAM6", "APLNR", "NOX4", "CCL18", "MS4A7"))

DotPlot(lymph, features = c("PECAM1", "CD34", "FCN3", "KDR", "SFTPD", "CA4", "LYZ", 
                            "GNG11", "RAMP2", "CLDN5", "RTKN2", "HEY1", "CEACAM6", 
                            "APLNR", "NOX4", "CCL18", "MS4A7")) + angle

## Cluster 2 - Fetal-enriched transitional AT2 (Subcluster)

# Temp object
obj <- subset(lymph, subset = lymph_brl_0.2 == 2)

FeaturePlot(obj, features = c("NKX2-1", "MGST1", "LAMP3", "KRT18", "RNASE1", "EPCAM", 
                              "KRT8", "PGC", "MMP7", "HYOU1", "PDIA4", "HSPA5", 
                              "SFTPC", "DMBT1", "XBP1", "NAPSA", "PDIA3", "SOD2", 
                              "SEC11C", "HMGA1"))

DotPlot(lymph, features = c(other_features, epithelial_features, "SOD2", "SEC11C")) + angle

## Cluster 3 - B/Plasma/Monocytes (Subcluster)

# Temp object
obj <- subset(lymph, subset = lymph_brl_0.2 == 3)

FeaturePlot(obj, features = c("MS4A1", "BANK1", "TNFRSF13C", "PTPRC", "CD79A", 
                              "CD19", "TCL1A", "ATP2A3", "CXCR4", "JCHAIN", 
                              "CD52", "CD79B", "LTB", "S100A8", "S100A9"))

DotPlot(lymph, features = c("MS4A1", "BANK1", "TNFRSF13C", "PTPRC", "CD79A", "CD19", 
                            "TCL1A", "ATP2A3", "CXCR4", "JCHAIN", "CD52", "CD79B", 
                            "LTB")) + angle

## Cluster 4 - T/NK & NKT/Treg (Subcluster)

# Temp object
obj <- subset(lymph, subset = lymph_brl_0.2 == 4)

FeaturePlot(obj, features = c("CD247", "IL7R", "TRAC", "KLRG1", "KLRB1", "CD69", 
                              "CXCR4", "LCK", "CD3E", "GZMA", "TGFB1", "NKG7", 
                              "CD3D", "KLRC1", "CD3G", "CD8A", "CTLA4", "FOXP3")) 

DotPlot(lymph, features = c("CD247", "IL7R", "TRAC", "KLRG1", "KLRB1", "CD69", 
                            "CXCR4", "LCK", "CD3E", "GZMA", "TGFB1", "NKG7", "CD3D", 
                            "KLRC1", "CD3G", "CD8A", "CTLA4", "FOXP3")) + angle

### Subclustering clusters 0, 1, 2, 3, and 4

## Cluster 0

lymph_subclustering <- readRDS("/scratch/smallapragada/rapids_pipeline/run2_lymph_subclustering_2024_08_15.rds")

## UMAP visualization of clusters in subclustered lymphoid object
DimPlot(lymph_subclustering,
        reduction = "umap",
        group.by = "lymph_sc_0.2_c0_0.1", 
        raster = T,
        label = T, 
        cols = randomcoloR::distinctColorPalette(20)) 

## UMAP visualization of phenotypes in subclustered lymphoid object
DimPlot(lymph_subclustering,
        reduction = "umap",
        group.by = "Phenotype", 
        raster = T,
        label = F, 
        cols = randomcoloR::distinctColorPalette(5)) 

## A summary of how many cells per cluster belong to a specific disease phenotype
table(lymph_subclustering$lymph_sc_0.2_c0_0.2, lymph_subclustering$Phenotype)

## Looking at the nCount_RNA and the nFeature_RNA distribution across clusters to determine low-quality clusters
ggplot(lymph_subclustering@meta.data, aes(factor(lymph_subclustering@meta.data$lymph_sc_0.2_c0_0.2), 
                                          (lymph_subclustering@meta.data$nCount_RNA))) + geom_violin()
ggplot(lymph_subclustering@meta.data, aes(factor(lymph_subclustering@meta.data$lymph_sc_0.2_c0_0.2), 
                                          (lymph_subclustering@meta.data$nFeature_RNA))) + geom_violin()

DotPlot(lymph_subclustering, group.by = "lymph_sc_0.2_c0_0.2", features = c(other_features, immune_features)) + angle

## Temp object to visualize clusters of interest
obj <- subset(lymph_subclustering, subset = lymph_sc_0.2_c0_0.2 == "0,0" | lymph_sc_0.2_c0_0.2 == "0,1" | 
                lymph_sc_0.2_c0_0.2 == "0,2" | lymph_sc_0.2_c0_0.2 == "0,3")

## Top markers for temp object
Idents(obj) = "lymph_sc_0.2_c0_0.2"
markers <- FindAllMarkers(obj)
markers %>%
  mutate(pct.diff = abs(pct.1 - pct.2)) %>%
  filter(avg_log2FC >= 0) %>%
  View()

FeaturePlot(obj, features = c("COL1A2", "COL3A1", "MEG3", "LUM", "FN1", "COL1A1", 
                              "DCN", "PDGFRA", "POSTN", "PDGFRB", "SNAI2", "HAS2", 
                              "ACTA2", "WNT2", "CSPG4", "FGF10", "CD4", "PTGDS"))

DotPlot(obj, features = c("COL1A2", "COL3A1", "MEG3", "LUM", "FN1", "COL1A1", 
                          "DCN", "PDGFRA", "POSTN", "PDGFRB", "SNAI2", "HAS2", 
                          "ACTA2", "WNT2", "CSPG4", "FGF10", "CD4", "PTGDS")) + angle

# Plot to spatially visualize where certain cells are located
samp <- subset(lymph_subclustering, subset = Sample == "s1_PDL002")

DimPlot(samp,
        reduction = "SP",
        group.by = "lymph_sc_0.2_c0_0.2", 
        raster = T,
        label = F,
        cols = c(`0,0` = "red", `0,1` = "green", `0,2` = "blue"), na.value = "grey90") + coord_equal()

# 0,0 - SMC
# 0,1 - SMC
# 0,2 - SMC
# 0,3 - Proliferating FB

## Cluster 1

## UMAP visualization of clusters in subclustered lymphoid object
DimPlot(lymph_subclustering,
        reduction = "umap",
        group.by = "lymph_sc_0.2_c1_0.1", 
        raster = T,
        label = T, 
        cols = randomcoloR::distinctColorPalette(20)) 

## UMAP visualization of phenotypes in subclustered lymphoid object
DimPlot(lymph_subclustering,
        reduction = "umap",
        group.by = "Phenotype", 
        raster = T,
        label = F, 
        cols = randomcoloR::distinctColorPalette(5)) 

## A summary of how many cells per cluster belong to a specific disease phenotype
table(lymph_subclustering$lymph_sc_0.2_c1_0.2, lymph_subclustering$Phenotype)

## Looking at the nCount_RNA and the nFeature_RNA distribution across clusters to determine low-quality clusters
ggplot(lymph_subclustering@meta.data, aes(factor(lymph_subclustering@meta.data$lymph_sc_0.2_c1_0.2), 
                                          (lymph_subclustering@meta.data$nCount_RNA))) + geom_violin()
ggplot(lymph_subclustering@meta.data, aes(factor(lymph_subclustering@meta.data$lymph_sc_0.2_c1_0.2), 
                                          (lymph_subclustering@meta.data$nFeature_RNA))) + geom_violin()

DotPlot(obj, group.by = "lymph_sc_0.2_c1_0.2", features = c(other_features, immune_features)) + angle

## Temp object to visualize clusters of interest
obj <- subset(lymph_subclustering, subset = lymph_sc_0.2_c1_0.2 == "1,0" | lymph_sc_0.2_c1_0.2 == "1,1" | 
                lymph_sc_0.2_c1_0.2 == "1,2" | lymph_sc_0.2_c1_0.2 == "1,3")

## Top markers
Idents(obj) = "lymph_sc_0.2_c1_0.2"
markers <- FindAllMarkers(obj)
markers %>%
  mutate(pct.diff = abs(pct.1 - pct.2)) %>%
  filter(avg_log2FC >= 0) %>%
  View()

FeaturePlot(obj, features = c("CTLA4", "FOXP3", "IL2RA", "S100A8", "S100A9", "S100A12", 
                              "SPP1", "MARCO", "FCER1G", "SNCA", "SLC7A11", "BCL2L11", 
                              "SFTPC", "SFTPD", "RTKN2", "CD34", "KDR", "CA4", "GNG11", 
                              "CLDN5", "BMPR2", "APLNR", "RAMP2", "ACKR1"))

DotPlot(obj, features = c("CTLA4", "FOXP3", "IL2RA", "S100A8", "S100A9", "S100A12", 
                          "SPP1", "MARCO", "FCER1G", "SNCA", "SLC7A11", "BCL2L11", 
                          "SFTPC", "SFTPD", "RTKN2", "CD34", "KDR", "CA4", "GNG11", 
                          "CLDN5", "BMPR2", "APLNR", "RAMP2", "ACKR1")) + angle

# 1,0 - Treg
# 1,1 - Monocytes
# 1,2 - AT2
# 1,3 - Capillaries

## Cluster 2

## UMAP visualization of clusters in subclustered lymphoid object
DimPlot(lymph_subclustering,
        reduction = "umap",
        group.by = "lymph_sc_0.2_c2_0.4", 
        raster = T,
        label = T, 
        cols = randomcoloR::distinctColorPalette(20)) 

## UMAP visualization of phenotypes in subclustered lymphoid object
DimPlot(lymph_subclustering,
        reduction = "umap",
        group.by = "Phenotype", 
        raster = T,
        label = F, 
        cols = randomcoloR::distinctColorPalette(5)) 

## A summary of how many cells per cluster belong to a specific disease phenotype
table(lymph_subclustering$lymph_sc_0.2_c2_0.4, lymph_subclustering$Phenotype)

## Looking at the nCount_RNA and the nFeature_RNA distribution across clusters to determine low-quality clusters
ggplot(lymph_subclustering@meta.data, aes(factor(lymph_subclustering@meta.data$lymph_sc_0.2_c2_0.4), 
                                          (lymph_subclustering@meta.data$nCount_RNA))) + geom_violin()
ggplot(lymph_subclustering@meta.data, aes(factor(lymph_subclustering@meta.data$lymph_sc_0.2_c2_0.4), 
                                          (lymph_subclustering@meta.data$nFeature_RNA))) + geom_violin()

## Temp object to visualize clusters of interest
obj <- subset(lymph_subclustering, subset = lymph_sc_0.2_c2_0.4 == "2,0" | lymph_sc_0.2_c2_0.4 == "2,1" | 
                lymph_sc_0.2_c2_0.4 == "2,2" | lymph_sc_0.2_c2_0.4 == "2,3" | lymph_sc_0.2_c2_0.4 == "2,4" |
                lymph_sc_0.2_c2_0.4 == "2,5" | lymph_sc_0.2_c2_0.4 == "2,6")

DotPlot(obj, group.by = "lymph_sc_0.2_c2_0.4", features = c(epithelial_features, immune_features)) + angle

## Top markers
Idents(obj) = "lymph_sc_0.2_c2_0.4"
lymph_markers <- FindAllMarkers(obj)
lymph_markers %>%
  mutate(pct.diff = abs(pct.1 - pct.2)) %>%
  filter(avg_log2FC >= 0) %>%
  View()

FeaturePlot(obj, features = c("DMBT1", "DUOX1", "SFTPC", "S100A8", 
                              "S100A9", "S100A12", "COL3A1", "COL1A2", "PDGFRA", "DCN",
                              "ICAM1", "MSLN", "HMGA1", "SCGB3A2", "AKR1C1", "KDR", "CA4",
                              "CLDN5", "TOP2A", "MKI67"))

DotPlot(obj, features = c("DMBT1", "DUOX1", "SFTPC", "SAA2", "LTF", "ERLEC1", "SPP1", "S100A8", 
                          "S100A9", "S100A12", "LYZ", "FCER1G", "COL3A1", "COL1A2", "PDGFRA", "DCN", "LUM",
                          "PDGFRB", "HAS2", "SNAI2", "FGF7", "ICAM1", "MSLN", "HMGA1", "SCGB3A2", "AKR1C1", "KDR", "CA4",
                          "ACKR1", "CD34", "GNG11", "KIT", "CLDN5", "TOP2A", "MKI67")) + angle

# 2,0 - Fetal-enriched AT2
# 2,1 - Fetal-enriched AT2
# 2,2 - Fetal FB
# 2,3 - Low quality
# 2,4 - Fetal trans AT2
# 2,5 - Fetal-enriched AT2
# 2,6 - Fetal proliferating at2

## Cluster 3

## UMAP visualization of clusters in subclustered lymphoid object
DimPlot(lymph_subclustering,
        reduction = "umap",
        group.by = "lymph_sc_0.2_c3_0.2", 
        raster = T,
        label = T, 
        cols = randomcoloR::distinctColorPalette(20)) 

## UMAP visualization of phenotypes in subclustered lymphoid object
DimPlot(lymph_subclustering,
        reduction = "umap",
        group.by = "Phenotype", 
        raster = T,
        label = F, 
        cols = randomcoloR::distinctColorPalette(5)) 

## A summary of how many cells per cluster belong to a specific disease phenotype
table(lymph_subclustering$lymph_sc_0.2_c3_0.2, lymph_subclustering$Phenotype)

## Looking at the nCount_RNA and the nFeature_RNA distribution across clusters to determine low-quality clusters
ggplot(lymph_subclustering@meta.data, aes(factor(lymph_subclustering@meta.data$lymph_sc_0.2_c3_0.2), 
                                          (lymph_subclustering@meta.data$nCount_RNA))) + geom_violin()
ggplot(lymph_subclustering@meta.data, aes(factor(lymph_subclustering@meta.data$lymph_sc_0.2_c3_0.2), 
                                          (lymph_subclustering@meta.data$nFeature_RNA))) + geom_violin()

## Temp object to visualize clusters of interest
obj <- subset(lymph_subclustering, subset = lymph_sc_0.2_c3_0.2 == "3,0" | lymph_sc_0.2_c3_0.2 == "3,1" | 
                lymph_sc_0.2_c3_0.2 == "3,2" | lymph_sc_0.2_c3_0.2 == "3,3")

DotPlot(obj, group.by = "lymph_sc_0.2_c3_0.2", features = c(other_features, immune_features)) + angle

## Top markers
Idents(obj) = "lymph_sc_0.2_c3_0.2"
lymph_markers <- FindAllMarkers(obj)
lymph_markers %>%
  mutate(pct.diff = abs(pct.1 - pct.2)) %>%
  filter(avg_log2FC >= 0) %>%
  View()

FeaturePlot(obj, features = c("MEG3", "APLNR", "CA4", "CD247", "APLN", "MS4A7", "BANK1", "TNFRSF13C", 
                              "CD79A", "CD19", "TCL1A", "CD52", "S100A8", "S100A9", "S100A12", "LYZ", "FCER1G",
                              "ACKR1", "LTF", "MCEMP1", "RETN", "JCHAIN", "XBP1", "SELENOS"))

DotPlot(obj, features = c("MEG3", "APLNR", "CA4", "CD247", "APLN", "MS4A7", "BANK1", "TNFRSF13C", 
                          "CD79A", "CD19", "TCL1A", "CD52", "S100A8", "S100A9", "S100A12", "LYZ", "FCER1G",
                          "ACKR1", "LTF", "MCEMP1", "RETN", "JCHAIN", "XBP1", "SELENOS")) + angle

# 3,0 - B
# 3,1 - B
# 3,2 - Monocytes
# 3,3 - Plasma

## Cluster 4

## UMAP visualization of clusters in subclustered lymphoid object
DimPlot(lymph_subclustering,
        reduction = "umap",
        group.by = "lymph_sc_0.2_c4_0.4", 
        raster = T,
        label = T, 
        cols = randomcoloR::distinctColorPalette(20)) 

## UMAP visualization of phenotypes in subclustered lymphoid object
DimPlot(lymph_subclustering,
        reduction = "umap",
        group.by = "Phenotype", 
        raster = T,
        label = F, 
        cols = randomcoloR::distinctColorPalette(5)) 

## A summary of how many cells per cluster belong to a specific disease phenotype
table(lymph_subclustering$lymph_sc_0.2_c4_0.4, lymph_subclustering$Phenotype)

## Looking at the nCount_RNA and the nFeature_RNA distribution across clusters to determine low-quality clusters
ggplot(lymph_subclustering@meta.data, aes(factor(lymph_subclustering@meta.data$lymph_sc_0.2_c4_0.4), 
                                          (lymph_subclustering@meta.data$nCount_RNA))) + geom_violin()
ggplot(lymph_subclustering@meta.data, aes(factor(lymph_subclustering@meta.data$lymph_sc_0.2_c4_0.4), 
                                          (lymph_subclustering@meta.data$nFeature_RNA))) + geom_violin()

## Temp object to visualize clusters of interest
obj <- subset(lymph_subclustering, subset = lymph_sc_0.2_c4_0.4 == "4,0" | lymph_sc_0.2_c4_0.4 == "4,1" | 
                lymph_sc_0.2_c4_0.4 == "4,2" | lymph_sc_0.2_c4_0.4 == "4,3" | lymph_sc_0.2_c4_0.4 == "4,4" | 
                lymph_sc_0.2_c4_0.4 == "4,5")

DotPlot(obj, group.by = "lymph_sc_0.2_c4_0.4", features = c(other_features, immune_features, "CTLA4", "FOXP3")) + angle

## Top markers
Idents(obj) = "lymph_sc_0.2_c4_0.4"
lymph_markers <- FindAllMarkers(obj)
lymph_markers %>%
  mutate(pct.diff = abs(pct.1 - pct.2)) %>%
  filter(avg_log2FC >= 0) %>%
  View()

FeaturePlot(obj, features = c("KLRC1", "GZMA", "NKG7", "KLRB1", "CD247", "GNLY", "CTLA4", "FOXP3", 
                              "IL2RA", "LEF1", "TRAC", "CCR7", "IL7R", "CD3E", "CD3D", "CD8A", "FCER1G",
                              "FCGR3A", "HAVCR2", "GPR183", "CD69", "LUM", "DCN", "COL3A1", "MKI67", "TOP2A"))

DotPlot(obj, features = c("KLRC1", "GZMA", "NKG7", "KLRB1", "CD247", "GNLY", "CTLA4", "FOXP3", 
                          "IL2RA", "LEF1", "TRAC", "CCR7", "IL7R", "CD3E", "CD3D", "CD8A", "FCER1G",
                          "FCGR3A", "HAVCR2", "GPR183", "CD69", "LUM", "DCN", "COL3A1", "MKI67", "TOP2A")) + angle

## Adding all CT labels
lymph_subclustering$CT <- ""
lymph_subclustering$CT[lymph_subclustering$lymph_sc_0.2_c0_0.2 == "0,0" | lymph_subclustering$lymph_sc_0.2_c0_0.2 == "0,1" |
                         lymph_subclustering$lymph_sc_0.2_c0_0.2 == "0,2"] <- "SMC"
lymph_subclustering$CT[lymph_subclustering$lymph_sc_0.2_c0_0.2 == "0,3"] <- "Proliferating fibroblasts"
lymph_subclustering$CT[lymph_subclustering$lymph_sc_0.2_c1_0.2 == "1,0" | lymph_subclustering$lymph_sc_0.2_c4_0.4 == "4,1"] <- "Treg"
lymph_subclustering$CT[lymph_subclustering$lymph_sc_0.2_c1_0.2 == "1,1" | lymph_subclustering$lymph_sc_0.2_c3_0.2 == "3,2"] <- "Monocytes"
lymph_subclustering$CT[lymph_subclustering$lymph_sc_0.2_c1_0.2 == "1,2"] <- "AT2"
lymph_subclustering$CT[lymph_subclustering$lymph_sc_0.2_c1_0.2 == "1,3"] <- "Capillaries"
lymph_subclustering$CT[lymph_subclustering$lymph_sc_0.2_c2_0.4 == "2,0"] <- "Fetal-enriched AT2"
lymph_subclustering$CT[lymph_subclustering$lymph_sc_0.2_c2_0.4 == "2,1"] <- "Fetal-enriched AT2"
lymph_subclustering$CT[lymph_subclustering$lymph_sc_0.2_c2_0.4 == "2,2"] <- "Fetal-enriched fibroblasts"
lymph_subclustering$CT[lymph_subclustering$lymph_sc_0.2_c2_0.4 == "2,3"] <- "Low quality"
lymph_subclustering$CT[lymph_subclustering$lymph_sc_0.2_c2_0.4 == "2,4"] <- "Fetal-enriched transitional AT2"
lymph_subclustering$CT[lymph_subclustering$lymph_sc_0.2_c2_0.4 == "2,5"] <- "Fetal-enriched AT2"
lymph_subclustering$CT[lymph_subclustering$lymph_sc_0.2_c2_0.4 == "2,6"] <- "Proliferating Fetal-enriched AT2"
lymph_subclustering$CT[lymph_subclustering$lymph_sc_0.2_c3_0.2 == "3,0" | lymph_subclustering$lymph_sc_0.2_c3_0.2 == "3,1"] <- "B"
lymph_subclustering$CT[lymph_subclustering$lymph_sc_0.2_c3_0.2 == "3,3"] <- "Plasma"
lymph_subclustering$CT[lymph_subclustering$lymph_sc_0.2_c4_0.4 == "4,0"] <- "NK & NKT"
lymph_subclustering$CT[lymph_subclustering$lymph_sc_0.2_c4_0.4 == "4,2"] <- "T"
lymph_subclustering$CT[lymph_subclustering$lymph_sc_0.2_c4_0.4 == "4,3"] <- "cDC"
lymph_subclustering$CT[lymph_subclustering$lymph_sc_0.2_c4_0.4 == "4,4"] <- "T"
lymph_subclustering$CT[lymph_subclustering$lymph_sc_0.2_c4_0.4 == "4,5"] <- "Proliferating lymphoid"

## Save lymphoid object
saveRDS(lymph_subclustering, "/scratch/smallapragada/run2/final_lymph_CT_2024_08_15.rds")

#### CELL TYPE ANNOTATION - MYELOID ----

## Reading in myeloid object
mye <- readRDS("/scratch/smallapragada/rapids_pipeline/run2_mye_2024_08_15.rds")

## Looking at the nCount_RNA and the nFeature_RNA distribution across clusters to determine low-quality clusters
ggplot(mye@meta.data, aes(factor(mye@meta.data$mye_brl_0.4), (mye@meta.data$nCount_RNA))) + geom_violin()
ggplot(mye@meta.data, aes(factor(mye@meta.data$mye_brl_0.4), (mye@meta.data$nFeature_RNA))) + geom_violin()

## UMAP visualization of clusters in myeloid object
DimPlot(mye,
        reduction = "umap",
        group.by = "mye_brl_0.4", 
        raster = T,
        label = T, 
        cols = randomcoloR::distinctColorPalette(20)) 

## UMAP visualization of phenotype in myeloid object
DimPlot(mye,
        reduction = "umap",
        group.by = "Phenotype", 
        raster = T,
        label = F, 
        cols = randomcoloR::distinctColorPalette(5)) 

## A summary of how many cells per cluster belong to a specific disease phenotype
table(mye$mye_brl_0.4, mye$Phenotype)

## Top markers for myeloid
Idents(mye) = "mye_brl_0.4"
mye_markers <- FindAllMarkers(mye)
mye_markers %>%
  mutate(pct.diff = abs(pct.1 - pct.2)) %>%
  filter(avg_log2FC >= 0) %>%
  View()

# Myeloid dotplot
DotPlot(mye, group.by = "mye_brl_0.4", features = c(other_features, immune_features)) + angle

## Cluster 0 - Monocytes

# Temp object
obj <- subset(mye, subset = mye_brl_0.4 == 0)

FeaturePlot(obj, features = c("S100A8", "S100A9", "S100A12", "FCER1G", "SLC25A37", "LYZ", "MCEMP1", "SOD2", "LTF", "ITGAM", 
                              "ITGAX", "RETN", "CXCR4", "RHOA"))

DotPlot(mye, features = c("S100A8", "S100A9", "S100A12", "FCER1G", "SLC25A37", "LYZ", "MCEMP1", "SOD2", "LTF", "ITGAM", 
                          "ITGAX", "RETN", "CXCR4", "RHOA", "SPP1")) + angle

## Cluster 1 - Proliferating Meg-Ery

# Temp object
obj <- subset(mye, subset = mye_brl_0.4 == 1)

FeaturePlot(obj, features = c("SLC2A1", "SLC25A37", "BCL2L1", "SNCA", "MKI67", 
                              "TOP2A", "HMGA1", "PCNA", "GCLM", "TOP1", "ATF4"))

DotPlot(mye, features = c("SLC2A1", "SLC25A37", "BCL2L1", "SNCA", "MKI67", "TOP2A", 
                          "HMGA1", "PCNA", "GCLM", "TOP1", "ATF4")) + angle

## Cluster 2 - FB

# Temp object
obj <- subset(mye_subclustering, subset = mye_brl_0.4 == 2)

FeaturePlot(obj, features = c("MEG3", "COL1A2", "SNAI2", "LUM", "FN1", "COL3A1", 
                              "COL1A1", "YAP1", "POSTN", "DCN", "WNT2", "HAS2", 
                              "VEGFA", "PLIN2", "PDGFRA", "PDGFRB", "SOD2"))

DotPlot(mye, features = c(mesenchymal_features)) + angle

## Cluster 3 - Transitonal AT2/Mast (recluster)

# Temp object
obj <- subset(mye, subset = mye_brl_0.4 == 3)

FeaturePlot(obj, features = c("CPA3", "TPSAB1", "KIT", "KLRG1", "LAMP3", "NAPSA", 
                              "SFTPC", "KRT8", "NKX2-1", "SFTPD", "PTGDS"))

DotPlot(mye, features = c("CPA3", "TPSAB1", "KIT", "KLRG1", "LAMP3", "NAPSA", 
                          "SFTPC", "KRT8", "NKX2-1", "SFTPD", "PTGDS")) + angle

## Cluster 4 - Alveolar macrophages (Subcluster)

# Temp object
obj <- subset(mye, subset = mye_brl_0.4 == 4)

FeaturePlot(obj, features = c("HLA-DRA", "CD68", "MS4A7", "LYZ", "MRC1", "PLIN2", 
                              "CCL18", "ITGAX", "FCER1G", "HLA-DQA1", "CD44", "PPARG", 
                              "MARCO", "FCGR3A", "C1QC", "FABP4"))

DotPlot(mye, features = c("HLA-DRA", "CD68", "MS4A7", "LYZ", "MRC1", "CCL18", "ITGAX", 
                          "FCER1G", "HLA-DQA1", "CD44", "PPARG", "MARCO", "FCGR3A", 
                          "C1QC", "FABP4", "PLIN2")) + angle

## Cluster 5 - Proliferating endothelial/low quality (Subcluster)

# Temp object
obj <- subset(mye, subset = mye_brl_0.4 == 5)

FeaturePlot(obj, features = c("TOP2A", "MKI67", "EPAS1", "CD34", "KDR", "RACGAP1", "CCNB2", "APLNR", "KIT", "GNG11", 
                              "WWTR1", "ZEB1"))

DotPlot(mye, features = c(other_features, endothelial_features)) + angle

## Cluster 6 - pDC/monocytes

# Temp object
obj <- subset(mye, subset = mye_brl_0.4 == 6)

FeaturePlot(obj, features = c("FCER1G", "AIF1", "PTPRC", "LYZ", "MS4A7", "CD68", 
                              "CD14", "CD44", "CST3", "S100A9", "S100A8", "FCGR3A", 
                              "FCN1", "GPR183", "MRC1", "XBP1"))

DotPlot(mye, features = c("FCER1G", "AIF1", "PTPRC", "LYZ", "MS4A7", "CD68", 
                          "CD14", "CD44", "CST3", "S100A9", "S100A8", "FCGR3A", 
                          "FCN1", "GPR183", "MRC1", "XBP1")) + angle

## Cluster 7 - Proliferating FB

# Temp object
obj <- subset(mye, subset = mye_brl_0.4 == 7)

FeaturePlot(obj, features = c("TOP2A", "MKI67", "LUM", "DCN", "COL3A1", "RTKN2", "POSTN", "MEG3", "PCNA", "YAP1", 
                              "FN1", "PDGFRA"))

DotPlot(mye, features = c("TOP2A", "MKI67", "LUM", "DCN", "COL3A1", "RTKN2", "POSTN", "MEG3", "PCNA", "YAP1", 
                          "FN1", "PDGFRA")) + angle

## Cluster 8 - Unsure (recluster)

# Temp object
obj <- subset(mye, subset = mye_brl_0.4 == 8)

FeaturePlot(obj, features = c("CCN2", "TNFRSF9", "CCR7", "LAMP3", "NFKB1", "GPR183", 
                              "IL7R", "CD86", "BCL2L11", "SELENOS", "CXCR4"))

DotPlot(mye, features = c("CCN2", "TNFRSF9", "CCR7", "LAMP3", "NFKB1", "GPR183", 
                          "IL7R", "CD86", "BCL2L11", "SELENOS", "CXCR4")) + angle

### Subclustering clusters 4, 5, and 6

mye_subclustering <- readRDS("/scratch/smallapragada/rapids_pipeline/run2_mye_subclustering_2024_08_15.rds")

## Cluster 4 

## UMAP visualization of clusters in subclustered myeloid object
DimPlot(mye_subclustering,
        reduction = "umap",
        group.by = "mye_sc_0.4_c4_0.2", 
        raster = T,
        label = T, 
        cols = randomcoloR::distinctColorPalette(20)) 

## UMAP visualization of phenotypes in subclustered myeloid object
DimPlot(mye_subclustering,
        reduction = "umap",
        group.by = "Phenotype", 
        raster = T,
        label = F, 
        cols = randomcoloR::distinctColorPalette(5)) 

## Looking at the nCount_RNA and the nFeature_RNA distribution across clusters to determine low-quality clusters
ggplot(mye_subclustering@meta.data, aes(factor(mye_subclustering@meta.data$mye_sc_0.4_c4_0.1), 
                                        (mye_subclustering@meta.data$nCount_RNA))) + geom_violin()
ggplot(mye_subclustering@meta.data, aes(factor(mye_subclustering@meta.data$mye_sc_0.4_c4_0.1), 
                                        (mye_subclustering@meta.data$nFeature_RNA))) + geom_violin()

## A summary of how many cells per cluster belong to a specific disease phenotype
table(mye_subclustering$mye_sc_0.4_c4_0.1, mye_subclustering$Phenotype)

DotPlot(mye_subclustering, group.by = "mye_sc_0.4_c4_0.1", features = c(other_features, immune_features)) + angle

## Temp object to visualize clusters of interest
obj <- subset(mye_subclustering, subset = mye_sc_0.4_c4_0.1 == "4,0" | mye_sc_0.4_c4_0.1 == "4,1")

## Top markers
Idents(obj) = "mye_sc_0.4_c4_0.1"
markers <- FindAllMarkers(obj)
markers %>%
  mutate(pct.diff = abs(pct.1 - pct.2)) %>%
  filter(avg_log2FC >= 0) %>%
  View()

FeaturePlot(obj, features = c("MARCO", "CD68", "MS4A7", "MRC1", "PLIN2", "CD44", 
                              "FCER1G", "CCL18", "PKM", "LYZ", "FCGR3A", "MCEMP1"))

DotPlot(obj, features = c("MARCO", "CD68", "MS4A7", "MRC1", "PLIN2", "CD44", 
                          "FCER1G", "CCL18", "PKM", "LYZ", "FCGR3A", "MCEMP1")) + angle

# 4,0 - Alveolar macrophages
# 4,1 - Low quality

## Cluster 5

## UMAP visualization of clusters in subclustered myeloid object
DimPlot(mye_subclustering,
        reduction = "umap",
        group.by = "mye_sc_0.4_c5_0.1", 
        raster = T,
        label = T, 
        cols = randomcoloR::distinctColorPalette(20)) 

## UMAP visualization of phenotypes in subclustered myeloid object
DimPlot(mye_subclustering,
        reduction = "umap",
        group.by = "Phenotype", 
        raster = T,
        label = F, 
        cols = randomcoloR::distinctColorPalette(5)) 

## A summary of how many cells per cluster belong to a specific disease phenotype
table(mye_subclustering$mye_sc_0.4_c5_0.1, mye_subclustering$Phenotype)

## Looking at the nCount_RNA and the nFeature_RNA distribution across clusters to determine low-quality clusters
ggplot(mye_subclustering@meta.data, aes(factor(mye_subclustering@meta.data$mye_sc_0.4_c5_0.1), 
                                        (mye_subclustering@meta.data$nCount_RNA))) + geom_violin()
ggplot(mye_subclustering@meta.data, aes(factor(mye_subclustering@meta.data$mye_sc_0.4_c5_0.1), 
                                        (mye_subclustering@meta.data$nFeature_RNA))) + geom_violin()

DotPlot(mye_subclustering, group.by = "mye_sc_0.4_c5_0.1", features = c(other_features, immune_features)) + angle

## Temp object to visualize clusters of interest
obj <- subset(mye_subclustering, subset = mye_sc_0.4_c5_0.1 == "5,0" | mye_sc_0.4_c5_0.1 == "5,1")

## Top markers
Idents(obj) = "mye_sc_0.4_c5_0.1"
markers <- FindAllMarkers(obj)
markers %>%
  mutate(pct.diff = abs(pct.1 - pct.2)) %>%
  filter(avg_log2FC >= 0) %>%
  View()

FeaturePlot(obj, features = c("PECAM1", "KDR", "GNG11", "APLNR", "CD34", "KIT", 
                              "CCL21", "CA4", "TOP2A", "MKI67"))

DotPlot(obj, features = c("PECAM1", "KDR", "GNG11", "APLNR", "CD34", "KIT", 
                          "CCL21", "CA4", "TOP2A", "MKI67")) + angle

# 5,0 - Proliferating endothelial
# 5,1 - Low quality

## Cluster 6 

## UMAP visualization of clusters in subclustered myeloid object
DimPlot(mye_subclustering,
        reduction = "umap",
        group.by = "mye_sc_0.4_c6_0.1", 
        raster = T,
        label = T, 
        cols = randomcoloR::distinctColorPalette(20)) 

## UMAP visualization of phenotype in subclustered myeloid object
DimPlot(mye_subclustering,
        reduction = "umap",
        group.by = "Phenotype", 
        raster = T,
        label = F, 
        cols = randomcoloR::distinctColorPalette(5)) 

## A summary of how many cells per cluster belong to a specific disease phenotype
table(mye_subclustering$mye_sc_0.4_c6_0.1, mye_subclustering$Phenotype)

## Looking at the nCount_RNA and the nFeature_RNA distribution across clusters to determine low-quality clusters
ggplot(mye_subclustering@meta.data, aes(factor(mye_subclustering@meta.data$mye_sc_0.4_c6_0.1), 
                                        (mye_subclustering@meta.data$nCount_RNA))) + geom_violin()
ggplot(mye_subclustering@meta.data, aes(factor(mye_subclustering@meta.data$mye_sc_0.4_c6_0.1), 
                                        (mye_subclustering@meta.data$nFeature_RNA))) + geom_violin()

DotPlot(mye_subclustering, group.by = "mye_sc_0.4_c6_0.1", features = c(other_features, immune_features)) + angle

## Temp object to visualize clusters of interest
obj <- subset(mye_subclustering, subset = mye_sc_0.4_c6_0.1 == "6,0" | mye_sc_0.4_c6_0.1 == "6,1" | mye_sc_0.4_c6_0.1 == "6,2")

## Top markers
Idents(obj) = "mye_sc_0.4_c6_0.1"
markers <- FindAllMarkers(obj)
markers %>%
  mutate(pct.diff = abs(pct.1 - pct.2)) %>%
  filter(avg_log2FC >= 0) %>%
  View()

FeaturePlot(obj, features = c("HLA-DQA1", "CD1C", "HLA-DQB1", "FCER1A", "FCER1G", 
                              "CST3", "BCL2L11", "LUM", "COL1A2", "COL3A1", "CD69", 
                              "CD1A", "MRC1", "C1QC", "GPR183", "SPP1", "SLC1A3", 
                              "MS4A7", "HMOX1", "FCGBP", "S100A8", "S100A9", "S100A12", 
                              "MARCO"))

DotPlot(obj, features = c("HLA-DQA1", "CD1C", "HLA-DQB1", "FCER1A", "FCER1G", 
                          "CST3", "BCL2L11", "LUM", "COL1A2", "COL3A1", "CD69", 
                          "CD1A", "MRC1", "C1QC", "GPR183", "SPP1", "SLC1A3", 
                          "MS4A7", "HMOX1", "FCGBP", "S100A8", "S100A9", "S100A12", 
                          "MARCO")) + angle

# 6,0 - cDC
# 6,1 - pDC
# 6,2 - Monocytes

### Reclustering clusters 3 and 8

cluster3 <- subset(mye, subset = mye_brl_0.4 == 3)
cluster8 <- subset(mye, subset = mye_brl_0.4 == 8)

saveRDS(cluster3, "/scratch/smallapragada/run2/run2_mye_reclustering_c3_2024_08_15.rds")
saveRDS(cluster8, "/scratch/smallapragada/run2/run2_mye_reclustering_c8_2024_08_15.rds")

## Cluster 3 

cluster3 <- readRDS("/scratch/smallapragada/rapids_pipeline/run2_mye_reclustering_c3_2024_08_15.rds")

## UMAP visualization of clusters in reclustered myeloid object
DimPlot(cluster3,
        reduction = "umap",
        group.by = "mye_c3_brl_0.3", 
        raster = T,
        label = T, 
        cols = randomcoloR::distinctColorPalette(20)) 

## UMAP visualization of phenotype in reclustered myeloid object
DimPlot(cluster3,
        reduction = "umap",
        group.by = "Phenotype", 
        raster = T,
        label = F, 
        cols = randomcoloR::distinctColorPalette(5)) 

## A summary of how many cells per cluster belong to a specific disease phenotype
table(cluster3$mye_c3_brl_0.3, cluster3$Phenotype)

## Looking at the nCount_RNA and the nFeature_RNA distribution across clusters to determine low-quality clusters
ggplot(cluster3@meta.data, aes(factor(cluster3@meta.data$mye_c3_brl_0.3), 
                               (cluster3@meta.data$nCount_RNA))) + geom_violin()
ggplot(cluster3@meta.data, aes(factor(cluster3@meta.data$mye_c3_brl_0.3), 
                               (cluster3@meta.data$nFeature_RNA))) + geom_violin()

DotPlot(cluster3, group.by = "mye_c3_brl_0.3", features = c(other_features, immune_features)) + angle

## Top markers
Idents(cluster3) = "mye_c3_brl_0.3"
markers <- FindAllMarkers(cluster3)
markers %>%
  mutate(pct.diff = abs(pct.1 - pct.2)) %>%
  filter(avg_log2FC >= 0) %>%
  View()

FeaturePlot(obj, features = c("LAMP3", "SFTPC", "NAPSA", "MGST1", "NKX2-1", "RNASE1", 
                              "KRT8", "KRT18", "SCGB3A2"))            ## Cluster 0

FeaturePlot(obj, features = c("NKX2-1", "SFTPC", "NAPSA", "SOX4", "SOX9", "RNASE1", 
                              "KRT8", "KRT18", "RTKN2", "TOP2A", "MKI67"))        ## Cluster 1

FeaturePlot(obj, features = c("SOD2", "LUM", "HIF1A", "DCN", "PDGFRA", "PLIN2", 
                              "HAS2", "COL1A2", "CCL2"))              ## Cluster 2

FeaturePlot(obj, features = c("COL1A2", "COL3A1", "MEG3", "COL1A1", "ACTA2", "FN1", 
                              "SPARCL1", "POSTN", "PDGFRB"))          ## Cluster 3

FeaturePlot(obj, features = c("POSTN", "SFRP2", "CTHRC1", "LUM", "DCN", "PTGDS", 
                              "MEG3", "FGF7", "FN1"))                 ## Cluster 4

FeaturePlot(obj, features = c("KIT", "TPSAB1", "ITGAX", "CPA3", "PTPRC", "FCER1G", 
                              "CD44", "CD69", "PLIN2"))               ## Cluster 5

FeaturePlot(obj, features = c("MSLN", "ICAM1", "NKX2-1", "KRT18", "VEGFA", "KRT8", 
                              "AGER", "RNASE1", "RTKN2"))             ## Cluster 6

FeaturePlot(obj, features = c("PGC", "KRT8", "KRT18", "SFTPD", "RTKN2", "NAPSA", 
                              "SCGB3A2", "MS4A7", "NOX4"))            ## Cluster 7

FeaturePlot(obj, features = c("PECAM1", "KLRG1", "CD34", "FCN3", "KDR", "IL7R", 
                              "CA4", "SPARCL1", "BMPR2", "RAMP2", "GNG11", "CLDN5"))          ## Cluster 8

FeaturePlot(obj, features = c("PTGDS", "CCN2", "DCN", "PDGFRA", "JCHAIN", 
                              "HERPUD1", "TPSAB1", "XBP1", "GSR"))    ## Cluster 9

FeaturePlot(obj, features = c("PTPRC", "CD8A", "PIM2", "IL7R", "TGFB1", "KLRB1", 
                              "CD247", "TRAC", "GZMA"))               ## Cluster 10

DotPlot(obj, features = c("LAMP3", "SFTPC", "NAPSA", "MGST1", "NKX2-1", "RNASE1", 
                          "KRT8", "KRT18", "SCGB3A2", "SOX4", "SOX9", "RTKN2", 
                          "TOP2A", "MKI67", "SOD2", "LUM", "HIF1A", "DCN", "PDGFRA", 
                          "PLIN2", "HAS2", "COL1A2", "CCL2", "COL3A1", "MEG3", 
                          "COL1A1", "ACTA2", "FN1", "SPARCL1", "POSTN", "PDGFRB", 
                          "SFRP2", "CTHRC1", "PTGDS", "FGF7", "KIT", "TPSAB1", 
                          "ITGAX", "CPA3", "PTPRC", "FCER1G", "CD44", "CD69", 
                          "MSLN", "ICAM1", "VEGFA", "AGER", "PGC", "SFTPD", "MS4A7", 
                          "NOX4", "PECAM1", "KLRG1", "CD34", "FCN3", "KDR", "IL7R", 
                          "CA4", "BMPR2", "RAMP2", "GNG11", "CLDN5", "CCN2", "JCHAIN", 
                          "HERPUD1", "XBP1", "GSR", "CD8A", "PIM2", "TGFB1", "KLRB1", 
                          "CD247", "TRAC", "GZMA")) + angle

# Cluster 0 - AT2
# Cluster 1 - Fetal-enriched AT2
# Cluster 2 - Fetal-enriched SOD2+ myoFB
# Cluster 3 - SMC
# Cluster 4 - Activated fibroblasts
# Cluster 5 - Mast
# Cluster 6 - AT1
# Cluster 7 - Transitional AT2
# Cluster 8 - Capillaries
# Cluster 9 - Activated fibroblasts
# Cluster 10 - T

## Cluster 8 

cluster8 <- readRDS("/scratch/smallapragada/rapids_pipeline/run2_mye_reclustering_c8_2024_08_15.rds")

## UMAP visualization of clusters in reclustered myeloid object
DimPlot(cluster8,
        reduction = "umap",
        group.by = "mye_c8_brl_0.3", 
        raster = T,
        label = T, 
        cols = randomcoloR::distinctColorPalette(20)) 

## UMAP visualization of phenotype in reclustered myeloid object
DimPlot(cluster8,
        reduction = "umap",
        group.by = "Phenotype", 
        raster = T,
        label = F, 
        cols = randomcoloR::distinctColorPalette(5))

## A summary of how many cells per cluster belong to a specific disease phenotype
table(cluster8$mye_c8_brl_0.3, cluster8$Phenotype)

## Looking at the nCount_RNA and the nFeature_RNA distribution across clusters to determine low-quality clusters
ggplot(cluster8@meta.data, aes(factor(cluster8@meta.data$mye_c8_brl_0.3), 
                               (cluster8@meta.data$nCount_RNA))) + geom_violin()
ggplot(cluster8@meta.data, aes(factor(cluster8@meta.data$mye_c8_brl_0.3), 
                               (cluster8@meta.data$nFeature_RNA))) + geom_violin()

DotPlot(cluster8, group.by = "mye_c8_brl_0.3", features = c(other_features, immune_features, "CTLA4", "FOXP3")) + angle

## Top markers
Idents(cluster8) = "mye_c8_brl_0.3"
markers <- FindAllMarkers(cluster8)
markers %>%
  mutate(pct.diff = abs(pct.1 - pct.2)) %>%
  filter(avg_log2FC >= 0) %>%
  View()

FeaturePlot(cluster8, features = c("MEG3", "POSTN", "SNAI2", "COL1A2", "FN1", 
                                   "COL3A1", "HAS2", "PLIN2", "PDGFRA", "SOD2", 
                                   "VEGFA"))          ## Cluster 0

FeaturePlot(obj, features = c("MEG3", "SOD2", "LAMP3", "GPR183", "CCR7", "NFKB1", 
                              "IL7R", "CD86", "CST3"))        ## Cluster 1

FeaturePlot(obj, features = c("S100A9", "MARCO", "S100A8", "CD68", "PLIN2", "CD14", 
                              "MCEMP1", "CCL18", "CPA3"))     ## Cluster 2

# Cluster 0 - FB
# Cluster 1 - pDC
# Cluster 2 - Low quality

## Adding cell type labels to each of the reclustered objects
cluster3$CT <- ""
cluster3$CT[cluster3$mye_c3_brl_0.3 == 0] <- "AT2"
cluster3$CT[cluster3$mye_c3_brl_0.3 == 1] <- "Fetal-enriched AT2"
cluster3$CT[cluster3$mye_c3_brl_0.3 == 2] <- "Fetal-enriched SOD2+ myoFB"
cluster3$CT[cluster3$mye_c3_brl_0.3 == 3] <- "SMC"
cluster3$CT[cluster3$mye_c3_brl_0.3 == 4 | cluster3$mye_c3_brl_0.3 == 9] <- "Activated fibroblasts"
cluster3$CT[cluster3$mye_c3_brl_0.3 == 5] <- "Mast"
cluster3$CT[cluster3$mye_c3_brl_0.3 == 6] <- "AT1"
cluster3$CT[cluster3$mye_c3_brl_0.3 == 7] <- "Transitional AT2"
cluster3$CT[cluster3$mye_c3_brl_0.3 == 8] <- "Capillaries"
cluster3$CT[cluster3$mye_c3_brl_0.3 == 10] <- "T"

cluster8$CT <- ""
cluster8$CT[cluster8$mye_c8_brl_0.3 == 0] <- "Fibroblasts"
cluster8$CT[cluster8$mye_c8_brl_0.3 == 1] <- "pDC"
cluster8$CT[cluster8$mye_c8_brl_0.3 == 2] <- "Low quality"

## Combining objects together and merging into one
firstpass <- merge(x = cluster3, y = cluster8)

## Adding labels from the merged reclustered object to the whole airway object + all other CT labels 
mye_subclustering$CT <- ""
mye_subclustering$CT <- firstpass$CT
mye_subclustering$CT[mye_subclustering$mye_brl_0.4 == 0] <- "Monocytes"
mye_subclustering$CT[mye_subclustering$mye_brl_0.4 == 1] <- "Proliferating Meg-Ery"
mye_subclustering$CT[mye_subclustering$mye_brl_0.4 == 2] <- "Fibroblasts"
mye_subclustering$CT[mye_subclustering$mye_sc_0.4_c4_0.1 == "4,0"] <- "Alveolar macrophages"
mye_subclustering$CT[mye_subclustering$mye_sc_0.4_c4_0.1 == "4,1"] <- "Low quality"
mye_subclustering$CT[mye_subclustering$mye_sc_0.4_c5_0.1 == "5,0"] <- "Proliferating endothelial"
mye_subclustering$CT[mye_subclustering$mye_sc_0.4_c5_0.1 == "5,1"] <- "Low quality"
mye_subclustering$CT[mye_subclustering$mye_sc_0.4_c6_0.1 == "6,0"] <- "cDC"
mye_subclustering$CT[mye_subclustering$mye_sc_0.4_c6_0.1 == "6,1"] <- "pDC"
mye_subclustering$CT[mye_subclustering$mye_sc_0.4_c6_0.1 == "6,2"] <- "Monocytes"
mye_subclustering$CT[mye_subclustering$mye_brl_0.4 == 7] <- "Proliferating fibroblasts"

## Save myeloid object
saveRDS(mye_subclustering, "/scratch/smallapragada/run2/final_mye_CT_2024_08_15.rds")

#### COMBINING SUBLINEAGE OBJECTS AND INTEGRATING CELL IDS TO ENTIRE OBJECT ----

## Matching cell ids from combined object to clustered object
run2$CT <- ""
run2$CT[run2$leiden_brl_0.5 == 3] <- "Low quality"
run2$CT <- airway_subclustering$CT
run2$CT <- alv$CT
run2$CT <- endo$CT
run2$CT <- mes_subclustering$CT
run2$CT <- lymph_subclustering$CT
run2$CT <- mye_subclustering$CT

## Removing low quality cells from dataset
run2_final <- subset(run2, subset = CT != "Low quality")

## Plotting all cell types so the labels don't overlap each other
plot <- DimPlot(run2_final,
                reduction = "umap",
                group.by = "CT", 
                raster = T,
                label = F, 
                cols = randomcoloR::distinctColorPalette(49)) + NoLegend()

LabelClusters(plot, id = "CT", repel = T, max.overlaps = Inf)

## Save annotated labels for whole dataset
saveRDS(run2_final, "/scratch/smallapragada/run2/final_CT_2024_08_15.rds")

#### CLEANING UP ANNOTATIONS AND FINALIZING LABELS ----

## Transferring current CT annotations to a new column (CT_final) to go more granular
run2_final$CT_final <- run2_final$CT

### Reclustering endothelial cell types & neutrophils to separate them out better

## Subsetting capillaries, venous, lymphatic endothelial, and neutrophil cells & saving for RAPIDS processing
endo_cleanup <- subset(run2_final, subset = CT == "Capillaries" | CT == "Venous" | 
                         CT == "Neutrophils" | CT == "Lymphatic endothelial")
saveRDS(endo_cleanup, "/scratch/smallapragada/run2/run2_endo_CTs_cleanup_2024_08_15.rds")

## Loading in reclustered object post-RAPIDS
endo_cleanup <- readRDS("/scratch/smallapragada/rapids_pipeline/run2_endo_CTs_cleanup_2024_08_15.rds")

## UMAP visualization of clusters in reclustered object
DimPlot(endo_cleanup,
        reduction = "umap",
        group.by = "endo_cleanup_brl_0.5", 
        raster = T,
        label = T, 
        cols = randomcoloR::distinctColorPalette(20)) 

## UMAP visualization of phenotype in reclustered object
DimPlot(endo_cleanup,
        reduction = "umap",
        group.by = "Phenotype", 
        raster = T,
        label = F, 
        cols = randomcoloR::distinctColorPalette(5)) 

## A summary of how many cells per cluster belong to a specific disease phenotype
table(endo_cleanup$endo_cleanup_brl_0.5, endo_cleanup$Phenotype)

## Looking at the nCount_RNA and the nFeature_RNA distribution across clusters to determine low-quality clusters
ggplot(endo_cleanup@meta.data, aes(factor(endo_cleanup@meta.data$endo_cleanup_brl_0.5), 
                                   (endo_cleanup@meta.data$nCount_RNA))) + geom_violin()
ggplot(endo_cleanup@meta.data, aes(factor(endo_cleanup@meta.data$endo_cleanup_brl_0.5), 
                                   (endo_cleanup@meta.data$nFeature_RNA))) + geom_violin()

## Top markers
Idents(endo_cleanup) = "endo_cleanup_brl_0.5"
markers <- FindAllMarkers(endo_cleanup)
markers %>%
  mutate(pct.diff = abs(pct.1 - pct.2)) %>%
  filter(avg_log2FC >= 0) %>%
  View()

DotPlot(endo_cleanup, features = c(endothelial_features, "ELANE", 
                                   mesenchymal_features, "MKI67", "TOP2A")) + angle

## Cluster 0 - Low quality
obj <- subset(endo_cleanup, subset = endo_cleanup_brl_0.5 == 0)

FeaturePlot(obj, features = c("PGC", "NAPSA", "SFTPD", "SFTPC", "KRT8", 
                              "SCGB3A2", "MGST1", "DUOX1"))

DotPlot(obj, features = c("PGC", "NAPSA", "SFTPD", "SFTPC", "KRT8", 
                          "SCGB3A2", "MGST1", "DUOX1")) + angle

## Cluster 1 - Venous/immune (Subcluster)
obj <- subset(endo_cleanup, subset = endo_cleanup_brl_0.5 == 1)

FeaturePlot(obj, features = c("SOD2", "CCL2", "ICAM1", "ACKR1", "HIF1A", "IRF1", 
                              "CCN2", "BMPR2"))

DotPlot(endo_cleanup, features = c("SOD2", "CCL2", "ICAM1", "ACKR1", "HIF1A", "IRF1", 
                                   "CCN2", "BMPR2")) + angle

## Cluster 2 - Capillaries
obj <- subset(endo_cleanup, subset = endo_cleanup_brl_0.5 == 2)

FeaturePlot(obj, features = c("APLNR", "POSTN", "SOX4", "KIT", "CA4", "MEG3", "CTNNB1",
                              "FN1", "PTGDS", "HEY1"))

DotPlot(obj, features = c("TPSAB1", "CPA3", "PTPRC"))

## Cluster 3 - Mast
obj <- subset(endo_cleanup, subset = endo_cleanup_brl_0.5 == 3)

FeaturePlot(obj, features = c("TPSAB1", "CPA3", "LUM"))

DotPlot(obj, features = c("TPSAB1", "CPA3", "LUM"))

## Cluster 4 - FB

obj <- subset(endo_cleanup, subset = endo_cleanup_brl_0.5 == 4)

FeaturePlot(obj, features = c("COL1A2", "COL1A1", "MEG3", "SPARCL1", "FN1", 
                              "COL3A1", "ACTA2", "POSTN", "PDIA4", "PDIA3", "PDGFRA",
                              "PLIN2", "WNT5A"))

## Cluster 5 - Arterial/FB (Subcluster)
obj <- subset(endo_cleanup, subset = endo_cleanup_brl_0.5 == 5)

FeaturePlot(obj, features = c("SLC7A11", "BMPR2", "FCN3", "ATF4", "GCLM", "XBP1", 
                              "HMGA1", "CCN2", "IL4R", "CA4", "COL1A2", "COL3A1"))

## Cluster 6 - Neutrophils
obj <- subset(endo_cleanup, subset = endo_cleanup_brl_0.5 == 6)

FeaturePlot(obj, features = c("ELANE", "LYZ", "NUCB2"))

FeaturePlot(obj, features = c("PTPRC", "KLRB1", "CD247", "PGC", "GZMA", "NKG7", 
                              "KLRC1", "GZMB", "FCGR3A", "CD3E"))

## Cluster 7 - Low quality
obj <- subset(endo_cleanup, subset = endo_cleanup_brl_0.5 == 7)

FeaturePlot(obj, features = c("PTPRC", "KLRB1", "CD247", "PGC", "GZMA", "NKG7", 
                              "KLRC1", "GZMB", "FCGR3A", "CD3E"))

## Cluster 8 - Venous
obj <- subset(endo_cleanup, subset = endo_cleanup_brl_0.5 == 8)

FeaturePlot(obj, features = c("COL15A1", "PLVAP", "COL3A1", "COL1A1", "CCN2", "SPARCL1", 
                              "COL1A2", "KDR", "LUM", "POSTN", "ACKR1", "APLNR"))

## Cluster 9 - Lymphatic endothelial
obj <- subset(endo_cleanup, subset = endo_cleanup_brl_0.5 == 9)

FeaturePlot(obj, features = c("CCL21", "MRC1"))

## Cluster 10 - Arterial
obj <- subset(endo_cleanup, subset = endo_cleanup_brl_0.5 == 10)

FeaturePlot(obj, features = c("IL7R", "RNASE1", "CD34", "HEY1", "BMPR2", 
                              "ITGAE", "RAMP2", "GNG11", "CA4"))

### Subclustering of clusters 1 and 5

## Loading in reclustered object post-RAPIDS
endo_sc_cleanup <- readRDS("/scratch/smallapragada/rapids_pipeline/run2_endo_CTs_cleanup_sc_2024_08_15.rds")

## Cluster 1

## UMAP visualization of clusters in reclustered object
DimPlot(endo_sc_cleanup,
        reduction = "umap",
        group.by = "endo_cleanup_sc_0.5_c1_0.2", 
        raster = T,
        label = T, 
        cols = randomcoloR::distinctColorPalette(20)) 

## UMAP visualization of phenotype in reclustered object
DimPlot(endo_sc_cleanup,
        reduction = "umap",
        group.by = "Phenotype", 
        raster = T,
        label = F, 
        cols = randomcoloR::distinctColorPalette(5)) 

## A summary of how many cells per cluster belong to a specific disease phenotype
table(endo_sc_cleanup$endo_cleanup_sc_0.5_c1_0.2, endo_sc_cleanup$Phenotype)

## Looking at the nCount_RNA and the nFeature_RNA distribution across clusters to determine low-quality clusters
ggplot(endo_sc_cleanup@meta.data, aes(factor(endo_sc_cleanup@meta.data$endo_cleanup_sc_0.5_c1_0.2), 
                                      (endo_sc_cleanup@meta.data$nCount_RNA))) + geom_violin()
ggplot(endo_sc_cleanup@meta.data, aes(factor(endo_sc_cleanup@meta.data$endo_cleanup_sc_0.5_c1_0.2), 
                                      (endo_sc_cleanup@meta.data$nFeature_RNA))) + geom_violin()

## Top markers
Idents(endo_sc_cleanup) = "endo_cleanup_sc_0.5_c1_0.2"
markers <- FindAllMarkers(endo_sc_cleanup)
markers %>%
  mutate(pct.diff = abs(pct.1 - pct.2)) %>%
  filter(avg_log2FC >= 0) %>%
  View()

DotPlot(endo_sc_cleanup, features = c(endothelial_features, "ELANE", 
                                      mesenchymal_features, "MKI67", "TOP2A")) + angle

## Temp object for the specific clusters of interest
obj <- subset(endo_sc_cleanup, subset = endo_cleanup_sc_0.5_c1_0.2 == "1,0" | 
                endo_cleanup_sc_0.5_c1_0.2 == "1,1")

FeaturePlot(obj, features = c(endothelial_features))

DotPlot(obj, features = c(endothelial_features, mesenchymal_features)) + angle

# 1,0 = Venous
# 1,1 = Capillaries

## Cluster 5

## UMAP visualization of clusters in reclustered object
DimPlot(endo_sc_cleanup,
        reduction = "umap",
        group.by = "endo_cleanup_sc_0.5_c5_0.2", 
        raster = T,
        label = T, 
        cols = randomcoloR::distinctColorPalette(20)) 

## UMAP visualization of phenotype in reclustered object
DimPlot(endo_sc_cleanup,
        reduction = "umap",
        group.by = "Phenotype", 
        raster = T,
        label = F, 
        cols = randomcoloR::distinctColorPalette(5)) 

## A summary of how many cells per cluster belong to a specific disease phenotype
table(endo_sc_cleanup$endo_cleanup_sc_0.5_c5_0.2, endo_sc_cleanup$Phenotype)

## Looking at the nCount_RNA and the nFeature_RNA distribution across clusters to determine low-quality clusters
ggplot(endo_sc_cleanup@meta.data, aes(factor(endo_sc_cleanup@meta.data$endo_cleanup_sc_0.5_c5_0.2), 
                                      (endo_sc_cleanup@meta.data$nCount_RNA))) + geom_violin()
ggplot(endo_sc_cleanup@meta.data, aes(factor(endo_sc_cleanup@meta.data$endo_cleanup_sc_0.5_c5_0.2), 
                                      (endo_sc_cleanup@meta.data$nFeature_RNA))) + geom_violin()

## Temp object for the specific clusters of interest
obj <- subset(endo_sc_cleanup, subset = endo_cleanup_sc_0.5_c5_0.2 == "5,0" | 
                endo_cleanup_sc_0.5_c5_0.2 == "5,1" | endo_cleanup_sc_0.5_c5_0.2 == "5,2")

## Top markers
Idents(obj) = "endo_cleanup_sc_0.5_c5_0.2"
markers <- FindAllMarkers(obj)
markers %>%
  mutate(pct.diff = abs(pct.1 - pct.2)) %>%
  filter(avg_log2FC >= 0) %>%
  View()

FeaturePlot(obj, features = c("HEY1", "ACKR1"))

DotPlot(obj, features = c(epithelial_features)) + angle

# 5,0 = Venous
# 5,1 & 5,2 = AKR1C1+/1C2+

## Adding cell type labels to reclustered object
endo_sc_cleanup$CT_final[endo_sc_cleanup$endo_cleanup_brl_0.5 == 0] <- "Low quality"
endo_sc_cleanup$CT_final[endo_sc_cleanup$endo_cleanup_sc_0.5_c1_0.2 == "1,0"] <- "Venous"
endo_sc_cleanup$CT_final[endo_sc_cleanup$endo_cleanup_sc_0.5_c1_0.2 == "1,1"] <- "Capillaries"
endo_sc_cleanup$CT_final[endo_sc_cleanup$endo_cleanup_brl_0.5 == 2] <- "Capillaries"
endo_sc_cleanup$CT_final[endo_sc_cleanup$endo_cleanup_brl_0.5 == 3] <- "Mast"
endo_sc_cleanup$CT_final[endo_sc_cleanup$endo_cleanup_brl_0.5 == 4] <- "HAS2+ fibroblasts"
endo_sc_cleanup$CT_final[endo_sc_cleanup$endo_cleanup_sc_0.5_c5_0.2 == "5,0"] <- "Venous"
endo_sc_cleanup$CT_final[endo_sc_cleanup$endo_cleanup_sc_0.5_c5_0.2 == "5,1"] <- "AKR1C1+ & AKR1C2+"
endo_sc_cleanup$CT_final[endo_sc_cleanup$endo_cleanup_sc_0.5_c5_0.2 == "5,2"] <- "AKR1C1+ & AKR1C2+"
endo_sc_cleanup$CT_final[endo_sc_cleanup$endo_cleanup_brl_0.5 == 6] <- "Neutrophils"
endo_sc_cleanup$CT_final[endo_sc_cleanup$endo_cleanup_brl_0.5 == 7] <- "Low quality"
endo_sc_cleanup$CT_final[endo_sc_cleanup$endo_cleanup_brl_0.5 == 8] <- "Venous"
endo_sc_cleanup$CT_final[endo_sc_cleanup$endo_cleanup_brl_0.5 == 9] <- "Lymphatic endothelial"
endo_sc_cleanup$CT_final[endo_sc_cleanup$endo_cleanup_brl_0.5 == 10] <- "Arterial"

## Adding new labels to original object
run2_final$CT_final <- endo_sc_cleanup$CT_final

DotPlot(run2_final, group.by = "CT_final", features = c(other_features, endothelial_features, "ELANE")) + angle

DimPlot(endo_sc_cleanup,
        reduction = "umap",
        group.by = "CT_final", 
        raster = T,
        label = T, 
        cols = randomcoloR::distinctColorPalette(10)) 

### Reclustering transitional AT2/mast population & proliferating AT2 cells

## Subsetting all AT2 cell populations & saving for RAPIDS processing
at2_cleanup <- subset(run2_final, subset = CT == "Fetal-enriched AT2" | 
                        CT == "Fetal-enriched transitional AT2" |
                        CT == "Proliferating AT2" | CT == "Transitional AT2" | 
                        CT == "Proliferating Fetal-enriched AT2" | CT == "AT2")
saveRDS(at2_cleanup, "/scratch/smallapragada/run2/run2_at2_CTs_cleanup_2024_08_15.rds")

## Loading in reclustered object post-RAPIDS
at2_cleanup <- readRDS("/scratch/smallapragada/rapids_pipeline/run2_at2_CTs_cleanup_2024_08_15.rds")

## UMAP visualization of clusters in reclustered object
DimPlot(at2_cleanup,
        reduction = "umap",
        group.by = "at2_cleanup_brl_0.5", 
        raster = T,
        label = T, 
        cols = randomcoloR::distinctColorPalette(10)) 

## UMAP visualization of phenotype in reclustered object
DimPlot(at2_cleanup,
        reduction = "umap",
        group.by = "Phenotype", 
        raster = T,
        label = F, 
        cols = randomcoloR::distinctColorPalette(5)) 

## A summary of how many cells per cluster belong to a specific disease phenotype
table(at2_cleanup$at2_cleanup_brl_0.5, at2_cleanup$Phenotype)

## Looking at the nCount_RNA and the nFeature_RNA distribution across clusters to determine low-quality clusters
ggplot(at2_cleanup@meta.data, aes(factor(at2_cleanup@meta.data$at2_cleanup_brl_0.5), 
                                  (at2_cleanup@meta.data$nCount_RNA))) + geom_violin()
ggplot(at2_cleanup@meta.data, aes(factor(at2_cleanup@meta.data$at2_cleanup_brl_0.5), 
                                  (at2_cleanup@meta.data$nFeature_RNA))) + geom_violin()

## Top markers
Idents(at2_cleanup) = "at2_cleanup_brl_0.5"
markers <- FindAllMarkers(at2_cleanup)
markers %>%
  mutate(pct.diff = abs(pct.1 - pct.2)) %>%
  filter(avg_log2FC >= 0) %>%
  View()

DotPlot(at2_cleanup, features = c(other_features, epithelial_features)) + angle

## Cluster 0 - Low quality

obj <- subset(at2_cleanup, subset = at2_cleanup_brl_0.5 == 0)

FeaturePlot(obj, features = c("PTPRC", "LYZ", "KLRC1", "KLRB1", "MS4A7", "CD2", 
                              "CD3E", "CD28", "NKG7", "MS4A1", "S100A8", "S100A9"))

## Cluster 1 - Low quality/AT1/Secretory 3A2+ (Subcluster)

obj <- subset(at2_cleanup, subset = at2_cleanup_brl_0.5 == 1)

FeaturePlot(obj, features = c("SCGB3A2", "RTKN2", "AGER", "KRT15", "SLC2A1", "GDF15", 
                              "KRT17", "KRT14", "AGR3", "SCGB1A1", "FGF2"))

## Cluster 2 - Fetal-enriched AT2

obj <- subset(at2_cleanup, subset = at2_cleanup_brl_0.5 == 2)

FeaturePlot(obj, features = c("HES1", "SOX4", "SOX9", "VEGFA", "CD44", "NKX2-1", 
                              "EPCAM", "CTNNB1", "TGFB2", "NAPSA", "SFTPC", "KRT8"))

## Cluster 3 - Fetal-enriched proliferating alveolar

obj <- subset(at2_cleanup, subset = at2_cleanup_brl_0.5 == 3)

FeaturePlot(obj, features = c("TOP2A", "MKI67", "RTKN2", "CCNB2", "LAMP3", "HIST1H1C", 
                              "HMGA1", "PCNA", "RACGAP1", "NKX2-1", "AGER", "SFTPC",
                              "SPRY2", "YAP1", "FCER1G", "PKM", "WNT5A", "ATF"))

## Cluster 4 - Mast

obj <- subset(at2_cleanup, subset = at2_cleanup_brl_0.5 == 4)

FeaturePlot(obj, features = c("CPA3", "TPSAB1", "LYZ", "CCL2", "S100A8", "S100A9", 
                              "MS4A7", "JCHAIN"))

DotPlot(obj, features = c("CPA3", "TPSAB1", "LYZ", "CCL2", "S100A8", "S100A9", 
                          "MS4A7", "JCHAIN")) + angle

## Cluster 5 - Basal

obj <- subset(at2_cleanup, subset = at2_cleanup_brl_0.5 == 5)

FeaturePlot(obj, features = c("KRT15", "RTKN2", "KRT8", "KRT5", "KRT17", 
                              "KRT14", "AGER", "ICAM1", "TP63"))

## Cluster 6 - Fetal FB

obj <- subset(at2_cleanup, subset = at2_cleanup_brl_0.5 == 6)

FeaturePlot(obj, features = c("MEG3", "PDGFRA", "COL1A2", "COL3A1", "SOD2", "FN1", 
                              "COL1A1", "LUM", "DCN", "POSTN", "PLIN2", "HAS2", 
                              "PDGFRB", "ACTA2", "CSPG4", "WNT5A"))

## Cluster 7 - Fetal transitional AT2

obj <- subset(at2_cleanup, subset = at2_cleanup_brl_0.5 == 7)

FeaturePlot(obj, features = c("MMP7", "HYOU1", "MGST1", "KRT18", "PDIA4", "DMBT1", 
                              "KRT8", "CPA3", "HSPA5", "RNASE1", "LAMP3", 
                              "RNASE1", "SOD2", "SFTPC"))

## Cluster 8 - Monocytes

obj <- subset(at2_cleanup, subset = at2_cleanup_brl_0.5 == 8)

FeaturePlot(obj, features = c("S100A8", "S100A9", "S100A12", "SPP1", "MARCO"))

## Cluster 9 - AT2

obj <- subset(at2_cleanup, subset = at2_cleanup_brl_0.5 == 9)

FeaturePlot(obj, features = c("MSLN", "DUOX1", "CEACAM6", "MGST1", "AGER", "ITGA3", "ICAM1", 
                              "SFTPD", "NAPSA", "CD34", "ILR7", "CA4", "LAMP3", "SCGB3A2", "KRT8", "SFTPC"))

### Subclustering of cluster 1 

## Loading in reclustered object post-RAPIDS
at2_sc_cleanup <- readRDS("/scratch/smallapragada/rapids_pipeline/run2_at2_CTs_cleanup_sc_2024_08_15.rds")

## Cluster 1

## UMAP visualization of clusters in reclustered object
DimPlot(at2_sc_cleanup,
        reduction = "umap",
        group.by = "at2_cleanup_brl_0.5_c1_0.2", 
        raster = T,
        label = T, 
        cols = randomcoloR::distinctColorPalette(20)) 

## UMAP visualization of phenotype in reclustered object
DimPlot(at2_sc_cleanup,
        reduction = "umap",
        group.by = "Phenotype", 
        raster = T,
        label = F, 
        cols = randomcoloR::distinctColorPalette(5)) 

## A summary of how many cells per cluster belong to a specific disease phenotype
table(at2_sc_cleanup$at2_cleanup_brl_0.5_c1_0.2, at2_sc_cleanup$Phenotype)

## Looking at the nCount_RNA and the nFeature_RNA distribution across clusters to determine low-quality clusters
ggplot(at2_sc_cleanup@meta.data, aes(factor(at2_sc_cleanup@meta.data$at2_cleanup_brl_0.5_c1_0.2), 
                                     (at2_sc_cleanup@meta.data$nCount_RNA))) + geom_violin()
ggplot(at2_sc_cleanup@meta.data, aes(factor(at2_sc_cleanup@meta.data$at2_cleanup_brl_0.5_c1_0.2), 
                                     (at2_sc_cleanup@meta.data$nFeature_RNA))) + geom_violin()

## Temp object for the specific clusters of interest
obj <- subset(at2_sc_cleanup, subset = at2_cleanup_brl_0.5_c1_0.2 == "1,0" | 
                at2_cleanup_brl_0.5_c1_0.2 == "1,1" | 
                at2_cleanup_brl_0.5_c1_0.2 == "1,2")

## Top markers
Idents(obj) = "at2_cleanup_brl_0.5_c1_0.2"
markers <- FindAllMarkers(obj)
markers %>%
  mutate(pct.diff = abs(pct.1 - pct.2)) %>%
  filter(avg_log2FC >= 0) %>%
  View()

FeaturePlot(obj, features = c("VEGFA", "RTKN2", "SCGB3A2", "AGER", "MSLN", 
                              "ICAM1", "SFTPD", "SFTPC", "MKI67", "PGC", 
                              "LAMP3", "SOD2", "KRT8", "MGST1", "MMP7", "HSPA5"))

DotPlot(obj, features = c("VEGFA", "RTKN2", "SCGB3A2", "AGER", "MSLN", 
                          "ICAM1", "SFTPC", "MKI67", "SFTPD", "PGC", 
                          "LAMP3", "SOD2", "KRT8", "MGST1", "MMP7", "HSPA5")) + angle

# 1,0 - AT1
# 1,1 - Low quality
# 1,2 - Transitional AT2

## Adding cell type labels to reclustered object
at2_sc_cleanup$CT_final[at2_sc_cleanup$at2_cleanup_brl_0.5 == 0] <- "Low quality"
at2_sc_cleanup$CT_final[at2_sc_cleanup$at2_cleanup_brl_0.5_c1_0.2 == "1,0"] <- "AT1"
at2_sc_cleanup$CT_final[at2_sc_cleanup$at2_cleanup_brl_0.5_c1_0.2 == "1,1"] <- "Low quality"
at2_sc_cleanup$CT_final[at2_sc_cleanup$at2_cleanup_brl_0.5_c1_0.2 == "1,2"] <- "Transitional AT2"
at2_sc_cleanup$CT_final[at2_sc_cleanup$at2_cleanup_brl_0.5 == 2] <- "Fetal-enriched AT2"
at2_sc_cleanup$CT_final[at2_sc_cleanup$at2_cleanup_brl_0.5 == 3] <- "Proliferating Fetal-enriched alveolar"
at2_sc_cleanup$CT_final[at2_sc_cleanup$at2_cleanup_brl_0.5 == 4] <- "Mast"
at2_sc_cleanup$CT_final[at2_sc_cleanup$at2_cleanup_brl_0.5 == 5] <- "Basal"
at2_sc_cleanup$CT_final[at2_sc_cleanup$at2_cleanup_brl_0.5 == 6] <- "Fetal-enriched fibroblasts"
at2_sc_cleanup$CT_final[at2_sc_cleanup$at2_cleanup_brl_0.5 == 7] <- "Fetal-enriched transitional AT2"
at2_sc_cleanup$CT_final[at2_sc_cleanup$at2_cleanup_brl_0.5 == 8] <- "Monocytes"
at2_sc_cleanup$CT_final[at2_sc_cleanup$at2_cleanup_brl_0.5 == 9] <- "AT2"

## Adding new labels to original object
run2_final$CT_final <- at2_sc_cleanup$CT_final

DotPlot(run2_final, group.by = "CT_final", features = c(other_features, epithelial_features)) + angle

### Reclustering airway cell types to separate them out better

## Subsetting Secretory cell types and proliferating airway cells & saving for RAPIDS processing
airway_cleanup <- subset(run2_final, subset = CT == "Secretory 3A2+ & 1A1+" | 
                           CT == "Multiciliated" | CT == "Secretory MUC5B+")

saveRDS(airway_cleanup, "/scratch/smallapragada/run2/run2_airway_CTs_cleanup_2024_08_15.rds")

## Loading in reclustered object post-RAPIDS
airway_cleanup <- readRDS("/scratch/smallapragada/rapids_pipeline/run2_airway_CTs_cleanup_2024_08_15.rds")

## UMAP visualization of clusters in reclustered object
DimPlot(airway_cleanup,
        reduction = "umap",
        group.by = "airway_cleanup_brl_0.3", 
        raster = T,
        label = T, 
        cols = randomcoloR::distinctColorPalette(15)) 

## UMAP visualization of phenotype in reclustered object
DimPlot(airway_cleanup,
        reduction = "umap",
        group.by = "Phenotype", 
        raster = T,
        label = F, 
        cols = randomcoloR::distinctColorPalette(5)) 

## A summary of how many cells per cluster belong to a specific disease phenotype
table(airway_cleanup$airway_cleanup_brl_0.3, airway_cleanup$Phenotype)

## Top markers
Idents(airway_cleanup) = "airway_cleanup_brl_0.3"
markers <- FindAllMarkers(airway_cleanup)
markers %>%
  mutate(pct.diff = abs(pct.1 - pct.2)) %>%
  filter(avg_log2FC >= 0) %>%
  View()

DotPlot(airway_cleanup, features = c(other_features, epithelial_features)) + angle

## Cluster 0 - Secretory MUC5B+/Multiciliated (Subcluster)

obj <- subset(airway_cleanup, subset = airway_cleanup_brl_0.3 == 0)

FeaturePlot(obj, features = c("CCN2", "XBP1", "WFDC2", "KRT8", "DUOX1", "NUCB2", 
                              "BPIFA1", "CDH26", "AGR3", "FOXJ1", "C20orf85", 
                              "MUC5B"))

DotPlot(obj, features = c("CCN2", "XBP1", "WFDC2", "KRT8", "DUOX1", "NUCB2", 
                          "BPIFA1", "CDH26", "AGR3", "FOXJ1", "C20orf85", 
                          "MUC5B")) + angle

## Cluster 1 - Monocyte

obj <- subset(airway_cleanup, subset = airway_cleanup_brl_0.3 == 1)

FeaturePlot(obj, features = c("PTPRC", "ITGAX", "S100A8", "S100A9", "FCGR3A", 
                              "SLC25A37", "FCER1G", "MCEMP1", "AIF1"))

DotPlot(obj, features = c("PTPRC", "ITGAX", "S100A8", "S100A9", "FCGR3A", 
                          "SLC25A37", "FCER1G", "MCEMP1", "AIF1")) + angle

## Cluster 2 - FB/Multiciliated (Subcluster)

obj <- subset(airway_cleanup, subset = airway_cleanup_brl_0.3 == 2)

FeaturePlot(obj, features = c("SPARCL1", "FN1", "DCN", "LUM", "COL1A2", "C20orf85", 
                              "AGR3", "SCGB3A2", "MUC5B", "SCGB1A1"))

DotPlot(obj, features = c("SPARCL1", "FN1", "DCN", "LUM", "COL1A2", "C20orf85", 
                          "AGR3", "SCGB3A2", "MUC5B", "SCGB1A1")) + angle

## Cluster 3 - Multiciliated

obj <- subset(airway_cleanup, subset = airway_cleanup_brl_0.3 == 3)

FeaturePlot(obj, features = c("C20orf85", "AGR3", "FOXJ1", "CCNA1", "CD4", "LMAN1", 
                              "SEC11C", "ITGB1", "MUC5B", "SCGB3A2"))

DotPlot(obj, features = c("C20orf85", "AGR3", "FOXJ1", "CCNA1", "CD4", "LMAN1", 
                          "SEC11C", "ITGB1", "MUC5B", "SCGB3A2")) + angle

## Cluster 4 - Multiciliated/Secretory MUC5B+ (Subcluster)

obj <- subset(airway_cleanup, subset = airway_cleanup_brl_0.3 == 4)

FeaturePlot(obj, features = c("STAT1", "SAA2", "SOD2", "MUC5B", "WFDC2", "IFIT3", 
                              "KRT8", "SPCS3", "SPCS2", "OAS2", "CXCL9", "AGR3", 
                              "SCGB3A2", "SCGB1A1"))

DotPlot(obj, features = c("STAT1", "SAA2", "SOD2", "MUC5B", "WFDC2", "IFIT3", 
                          "KRT8", "SPCS3", "SPCS2", "OAS2", "CXCL9", "AGR3", 
                          "SCGB3A2", "SCGB1A1")) + angle

## Cluster 5 - Multiciliated

obj <- subset(airway_cleanup, subset = airway_cleanup_brl_0.3 == 5)

FeaturePlot(obj, features = c("RNASE1", "CEACAM6", "TOP2A", "VEGFA", "MKI67", 
                              "SCGB3A2", "KRT8", "IDH1", "VIM", "MUC5B", "SCGB1A1",
                              "C20orf85"))
DotPlot(obj, features = c("RNASE1", "CEACAM6", "TOP2A", "VEGFA", "MKI67", 
                          "SCGB3A2", "KRT8", "IDH1", "VIM", "MUC5B", "SCGB1A1",
                          "C20orf85")) + angle

### Subclustering of clusters 0, 2, and 4

## Loading in reclustered object post-RAPIDS
airway_sc_cleanup <- readRDS("/scratch/smallapragada/rapids_pipeline/run2_airway_CTs_cleanup_sc_2024_08_15.rds")

## Cluster 0

## UMAP visualization of clusters in reclustered object
DimPlot(airway_sc_cleanup,
        reduction = "umap",
        group.by = "airway_cleanup_brl_0.3_c0_0.2", 
        raster = T,
        label = T, 
        cols = randomcoloR::distinctColorPalette(20)) 

## UMAP visualization of phenotype in reclustered object
DimPlot(airway_sc_cleanup,
        reduction = "umap",
        group.by = "Phenotype", 
        raster = T,
        label = F, 
        cols = randomcoloR::distinctColorPalette(5)) 

## A summary of how many cells per cluster belong to a specific disease phenotype
table(airway_sc_cleanup$airway_cleanup_brl_0.3_c0_0.2, airway_sc_cleanup$Phenotype)

## Looking at the nCount_RNA and the nFeature_RNA distribution across clusters to determine low-quality clusters
ggplot(airway_sc_cleanup@meta.data, aes(factor(airway_sc_cleanup@meta.data$airway_cleanup_brl_0.3_c0_0.2), 
                                        (airway_sc_cleanup@meta.data$nCount_RNA))) + geom_violin()
ggplot(airway_sc_cleanup@meta.data, aes(factor(airway_sc_cleanup@meta.data$airway_cleanup_brl_0.3_c0_0.2), 
                                        (airway_sc_cleanup@meta.data$nFeature_RNA))) + geom_violin()

## Temp object for the specific clusters of interest
obj <- subset(airway_sc_cleanup, subset = airway_cleanup_brl_0.3_c0_0.2 == "0,0" | 
                airway_cleanup_brl_0.3_c0_0.2 == "0,1" |
                airway_cleanup_brl_0.3_c0_0.2 == "0,2")

## Top markers
Idents(obj) = "airway_cleanup_brl_0.3_c0_0.2"
markers <- FindAllMarkers(obj)
markers %>%
  mutate(pct.diff = abs(pct.1 - pct.2)) %>%
  filter(avg_log2FC >= 0) %>%
  View()

FeaturePlot(obj, features = c("SCGB1A1", "SCGB3A2", "MUC5B", "SFTPD", 
                              "MUC5AC", "AGR3", "C20orf85"))

DotPlot(obj, features = c("SCGB1A1", "SCGB3A2", "MUC5B", "SFTPD", 
                          "MUC5AC", "AGR3", "C20orf85")) + angle

# 0,0 - MUC5B+
# 0,1 - Multiciliated
# 0,2 - Multiciliated

## Cluster 2

## UMAP visualization of clusters in reclustered object
DimPlot(airway_sc_cleanup,
        reduction = "umap",
        group.by = "airway_cleanup_brl_0.3_c2_0.2", 
        raster = T,
        label = T, 
        cols = randomcoloR::distinctColorPalette(20)) 

## UMAP visualization of phenotype in reclustered object
DimPlot(airway_sc_cleanup,
        reduction = "umap",
        group.by = "Phenotype", 
        raster = T,
        label = F, 
        cols = randomcoloR::distinctColorPalette(5)) 

## A summary of how many cells per cluster belong to a specific disease phenotype
table(airway_sc_cleanup$airway_cleanup_brl_0.3_c2_0.2, airway_sc_cleanup$Phenotype)

## Looking at the nCount_RNA and the nFeature_RNA distribution across clusters to determine low-quality clusters
ggplot(airway_sc_cleanup@meta.data, aes(factor(airway_sc_cleanup@meta.data$airway_cleanup_brl_0.3_c2_0.2), 
                                        (airway_sc_cleanup@meta.data$nCount_RNA))) + geom_violin()
ggplot(airway_sc_cleanup@meta.data, aes(factor(airway_sc_cleanup@meta.data$airway_cleanup_brl_0.3_c2_0.2), 
                                        (airway_sc_cleanup@meta.data$nFeature_RNA))) + geom_violin()

## Temp object for the specific clusters of interest
obj <- subset(airway_sc_cleanup, subset = airway_cleanup_brl_0.3_c2_0.2 == "2,0" | 
                airway_cleanup_brl_0.3_c2_0.2 == "2,1" |
                airway_cleanup_brl_0.3_c2_0.2 == "2,2")

## Top markers
Idents(obj) = "airway_cleanup_brl_0.3_c2_0.2"
markers <- FindAllMarkers(obj)
markers %>%
  mutate(pct.diff = abs(pct.1 - pct.2)) %>%
  filter(avg_log2FC >= 0) %>%
  View()

# Plot to spatially visualize where certain cells are located
samp <- subset(airway_sc_cleanup, subset = Sample == "PDL002")

DimPlot(samp,
        reduction = "SP",
        group.by = "airway_cleanup_brl_0.3_c2_0.2", 
        raster = T,
        label = F,
        cols = c(`2,1` = "red"), na.value = "grey90") + coord_equal()

FeaturePlot(obj, features = c("SCGB1A1", "SCGB3A2", "SPARCL1", "VEGFA", "SFTPD", 
                              "COL1A1", "COL1A2", "ACTA2", "WNT5A", "PDGFRA"))

DotPlot(obj, features = c("SCGB1A1", "SCGB3A2", "SPARCL1", "VEGFA", "SFTPD", 
                          "COL1A1", "COL1A2", "ACTA2", "WNT5A", "PDGFRA")) + angle

# 2,0 - Multiciliated
# 2,1 - FB
# 2,2 - Secretory 3A2/1A1

## Cluster 4

## UMAP visualization of clusters in reclustered object
DimPlot(airway_sc_cleanup,
        reduction = "umap",
        group.by = "airway_cleanup_brl_0.3_c4_0.2", 
        raster = T,
        label = T, 
        cols = randomcoloR::distinctColorPalette(20)) 

## UMAP visualization of phenotype in reclustered object
DimPlot(airway_sc_cleanup,
        reduction = "umap",
        group.by = "Phenotype", 
        raster = T,
        label = F, 
        cols = randomcoloR::distinctColorPalette(5)) 

## A summary of how many cells per cluster belong to a specific disease phenotype
table(airway_sc_cleanup$airway_cleanup_brl_0.3_c4_0.2, airway_sc_cleanup$Phenotype)

## Looking at the nCount_RNA and the nFeature_RNA distribution across clusters to determine low-quality clusters
ggplot(airway_sc_cleanup@meta.data, aes(factor(airway_sc_cleanup@meta.data$airway_cleanup_brl_0.3_c4_0.2), 
                                        (airway_sc_cleanup@meta.data$nCount_RNA))) + geom_violin()
ggplot(airway_sc_cleanup@meta.data, aes(factor(airway_sc_cleanup@meta.data$airway_cleanup_brl_0.3_c4_0.2), 
                                        (airway_sc_cleanup@meta.data$nFeature_RNA))) + geom_violin()

## Temp object for the specific clusters of interest
obj <- subset(airway_sc_cleanup, subset = airway_cleanup_brl_0.3_c4_0.2 == "4,0" | 
                airway_cleanup_brl_0.3_c4_0.2 == "4,1")

## Top markers
Idents(obj) = "airway_cleanup_brl_0.3_c4_0.2"
markers <- FindAllMarkers(obj)
markers %>%
  mutate(pct.diff = abs(pct.1 - pct.2)) %>%
  filter(avg_log2FC >= 0) %>%
  View()

FeaturePlot(obj, features = c("SCGB1A1", "SCGB3A2", "MUC5B", "SFTPD", 
                              "MUC5AC", "AGR3", "C20orf85"))

DotPlot(obj, features = c("SCGB1A1", "SCGB3A2", "MUC5B", "SFTPD", 
                          "MUC5AC", "AGR3", "C20orf85")) + angle

# 4,0 - Multiciliated
# 4,1 - Secretory MUC5B+

## Adding cell type labels to reclustered object
airway_sc_cleanup$CT_final <- ""
airway_sc_cleanup$CT_final[airway_sc_cleanup$airway_cleanup_brl_0.3_c0_0.2 == "0,0"] <- "Secretory MUC5B+"
airway_sc_cleanup$CT_final[airway_sc_cleanup$airway_cleanup_brl_0.3_c0_0.2 == "0,1"] <- "Multiciliated"
airway_sc_cleanup$CT_final[airway_sc_cleanup$airway_cleanup_brl_0.3_c0_0.2 == "0,2"] <- "Multiciliated"
airway_sc_cleanup$CT_final[airway_sc_cleanup$airway_cleanup_brl_0.3_c0_0.2 == 1] <- "Monocytes"
airway_sc_cleanup$CT_final[airway_sc_cleanup$airway_cleanup_brl_0.3_c2_0.2 == "2,0"] <- "Multiciliated"
airway_sc_cleanup$CT_final[airway_sc_cleanup$airway_cleanup_brl_0.3_c2_0.2 == "2,1"] <- "Capillaries"
airway_sc_cleanup$CT_final[airway_sc_cleanup$airway_cleanup_brl_0.3_c2_0.2 == "2,2"] <- "Secretory 3A2+ & 1A1+"
airway_sc_cleanup$CT_final[airway_sc_cleanup$airway_cleanup_brl_0.3_c2_0.2 == 3] <- "Multiciliated"
airway_sc_cleanup$CT_final[airway_sc_cleanup$airway_cleanup_brl_0.3_c4_0.2 == "4,0"] <- "Multiciliated"
airway_sc_cleanup$CT_final[airway_sc_cleanup$airway_cleanup_brl_0.3_c4_0.2 == "4,1"] <- "Secretory MUC5B+"
airway_sc_cleanup$CT_final[airway_sc_cleanup$airway_cleanup_brl_0.3_c4_0.2 == 5] <- "Multiciliated"

## Adding new labels to original object
run2_final$CT_final <- airway_sc_cleanup$CT_final

DotPlot(run2_final, group.by = "CT_final", features = c(epithelial_features, "MKI67", "TOP2A")) + angle

### Reclustering some myeloid CTs to separate them out better

## Subsetting plasma, venous, pDC, monocytes, and cDC cells & saving for RAPIDS processing
mye_cleanup <- subset(run2_final, subset = CT == "Plasma" | CT == "pDC" | 
                        CT == "cDC" | CT == "Monocytes")
saveRDS(mye_cleanup, "/scratch/smallapragada/run2/run2_mye_CTs_cleanup_2024_08_15.rds")

## Loading in reclustered object post-RAPIDS
mye_cleanup <- readRDS("/scratch/smallapragada/rapids_pipeline/run2_mye_CTs_cleanup_2024_08_15.rds")

## UMAP visualization of clusters in reclustered object
DimPlot(mye_cleanup,
        reduction = "umap",
        group.by = "mye_cleanup_brl_0.4", 
        raster = T,
        label = T, 
        cols = randomcoloR::distinctColorPalette(15)) 

## UMAP visualization of phenotype in reclustered object
DimPlot(mye_cleanup,
        reduction = "umap",
        group.by = "Phenotype", 
        raster = T,
        label = F, 
        cols = randomcoloR::distinctColorPalette(5)) 

## A summary of how many cells per cluster belong to a specific disease phenotype
table(mye_cleanup$mye_cleanup_brl_0.4, mye_cleanup$Phenotype)

## Looking at the nCount_RNA and the nFeature_RNA distribution across clusters to determine low-quality clusters
ggplot(mye_cleanup@meta.data, aes(factor(mye_cleanup@meta.data$mye_cleanup_brl_0.4), 
                                  (mye_cleanup@meta.data$nCount_RNA))) + geom_violin()
ggplot(mye_cleanup@meta.data, aes(factor(mye_cleanup@meta.data$mye_cleanup_brl_0.4), 
                                  (mye_cleanup@meta.data$nFeature_RNA))) + geom_violin()

## Top markers
Idents(mye_cleanup) = "mye_cleanup_brl_0.4"
markers <- FindAllMarkers(mye_cleanup)
markers %>%
  mutate(pct.diff = abs(pct.1 - pct.2)) %>%
  filter(avg_log2FC >= 0) %>%
  View()

DotPlot(mye_cleanup, features = c(other_features, immune_features)) + angle

## Cluster 0 - Venous

obj <- subset(mye_cleanup, subset = mye_cleanup_brl_0.4 == 0)

FeaturePlot(obj, features = c("ACKR1", "CCN2", "ICAM1", "CCL2", "WWTR1", "SPARCL1", 
                              "SOD2", "KDR", "BMPR2", "ITGB1", "ITGAV", "POSTN"))

DotPlot(obj, features = c("ACKR1", "CCN2", "ICAM1", "CCL2", "WWTR1", "SPARCL1", 
                          "SOD2", "KDR", "BMPR2", "ITGB1", "ITGAV", "POSTN"))

## Cluster 1 - Monocytes

obj <- subset(mye_cleanup, subset = mye_cleanup_brl_0.4 == 1)

FeaturePlot(obj, features = c("FCN1", "LYZ", "AIF1", "PTPRC", "CD68", "FCER1G", 
                              "CD14", "MS4A7", "CD44", "S100A8", "S100A9", "FCGR3A", 
                              "SPP1")) 

DotPlot(obj, features = c("FCN1", "LYZ", "AIF1", "PTPRC", "CD68", "FCER1G", 
                          "CD14", "MS4A7", "CD44", "S100A8", "S100A9", "FCGR3A", 
                          "SPP1")) + angle

## Cluster 2 - Transitional AT2

obj <- subset(mye_cleanup, subset = mye_cleanup_brl_0.4 == 2)

FeaturePlot(obj, features = c("NAPSA", "SFTPC", "SFTPD", "LAMP3", "KRT8", 
                              "KRT18", "SCGB3A2", "KRT5", "KRT6A"))

DotPlot(obj, features = c("NAPSA", "SFTPC", "SFTPD", "LAMP3", "KRT8", 
                          "KRT18", "SCGB3A2", "KRT5", "KRT6A")) + angle

## Cluster 3 - Monocytes

obj <- subset(mye_cleanup, subset = mye_cleanup_brl_0.4 == 3)

FeaturePlot(obj, features = c("FCGR3A", "SLC25A37", "ITGAX", "LPAR2", "CD14", 
                              "AIF1", "IL1B", "CD68", "BCL2L11", "S100A8", "S100A9"))

DotPlot(obj, features = c("FCGR3A", "SLC25A37", "ITGAX", "LPAR2", "CD14", 
                          "AIF1", "IL1B", "CD68", "BCL2L11", "S100A8", "S100A9")) + angle

## Cluster 4 - Proliferating monocytes/neutrophils (Subcluster)

obj <- subset(mye_cleanup, subset = mye_cleanup_brl_0.4 == 4)

FeaturePlot(obj, features = c("MKI67", "TOP2A", "RETN", "LYZ", "ELANE", "PKM", 
                              "S100A9", "S100A8", "CEACAM6", "CCNB2", "PCNA", 
                              "HIST1H1C"))

DotPlot(obj, features = c("MKI67", "TOP2A", "RETN", "LYZ", "ELANE", "PKM", 
                          "S100A9", "S100A8", "CEACAM6", "CCNB2", "PCNA", 
                          "HIST1H1C")) + angle

## Cluster 5 - Capillaries

obj <- subset(mye_cleanup, subset = mye_cleanup_brl_0.4 == 5)

FeaturePlot(obj, features = c("CD34", "KDR", "FCN3", "GNG11", "BMPR2", "CA4", 
                              "CLDN5", "APLNR", "IL7R", "RAMP2", "ACKR1", "HEY1"))

DotPlot(obj, features = c("CD34", "KDR", "FCN3", "GNG11", "BMPR2", "CA4", 
                          "CLDN5", "APLNR", "IL7R", "RAMP2", "ACKR1", "HEY1")) + angle


## Cluster 6 - SPP1+ macrophages/cDC (Subcluster)

obj <- subset(mye_cleanup, subset = mye_cleanup_brl_0.4 == 6)

FeaturePlot(obj, features = c("HLA-DRA", "MRC1", "MS4A7", "HLA-DQA1", "GPR183", 
                              "CD68", "C1QC", "HERPUD1", "HLA-DQB1", "CD14", 
                              "CD4", "FCER1G", "CD86", "SPP1"))

DotPlot(obj, features = c("HLA-DRA", "MRC1", "MS4A7", "HLA-DQA1", "GPR183", 
                          "CD68", "C1QC", "HERPUD1", "HLA-DQB1", "CD14", 
                          "CD4", "FCER1G", "CD86", "SPP1")) + angle

## Cluster 7 - NK & NKT

obj <- subset(mye_cleanup, subset = mye_cleanup_brl_0.4 == 7)

FeaturePlot(obj, features = c("CD247", "GZMA", "NKG7", "GNLY", "KLRC1", "KLRB1"))

DotPlot(obj, features = c("CD247", "GZMA", "NKG7", "GNLY", "KLRC1", "KLRB1"))

## Cluster 8 - Meg-ery

obj <- subset(mye_cleanup, subset = mye_cleanup_brl_0.4 == 8)

FeaturePlot(obj, features = c("ATP2A3", "MEG3", "TGFB1", "BCL2L1", "SNCA", 
                              "TBXA2R", "TOP1", "PECAM1", "COL1A2"))

DotPlot(obj, features = c("ATP2A3", "MEG3", "TGFB1", "BCL2L1", "SNCA", 
                          "TBXA2R", "TOP1", "PECAM1", "COL1A2")) + angle

## Cluster 9 - FB

obj <- subset(mye_cleanup, subset = mye_cleanup_brl_0.4 == 9)

FeaturePlot(obj, features = c("MEG3", "COL1A2", "PDGFRA", "COL3A1", "LUM", 
                              "PDGFRB", "CSPG4", "ACTA2", "PLIN2", "HAS2", 
                              "WNT5A"))

DotPlot(obj, features = c("MEG3", "COL1A2", "PDGFRA", "COL3A1", "LUM", 
                          "PDGFRB", "CSPG4", "ACTA2", "PLIN2", "HAS2", 
                          "WNT5A")) + angle

## Cluster 10 - pDC/B/Plasma

obj <- subset(mye_cleanup, subset = mye_cleanup_brl_0.4 == 10)

FeaturePlot(obj, features = c("JCHAIN", "XBP1", "GZMB", "PTGDS", "SELENOS", "LILRA4", "GPR183", 
                              "PIM2", "SEC11C", "LMAN1", "SSR3", "HERPUD1", "PRDX4", "UBE2J1", 
                              "S100A8", "S100A9", "FCER1G", "LYZ", "AIF1", "FCGR3A",
                              "MS4A1", "BANK1", "TNFRSF13C", "CD19"))

DotPlot(obj, features = c("JCHAIN", "XBP1", "GZMB", "PTGDS", "SELENOS", "LILRA4", "GPR183", 
                          "PIM2", "SEC11C", "LMAN1", "SSR3", "HERPUD1", "PRDX4", "UBE2J1", 
                          "S100A8", "S100A9", "FCER1G", "LYZ", "AIF1", "FCGR3A",
                          "MS4A1", "BANK1", "TNFRSF13C", "CD19")) + angle

### Subclustering of clusters 4, 6, and 10

## Loading in reclustered object post-RAPIDS
mye_sc_cleanup <- readRDS("/scratch/smallapragada/rapids_pipeline/run2_mye_CTs_cleanup_sc_2024_08_15.rds")

## Cluster 4

## UMAP visualization of clusters in reclustered object
DimPlot(mye_sc_cleanup,
        reduction = "umap",
        group.by = "mye_cleanup_brl_0.4_c4_0.2", 
        raster = T,
        label = T, 
        cols = randomcoloR::distinctColorPalette(15)) 

## UMAP visualization of phenotype in reclustered object
DimPlot(mye_sc_cleanup,
        reduction = "umap",
        group.by = "Phenotype", 
        raster = T,
        label = F, 
        cols = randomcoloR::distinctColorPalette(5)) 

## A summary of how many cells per cluster belong to a specific disease phenotype
table(mye_sc_cleanup$mye_cleanup_brl_0.4_c4_0.2, mye_sc_cleanup$Phenotype)

## Looking at the nCount_RNA and the nFeature_RNA distribution across clusters to determine low-quality clusters
ggplot(mye_sc_cleanup@meta.data, aes(factor(mye_sc_cleanup@meta.data$mye_cleanup_brl_0.4_c4_0.2), 
                                     (mye_sc_cleanup@meta.data$nCount_RNA))) + geom_violin()
ggplot(mye_sc_cleanup@meta.data, aes(factor(mye_sc_cleanup@meta.data$mye_cleanup_brl_0.4_c4_0.2), 
                                     (mye_sc_cleanup@meta.data$nFeature_RNA))) + geom_violin()

obj <- subset(mye_sc_cleanup, subset = mye_cleanup_brl_0.4_c4_0.2 == "4,0" | mye_cleanup_brl_0.4_c4_0.2 == "4,1")

## Top markers
Idents(obj) = "mye_cleanup_brl_0.4_c4_0.2"
markers <- FindAllMarkers(obj)
markers %>%
  mutate(pct.diff = abs(pct.1 - pct.2)) %>%
  filter(avg_log2FC >= 0) %>%
  View()

FeaturePlot(obj, features = c("MKI67", "TOP2A", "RETN", "LYZ", "ELANE", "PKM", 
                              "S100A9", "S100A8", "CEACAM6", "CCNB2", "PCNA", 
                              "HIST1H1C"))

DotPlot(obj, features = c("MKI67", "TOP2A", "RETN", "LYZ", "ELANE", "PKM", 
                          "S100A9", "S100A8", "CEACAM6", "CCNB2", "PCNA", 
                          "HIST1H1C")) + angle

# 4,0 - Proliferating monocytes
# 4,1 - Neutrophils

## Cluster 6 

## UMAP visualization of clusters in reclustered object
DimPlot(mye_sc_cleanup,
        reduction = "umap",
        group.by = "mye_cleanup_brl_0.4_c6_0.1", 
        raster = T,
        label = T, 
        cols = randomcoloR::distinctColorPalette(15)) 

## UMAP visualization of phenotype in reclustered object
DimPlot(mye_sc_cleanup,
        reduction = "umap",
        group.by = "Phenotype", 
        raster = T,
        label = F, 
        cols = randomcoloR::distinctColorPalette(5)) 

## A summary of how many cells per cluster belong to a specific disease phenotype
table(mye_sc_cleanup$mye_cleanup_brl_0.4_c6_0.1, mye_sc_cleanup$Phenotype)

## Looking at the nCount_RNA and the nFeature_RNA distribution across clusters to determine low-quality clusters
ggplot(mye_sc_cleanup@meta.data, aes(factor(mye_sc_cleanup@meta.data$mye_cleanup_brl_0.4_c6_0.1), 
                                     (mye_sc_cleanup@meta.data$nCount_RNA))) + geom_violin()
ggplot(mye_sc_cleanup@meta.data, aes(factor(mye_sc_cleanup@meta.data$mye_cleanup_brl_0.4_c6_0.1), 
                                     (mye_sc_cleanup@meta.data$nFeature_RNA))) + geom_violin()

obj <- subset(mye_sc_cleanup, subset = mye_cleanup_brl_0.4_c6_0.1 == "6,0" | mye_cleanup_brl_0.4_c6_0.1 == "6,1")

## Top markers
Idents(obj) = "mye_cleanup_brl_0.4_c4_0.2"
markers <- FindAllMarkers(obj)
markers %>%
  mutate(pct.diff = abs(pct.1 - pct.2)) %>%
  filter(avg_log2FC >= 0) %>%
  View()

FeaturePlot(obj, features = c("HLA-DRA", "MRC1", "MS4A7", "HLA-DQA1", "GPR183", 
                              "CD68", "C1QC", "HERPUD1", "HLA-DQB1",
                              "CD14", "CD4", "FCER1G", "CD86", "SPP1"))

DotPlot(obj, features = c("HLA-DRA", "MRC1", "MS4A7", "HLA-DQA1", "GPR183", 
                          "CD68", "C1QC", "HERPUD1", "HLA-DQB1", "CD14", 
                          "CD4", "FCER1G", "CD86", "SPP1")) + angle

# 6,0 - SPP1+ macrophages
# 6,1 - cDC

## Cluster 10

## UMAP visualization of clusters in reclustered object
DimPlot(mye_sc_cleanup,
        reduction = "umap",
        group.by = "mye_cleanup_brl_0.4_c10_0.3", 
        raster = T,
        label = T, 
        cols = randomcoloR::distinctColorPalette(15)) 

## UMAP visualization of phenotype in reclustered object
DimPlot(mye_sc_cleanup,
        reduction = "umap",
        group.by = "Phenotype", 
        raster = T,
        label = F, 
        cols = randomcoloR::distinctColorPalette(5)) 

## A summary of how many cells per cluster belong to a specific disease phenotype
table(mye_sc_cleanup$mye_cleanup_brl_0.4_c10_0.3, mye_sc_cleanup$Phenotype)

## Looking at the nCount_RNA and the nFeature_RNA distribution across clusters to determine low-quality clusters
ggplot(mye_sc_cleanup@meta.data, aes(factor(mye_sc_cleanup@meta.data$mye_cleanup_brl_0.4_c10_0.3), 
                                     (mye_sc_cleanup@meta.data$nCount_RNA))) + geom_violin()
ggplot(mye_sc_cleanup@meta.data, aes(factor(mye_sc_cleanup@meta.data$mye_cleanup_brl_0.4_c10_0.3), 
                                     (mye_sc_cleanup@meta.data$nFeature_RNA))) + geom_violin()

obj <- subset(mye_sc_cleanup, subset = mye_cleanup_brl_0.4_c10_0.3 == "10,0" | 
                mye_cleanup_brl_0.4_c10_0.3 == "10,1" | 
                mye_cleanup_brl_0.4_c10_0.3 == "10,2" | 
                mye_cleanup_brl_0.4_c10_0.3 == "10,3")

## Top markers
Idents(obj) = "mye_cleanup_brl_0.4_c10_0.3"
markers <- FindAllMarkers(obj)
markers %>%
  mutate(pct.diff = abs(pct.1 - pct.2)) %>%
  filter(avg_log2FC >= 0) %>%
  View()

FeaturePlot(obj, features = c("JCHAIN", "XBP1", "GZMB", "PTGDS", "SELENOS", "LILRA4", "GPR183", 
                              "PIM2", "SEC11C", "LMAN1", "SSR3", "HERPUD1", "PRDX4", "UBE2J1", 
                              "S100A8", "S100A9", "FCER1G", "LYZ", "AIF1", "FCGR3A",
                              "MS4A1", "BANK1", "TNFRSF13C", "CD19"))

FeaturePlot(obj, features = c("JCHAIN", "XBP1", "GZMB", "PTGDS", "SELENOS", "LILRA4", "GPR183", 
                              "PIM2", "SEC11C", "LMAN1", "SSR3", "HERPUD1", "PRDX4", "UBE2J1", 
                              "S100A8", "S100A9", "FCER1G", "LYZ", "AIF1", "FCGR3A",
                              "MS4A1", "BANK1", "TNFRSF13C", "CD19")) + angle

# 10,0 - pDC
# 10,1 - Plasma
# 10,2 - Monocytes
# 10,3 - B

## Adding cell type labels to reclustered object
mye_sc_cleanup$CT_final[mye_sc_cleanup$mye_cleanup_brl_0.4 == 0] <- "Venous"
mye_sc_cleanup$CT_final[mye_sc_cleanup$mye_cleanup_brl_0.4 == 1] <- "Monocytes"
mye_sc_cleanup$CT_final[mye_sc_cleanup$mye_cleanup_brl_0.4 == 2] <- "Transitional AT2"
mye_sc_cleanup$CT_final[mye_sc_cleanup$mye_cleanup_brl_0.4 == 3] <- "Monocytes"
mye_sc_cleanup$CT_final[mye_sc_cleanup$mye_cleanup_brl_0.4_c4_0.2 == "4,0"] <- "Proliferating monocytes"
mye_sc_cleanup$CT_final[mye_sc_cleanup$mye_cleanup_brl_0.4_c4_0.2 == "4,1"] <- "Neutrophils"
mye_sc_cleanup$CT_final[mye_sc_cleanup$mye_cleanup_brl_0.4 == 5] <- "Capillaries"
mye_sc_cleanup$CT_final[mye_sc_cleanup$mye_cleanup_brl_0.4_c6_0.1 == "6,0"] <- "cDC"
mye_sc_cleanup$CT_final[mye_sc_cleanup$mye_cleanup_brl_0.4_c6_0.1 == "6,1"] <- "SPP1+ macrophages"
mye_sc_cleanup$CT_final[mye_sc_cleanup$mye_cleanup_brl_0.4 == 7] <- "NK & NKT"
mye_sc_cleanup$CT_final[mye_sc_cleanup$mye_cleanup_brl_0.4 == 8] <- "Meg-Ery"
mye_sc_cleanup$CT_final[mye_sc_cleanup$mye_cleanup_brl_0.4 == 9] <- "Capillaries"
mye_sc_cleanup$CT_final[mye_sc_cleanup$mye_cleanup_brl_0.4_c10_0.3 == "10,0"] <- "pDC"
mye_sc_cleanup$CT_final[mye_sc_cleanup$mye_cleanup_brl_0.4_c10_0.3 == "10,1"] <- "Plasma"
mye_sc_cleanup$CT_final[mye_sc_cleanup$mye_cleanup_brl_0.4_c10_0.3 == "10,2"] <- "Monocytes"
mye_sc_cleanup$CT_final[mye_sc_cleanup$mye_cleanup_brl_0.4_c10_0.3 == "10,3"] <- "B"

## Adding new labels to original object
run2_final$CT_final <- mye_sc_cleanup$CT_final

DotPlot(mye_sc_cleanup, group.by = "CT_final", features = c(immune_features, "SPP1")) + angle

### Reclustering alveolar macrophages to pull out SPP1+ macrophages

## Subsetting alveolar macrophages & saving for RAPIDS processing
macs_cleanup <- subset(run2_final, subset = CT == "Alveolar macrophages")
saveRDS(macs_cleanup, "/scratch/smallapragada/run2/run2_macs_cleanup_2024_08_15.rds")

## Loading in reclustered object post-RAPIDS
macs_cleanup <- readRDS("/scratch/smallapragada/rapids_pipeline/run2_macs_CTs_cleanup_2024_08_15.rds")

## UMAP visualization of clusters in reclustered object
DimPlot(macs_cleanup,
        reduction = "umap",
        group.by = "macs_cleanup_brl_0.5", 
        raster = T,
        label = T, 
        cols = randomcoloR::distinctColorPalette(15)) 

## UMAP visualization of phenotype in reclustered object
DimPlot(macs_cleanup,
        reduction = "umap",
        group.by = "Phenotype", 
        raster = T,
        label = T, 
        cols = randomcoloR::distinctColorPalette(15)) 

## A summary of how many cells per cluster belong to a specific disease phenotype
table(macs_cleanup$macs_cleanup_brl_0.5, macs_cleanup$Phenotype)

## Looking at the nCount_RNA and the nFeature_RNA distribution across clusters to determine low-quality clusters
ggplot(macs_cleanup@meta.data, aes(factor(macs_cleanup@meta.data$macs_cleanup_brl_0.5), 
                                   (macs_cleanup@meta.data$nCount_RNA))) + geom_violin()
ggplot(macs_cleanup@meta.data, aes(factor(macs_cleanup@meta.data$macs_cleanup_brl_0.5), 
                                   (macs_cleanup@meta.data$nFeature_RNA))) + geom_violin()

DotPlot(macs_cleanup, features = c(other_features, immune_features)) + angle

## Top markers
Idents(macs_cleanup) = "macs_cleanup_brl_0.5"
markers <- FindAllMarkers(macs_cleanup)
markers %>%
  mutate(pct.diff = abs(pct.1 - pct.2)) %>%
  filter(avg_log2FC >= 0) %>%
  View()

## Cluster 0 - Plasma

FeaturePlot(macs_cleanup, features = c("CXCR5", "TERT", "SPRY2", "CCNB2"))

DotPlot(macs_cleanup, features = c("CXCR5", "TERT", "SPRY2", "CCNB2"))

## Cluster 1 - SPP1+ macrophages

FeaturePlot(macs_cleanup, features = c("SPP1", "CCL2", "RNASE1", "HIF1A"))

DotPlot(macs_cleanup, features = c("SPP1", "CCL2", "RNASE1", "HIF1A"))

## Cluster 2 - Fibroblasts

FeaturePlot(macs_cleanup, features = c("COL1A2", "COL3A1", "POSTN", "LUM", 
                                       "DCN", "ACTA2", "WNT5A", "HAS2", "PLIN2"))

DotPlot(macs_cleanup, features = c("COL1A2", "COL3A1", "POSTN", "LUM", 
                                   "DCN", "ACTA2", "WNT5A", "HAS2", "PLIN2")) + angle

## Cluster 3 - Proliferating macrophages

FeaturePlot(macs_cleanup, features = c("MKI67", "TOP2A", "FABP4", "MARCO", 
                                       "S100A8", "S100A9"))

DotPlot(macs_cleanup, features = c("MKI67", "TOP2A", "FABP4", "MARCO", 
                                   "S100A8", "S100A9")) + angle

## Cluster 4 - Alveolar macrophages

FeaturePlot(macs_cleanup, features = c("HLA-DQA1", "ICAM1", "SOD2", "HLA-DQB1", 
                                       "FABP4", "CCL18"))

DotPlot(macs_cleanup, features = c("HLA-DQA1", "ICAM1", "SOD2", "HLA-DQB1", 
                                   "FABP4", "CCL18")) + angle

## Cluster 5 - Alveolar macrophages

FeaturePlot(macs_cleanup, features = c("MCEMP1", "CCL18", "FCGR3A", "MARCO", 
                                       "CD52", "PPARG"))

DotPlot(macs_cleanup, features = c("MCEMP1", "CCL18", "FCGR3A", "MARCO", 
                                   "CD52", "PPARG")) + angle

## Cluster 6 - Monocytes

FeaturePlot(macs_cleanup, features = c("IDH1", "SLC1A3", "PLIN2", "ITGAX", 
                                       "SLC7A11", "CD14", "S100A9"))

DotPlot(macs_cleanup, features = c("IDH1", "SLC1A3", "PLIN2", "ITGAX", 
                                   "SLC7A11", "CD14", "S100A9")) + angle

## Cluster 7 - Alveolar macrophages

FeaturePlot(macs_cleanup, features = c("FABP4", "PECAM1", "MCEMP1", "RETN", 
                                       "CST3", "AIF1"))

DotPlot(macs_cleanup, features = c("FABP4", "PECAM1", "MCEMP1", "RETN", 
                                   "CST3", "AIF1")) + angle

## Adding cell type labels to reclustered object
macs_cleanup$CT_final[macs_cleanup$macs_cleanup_brl_0.5 == 0] <- "Plasma"
macs_cleanup$CT_final[macs_cleanup$macs_cleanup_brl_0.5 == 1] <- "SPP1+ macrophages"
macs_cleanup$CT_final[macs_cleanup$macs_cleanup_brl_0.5 == 2] <- "Capillaries"
macs_cleanup$CT_final[macs_cleanup$macs_cleanup_brl_0.5 == 3] <- "Proliferating macrophages"
macs_cleanup$CT_final[macs_cleanup$macs_cleanup_brl_0.5 == 4] <- "Alveolar macrophages"
macs_cleanup$CT_final[macs_cleanup$macs_cleanup_brl_0.5 == 5] <- "Alveolar macrophages"
macs_cleanup$CT_final[macs_cleanup$macs_cleanup_brl_0.5 == 6] <- "Monocytes"
macs_cleanup$CT_final[macs_cleanup$macs_cleanup_brl_0.5 == 7] <- "Alveolar macrophages"

## Adding new labels to original object
run2_final$CT_final <- macs_cleanup$CT_final

### Reclustering fibroblasts to separate them out better

## Subsetting all fibroblast cell types & saving for RAPIDS processing
fib_cleanup <- subset(run2_final, subset = CT == "Activated fibroblasts" | 
                        CT == "Adventitial fibroblasts" |
                        CT == "Fetal-enriched fibroblasts" | 
                        CT == "Fetal-enriched SOD2+ myoFB" | 
                        CT == "Fibroblasts" | CT == "PLIN2+ fibroblasts" | 
                        CT == "Proliferating fibroblasts")
saveRDS(fib_cleanup, "/scratch/smallapragada/run2/run2_fib_CTs_cleanup_2024_08_15.rds")

## Loading in reclustered object post-RAPIDS
fib_cleanup <- readRDS("/scratch/smallapragada/rapids_pipeline/run2_fib_CTs_cleanup_2024_08_15.rds") 

## UMAP visualization of clusters in reclustered object
DimPlot(fib_cleanup,
        reduction = "umap",
        group.by = "fib_cleanup_brl_0.4", 
        raster = T,
        label = T, 
        cols = randomcoloR::distinctColorPalette(15)) 

## UMAP visualization of phenotype in reclustered object
DimPlot(fib_cleanup,
        reduction = "umap",
        group.by = "Phenotype", 
        raster = T,
        label = F, 
        cols = randomcoloR::distinctColorPalette(5)) 

## A summary of how many cells per cluster belong to a specific disease phenotype
table(fib_cleanup$fib_cleanup_brl_0.4, fib_cleanup$Phenotype)

## Looking at the nCount_RNA and the nFeature_RNA distribution across clusters to determine low-quality clusters
ggplot(fib_cleanup@meta.data, aes(factor(fib_cleanup@meta.data$fib_cleanup_brl_0.4), 
                                  (fib_cleanup@meta.data$nCount_RNA))) + geom_violin()
ggplot(fib_cleanup@meta.data, aes(factor(fib_cleanup@meta.data$fib_cleanup_brl_0.4), 
                                  (fib_cleanup@meta.data$nFeature_RNA))) + geom_violin()

DotPlot(fib_cleanup, features = c(other_features, mesenchymal_features)) + angle

## Top markers
Idents(fib_cleanup) = "fib_cleanup_brl_0.4"
markers <- FindAllMarkers(fib_cleanup)
markers %>%
  mutate(pct.diff = abs(pct.1 - pct.2)) %>%
  filter(avg_log2FC >= 0) %>%
  View()

## Cluster 0 - Activated FB

obj <- subset(fib_cleanup, subset = fib_cleanup_brl_0.4 == 0)

FeaturePlot(obj, features = c("PTGDS", "SFRP2", "CTHRC1", "DCN", "LUM"))

DotPlot(obj, features = c("PTGDS", "SFRP2", "CTHRC1", "DCN", "LUM"))

## Cluster 1 - HAS2+ FB/myoFB (Subcluster)

obj <- subset(fib_cleanup, subset = fib_cleanup_brl_0.4 == 1)

FeaturePlot(obj, features = c("POSTN", "DCN", "ACTA2", "PDIA3", "LUM", "WNT5A", 
                              "FGF18", "COL3A1", "COL1A2", "PLIN2", "HAS2", 
                              "PDGFRA", "CXCL14"))

DotPlot(obj, features = c("POSTN", "DCN", "ACTA2", "PDIA3", "LUM", "WNT5A", 
                          "FGF18", "COL3A1", "COL1A2", "PLIN2", "HAS2", 
                          "PDGFRA", "CXCL14")) + angle

# Plot to spatially visualize where certain cells are located
samp <- subset(fib_cleanup, subset = Sample == "s1_PDL001")

DimPlot(samp,
        reduction = "SP",
        group.by = "fib_cleanup_brl_0.4", 
        raster = T,
        label = F,
        cols = c(`1` = "green"), na.value = "red") + coord_equal()

## Cluster 2 - HAS2+ FB

obj <- subset(fib_cleanup, subset = fib_cleanup_brl_0.4 == 2)

FeaturePlot(obj, features = c("SNAI2", "WNT2", "MEG3", "HAS2", "ITM2C", "FN1", 
                              "FGF10", "CD4", "VEGFA", "BMP4", "SNCA", "BCL2L11"))

DotPlot(obj, features = c("SNAI2", "WNT2", "MEG3", "HAS2", "ITM2C", "FN1", 
                          "FGF10", "CD4", "VEGFA", "BMP4", "SNCA", "BCL2L11")) + angle

## Cluster 3 - Transitional AT2

obj <- subset(fib_cleanup, subset = fib_cleanup_brl_0.4 == 3)

FeaturePlot(obj, features = c("LAMP3", "SFTPC", "NAPSA", "MSLN", "SCGB3A2", 
                              "S100A8", "KRT8"))

DotPlot(obj, features = c("LAMP3", "SFTPC", "NAPSA", "MSLN", "SCGB3A2", 
                          "S100A8", "KRT8")) 

## Cluster 4 - Fetal-enriched FB

obj <- subset(fib_cleanup, subset = fib_cleanup_brl_0.4 == 4)

FeaturePlot(obj, features = c("CCN2", "STAT1", "CXCL9", "SOD2", "CXCL14", "IFIT3", 
                              "ZEB1", "WNT5A", "PDGFRA", "ACTA2"))

DotPlot(obj, features = c("CCN2", "STAT1", "CXCL9", "SOD2", "CXCL14", "IFIT3", 
                          "ZEB1", "WNT5A", "PDGFRA", "ACTA2")) + angle

## Cluster 5 - Pericytes

obj <- subset(fib_cleanup, subset = fib_cleanup_brl_0.4 == 5)

FeaturePlot(obj, features = c("CD34", "PDGFRB", "ACTA2", "CSPG4", "CA4", "APLNR", 
                              "CLDN5", "GNG11", "PDGFRA", "WNT5A"))

DotPlot(obj, features = c("CD34", "PDGFRB", "ACTA2", "CSPG4", "CA4", "APLNR", 
                          "CLDN5", "GNG11", "PDGFRA", "WNT5A")) + angle

## Cluster 6 - Mast

obj <- subset(fib_cleanup, subset = fib_cleanup_brl_0.4 == 6)

FeaturePlot(obj, features = c("CPA3", "TPSAB1"))

DotPlot(obj, features = c("CPA3", "TPSAB1"))

## Cluster 7 - Fetal-enriched FB/Adventitial FB (Subcluster)

obj <- subset(fib_cleanup, subset = fib_cleanup_brl_0.4 == 7)

FeaturePlot(obj, features = c("MFAP5", "PI16", "DCN", "LUM", "COL2A1", "COL3A1", 
                              "SOX9", "SOD2"))

DotPlot(obj, features = c("MFAP5", "PI16", "DCN", "LUM", "COL2A1", "COL3A1", 
                          "SOX9", "SOD2")) + angle

## Cluster 8 - HAS2+ FB

obj <- subset(fib_cleanup, subset = fib_cleanup_brl_0.4 == 8)

FeaturePlot(obj, features = c("TOP2A", "MKI67", "HAS2", "LUM", "DCN", "SOX4", "FAP"))

DotPlot(obj, features = c("TOP2A", "MKI67", "HAS2", "LUM", "DCN", "SOX4", "FAP"))

## Cluster 9 - HAS2+ FB

obj <- subset(fib_cleanup, subset = fib_cleanup_brl_0.4 == 9)

FeaturePlot(obj, features = c("CCL2", "CD44", "S100A8", "LUM", "DCN", "LYZ", 
                              "S100A9", "FAS", "PLIN2", "HAS2", "HAS1"))

DotPlot(obj, features = c("CCL2", "CD44", "S100A8", "LUM", "DCN", "LYZ", 
                          "S100A9", "FAS", "PLIN2", "HAS2", "HAS1")) + angle

## Cluster 10 - SPP1+ macrophages

obj <- subset(fib_cleanup, subset = fib_cleanup_brl_0.4 == 10)

FeaturePlot(obj, features = c("FCER1G", "SPP1", "MS4A7", "MRC1", "AIF1", "GPR183", 
                              "C1QC", "S100A8", "MARCO"))

DotPlot(obj, features = c("FCER1G", "SPP1", "MS4A7", "MRC1", "AIF1", "GPR183", 
                          "C1QC", "S100A8", "MARCO")) + angle

## Cluster 11 - Fetal-enriched transitional AT2

obj <- subset(fib_cleanup, subset = fib_cleanup_brl_0.4 == 11)

FeaturePlot(obj, features = c("LAMP3", "DMBT1", "SFTPC", "KRT8", "NAPSA", "SCGB3A2"))

DotPlot(obj, features = c("LAMP3", "DMBT1", "SFTPC", "KRT8", "NAPSA", "SCGB3A2"))

## Cluster 12 - Fetal-enriched FB

obj <- subset(fib_cleanup, subset = fib_cleanup_brl_0.4 == 12)

FeaturePlot(obj, features = c("PDGFRA", "SOD2", "CD44", "HIF1A", "AXIN2", "HSPA5", 
                              "SPCS2", "LUM", "DCN", "WNT5A", "ACTA2"))

DotPlot(obj, features = c("PDGFRA", "SOD2", "CD44", "HIF1A", "AXIN2", "HSPA5", 
                          "SPCS2", "LUM", "DCN", "WNT5A", "ACTA2")) + angle

## Cluster 13 - Capillaries

obj <- subset(fib_cleanup, subset = fib_cleanup_brl_0.4 == 13)

FeaturePlot(obj, features = c("KDR", "CD34", "CA4", "KIT", "GNG11", "ACKR1", "HEY1"))

DotPlot(obj, features = c("KDR", "CD34", "CA4", "KIT", "GNG11", "ACKR1", "HEY1"))

### Subclustering of clusters 1 and 7

## Cluster 1

## Loading in reclustered object post-RAPIDS
fib_sc_cleanup <- readRDS("/scratch/smallapragada/rapids_pipeline/run2_fib_CTs_cleanup_sc_2024_08_15.rds") 

## UMAP visualization of clusters in reclustered object
DimPlot(fib_sc_cleanup,
        reduction = "umap",
        group.by = "fib_cleanup_brl_0.4_c1_0.3", 
        raster = T,
        label = T, 
        cols = randomcoloR::distinctColorPalette(15)) 

## UMAP visualization of phenotype in reclustered object
DimPlot(fib_sc_cleanup,
        reduction = "umap",
        group.by = "Phenotype", 
        raster = T,
        label = F, 
        cols = randomcoloR::distinctColorPalette(5)) 

## A summary of how many cells per cluster belong to a specific disease phenotype
table(fib_sc_cleanup$fib_cleanup_brl_0.4_c1_0.3, fib_sc_cleanup$Phenotype)

## Looking at the nCount_RNA and the nFeature_RNA distribution across clusters to determine low-quality clusters
ggplot(fib_sc_cleanup@meta.data, aes(factor(fib_sc_cleanup@meta.data$fib_cleanup_brl_0.4_c1_0.3), 
                                     (fib_sc_cleanup@meta.data$nCount_RNA))) + geom_violin()
ggplot(fib_sc_cleanup@meta.data, aes(factor(fib_sc_cleanup@meta.data$fib_cleanup_brl_0.4_c1_0.3), 
                                     (fib_sc_cleanup@meta.data$nFeature_RNA))) + geom_violin()

## Temp object for the specific clusters of interest
obj <- subset(fib_sc_cleanup, subset = fib_cleanup_brl_0.4_c1_0.3 == "1,0" | fib_cleanup_brl_0.4_c1_0.3 == "1,1")

## Top markers
Idents(obj) = "fib_cleanup_brl_0.4_c1_0.3"
markers <- FindAllMarkers(obj)
markers %>%
  mutate(pct.diff = abs(pct.1 - pct.2)) %>%
  filter(avg_log2FC >= 0) %>%
  View()

FeaturePlot(obj, features = c("POSTN", "DCN", "ACTA2", "PDIA3", "LUM", "WNT5A", 
                              "FGF18", "COL3A1", "COL1A2", "PLIN2", "HAS2", 
                              "PDGFRA", "CXCL14"))

DotPlot(obj, features = c("POSTN", "DCN", "ACTA2", "PDIA3", "LUM", "WNT5A", 
                          "FGF18", "COL3A1", "COL1A2", "PLIN2", "HAS2", 
                          "PDGFRA", "CXCL14")) + angle

# Cluster 1,0 - HAS2+ FB
# Cluster 1,1 - MyoFB

## Cluster 7

## UMAP visualization of clusters in reclustered object
DimPlot(fib_sc_cleanup,
        reduction = "umap",
        group.by = "fib_cleanup_brl_0.4_c7_0.1", 
        raster = T,
        label = T, 
        cols = randomcoloR::distinctColorPalette(15)) 

## UMAP visualization of phenotype in reclustered object
DimPlot(fib_sc_cleanup,
        reduction = "umap",
        group.by = "Phenotype", 
        raster = T,
        label = F, 
        cols = randomcoloR::distinctColorPalette(5)) 

## A summary of how many cells per cluster belong to a specific disease phenotype
table(fib_sc_cleanup$fib_cleanup_brl_0.4_c7_0.1, fib_sc_cleanup$Phenotype)

## Looking at the nCount_RNA and the nFeature_RNA distribution across clusters to determine low-quality clusters
ggplot(fib_sc_cleanup@meta.data, aes(factor(fib_sc_cleanup@meta.data$fib_cleanup_brl_0.4_c7_0.1), 
                                     (fib_sc_cleanup@meta.data$nCount_RNA))) + geom_violin()
ggplot(fib_sc_cleanup@meta.data, aes(factor(fib_sc_cleanup@meta.data$fib_cleanup_brl_0.4_c7_0.1), 
                                     (fib_sc_cleanup@meta.data$nFeature_RNA))) + geom_violin()

## Temp object for the specific clusters of interest
obj <- subset(fib_sc_cleanup, subset = fib_cleanup_brl_0.4_c7_0.1 == "7,0" | fib_cleanup_brl_0.4_c7_0.1 == "7,1")

## Top markers
Idents(obj) = "fib_cleanup_brl_0.4_c7_0.1"
markers <- FindAllMarkers(obj)
markers %>%
  mutate(pct.diff = abs(pct.1 - pct.2)) %>%
  filter(avg_log2FC >= 0) %>%
  View()

FeaturePlot(obj, features = c("MFAP5", "PI16", "DCN", "LUM", "COL2A1", "COL3A1", 
                              "SOX9", "SOD2"))

DotPlot(obj, features = c("MFAP5", "PI16", "DCN", "LUM", "COL2A1", "COL3A1", 
                          "SOX9", "SOD2")) + angle

# 7,0 - Fetal-enriched FB
# 7,1 - Adventitial FB

## Adding cell type labels to reclustered object
fib_sc_cleanup$CT_final[fib_sc_cleanup$fib_cleanup_brl_0.4 == 0] <- "Activated fibroblasts"
fib_sc_cleanup$CT_final[fib_sc_cleanup$fib_cleanup_brl_0.4_c1_0.3 == "1,0"] <- "HAS2+ fibroblasts"
fib_sc_cleanup$CT_final[fib_sc_cleanup$fib_cleanup_brl_0.4_c1_0.3 == "1,1"] <- "MyoFB"
fib_sc_cleanup$CT_final[fib_sc_cleanup$fib_cleanup_brl_0.4 == 2] <- "HAS2+ fibroblasts"
fib_sc_cleanup$CT_final[fib_sc_cleanup$fib_cleanup_brl_0.4 == 3] <- "Transitional AT2"
fib_sc_cleanup$CT_final[fib_sc_cleanup$fib_cleanup_brl_0.4 == 4] <- "Fetal-enriched fibroblasts"
fib_sc_cleanup$CT_final[fib_sc_cleanup$fib_cleanup_brl_0.4 == 5] <- "Pericytes"
fib_sc_cleanup$CT_final[fib_sc_cleanup$fib_cleanup_brl_0.4 == 6] <- "Mast"
fib_sc_cleanup$CT_final[fib_sc_cleanup$fib_cleanup_brl_0.4_c7_0.1 == "7,0"] <- "Fetal-enriched fibroblasts"
fib_sc_cleanup$CT_final[fib_sc_cleanup$fib_cleanup_brl_0.4_c7_0.1 == "7,1"] <- "Adventitial fibroblasts"
fib_sc_cleanup$CT_final[fib_sc_cleanup$fib_cleanup_brl_0.4 == 8] <- "HAS2+ fibroblasts"
fib_sc_cleanup$CT_final[fib_sc_cleanup$fib_cleanup_brl_0.4 == 9] <- "HAS2+ fibroblasts"
fib_sc_cleanup$CT_final[fib_sc_cleanup$fib_cleanup_brl_0.4 == 10] <- "SPP1+ macrophages"
fib_sc_cleanup$CT_final[fib_sc_cleanup$fib_cleanup_brl_0.4 == 11] <- "Fetal-enriched transitional AT2"
fib_sc_cleanup$CT_final[fib_sc_cleanup$fib_cleanup_brl_0.4 == 12] <- "Fetal-enriched fibroblasts"
fib_sc_cleanup$CT_final[fib_sc_cleanup$fib_cleanup_brl_0.4 == 13] <- "Capillaries"

## Adding new labels to original object
run2_final$CT_final <- fib_sc_cleanup$CT_final

## Filtering out low quality cells
run2_final <- subset(run2_final, subset = CT_final != "Low quality")

## Updated CT label visualization

plot <- DimPlot(run2_final,
                reduction = "umap",
                group.by = "CT_final",
                raster = T,
                label = F, 
                cols = randomcoloR::distinctColorPalette(50)) + NoLegend()

LabelClusters(plot, id = "CT_final", repel = T, max.overlaps = Inf) 

## DotPlot of marker genes 
DotPlot(run2_final, group.by = "CT_final", features = c("SFRP2", "MFAP5", "AKR1C1", "MARCO", "HEY1", "AGER", "SFTPC",
                                                        "BANK1", "KRT5", "CA4", "FCGR3A", "SOX9", "COL1A2", 
                                                        "SOD2", "HAS2", "CCL21", "CPA3", 
                                                        "SLC25A37", "S100A8", "C20orf85", "PDGFRA", "ELANE", "NKG7", "JCHAIN", 
                                                        "PDGFRB", "PIM2", "CALCA", "MKI67", "SCGB3A2", "PECAM1", "PTPRC", "SCGB1A1",
                                                        "MUC5B", "ACTA2", "SPP1", "TRAC", "KRT8", "CTLA4", "ACKR1")) + angle

## Adding lineage labels
run2_final$Lineage <- ""
run2_final$Lineage[run2_final$CT_final %in% c("AKR1C1+ & AKR1C2+", "AT1", "AT2", "Basal", 
                                              "Fetal-enriched AT2", "Fetal-enriched transitional AT2", 
                                              "Multiciliated", "PNEC", "Proliferating airway", 
                                              "Proliferating basal", "Proliferating Fetal-enriched alveolar", 
                                              "Secretory 3A2+ & 1A1+", "Secretory MUC5B+", 
                                              "Transitional AT2")] <- "Epithelial"
run2_final$Lineage[run2_final$CT_final %in% c("Arterial", "Lymphatic endothelial", 
                                              "Proliferating endothelial", "Venous", 
                                              "Capillaries")] <- "Endothelial"
run2_final$Lineage[run2_final$CT_final %in% c("Activated fibroblasts", "Adventitial fibroblasts", 
                                              "Fetal-enriched fibroblasts", 
                                              "HAS2+ fibroblasts", 
                                              "MyoFB", "Pericytes", "SMC")] <- "Mesenchymal"
run2_final$Lineage[run2_final$CT_final %in% c("B", "cDC", "Mast", "Meg-Ery", "Monocytes", 
                                              "Neutrophils", "NK & NKT", "pDC", "Plasma",
                                              "Proliferating lymphoid", "Proliferating macrophages", 
                                              "Proliferating Meg-Ery", "Proliferating monocytes",
                                              "SPP1+ macrophages", "Alveolar macrophages",
                                              "T", "Treg")] <- "Immune"

## Adding sublineage labels
run2_final$Sublineage <- ""
run2_final$Sublineage[run2_final$CT_final %in% c("Alveolar macrophages", "Proliferating macrophages", "Meg-Ery", "Mast", "Proliferating Meg-Ery",
                                                 "Monocytes", "Proliferating monocytes", "Neutrophils", "SPP1+ macrophages")] <- "Myeloid"
run2_final$Sublineage[run2_final$CT_final %in% c("pDC", "NK & NKT", "Treg", "Plasma", "cDC", "T", "B", "Proliferating lymphoid")] <- "Lymphoid"
run2_final$Sublineage[run2_final$CT_final %in% c("Activated fibroblasts", "MyoFB", "HAS2+ fibroblasts", "Adventitial fibroblasts", 
                                                 "Fetal-enriched fibroblasts", "SMC", "Pericytes")] <- "Mesenchymal"
run2_final$Sublineage[run2_final$CT_final %in% c("Arterial", "Proliferating endothelial", "Capillaries", "Venous", "Lymphatic endothelial")] <- "Endothelial"
run2_final$Sublineage[run2_final$CT_final %in% c("AT1", "AT2", "Transitional AT2", "Fetal-enriched AT2", "Fetal-enriched transitional AT2", 
                                                 "Proliferating Fetal-enriched alveolar")] <- "Alveolar"
run2_final$Sublineage[run2_final$CT_final %in% c("AKR1C1+ & AKR1C2+", "Proliferating airway", "Secretory 3A2+ & 1A1+", "PNEC", "Multiciliated",
                                                 "Basal", "Proliferating basal", "Secretory MUC5B+")] <- "Airway"

## Adding final metadata into object
run2_final$Gestational_age_weeks <- ""
run2_final$Gestational_age_weeks[run2_final$Sample_combined %in% "PDL001"] <- "16"
run2_final$Gestational_age_weeks[run2_final$Sample_combined %in% "PDL002"] <- "18"
run2_final$Gestational_age_weeks[run2_final$Sample_combined %in% "PDL003"] <- "21"
run2_final$Gestational_age_weeks[run2_final$Sample_combined %in% "PDL004"] <- "23"
run2_final$Gestational_age_weeks[run2_final$Sample_combined %in% c("PDL006", "PDL005")] <- "24"
run2_final$Gestational_age_weeks[run2_final$Sample_combined %in% "PDL007"] <- "26"
run2_final$Gestational_age_weeks[run2_final$Sample_combined %in% "PDL008"] <- "27"
run2_final$Gestational_age_weeks[run2_final$Sample_combined %in% "PDL009"] <- "28"
run2_final$Gestational_age_weeks[run2_final$Sample_combined %in% "PDL014"] <- "29"
run2_final$Gestational_age_weeks[run2_final$Sample_combined %in% "PDL015"] <- "33"
run2_final$Gestational_age_weeks[run2_final$Sample_combined %in% "PDL010"] <- "34"
run2_final$Gestational_age_weeks[run2_final$Sample_combined %in% c("PDL011", "PDL016")] <- "35"
run2_final$Gestational_age_weeks[run2_final$Sample_combined %in% "PDL013"] <- "36"
run2_final$Gestational_age_weeks[run2_final$Sample_combined %in% "PDL012"] <- "40"
run2_final$Gestational_age_weeks[run2_final$Sample_combined %in% "PDL017"] <- "28"

run2_final$Gestational_age_weeks_scaled <- ""
run2_final$Gestational_age_weeks_scaled[run2_final$Sample_combined %in% "PDL001"] <- "-4"
run2_final$Gestational_age_weeks_scaled[run2_final$Sample_combined %in% "PDL002"] <- "-2"
run2_final$Gestational_age_weeks_scaled[run2_final$Sample_combined %in% "PDL003"] <- "1"
run2_final$Gestational_age_weeks_scaled[run2_final$Sample_combined %in% "PDL004"] <- "3"
run2_final$Gestational_age_weeks_scaled[run2_final$Sample_combined %in% c("PDL006", "PDL005")] <- "4"
run2_final$Gestational_age_weeks_scaled[run2_final$Sample_combined %in% "PDL007"] <- "6"
run2_final$Gestational_age_weeks_scaled[run2_final$Sample_combined %in% "PDL008"] <- "7"
run2_final$Gestational_age_weeks_scaled[run2_final$Sample_combined %in% "PDL009"] <- "8"
run2_final$Gestational_age_weeks_scaled[run2_final$Sample_combined %in% "PDL014"] <- "9"
run2_final$Gestational_age_weeks_scaled[run2_final$Sample_combined %in% "PDL015"] <- "13"
run2_final$Gestational_age_weeks_scaled[run2_final$Sample_combined %in% "PDL010"] <- "14"
run2_final$Gestational_age_weeks_scaled[run2_final$Sample_combined %in% c("PDL011", "PDL016")] <- "15"
run2_final$Gestational_age_weeks_scaled[run2_final$Sample_combined %in% "PDL013"] <- "16"
run2_final$Gestational_age_weeks_scaled[run2_final$Sample_combined %in% "PDL012"] <- "20"
run2_final$Gestational_age_weeks_scaled[run2_final$Sample_combined %in% "PDL017"] <- "8"

run2_final$Life_span_weeks <- ""
run2_final$Life_span_weeks[run2_final$Sample_combined %in% "PDL001"] <- "0"
run2_final$Life_span_weeks[run2_final$Sample_combined %in% "PDL002"] <- "0"
run2_final$Life_span_weeks[run2_final$Sample_combined %in% "PDL003"] <- "0"
run2_final$Life_span_weeks[run2_final$Sample_combined %in% "PDL004"] <- "7"
run2_final$Life_span_weeks[run2_final$Sample_combined %in% "PDL006"] <- "0"
run2_final$Life_span_weeks[run2_final$Sample_combined %in% "PDL005"] <- "0"
run2_final$Life_span_weeks[run2_final$Sample_combined %in% "PDL007"] <- "2"
run2_final$Life_span_weeks[run2_final$Sample_combined %in% "PDL008"] <- "16"
run2_final$Life_span_weeks[run2_final$Sample_combined %in% "PDL009"] <- "3"
run2_final$Life_span_weeks[run2_final$Sample_combined %in% "PDL014"] <- "0"
run2_final$Life_span_weeks[run2_final$Sample_combined %in% "PDL015"] <- "0"
run2_final$Life_span_weeks[run2_final$Sample_combined %in% "PDL010"] <- "4"
run2_final$Life_span_weeks[run2_final$Sample_combined %in% "PDL011"] <- "0"
run2_final$Life_span_weeks[run2_final$Sample_combined %in% "PDL016"] <- "28"
run2_final$Life_span_weeks[run2_final$Sample_combined %in% "PDL013"] <- "0"
run2_final$Life_span_weeks[run2_final$Sample_combined %in% "PDL012"] <- "0"
run2_final$Life_span_weeks[run2_final$Sample_combined %in% "PDL017"] <- "1560"

## Adjusting certain cell type labels 
run2_final$CT_final <- gsub("Fetal-enriched", "Immature", run2_final$CT_final)

run2_final$CT_final <- gsub("HAS2\\+ fibroblasts", "Fibroblasts", run2_final$CT_final)

run2_final$CT_final <- gsub("Immature transitional AT2", "Transitional AT2", run2_final$CT_final)

run2_final$CT_final <- gsub("Proliferating macrophages", "Proliferating monocytes", run2_final$CT_final)

## Removing the "s" naming convention before sample name 
run2_final$Sample_combined <- gsub("s1_|s2_", "", run2_final$Sample)

## Organizing columns and removing extraneous information in final object
run2_final$X <- NULL
run2_final$ident <- NULL

run2_final@meta.data <- run2_final@meta.data %>%
  relocate(adj_x_centroid, .after = y_centroid) %>%
  relocate(adj_y_centroid, .after = adj_x_centroid) %>%
  relocate(Sample_combined, .after = Sample) %>%
  relocate(Gestational_age_weeks, .after = Phenotype) %>%
  relocate(Gestational_age_weeks_scaled, .after = Gestational_age_weeks) %>%
  relocate(Life_span_weeks, .after = Gestational_age_weeks_scaled)

## Adjusting some coordinates of different samples to line them up correctly

# s2_PDL007 & s1_PDL008
run2_final@meta.data <- run2_final@meta.data %>%
  mutate(adj_x_centroid = case_when(run2_final$Sample == "s1_PDL008" ~ (-run2_final$x_centroid),
                                    TRUE ~ run2_final$adj_x_centroid))

run2_final@meta.data <- run2_final@meta.data %>%
  mutate(adj_y_centroid = case_when(run2_final$Sample == "s2_PDL007" ~ (-run2_final$y_centroid + 25000),
                                    run2_final$Sample == "s1_PDL008" ~ (-run2_final$y_centroid + 14000),
                                    TRUE ~ run2_final$adj_y_centroid))

# Add in spatial information as dimension reduction objects
position_xy <- cbind(run2_final$adj_x_centroid, run2_final$adj_y_centroid)
row.names(position_xy) <- row.names(run2_final@meta.data)
colnames(position_xy) <- c("SP_1", "SP_2")
run2_final[["sp"]] <- CreateDimReducObject(embeddings = position_xy, key = "SP_",
                                           assay = DefaultAssay(run2_final))

# Spatial visualization of all samples from both slides
DimPlot(run2_final, 
        reduction = "sp",
        group.by = "Sample",
        raster = T,
        label = T,
        cols = randomcoloR::distinctColorPalette(34)) +
  NoLegend()

## Saving all CT labels
saveRDS(run2_final, "/scratch/smallapragada/run2/final_complete_CT_2024_09_30.rds")

## Adjusting coordinates for the last time

# s1_PDL007 & s2_PDL007
run2_final@meta.data <- run2_final@meta.data %>%
  mutate(adj_y_centroid = case_when(run2_final$Sample == "s1_PDL007" ~ (-run2_final$y_centroid + 14000),
                                    run2_final$Sample == "s2_PDL007" ~ (run2_final$y_centroid),
                                    TRUE ~ run2_final$adj_y_centroid))

# Add in spatial information as dimension reduction objects
position_xy <- cbind(run2_final$adj_x_centroid, run2_final$adj_y_centroid)
row.names(position_xy) <- row.names(run2_final@meta.data)
colnames(position_xy) <- c("SP_1", "SP_2")
run2_final[["sp"]] <- CreateDimReducObject(embeddings = position_xy, key = "SP_",
                                           assay = DefaultAssay(run2_final))

# Spatial visualization of all samples from both slides
DimPlot(run2_final, 
        reduction = "sp",
        group.by = "Sample",
        raster = T,
        label = T,
        cols = randomcoloR::distinctColorPalette(34)) +
  NoLegend()

## Re-adjudicating immature fibroblasts (grouping them into fibroblasts)
run2_final$CT_final[run2_final$CT_final %in% "Immature fibroblasts"] <- "Fibroblasts"

# Dev_stages and rare disease
run2_final@meta.data <- run2_final@meta.data %>%
  mutate(dev_stages = case_when(run2_final$Sample_combined %in% c("PDL001", "PDL002") ~ "Early canalicular",
                                run2_final$Sample_combined %in% c("PDL003", "PDL004", "PDL005", "PDL006") ~ "Late canalicular",
                                run2_final$Sample_combined %in% c("PDL007", "PDL008", "PDL009", "PDL010",  "PDL011") ~ "Saccular",
                                run2_final$Sample_combined %in% c("PDL014") ~ "Rare disease (CHAOS)",
                                run2_final$Sample_combined %in% c("PDL015") ~ "Rare disease (PH + Kidneys)",
                                run2_final$Sample_combined %in% c("PDL016") ~ "Rare disease (Infant BPD)",
                                run2_final$Sample_combined %in% c("PDL017") ~ "Rare disease (Adult BPD)",
                                TRUE ~ "Alveolar"))

## Saving all final CT labels
saveRDS(run2_final, "/scratch/smallapragada/run2/final_complete_CT_2025_03_26.rds")

# Adding DX score
run2_final@meta.data <- run2_final@meta.data %>%
  mutate(dx_score = case_when(run2_final$Sample_combined %in% c("PDL006", "PDL011", "PDL012", "PDL013") ~ 0,
                              run2_final$Sample_combined %in% c("PDL003", "PDL004", "PDL007") ~ 1,
                              run2_final$Sample_combined %in% c("PDL008", "PDL009") ~ 2,
                              run2_final$Sample_combined %in% c("PDL010", "PDL005") ~ 3,
                              TRUE ~ NA))

run2_final@meta.data <- run2_final@meta.data %>%
  relocate(dx_score, .after = Life_span_weeks) %>%
  relocate(dev_stages, .after = dx_score)

## FINAL FINAL OBJECT
saveRDS(run2_final, "/scratch/smallapragada/run2/final_complete_CT_updated_ids_2025_03_26.rds")
