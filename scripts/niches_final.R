#### LOAD PACKAGES AND SET ENVIRONMENT ----
library(SeuratObject)
library(tidyverse)
library(gplots)
library(rlist)
library(Seurat, lib.loc = "/usr/local/lib/R/site-library")

# Set seed & working directory
set.seed(0712)
work_dir <- setwd("/scratch/smallapragada/run2")

# Load data in
run2_final <- readRDS("/scratch/smallapragada/run2/final_complete_CT_updated_ids_2025_03_26.rds")

#### CELL-BASED NICHES ----

## Create FOV
coords <- run2_final@meta.data[, c("adj_x_centroid", "adj_y_centroid")]

coords$adj_x_centroid <- as.numeric(coords$adj_x_centroid)
coords$adj_y_centroid <- as.numeric(coords$adj_y_centroid)

run2_final[["fov"]]  <- CreateFOV(coords, assay = "RNA", type = "centroids")

## Niche analyses - n10 (all samples)

run2_final <- BuildNicheAssay(object = run2_final, fov = "fov", group.by = "CT_final",
                              niches.k = 11, neighbors.k = 10)

run2_final@meta.data <- run2_final@meta.data %>%
  mutate(CNiches = paste0("C", run2_final@meta.data$niches)) %>%
  select(!niches)

saveRDS(run2_final, "/scratch/smallapragada/run2/final_complete_cniches_allsamples_2025_03_26.rds")

## Plotting cell niches
DimPlot(run2_final,
        reduction = "SP",
        group.by = "CNiches",
        raster = T,
        label = F,
        cols = randomcoloR::distinctColorPalette(13))

#### TRANSCRIPT-BASED NICHES ----

## Pulling all .csv files per sample that contains the transcript niche information
transcript_files <- list.files(path = "/scratch/smallapragada/run2/updated_transcript_niche_files/", full.names = TRUE, recursive = FALSE)

## Parsing through files and fixing the format
transcripts <- lapply(transcript_files, function(XX)
  read.delim(XX, sep = ","))

## Binding all transcript files together into one
transcripts_all <- bind_rows(transcripts)

## Appending bound transcript file to seurat object
transcripts_gmm11 <- transcripts_all %>%
  select(cell_id, gmm11_5k_trained_hex_gmm)

## Create a named vector with cell IDs
tniche_values <- paste0("T", transcripts_gmm11$gmm11_5k_trained_hex_gmm)
names(tniche_values) <- transcripts_gmm11$cell_id 

## Assign values by matching on cell_id (rownames assumed to be cell IDs)
run2_final@meta.data$TNiche <- tniche_values[rownames(run2_final@meta.data)]

## Adding rownames back into object
rownames(run2_final@meta.data) <- run2_final@meta.data$cell_id

## Shift all the transcript niches up a number (they start with 0)
run2_final@meta.data <- run2_final@meta.data %>%
  mutate(TNiches = case_when(run2_final@meta.data$TNiche == "T0" ~ "T1", 
                             run2_final@meta.data$TNiche == "T1" ~ "T2",
                             run2_final@meta.data$TNiche == "T2" ~ "T3",
                             run2_final@meta.data$TNiche == "T3" ~ "T4",
                             run2_final@meta.data$TNiche == "T4" ~ "T5",
                             run2_final@meta.data$TNiche == "T5" ~ "T6",
                             run2_final@meta.data$TNiche == "T6" ~ "T7",
                             run2_final@meta.data$TNiche == "T7" ~ "T8",
                             run2_final@meta.data$TNiche == "T8" ~ "T9",
                             run2_final@meta.data$TNiche == "T9" ~ "T10",
                             run2_final@meta.data$TNiche == "T10" ~ "T11",
                             TRUE ~ NA_character_)) %>%
  select(!TNiche)

## Plotting transcript niches
DimPlot(run2_final,
        reduction = "SP",
        group.by = "TNiche",
        raster = T,
        label = F,
        cols = randomcoloR::distinctColorPalette(13))

## Save .RDS
saveRDS(run2_final, "/scratch/smallapragada/run2/final_complete_cniches_tniches_allsamples_2025_05_26.rds")
