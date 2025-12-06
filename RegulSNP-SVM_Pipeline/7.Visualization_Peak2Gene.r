#!/usr/bin/env Rscript
#
# 6.Visualization_Peak2Gene.r
#
# Purpose:
#   Integrate deltaSVM scores with peak annotations and Peak2Gene links.
#   Generate SNP → Peak → Gene mappings and visualization.
#
# Inputs:
#   results/svm_summary/deltaSVM_summary.csv
#   results/peaks/<celltype>/training_peaks.rds
#   results/peak2gene/<celltype>/Peak2Gene.csv   (standard ArchR output)
#
# Outputs:
#   results/peak2gene_summary/peak2gene_mapping.csv
#   results/peak2gene_summary/peak2gene_network.pdf
#
# ======================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(igraph)
  library(ggraph)
  library(readr)
  library(GenomicRanges)
})

# ---------------------------------------------------------------
# 1. Load configs
# ---------------------------------------------------------------
cfg_paths <- read_yaml("config/config_paths.yaml")
cfg_param <- read_yaml("config/config_parameters.yaml")

svm_summary_file <- file.path(cfg_paths$project_root, "results/svm_summary/deltaSVM_summary.csv")
peaks_dir        <- cfg_paths$input$training_peaks
p2g_dir          <- cfg_paths$results$peak2gene
outdir           <- file.path(cfg_paths$project_root, "results/peak2gene_summary")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

models <- cfg_param$models$celltypes

# ---------------------------------------------------------------
# 2. Load deltaSVM matrix
# ---------------------------------------------------------------
message("Loading deltaSVM summary matrix...")

svm_mat <- read_csv(svm_summary_file, show_col_types = FALSE)

if (nrow(svm_mat) == 0) stop("deltaSVM summary is empty.")

# Long format for merging
svm_long <- svm_mat %>% pivot_longer(-SNP_ID, names_to = "model", values_to = "deltaSVM")

# ---------------------------------------------------------------
# 3. Helper: Convert "chr_start_end" → GRanges
# ---------------------------------------------------------------
to_gr <- function(x) {
  df <- tidyr::separate(x, SNP_ID, into = c("chr", "start", "end"), sep = "_", convert = TRUE)
  makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)
}

# ---------------------------------------------------------------
# 4. Build SNP GRanges for overlap with peaks
# ---------------------------------------------------------------
svm_gr <- to_gr(svm_long)

# ---------------------------------------------------------------
# 5. Initialize final result storage
# ---------------------------------------------------------------
all_results <- list()

# ---------------------------------------------------------------
# 6. Loop through cell types (models)
# ---------------------------------------------------------------
for (model_name in models) {

  message("---------------------------------------------------------")
  message("Processing Peak2Gene for model: ", model_name)
  message("---------------------------------------------------------")

  # ------------------------------------------
  # Load model-specific peaks
  # ------------------------------------------
  peak_file <- file.path(peaks_dir, paste0(model_name, "_training_peaks.rds"))
  if (!file.exists(peak_file)) {
    warning("Missing peaks for ", model_name)
    next
  }

  peaks <- readRDS(peak_file)
  seqlevelsStyle(peaks) <- "UCSC"
  peaks$peakName <- paste0(seqnames(peaks), "_", start(peaks), "_", end(peaks))

  # ------------------------------------------
  # Load Peak2Gene file
  # ------------------------------------------
  p2g_file <- file.path(p2g_dir, paste0(model_name, "_Peak2Gene.csv"))
  if (!file.exists(p2g_file)) {
    warning("Missing Peak2Gene: ", p2g_file)
    next
  }

  p2g <- read_csv(p2g_file, show_col_types = FALSE) %>%
    rename(peakName = Peak, gene = Gene)

  # ------------------------------------------
  # 1) Overlap SNPs with peaks
  # ------------------------------------------
  hits <- findOverlaps(svm_gr, peaks)
  if (length(hits) == 0) {
    message("No SNPs overlap peaks for model: ", model_name)
    next
  }

  snp_overlaps <- cbind(
    as.data.frame(svm_gr[queryHits(hits)]),
    peakName = peaks$peakName[subjectHits(hits)]
  )

  # Attach model-specific deltaSVM
  snp_overlaps <- snp_overlaps %>%
    left_join(
      svm_long %>% filter(model == model_name),
      by = c("SNP_ID")
    )

  # ------------------------------------------
  # 2) Attach Peak2Gene links
  # ------------------------------------------
  snp_p2g <- snp_overlaps %>%
    left_join(p2g, by = "peakName")

  all_results[[model_name]] <- snp_p2g
}

# ---------------------------------------------------------------
# 7. Combine all SNP–Peak–Gene mappings
# ---------------------------------------------------------------
final_map <- bind_rows(all_results)
write_csv(final_map, file.path(outdir, "peak2gene_mapping.csv"))
message("✓ SNP → Peak → Gene mapping saved.")

# ---------------------------------------------------------------
# 8. Network visualization
# ---------------------------------------------------------------
message("Generating network visualization...")

# Build graph
edges <- final_map %>%
  select(SNP_ID, gene, deltaSVM) %>%
  mutate(weight = abs(deltaSVM)) %>%
  drop_na()

g <- graph_from
