#!/usr/bin/env Rscript
#
# 6.Visualization_SumSVMPrediction.r
#
# Purpose:
#   Summarize deltaSVM scores across models and generate QC/visualization.
#
# Inputs:
#   results/snp_scores/<model>.SNP.deltaSVM.txt
#   config/config_paths.yaml
#   config/config_parameters.yaml
#
# Output:
#   results/svm_summary/deltaSVM_summary.csv
#   results/svm_summary/deltaSVM_heatmap.pdf
#
# ======================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(pheatmap)
})

# ---------------------------------------------------------------
# 1. Load configs
# ---------------------------------------------------------------
cfg_paths <- read_yaml("config/config_paths.yaml")
cfg_param <- read_yaml("config/config_parameters.yaml")

snp_score_dir <- cfg_paths$results$snp_scores
outdir        <- file.path(cfg_paths$project_root, "results/svm_summary")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

models <- cfg_param$models$celltypes

# ---------------------------------------------------------------
# 2. Helper: Load one model’s deltaSVM file
# ---------------------------------------------------------------
load_delta_file <- function(model) {
  file <- file.path(snp_score_dir, paste0(model, ".SNP.deltaSVM.txt"))
  if (!file.exists(file)) return(NULL)
  
  df <- read.table(file, header = FALSE, stringsAsFactors = FALSE)
  colnames(df) <- c("SNP_ID", "REF", "ALT", "deltaSVM")
  df$model <- model
  return(df)
}

# ---------------------------------------------------------------
# 3. Load all deltaSVM outputs
# ---------------------------------------------------------------
message("Loading all deltaSVM predictions...")

all_scores <- lapply(models, load_delta_file)
all_scores <- bind_rows(all_scores)

if (nrow(all_scores) == 0) {
  stop("No deltaSVM files found. Check snp_scores directory.")
}

# ---------------------------------------------------------------
# 4. Summarize (wide format: SNP × model matrix)
# ---------------------------------------------------------------
delta_mat <- all_scores %>%
  select(SNP_ID, model, deltaSVM) %>%
  pivot_wider(names_from = model, values_from = deltaSVM)

write_csv(delta_mat, file.path(outdir, "deltaSVM_summary.csv"))

# ---------------------------------------------------------------
# 5. Heatmap visualization
# ---------------------------------------------------------------
message("Generating heatmap...")

mat <- as.matrix(delta_mat[,-1])
rownames(mat) <- delta_mat$SNP_ID

pdf(file.path(outdir, "deltaSVM_heatmap.pdf"), width = 8, height = 10)

pheatmap(
  mat,
  color = colorRampPalette(c("#2C3E50", "white", "#E74C3C"))(100),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize = 8,
  main = "deltaSVM Scores Across Models"
)

dev.off()

message("✓ deltaSVM summary & visualization complete.")
