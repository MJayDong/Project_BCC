#!/usr/bin/env Rscript
#
# 4.Visualization_PeakModel.r
#
# Purpose:
#   Visualize gkm-SVM model feature weights and annotate high-scoring peaks.
#
# Inputs:
#   config/config_paths.yaml
#   config/config_parameters.yaml
#   results/models/full_models_1000bp/*.model.txt
#   results/peaks/<celltype>/training_peaks.rds
#
# Outputs:
#   results/model_plots/<celltype>/*pdf
#
# ======================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(ggplot2)
  library(GenomicRanges)
  library(Biostrings)
  library(dplyr)
})

# ---------------------------------------------------------------
# 1. Load config
# ---------------------------------------------------------------
cfg_paths <- read_yaml("config/config_paths.yaml")
cfg_param <- read_yaml("config/config_parameters.yaml")

model_dir   <- cfg_paths$models$full_models
peak_dir    <- cfg_paths$input$training_peaks
outdir      <- file.path(cfg_paths$project_root, "results/model_plots")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------
# 2. Helper: Load gkm-SVM weight matrix
# ---------------------------------------------------------------
load_gkmsvm_weights <- function(model_file) {
  w <- read.table(model_file, header = FALSE, stringsAsFactors = FALSE)
  colnames(w) <- c("kmer", "weight")
  w <- w[order(-w$weight), ]
  return(w)
}

# ---------------------------------------------------------------
# 3. Helper: Plot Top-K kmers
# ---------------------------------------------------------------
plot_top_kmers <- function(weight_df, model_name, k = 30) {
  df <- head(weight_df, k)
  df$kmer <- factor(df$kmer, levels = df$kmer)

  ggplot(df, aes(x = kmer, y = weight)) +
    geom_bar(stat = "identity", fill = "#2874A6") +
    coord_flip() +
    theme_bw(base_size = 12) +
    labs(
      title = paste0("Top ", k, " gkm-SVM kmers: ", model_name),
      x = "k-mer",
      y = "Weight"
    )
}

# ---------------------------------------------------------------
# 4. Loop through cell-type models
# ---------------------------------------------------------------
models <- cfg_param$models$celltypes

for (model_name in models) {

  message("=========================================")
  message("Visualizing model: ", model_name)
  message("=========================================")

  # ---- 4.1 Identify model file ----
  model_file <- list.files(model_dir, pattern = paste0(model_name, ".*model.txt"), full.names = TRUE)
  if (length(model_file) == 0) {
    warning("No model file found for: ", model_name)
    next
  }
  model_file <- model_file[1]

  # ---- 4.2 Load weights ----
  w <- load_gkmsvm_weights(model_file)

  # ---- 4.3 Plot top-k k-mer weights ----
  p <- plot_top_kmers(w, model_name, k = cfg_param$visualization$top_k_kmers)

  pdf(file.path(outdir, paste0(model_name, "_topKmers.pdf")), width = 6, height = 8)
  print(p)
  dev.off()

  # ---- 4.4 Save weight table ----
  write.csv(w, file.path(outdir, paste0(model_name, "_weights.csv")), row.names = FALSE)
}

message("✓ Finished model visualization.")
