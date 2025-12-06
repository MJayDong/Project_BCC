#!/usr/bin/env Rscript
#
# ================================================================
# Script: generate_null_sequences.R
# Purpose:
#   Generate GC-matched background ("null") sequences for gkm-SVM.
#
# Input:
#   config/config_paths.yaml
#   config/config_parameters.yaml
#   results/peaks/<celltype>/training_peaks.rds
#
# Output:
#   results/null_sequences/<celltype>/null_sequences.fasta
#   results/null_sequences/<celltype>/null_regions.rds
#
# Core Methods:
#   1. Sample random genomic positions
#   2. Remove blacklist & true peak overlap
#   3. Match GC distribution to training peaks
#   4. Extract FASTA sequences
#
# Dependencies:
#   00_core/utils_randomSeq.R
#   BSgenome.Hsapiens.*
# ================================================================

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(Biostrings)
  library(parallel)
  library(yaml)
})

# ---------------------------------------------------------------
# 1. Load config
# ---------------------------------------------------------------
cfg_paths <- yaml::read_yaml("config/config_paths.yaml")
cfg_param <- yaml::read_yaml("config/config_parameters.yaml")

genome_id <- cfg_param$genome$build
peak_dir  <- cfg_paths$input$training_peaks
outdir    <- file.path(cfg_paths$project_root, "results/null_sequences")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------
# 2. Load genome
# ---------------------------------------------------------------
if (genome_id == "hg38") {
  library(BSgenome.Hsapiens.UCSC.hg38.masked)
  genome <- BSgenome.Hsapiens.UCSC.hg38.masked
} else if (genome_id == "hg19") {
  library(BSgenome.Hsapiens.UCSC.hg19.masked)
  genome <- BSgenome.Hsapiens.UCSC.hg19.masked
} else {
  stop("Unsupported genome build.")
}

# ---------------------------------------------------------------
# 3. Load helper functions
# ---------------------------------------------------------------
source("00_core/utils_randomSeq.R")   # getRandomPos, selectTargetSeqs
source("00_core/utils_gc.R")          # gcContent
source("00_core/utils_granges.R")     # trim_oob, trim_N_seqs

# ---------------------------------------------------------------
# 4. Load training peak set
# ---------------------------------------------------------------
peaks <- readRDS(file.path(peak_dir, "training_peaks.rds"))
seqlevelsStyle(peaks) <- "UCSC"

peak_width <- cfg_param$peaks$window_bp * 2 + 1
peaks <- resize(peaks, width = peak_width, fix="center")
peaks$GC <- gcContent(peaks, genome)

message("Training peaks loaded: ", length(peaks))

# ---------------------------------------------------------------
# 5. Build blacklist (from masked genome)
# ---------------------------------------------------------------
message("Constructing blacklist from masked genome...")

chroms <- seqlevels(peaks)
masked_gr <- GRanges()

for (chr in chroms) {
  masks <- Biostrings::masks(genome[[chr]])
  valid <- masks@NAMES[sapply(masks@nir_list, length) > 0]
  if (length(valid) == 0) next
  
  gr <- unlist(GRangesList(
    lapply(valid, function(m) GRanges(chr, as(masks[[m]], "IRanges")))
  ))
  masked_gr <- c(masked_gr, gr)
}

masked_gr <- reduce(masked_gr)
message("Blacklist regions: ", length(masked_gr))

# ---------------------------------------------------------------
# 6. Sample random genome windows
# ---------------------------------------------------------------
set.seed(cfg_param$random$seed)

n_null <- cfg_param$nullseq$n_per_celltype
message("Sampling initial random set (x5 oversampling)...")

random_raw <- getRandomPos(
  n = n_null * 5,
  genome = genome,
  use_chr = chroms,
  width = peak_width,
  blacklist_gr = masked_gr,
  non_overlapping = FALSE
)

random_raw$GC <- gcContent(random_raw, genome)

# ---------------------------------------------------------------
# 7. Select GC-matched null sequences
# ---------------------------------------------------------------
message("Selecting GC-matched null sequences...")

null_regions <- selectTargetSeqs(
  peaks,
  targets_gr = random_raw,
  blacklist_gr = masked_gr,
  nseqs = n_null,
  nbins = cfg_param$nullseq$gc_bins,
  bin_type = "size"
)

message("Final null sequences: ", length(null_regions))

# ---------------------------------------------------------------
# 8. Extract FASTA
# ---------------------------------------------------------------
null_fasta <- getSeq(genome, null_regions)
names(null_fasta) <- paste0("null_", seq_along(null_fasta))

fasta_file <- file.path(outdir, "null_sequences.fasta")
writeXStringSet(null_fasta, fasta_file)

saveRDS(null_regions, file.path(outdir, "null_regions.rds"))

message("✓ Null sequence generation complete.")
message("FASTA: ", fasta_file)
