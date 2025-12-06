#!/usr/bin/env Rscript
#
# ================================================================
# Script: 1.snp_seq_generator.R
# Purpose:
#   Generate REF/ALT SNP-centered FASTA sequences from a GRanges object
#   containing SNP coordinates + reference/alternative alleles.
#
# Input:
#   config/config_paths.yaml
#   config/config_parameters.yaml
#   preprocessing/finemap_snps_hg38.rds (or hg19)
#
# Output:
#   results/snp_sequences/<celltype>/ref_snp_seqs.*.fasta
#   results/snp_sequences/<celltype>/alt_snp_seqs.*.fasta
#   results/snp_sequences/<celltype>/<celltype>.rds
#
# Dependencies:
#   BSgenome.Hsapiens.*
#   Biostrings
#   GenomicRanges
# ================================================================

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(Biostrings)
  library(dplyr)
  library(yaml)
})

# ---------------------------------------------------------------
# 1. Load configuration
# ---------------------------------------------------------------
cfg_paths <- yaml::read_yaml("config/config_paths.yaml")
cfg_param <- yaml::read_yaml("config/config_parameters.yaml")

snps_rds  <- cfg_paths$input$snps_granges        # GRanges of SNPs
genome_id <- cfg_param$genome$build              # "hg19" or "hg38"

outdir    <- file.path(cfg_paths$project_root, "results/snp_sequences")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Load genome
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
# 2. Load SNP table (already standardized)
# ---------------------------------------------------------------
gr <- readRDS(snps_rds)

if (!all(c("REF","ALT") %in% colnames(mcols(gr)))) {
  stop("GRanges must contain REF and ALT columns.")
}

# Standardize chromosome format
seqlevelsStyle(gr) <- "UCSC"

# Add ±25 bp windows
window <- cfg_param$snp$window_bp   # e.g., 25
gr <- resize(gr, width = (window*2 + 1), fix="center")

# ---------------------------------------------------------------
# 3. Extract REF and ALT sequences
# ---------------------------------------------------------------
message("Extracting REF sequences...")
ref_seqs <- getSeq(genome, gr)

message("Extracting ALT sequences...")
alt_seqs <- ref_seqs

# Replace central base
mid <- ceiling(width(ref_seqs) / 2)
ref_nt <- as.character(gr$REF)
alt_nt <- as.character(gr$ALT)

ref_seqs <- replaceLetterAt(ref_seqs, at=mid, letter=ref_nt)
alt_seqs <- replaceLetterAt(alt_seqs, at=mid, letter=alt_nt)

names(ref_seqs) <- paste0(gr$rsid, "_REF_", as.character(seqnames(gr)))
names(alt_seqs) <- paste0(gr$rsid, "_ALT_", as.character(seqnames(gr)))

# ---------------------------------------------------------------
# 4. Write FASTA outputs
# ---------------------------------------------------------------
message("Saving output FASTA files...")

ref_file <- file.path(outdir, "ref_snp_seqs.fasta")
alt_file <- file.path(outdir, "alt_snp_seqs.fasta")

writeXStringSet(ref_seqs, filepath = ref_file)
writeXStringSet(alt_seqs, filepath = alt_file)

# Save processed GRanges
saveRDS(gr, file.path(outdir, "snp_windows.rds"))

message("✓ Finished generating SNP-centered sequences.")
message("REF FASTA: ", ref_file)
message("ALT FASTA: ", alt_file)
