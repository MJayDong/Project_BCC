#!/usr/bin/env bash
#
# 4.Running_gkmsvm_explain_snps.sh
# Score SNP sequences (REF/ALT) using trained gkm-SVM models
# and compute deltaSVM scores.
#
# Part of the RegulSNP-SVM pipeline.

set -euo pipefail

############################################
# User-configurable paths
############################################

# gkmpredict binary (from lsgkm)
GKMPREDICT_BIN="${GKMPREDICT_BIN:-gkmpredict}"

# Directory with FASTA files:
#   <model>.SNP.ref.fa
#   <model>.SNP.alt.fa
SNP_FASTA_DIR="${SNP_FASTA_DIR:-./snp_fastas}"

# Directory containing trained full models (*.model.txt)
MODEL_DIR="${MODEL_DIR:-./models/full_models_1000bp}"

# Output directory
OUTDIR="${OUTDIR:-./snp_scores}"
mkdir -p "${OUTDIR}"

############################################
# Models (same list as training stage)
############################################

models=(
  "CD8_eff_Post"
  "CD8_ex_Post"
  "CD8_mem_Pre"
  "CD8_mem_Post"
  "Naive_Pre"
  "Naive_Post"
  "Tfh_Pre"
  "Tfh_Post"
  "Th17_Pre"
  "Th17_Post"
  "Tregs_Pre"
  "Tregs_Post"
)

############################################
# Main loop
############################################

for model in "${models[@]}"; do
  echo "==============================="
  echo "Scoring SNPs for model: ${model}"
  echo "==============================="

  model_file=$(ls "${MODEL_DIR}/${model}"*.model.txt 2>/dev/null || true)
  if [[ -z "${model_file}" ]]; then
    echo "[WARN] No model found for ${model}, skipping"
    continue
  fi

  ref_fa="${SNP_FASTA_DIR}/${model}.SNP.ref.fa"
  alt_fa="${SNP_FASTA_DIR}/${model}.SNP.alt.fa"

  if [[ ! -f "${ref_fa}" ]]; then
    echo "[WARN] Missing REF fasta: ${ref_fa}, skipping"
    continue
  fi
  if [[ ! -f "${alt_fa}" ]]; then
    echo "[WARN] Missing ALT fasta: ${alt_fa}, skipping"
    continue
  fi

  out_ref="${OUTDIR}/${model}.SNP.ref.svm.txt"
  out_alt="${OUTDIR}/${model}.SNP.alt.svm.txt"
  out_delta="${OUTDIR}/${model}.SNP.deltaSVM.txt"

  ############################################
  # Score REF alleles
  ############################################
  if [[ -f "${out_ref}" ]]; then
    echo "  REF scores exist, skip"
  else
    echo "  Predicting REF alleles..."
    "${GKMPREDICT_BIN}" \
      "${ref_fa}" \
      "${model_file}" \
      "${out_ref}"
  fi

  ############################################
  # Score ALT alleles
  ############################################
  if [[ -f "${out_alt}" ]]; then
    echo "  ALT scores exist, skip"
  else
    echo "  Predicting ALT alleles..."
    "${GKMPREDICT_BIN}" \
      "${alt_fa}" \
      "${model_file}" \
      "${out_alt}"
  fi

  ############################################
  # Compute deltaSVM = ALT - REF
  ############################################
  if [[ -f "${out_delta}" ]]; then
    echo "  deltaSVM exists, skip"
  else
    echo "  Computing deltaSVM..."

    paste "${out_ref}" "${out_alt}" | \
      awk 'BEGIN{OFS="\t"} {print $1, $2, $4, $4-$2}' \
      > "${out_delta}"

    # Output columns:
    #   ID | REF_score | ALT_score | deltaSVM
  fi

  echo "  Done."
done

echo "All SNP scoring completed."
