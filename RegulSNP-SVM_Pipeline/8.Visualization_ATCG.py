#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
8.Visualization_ATCG.py

Purpose:
    Visualize nucleotide composition (ATCG) for training or SNP sequences.
    Useful for QC and comparing REF/ALT flanking regions.

Inputs:
    - FASTA file (training peaks, null sequences, or SNP FASTA)
    - config/config_paths.yaml

Outputs:
    - results/seq_visualization/<name>_nt_composition.pdf
    - results/seq_visualization/<name>_logo.pdf (if logomaker available)

Author:
    RegulSNP-SVM pipeline
"""

import os
import sys
import yaml
import argparse
import pandas as pd
from collections import Counter
from Bio import SeqIO
import matplotlib.pyplot as plt

# Optional: sequence logo
LOGOMAKER_AVAILABLE = False
try:
    import logomaker
    LOGOMAKER_AVAILABLE = True
except ImportError:
    pass


def load_config():
    """Load YAML config files."""
    with open("config/config_paths.yaml") as f:
        cfg_paths = yaml.safe_load(f)
    with open("config/config_parameters.yaml") as f:
        cfg_param = yaml.safe_load(f)
    return cfg_paths, cfg_param


def compute_nt_composition(fasta_file):
    """Count frequencies of A/T/C/G across all sequences."""
    counts = Counter()

    total = 0
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq).upper()
        for nt in seq:
            if nt in "ATCG":
                counts[nt] += 1
                total += 1

    df = pd.DataFrame({
        "nt": ["A", "T", "C", "G"],
        "count": [counts["A"], counts["T"], counts["C"], counts["G"]],
    })
    df["freq"] = df["count"] / df["count"].sum()

    return df


def plot_nt_composition(df, name, outdir):
    """Barplot of nucleotide composition."""
    plt.figure(figsize=(6,4))
    plt.bar(df["nt"], df["freq"], color=["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"])
    plt.title(f"Nucleotide Composition: {name}")
    plt.ylabel("Frequency")
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, f"{name}_nt_composition.pdf"))
    plt.close()


def compute_position_matrix(fasta_file):
    """Compute position frequency matrix for sequence logo."""
    seqs = [str(r.seq).upper() for r in SeqIO.parse(fasta_file, "fasta")]
    if len(seqs) == 0:
        raise ValueError("No sequences found in FASTA.")

    L = len(seqs[0])
    for s in seqs:
        if len(s) != L:
            raise ValueError("Sequences must have equal length for logos.")

    mat = pd.DataFrame(0, index=range(L), columns=list("ATCG"))
    for s in seqs:
        for i, nt in enumerate(s):
            if nt in "ATCG":
                mat.loc[i, nt] += 1

    mat = mat.div(mat.sum(axis=1), axis=0)  # convert to frequency
    return mat


def plot_logo(mat, name, outdir):
    """Generate sequence logo using logomaker."""
    if not LOGOMAKER_AVAILABLE:
        print("[WARN] logomaker not installed: skipping logo.")
        return

    plt.figure(figsize=(10,3))
    logomaker.Logo(mat, shade_below=.5, fade_below=.5)
    plt.title(f"Sequence Logo: {name}")
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, f"{name}_logo.pdf"))
    plt.close()


def main():
    parser = argparse.ArgumentParser(description="ATCG visualization tool")
    parser.add_argument("--fasta", required=True, help="Input FASTA file")
    parser.add_argument("--name", required=False, default="sequence", help="Name tag for output files")
    args = parser.parse_args()

    cfg_paths, cfg_param = load_config()
    outdir = os.path.join(cfg_paths["project_root"], "results/seq_visualization")
    os.makedirs(outdir, exist_ok=True)

    fasta_file = args.fasta
    name = args.name

    print("Computing nucleotide composition...")
    df = compute_nt_composition(fasta_file)
    plot_nt_composition(df, name, outdir)
    print("✓ Composition plot saved.")

    print("Computing sequence logo matrix...")
    try:
        mat = compute_position_matrix(fasta_file)
        plot_logo(mat, name, outdir)
        print("✓ Sequence logo saved.")
    except Exception as e:
        print("[WARN] logo skipped:", str(e))


if __name__ == "__main__":
    main()
