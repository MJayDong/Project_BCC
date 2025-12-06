# RegulSNP-SVM
A scalable and modular framework for predicting the regulatory impact of genetic variants using scATAC-seq–derived chromatin accessibility profiles and gkm-SVM models.

---

## Overview
RegulSNP-SVM provides an end-to-end workflow for:

- Extracting SNP-centered reference and alternative sequences
- Building GC-matched null sequences for negative training
- Training cell-type–specific gkm-SVM models using scATAC-seq peak sets
- Performing SNP effect prediction (ΔSVM)
- Visualizing peak models, SVM scores, regulatory maps, and peak–gene links

This framework is designed for single-cell chromatin datasets (e.g., ArchR projects) and supports millions of peaks across diverse immune and non-immune cell types. All components are fully modular and can be reused independently for other regulatory genomics projects.

---

## Key Features

### **1. SNP-centered sequence extraction**
- Unified processing for fine-mapping tables (Excel/CSV)
- Automatic coordinate standardization, duplicate removal, QC
- Supports liftover between hg19 and hg38
- Generates ±25 bp sequence windows for both REF/ALT alleles

### **2. GC-matched null sequence generation**
- Genome-wide random region sampling  
- Removal of blacklist and scATAC peaks  
- GC-distribution matching using adaptive binning  

Ensures unbiased negative sets for gkm-SVM training.

### **3. Cell-type–specific gkm-SVM training**
- scATAC peak selection based on marker peaks + top-scoring peaks  
- Optional Jaccard similarity analysis across cell types  
- Automatic FASTA writing and split into train/test  

### **4. SNP impact scoring (ΔSVM)**
- Predicts allele-specific accessibility differences  
- Supports millions of SNPs  
- Includes motif-centric interpretation tools  

### **5. Visualization modules**
- Peak model heatmaps  
- ΔSVM score distribution and annotation  
- Peak-to-gene regulatory maps  
- Base-resolution ATCG contribution plots  

---

## Requirements

### **R packages**
