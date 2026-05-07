# Computational Reanalysis of GSE228048
## SPP1/OPN-Associated Cell States in Arrhythmogenic Cardiomyopathy

**Author:** Gallage H.D.N. Ariyaratne 
**Institution:** Florida State University College of Medicine  
**Lab:** Stephen P. Chelko, Ph.D., FHRS, FACC - Laboratory  
**Committee Member:** Dr. Julia Wang  

## Overview
This repository contains the full computational analysis code supporting 
the qualifying exam essay. The analysis interrogates the role of 
SPP1/OPN-associated cell states in ACM using the published GSE228048 
snRNA-seq and CITE-seq dataset (PMID: 38949031), without a direct Spp1 
perturbation condition.

## Dataset
GSE228048 — available at NCBI GEO:  
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE228048

## Analysis Modules
| Script | Module | Description |
|--------|--------|-------------|
| 01_QC_preprocessing.R | II | Quality control and cell type validation |
| 02_module1_spp1_atlas.R | III | Global Spp1 expression and OPN receptor landscape |
| 03_module2_DAM_fibroblast.R | IV | DAM identity and fibroblast OPN responsiveness |
| 04_module3_cellchat_nichenet.R | V | Intercellular communication and ligand activity |
| 05_module4_pseudotime.R | VI | Monocle3 pseudotime trajectory analysis |
| 06_module5_regulon.R | VII | Transcription factor regulon activity scoring |
| 07_module6_CITE_seq_WNN.R | VIII | CITE-seq WNN multimodal integration |

## Software Requirements
- R 4.3.1
- Seurat 4.x
- Monocle3
- CellChat
- DESeq2
- ggplot2
- miloR

## HPC Environment
All analyses were performed on the FSU HPC cluster.  
Conda environment: acm_analysis  
Module: R/4.3.1  

## Statistical Methods Note
See `supplementary/statistical_methods_note.md` for full declaration 
of all statistical tests, correction methods, and fold change 
calculation procedures referenced in the essay.
