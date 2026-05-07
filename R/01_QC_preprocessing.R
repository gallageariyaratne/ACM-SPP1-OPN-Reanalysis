# ============================================================
# Module II: Quality Control and Preprocessing
# GSE228048 SPP1/OPN Reanalysis
# Author: Gallage H.D.N. Ariyaratne (Nipun)
# Chelko Lab, FSU College of Medicine
# May 2026
# ============================================================

library(Seurat)
library(ggplot2)
library(patchwork)

# ── Load data ─────────────────────────────────────────────
sn <- readRDS("data/GSE228048_sn_global.rds.gz")

# ── QC metrics ────────────────────────────────────────────
# nUMI: total counts per nucleus
# nFeature: number of genes detected per nucleus  
# percent.mt: mitochondrial transcript percentage

sn[["percent.mt"]] <- PercentageFeatureSet(sn, pattern = "^mt-")

# ── QC violin plots (Supplementary Figure 1) ──────────────
p1 <- VlnPlot(sn, features = "nCount_RNA", 
              group.by = "genotype_timepoint",
              pt.size = 0) + 
      labs(title = "nUMI (Counts per Nucleus)") +
      theme(legend.position = "none")

p2 <- VlnPlot(sn, features = "nFeature_RNA",
              group.by = "genotype_timepoint", 
              pt.size = 0) +
      labs(title = "nFeature (Genes per Nucleus)") +
      theme(legend.position = "none")

p3 <- VlnPlot(sn, features = "percent.mt",
              group.by = "genotype_timepoint",
              pt.size = 0) +
      labs(title = "Mitochondrial % per Nucleus") +
      theme(legend.position = "none")

p1 | p2 | p3
ggsave("figures/FigS1_QC_metrics.pdf", width = 14, height = 4)

# ── Filter nuclei ─────────────────────────────────────────
# Exclude nuclei with percent.mt > 2.5% (damaged/necrotic)
sn <- subset(sn, subset = percent.mt < 2.5)

# ── Validate published cell type annotations ───────────────
# Use published annotations directly — validated by UMAP separation
DimPlot(sn, group.by = "cell_type", label = TRUE) +
  labs(title = "Global Cell Type Annotations")
ggsave("figures/Fig4_CellType_UMAP.pdf", width = 8, height = 7)

# ── Save filtered object ──────────────────────────────────
saveRDS(sn, "results/sn_filtered.rds")
message("QC complete. Cells retained: ", ncol(sn))
