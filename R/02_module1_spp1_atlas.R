# ============================================================
# Module III: Global Spp1 Expression Atlas + OPN Receptor Landscape
# ============================================================

library(Seurat)
library(ggplot2)

sn <- readRDS("results/sn_filtered.rds")

# ── Figure 2: Spp1 dot plot across all cell types ─────────
DotPlot(sn, 
        features = "Spp1",
        group.by = "cell_type",
        split.by = "genotype_timepoint") +
  coord_flip() +
  labs(title = "Spp1 Expression Across Cell Types and Genotypes",
       subtitle = "dot size = % expressing, color = mean expression")
ggsave("figures/Fig2_Spp1_dotplot.pdf", width = 12, height = 6)

# ── Figure 3: Spp1 UMAP by genotype ──────────────────────
FeaturePlot(sn, 
            features = "Spp1",
            split.by = "genotype_timepoint",
            cols = c("lightgrey", "darkred"),
            order = TRUE)
ggsave("figures/Fig3_Spp1_UMAP.pdf", width = 18, height = 6)

# ── Supplementary Figure 4: OPN receptor landscape ────────
receptors <- c("Cd44", "Itgav", "Itgb3", "Itgb5")

for (gene in receptors) {
  DotPlot(sn,
          features = gene,
          group.by = "cell_type",
          split.by = "genotype_timepoint") +
    coord_flip() +
    labs(title = paste(gene, "Expression Across Cell Types"))
  ggsave(paste0("figures/FigS_receptor_", gene, ".pdf"),
         width = 12, height = 6)
}

# ── Supplementary Figure S1: Temporal dynamics ────────────
myeloid <- subset(sn, cell_type == "Myeloid")

# Mean Spp1 per condition and timepoint
avg_spp1 <- AverageExpression(myeloid, 
                               features = "Spp1",
                               group.by = "genotype_timepoint")
print(avg_spp1)

# Myeloid cell count per condition
table(myeloid$genotype_timepoint)

saveRDS(myeloid, "results/myeloid_obj.rds")
message("Module 1 complete")
