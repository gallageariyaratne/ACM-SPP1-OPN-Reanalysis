# ============================================================
# Module IV: DAM Identity + Fibroblast OPN Responsiveness
# ============================================================

library(Seurat)
library(ggplot2)

myeloid <- readRDS("results/myeloid_obj.rds")
sn      <- readRDS("results/sn_filtered.rds")

# ── Figure 7: Spp1 co-expression in Dsg2mut myeloid ───────
dsg2_myeloid <- subset(myeloid, genotype == "Dsg2mut")

# Pearson correlation of all genes with Spp1
spp1_expr <- GetAssayData(dsg2_myeloid, 
                           assay = "SCT", 
                           slot  = "data")["Spp1", ]

gene_matrix <- GetAssayData(dsg2_myeloid, 
                             assay = "SCT", 
                             slot  = "data")

correlations <- apply(gene_matrix, 1, function(x) {
  cor(x, spp1_expr, method = "pearson")
})

top_corr <- sort(correlations, decreasing = TRUE)[1:20]
write.csv(data.frame(gene = names(top_corr), 
                     pearson_r = top_corr),
          "results/tables/spp1_pearson_correlations.csv")

# ── Supplementary Figure S3: Spp1+ vs Spp1- comparison ───
myeloid$spp1_pos <- GetAssayData(myeloid, 
                                  assay = "SCT",
                                  slot  = "data")["Spp1", ] > 0

spp1_pos <- myeloid@meta.data[myeloid$spp1_pos == TRUE, ]
spp1_neg <- myeloid@meta.data[myeloid$spp1_pos == FALSE, ]

message("Spp1+ nuclei: ", sum(myeloid$spp1_pos))
message("Spp1- nuclei: ", sum(!myeloid$spp1_pos))

# Mean expression of DAM markers in each group
dam_genes <- c("Spp1", "Trem2", "Fn1", "Lgals3", "Gpnmb",
               "Timd4", "Lyve1", "Ccr2", "Adgre1", "Cd44")

avg_pos <- rowMeans(GetAssayData(
  subset(myeloid, spp1_pos == TRUE), 
  assay = "SCT", slot = "data")[dam_genes, ])

avg_neg <- rowMeans(GetAssayData(
  subset(myeloid, spp1_pos == FALSE),
  assay = "SCT", slot = "data")[dam_genes, ])

fold_changes <- avg_pos / (avg_neg + 1e-6)  # pseudocount for zero denom

dam_results <- data.frame(
  gene       = dam_genes,
  mean_pos   = round(avg_pos, 3),
  mean_neg   = round(avg_neg, 3),
  fold_change = round(fold_changes, 1)
)
write.csv(dam_results, 
          "results/tables/DAM_marker_fold_changes.csv",
          row.names = FALSE)
print(dam_results)

# ── Figure 8: Fibroblast compartment analysis ─────────────
fibroblasts <- subset(sn, cell_type == "Fibroblast")

# CD44-high vs CD44-low stratification
cd44_expr <- GetAssayData(fibroblasts, 
                           assay = "SCT", 
                           slot  = "data")["Cd44", ]

fibroblasts$cd44_group <- ifelse(
  cd44_expr > quantile(cd44_expr, 0.75),
  "CD44-high", "CD44-low"
)

# Myofibroblast score
activation_genes <- c("Postn", "Acta2", "Col1a1", "Col3a1")
fibroblasts <- AddModuleScore(fibroblasts,
                               features = list(activation_genes),
                               name     = "Myofibroblast")

# Compare scores between CD44-high and CD44-low
high_scores <- fibroblasts$Myofibroblast1[
  fibroblasts$cd44_group == "CD44-high"]
low_scores  <- fibroblasts$Myofibroblast1[
  fibroblasts$cd44_group == "CD44-low"]

mw_test <- wilcox.test(high_scores, low_scores)
message("Mann-Whitney p = ", round(mw_test$p.value, 4))

# Sign test: directional consistency across 6 conditions
conditions <- unique(fibroblasts$genotype_timepoint)
direction  <- sapply(conditions, function(cond) {
  sub <- subset(fibroblasts, genotype_timepoint == cond)
  mean(sub$Myofibroblast1[sub$cd44_group == "CD44-high"]) >
  mean(sub$Myofibroblast1[sub$cd44_group == "CD44-low"])
})
message("Consistent direction across conditions: ", sum(direction), "/6")
message("Sign test probability: ", 0.5^6)

saveRDS(fibroblasts, "results/fibroblast_obj.rds")
message("Module 2 complete")
