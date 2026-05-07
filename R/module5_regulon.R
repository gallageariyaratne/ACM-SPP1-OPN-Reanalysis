

library(Seurat)
library(ggplot2)

myeloid <- readRDS("results/myeloid_obj.rds")


gene_sets <- list(
  NFkB_Rela   = c("Nfkb1", "Rela", "Ikbkb", "Tnf", "Il6", "Cxcl10"),
  AP1_FosJun  = c("Fos", "Jun", "Junb", "Fosl1", "Fosl2", "Atf3"),
  DAM_Program = c("Spp1", "Trem2", "Fn1", "Lgals3", "Gpnmb", "Cd44"),
  TRM_Program = c("Timd4", "Lyve1", "Folr2", "Cx3cr1"),
  MMP_Remodel = c("Mmp12", "Mmp14", "Mmp9", "Mmp2"),
  Complement  = c("C1qa", "C1qb", "C1qc", "C3", "Cfb")
)

myeloid <- AddModuleScore(myeloid,
                           features = gene_sets,
                           name     = "Regulon")



myeloid$spp1_pos <- GetAssayData(myeloid,
                                  assay = "SCT",
                                  slot  = "data")["Spp1", ] > 0

spp1_pos_meta <- myeloid@meta.data[myeloid$spp1_pos == TRUE, ]
spp1_neg_meta <- myeloid@meta.data[myeloid$spp1_pos == FALSE, ]

message("Spp1+ n = ", nrow(spp1_pos_meta))
message("Spp1- n = ", nrow(spp1_neg_meta))


regulon_cols <- paste0("Regulon", 1:6)
program_names <- names(gene_sets)

results <- data.frame()
for (i in seq_along(regulon_cols)) {
  col      <- regulon_cols[i]
  pos_mean <- mean(spp1_pos_meta[[col]])
  neg_mean <- mean(spp1_neg_meta[[col]])
  fc       <- pos_mean / (neg_mean + 1e-6)
  
  results <- rbind(results, data.frame(
    program     = program_names[i],
    spp1_pos    = round(pos_mean, 3),
    spp1_neg    = round(neg_mean, 3),
    fold_change  = round(fc, 1)
  ))
}



p_values <- c()
for (col in regulon_cols) {
  test     <- wilcox.test(spp1_pos_meta[[col]],
                           spp1_neg_meta[[col]],
                           alternative = "two.sided")
  p_values <- c(p_values, test$p.value)
}

p_bonf              <- p.adjust(p_values, method = "bonferroni")
results$p_raw       <- round(p_values, 4)
results$p_bonf      <- round(p_bonf, 4)
results$significant <- ifelse(p_bonf < 0.001, "***",
                        ifelse(p_bonf < 0.05, "*", "ns"))

print(results)
write.csv(results,
          "results/tables/regulon_fold_changes.csv",
          row.names = FALSE)

avg_regulon <- AverageExpression(myeloid,
                                  features   = regulon_cols,
                                  group.by   = "genotype_timepoint")
write.csv(avg_regulon$SCT,
          "results/tables/regulon_by_genotype.csv")

message("Module 5 complete")
message("Results saved to results/tables/regulon_fold_changes.csv")
