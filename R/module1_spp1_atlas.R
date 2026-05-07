
library(Seurat)
library(ggplot2)

sn <- readRDS("results/sn_filtered.rds")

DotPlot(sn, 
        features = "Spp1",
        group.by = "cell_type",
        split.by = "genotype_timepoint") +
  coord_flip() +
  labs(title = "Spp1 Expression Across Cell Types and Genotypes",
       subtitle = "dot size = % expressing, color = mean expression")
ggsave("figures/Fig2_Spp1_dotplot.pdf", width = 12, height = 6)

FeaturePlot(sn, 
            features = "Spp1",
            split.by = "genotype_timepoint",
            cols = c("lightgrey", "darkred"),
            order = TRUE)
ggsave("figures/Fig3_Spp1_UMAP.pdf", width = 18, height = 6)


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


myeloid <- subset(sn, cell_type == "Myeloid")


avg_spp1 <- AverageExpression(myeloid, 
                               features = "Spp1",
                               group.by = "genotype_timepoint")
print(avg_spp1)

table(myeloid$genotype_timepoint)

saveRDS(myeloid, "results/myeloid_obj.rds")
message("Module 1 complete")
