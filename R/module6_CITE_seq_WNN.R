

library(Seurat)
library(ggplot2)

cite <- readRDS("data/GSE228048_CITE_global.rds.gz")



cite <- NormalizeData(cite,
                       normalization.method = "CLR",
                       margin               = 2,
                       assay                = "ADT")


cite <- FindMultiModalNeighbors(
  cite,
  reduction.list = list("pca", "apca"),
  dims.list      = list(1:30, 1:18),
  modality.weight.name = "RNA.weight"
)

cite <- RunUMAP(cite,
                nn.name = "weighted.nn",
                reduction.name = "wnn.umap",
                reduction.key  = "wnnUMAP_")

DimPlot(cite,
        reduction  = "wnn.umap",
        group.by   = "cell_type",
        label      = TRUE) +
  labs(title = "WNN UMAP: Cell Types")
ggsave("figures/Fig14a_WNN_celltypes.pdf", width = 7, height = 6)

FeaturePlot(cite,
            features   = "Spp1",
            reduction  = "wnn.umap") +
  labs(title = "WNN UMAP: Spp1 RNA")
ggsave("figures/Fig14b_WNN_Spp1RNA.pdf", width = 6, height = 6)

FeaturePlot(cite,
            features  = "CD44",
            reduction = "wnn.umap") +
  labs(title = "WNN UMAP: CD44 Protein (ADT)")
ggsave("figures/Fig14c_WNN_CD44protein.pdf", width = 6, height = 6)

cd44_by_celltype <- AverageExpression(
  cite,
  features = "CD44",
  group.by = c("cell_type", "genotype"),
  assay    = "ADT"
)
print(cd44_by_celltype$ADT)
write.csv(cd44_by_celltype$ADT,
          "results/tables/CD44_ADT_by_celltype_genotype.csv")



fibroblast_cd44_rna     <- mean(
  GetAssayData(subset(cite, cell_type == "Fibroblast"),
               assay = "SCT", slot = "data")["Cd44", ])

fibroblast_cd44_protein <- mean(
  GetAssayData(subset(cite, cell_type == "Fibroblast"),
               assay = "ADT", slot = "data")["CD44", ])

message("Fibroblast Cd44 RNA mean: ",    round(fibroblast_cd44_rna, 3))
message("Fibroblast CD44 protein mean: ", round(fibroblast_cd44_protein, 3))

message("Module 6 complete")
