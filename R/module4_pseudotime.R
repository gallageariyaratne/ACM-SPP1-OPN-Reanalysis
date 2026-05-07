
library(Seurat)
library(monocle3)
library(ggplot2)

myeloid <- readRDS("results/myeloid_obj.rds")


cds <- as.cell_data_set(myeloid)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)


wt_cells     <- colnames(myeloid)[myeloid$genotype == "WT"]
timd4_expr   <- GetAssayData(myeloid, 
                              assay = "SCT",
                              slot  = "data")["Timd4", wt_cells]
lyve1_expr   <- GetAssayData(myeloid,
                              assay = "SCT", 
                              slot  = "data")["Lyve1", wt_cells]

root_cells <- wt_cells[timd4_expr > 0 & lyve1_expr > 0]
message("Root cells identified: ", length(root_cells))

cds <- order_cells(cds, root_cells = root_cells)


plot_cells(cds,
           color_cells_by = "pseudotime",
           label_roots    = TRUE,
           label_leaves   = FALSE) +
  labs(title = "Myeloid Pseudotime (distance from TR-Mac root)")
ggsave("figures/Fig11_pseudotime_UMAP.pdf", width = 7, height = 6)

pseudotime_df <- data.frame(
  cell        = colnames(cds),
  pseudotime  = pseudotime(cds),
  spp1        = as.numeric(
    GetAssayData(myeloid, assay = "SCT", slot = "data")["Spp1", ]),
  genotype    = myeloid$genotype_timepoint
)


kw_test <- kruskal.test(spp1 ~ genotype, data = pseudotime_df)
message("Kruskal-Wallis p = ", kw_test$p.value)

write.csv(pseudotime_df,
          "results/tables/pseudotime_spp1_by_genotype.csv",
          row.names = FALSE)
message("Module 4 complete")
