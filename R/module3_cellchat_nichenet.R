
library(Seurat)
library(CellChat)
library(ggplot2)

sn <- readRDS("results/sn_filtered.rds")


conditions <- unique(sn$genotype_timepoint)
cellchat_list <- list()

for (cond in conditions) {
  sub <- subset(sn, genotype_timepoint == cond)
  
  cc <- createCellChat(object   = sub,
                        group.by = "cell_type")
  
  CellChatDB         <- CellChatDB.mouse
  cc@DB              <- CellChatDB
  cc                 <- subsetData(cc)
  cc                 <- identifyOverExpressedGenes(cc)
  cc                 <- identifyOverExpressedInteractions(cc)
  cc                 <- computeCommunProb(cc)
  cc                 <- computeCommunProbPathway(cc)
  cc                 <- aggregateNet(cc)
  
  cellchat_list[[cond]] <- cc
}


for (cond in conditions) {
  cc <- cellchat_list[[cond]]
  if ("SPP1" %in% cc@netP$pathways) {
    spp1_scores <- cc@netP$prob[,, "SPP1"]
    message(cond, " SPP1 myeloid->fibroblast: ",
            round(spp1_scores["Myeloid", "Fibroblast"], 4))
  } else {
    message(cond, ": SPP1 pathway not detected")
  }
}


myeloid <- readRDS("results/myeloid_obj.rds")
dsg2_mac <- subset(myeloid, genotype == "Dsg2mut")

macrophage_ligands <- c("Apoe", "C1qb", "C1qa", "Lgals3",
                         "Fn1", "Gpnmb", "Trem2", "Mmp12",
                         "Spp1", "Mmp9")

avg_expr <- rowMeans(GetAssayData(dsg2_mac, 
                                   assay = "SCT",
                                   slot  = "data")[macrophage_ligands, ])

ligand_ranking <- sort(avg_expr, decreasing = TRUE)
print(ligand_ranking)


fibroblasts   <- readRDS("results/fibroblast_obj.rds")
target_genes  <- c("Postn", "Acta2", "Col1a1", "Col3a1", "Fn1", "Cd44")

target_expr <- AverageExpression(fibroblasts,
                                  features   = target_genes,
                                  group.by   = "genotype_timepoint")
print(target_expr$SCT)
write.csv(target_expr$SCT,
          "results/tables/nichenet_target_genes.csv")

message("Module 3 complete")
