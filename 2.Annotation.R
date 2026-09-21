############################################################
# 2. INTEGRATION AND CELL TYPE ANNOTATION
#   - Integrates the 10 patients with Seurat CCA integration
#     to remove patient-level batch effects
#   - Clusters the cells and annotates each cluster
#     using canonical marker genes
#
# Input:
#   output/ccRCC_merged_10_patients.rds   (from 1.Data_Preprocessing.R)
#
# Output:
#   output/ccRCC_integrated_annotated.rds
#
# Next step: 3.InferCNV.R
############################################################

rm(list = ls()); gc()

library(Seurat)
library(dplyr)

options(future.globals.maxSize = 20 * 1024^3)

output_dir <- "output"


## ---------------------------------------------------------
## Load the merged object
## ---------------------------------------------------------
ccRCC_merged <- readRDS(file.path(output_dir, "ccRCC_merged_10_patients.rds"))


## ---------------------------------------------------------
## Integration across patients
## Each patient is normalised on its own, then Seurat finds
## shared "anchor" cells to align them into one space.
## ---------------------------------------------------------
patient_list <- SplitObject(ccRCC_merged, split.by = "Patient")
rm(ccRCC_merged); gc()

for (i in seq_along(patient_list)) {
  DefaultAssay(patient_list[[i]]) <- "RNA"
  patient_list[[i]] <- NormalizeData(patient_list[[i]])
  patient_list[[i]] <- FindVariableFeatures(patient_list[[i]], nfeatures = 2000)
}

anchors    <- FindIntegrationAnchors(object.list = patient_list, dims = 1:30)
integrated <- IntegrateData(anchorset = anchors)
rm(patient_list, anchors); gc()


## ---------------------------------------------------------
## Dimensional reduction and clustering
## ---------------------------------------------------------
DefaultAssay(integrated) <- "integrated"

integrated <- ScaleData(integrated)
integrated <- RunPCA(integrated)
integrated <- RunUMAP(integrated, dims = 1:30)
integrated <- FindNeighbors(integrated, dims = 1:30)
integrated <- FindClusters(integrated, resolution = 0.4)

## Quick look at the clusters
DimPlot(integrated, label = TRUE) + NoLegend()
DimPlot(integrated, group.by = "Patient")


## ---------------------------------------------------------
## Cell type annotation
## Labels were assigned by checking canonical markers per
## cluster at resolution 0.4. For example:
##   Tumor       CA9, NDUFA4L2
##   Macrophage  CD68, C1QA
##   Endothelial PECAM1, VWF
##   T cells     CD3D, CD4, CD8A
## Uncomment the lines below to re-check the markers yourself.
## ---------------------------------------------------------
# DefaultAssay(integrated) <- "RNA"
# markers <- FindAllMarkers(integrated, only.pos = TRUE)
# DotPlot(integrated, features = c("CA9", "NDUFA4L2", "CD68", "C1QA",
#                                  "PECAM1", "VWF", "CD3D", "CD8A")) + RotatedAxis()

cluster_labels <- c(
  "0"  = "Tumor",       "1"  = "Macro",      "2"  = "Tumor",
  "3"  = "vSMC",        "4"  = "Endo",       "5"  = "Endo",
  "6"  = "CD4+ T",      "7"  = "CD8+ T",     "8"  = "Mono",
  "9"  = "Endo",        "10" = "DC",         "11" = "Tumor",
  "12" = "Cycling",     "13" = "Tumor",      "14" = "B Cells",
  "15" = "Pericytes",   "16" = "CD4+ Treg",  "17" = "Fibroblast",
  "18" = "Mast",        "19" = "Podocytes",  "20" = "Endo"
)

integrated$CellType <- cluster_labels[as.character(integrated$integrated_snn_res.0.4)]

## Annotated UMAP
DimPlot(integrated, group.by = "CellType", label = TRUE, repel = TRUE)


## ---------------------------------------------------------
## Save
## ---------------------------------------------------------
saveRDS(integrated, file.path(output_dir, "ccRCC_integrated_annotated.rds"))
