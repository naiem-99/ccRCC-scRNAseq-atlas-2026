############################################################
# 3. COPY NUMBER INFERENCE WITH inferCNV
#   - Takes tumor cells and a set of immune reference cells
#   - Infers large-scale copy number changes in tumor cells
#     by comparing their expression along each chromosome
#     to the immune reference (assumed to be diploid)
#   - Adds a per-cell CNV score back to the Seurat object
#
# Reference cells: monocytes, dendritic cells, B cells
#
# Input:
#   output/ccRCC_integrated_annotated.rds   (from 2.Annotation.R)
#   data/hg38_gencode_v27.txt                (gene position file)
#
# Output:
#   output/infercnv/                          (heatmaps and HMM calls)
#   output/ccRCC_tumor_CNV.rds
############################################################

rm(list = ls()); gc()

library(Seurat)
library(dplyr)
library(infercnv)

options(scipen = 100)

data_dir     <- "data"
output_dir   <- "output"
infercnv_dir <- file.path(output_dir, "infercnv")

dir.create(infercnv_dir, showWarnings = FALSE, recursive = TRUE)


## ---------------------------------------------------------
## Load the annotated object
## ---------------------------------------------------------
integrated <- readRDS(file.path(output_dir, "ccRCC_integrated_annotated.rds"))


## ---------------------------------------------------------
## Keep tumor cells and the immune reference cells
## ---------------------------------------------------------
reference_types <- c("Mono", "DC", "B Cells")
keep_types      <- c(reference_types, "Tumor")

cnv_obj <- subset(integrated, subset = CellType %in% keep_types)
rm(integrated); gc()


## ---------------------------------------------------------
## Prepare inferCNV inputs
##   1. raw count matrix
##   2. cell annotation table (cell barcode -> cell type)
##   3. gene order file (gene -> chromosome position)
## ---------------------------------------------------------
counts_matrix <- GetAssayData(cnv_obj, assay = "RNA", slot = "counts")

cell_annotation <- data.frame(
  cell_type = as.character(cnv_obj$CellType),
  row.names = colnames(cnv_obj)
)

## Make sure every reference group is present before running
missing_refs <- setdiff(reference_types, unique(cell_annotation$cell_type))
if (length(missing_refs) > 0) {
  stop("Missing reference groups: ", paste(missing_refs, collapse = ", "))
}


## ---------------------------------------------------------
## Create and run inferCNV
## Subcluster mode lets inferCNV find subclones within tumor cells.
## The i6 HMM calls six copy number states per region.
## ---------------------------------------------------------
infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = counts_matrix,
  annotations_file  = cell_annotation,
  gene_order_file   = file.path(data_dir, "hg38_gencode_v27.txt"),
  ref_group_names   = reference_types
)
rm(counts_matrix); gc()

infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff            = 0.1,     # recommended value for 10x data
  out_dir           = infercnv_dir,
  cluster_by_groups = TRUE,
  analysis_mode     = "subclusters",
  denoise           = TRUE,
  HMM               = TRUE,
  HMM_type          = "i6",
  BayesMaxPNormal   = 0,
  output_format     = "pdf",
  num_threads       = 3
)

saveRDS(infercnv_obj, file.path(infercnv_dir, "infercnv_obj.rds"))


## ---------------------------------------------------------
## Add inferCNV results back to the Seurat object
## ---------------------------------------------------------
cnv_obj <- infercnv::add_to_seurat(
  infercnv_output_path = infercnv_dir,
  seurat_obj           = cnv_obj,
  top_n                = 10
)
rm(infercnv_obj); gc()
## ---------------------------------------------------------
## Keep tumor cells and compute a simple CNV burden score
## The score adds up the scaled proportion of each chromosome
## affected by a copy number change in that cell.
## ---------------------------------------------------------
tumor_cnv <- subset(cnv_obj, subset = CellType == "Tumor")
rm(cnv_obj); gc()

cnv_columns <- grep("^proportion_scaled_cnv_chr",
                    colnames(tumor_cnv@meta.data), value = TRUE)

tumor_cnv$CNV_score <- rowSums(tumor_cnv@meta.data[, cnv_columns], na.rm = TRUE)

## ---------------------------------------------------------
## CNV burden across patients
## ---------------------------------------------------------
VlnPlot(tumor_cnv, features = "CNV_score", group.by = "Patient", pt.size = 0)

## ---------------------------------------------------------
## Save
## ---------------------------------------------------------
saveRDS(tumor_cnv, file.path(output_dir, "ccRCC_tumor_CNV.rds"))
