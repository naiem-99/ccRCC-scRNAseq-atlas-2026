############################################################
# 1. DATA PREPROCESSING
#
#   - Loads two public ccRCC single-cell datasets from GEO
#       GSE304466: 3 ccRCC tumor samples
#       GSE159115: 7 ccRCC tumor samples (normal and chRCC removed)
#   - Runs quality control on each dataset
#   - Merges them into one 10-patient object
#
# Output:
#   output/ccRCC_merged_10_patients.rds
#
# Next step: 2.Annotation.R
############################################################

rm(list = ls()); gc()

library(Seurat)
library(dplyr)
library(scDblFinder)
library(BiocParallel)
library(SingleCellExperiment)

# Allow Seurat to handle large objects
options(future.globals.maxSize = 20 * 1024^3)


## ---------------------------------------------------------
## Folder setup
## Change these paths to match where your data is stored.
## ---------------------------------------------------------
data_dir   <- "data"      # raw downloaded GEO files
output_dir <- "output"    # all saved objects go here

dir.create(output_dir, showWarnings = FALSE)


############################################################
# PART A: GSE304466 (3 tumor samples)
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE304466
############################################################

## Load the 10x count matrices for each sample
ccRCC001_counts <- Read10X(file.path(data_dir, "GSE304466", "51"))
ccRCC003_counts <- Read10X(file.path(data_dir, "GSE304466", "53"))
ccRCC005_counts <- Read10X(file.path(data_dir, "GSE304466", "55"))

## Create Seurat objects
## Keep genes seen in at least 10 cells, and cells with at least 600 genes
ccRCC001 <- CreateSeuratObject(ccRCC001_counts, project = "ccRCC001",
                               min.cells = 10, min.features = 600)
ccRCC003 <- CreateSeuratObject(ccRCC003_counts, project = "ccRCC003",
                               min.cells = 10, min.features = 600)
ccRCC005 <- CreateSeuratObject(ccRCC005_counts, project = "ccRCC005",
                               min.cells = 10, min.features = 600)

rm(ccRCC001_counts, ccRCC003_counts, ccRCC005_counts); gc()

## Add the sample name in front of each cell barcode
## so barcodes stay unique after merging
ccRCC001 <- RenameCells(ccRCC001, add.cell.id = "ccRCC001")
ccRCC003 <- RenameCells(ccRCC003, add.cell.id = "ccRCC003")
ccRCC005 <- RenameCells(ccRCC005, add.cell.id = "ccRCC005")

## Merge the three samples
gse304466 <- merge(ccRCC001, y = list(ccRCC003, ccRCC005), project = "GSE304466")
rm(ccRCC001, ccRCC003, ccRCC005); gc()

## Mitochondrial gene percentage (high values mean damaged cells)
gse304466[["percent.mt"]] <- PercentageFeatureSet(gse304466, pattern = "^MT-")

## Doublet detection with scDblFinder, run separately per sample
sce <- as.SingleCellExperiment(gse304466)
sce <- scDblFinder(sce, samples = "orig.ident",
                   BPPARAM = SnowParam(workers = 3))

gse304466$doublet_class <- sce$scDblFinder.class
rm(sce); gc()

## Keep only singlets with less than 10% mitochondrial reads
gse304466 <- subset(gse304466,
                    subset = doublet_class == "singlet" & percent.mt < 10)

## Patient label used throughout the pipeline
gse304466$Patient <- gse304466$orig.ident


############################################################
# PART B: GSE159115 (7 ccRCC tumor samples)
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159115
############################################################

gse159115_dir <- file.path(data_dir, "GSE159115")

## Find all H5 count files in the folder
h5_files <- list.files(gse159115_dir,
                       pattern = "filtered_gene_bc_matrices_h5.h5$",
                       full.names = TRUE)

## Read each file into a Seurat object and record its sample name
sample_list <- list()

for (f in h5_files) {
  sample_id <- gsub("_filtered_gene_bc_matrices_h5.h5", "", basename(f))
  counts    <- Read10X_h5(f)
  obj       <- CreateSeuratObject(counts = counts)
  obj$sample <- sample_id
  sample_list[[sample_id]] <- obj
}

## Merge all samples, adding the sample name to each barcode
gse159115 <- merge(sample_list[[1]],
                   y = sample_list[-1],
                   add.cell.ids = names(sample_list),
                   project = "GSE159115")
rm(sample_list, counts, obj); gc()

## Trim barcodes so they match the cell names in the author annotation files
new_names <- sub("^.*?(SI_)", "\\1", colnames(gse159115))
gse159115 <- RenameCells(gse159115, new.names = new_names)

## Load the authors' cell annotations
anno_ccRCC  <- read.csv(file.path(gse159115_dir, "GSE159115_ccRCC_anno.csv"))
anno_chRCC  <- read.csv(file.path(gse159115_dir, "GSE159115_chRCC_anno.csv"))
anno_normal <- read.csv(file.path(gse159115_dir, "GSE159115_normal_anno.csv"))

anno_all <- rbind(anno_ccRCC, anno_chRCC, anno_normal)
rownames(anno_all) <- anno_all$cell

gse159115 <- AddMetaData(gse159115, metadata = anno_all)

## Drop cells the authors did not annotate
gse159115 <- gse159115[, !is.na(gse159115$anno)]

## Label each sample as normal kidney, chRCC, or ccRCC
normal_samples <- c("SI_18856", "SI_19704", "SI_21255",
                    "SI_21256", "SI_22369", "SI_22605")
chRCC_sample   <- "SI_21561"

gse159115$tissue <- "ccRCC"
gse159115$tissue[gse159115$sample %in% normal_samples] <- "Normal"
gse159115$tissue[gse159115$sample == chRCC_sample]     <- "chRCC"

## Keep ccRCC tumor samples only
gse159115 <- subset(gse159115, subset = tissue == "ccRCC")

## Patient label used throughout the pipeline
gse159115$Patient <- gse159115$sample


############################################################
# PART C: MERGE BOTH DATASETS (10 patients)
############################################################

## Keep only the columns both datasets share
keep_cols <- c("orig.ident", "nCount_RNA", "nFeature_RNA", "Patient")
gse304466@meta.data <- gse304466@meta.data[, keep_cols]
gse159115@meta.data <- gse159115@meta.data[, keep_cols]

ccRCC_merged <- merge(gse159115, y = gse304466, project = "ccRCC_10")
rm(gse304466, gse159115); gc()

## Recompute mitochondrial and ribosomal percentages on the merged object
ccRCC_merged[["percent.mt"]] <- PercentageFeatureSet(ccRCC_merged, pattern = "^MT-")
ccRCC_merged[["percent.rb"]] <- PercentageFeatureSet(ccRCC_merged, pattern = "^RP[SL]")

## Final mitochondrial filter applied to all 10 patients
ccRCC_merged <- subset(ccRCC_merged, subset = percent.mt < 10)


## Save
saveRDS(ccRCC_merged, file.path(output_dir, "ccRCC_merged_10_patients.rds"))
