#!/bin/env -S Rscript --vanilla
#SBATCH --cpus-per-task 4
#SBATCH --mem-per-cpu 20G


# Load needed libraries
library(Seurat)
library(dplyr) # Provide pipe operator %>%
library(harmony)
library(scDblFinder)
library(ggplot2)

# Read in Seurat object
obj <- readRDS("filtered_merged_unintegrated.rds")

# Split merged object it into separate sample layers
# We must split before we Join and we must Join after we split
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$sample)

# NormalizeData         Normalize gene expression to 10000 counts/cell, then for
#                         each gene expression value x, return natural log of x + 1.
#                         (Each value will be in the range of 0 to 9.21044)
# FindVariableFeatures  Pick the most variable genes.
# ScaleData             Shift and scale each gene's expression across all the
#                         data to have an average of 0 and a variance of 1.
#                         The default is to work only on variable features.
# RunPCA                Dimensionality reduction which summarizes expression
#                         values for thousands of "variable" genes in a few
#                         (e.g., 50) principle components, thus making
#                         subsequent calculations feasible.
# FindNeighbors         Calculate neighbor graphs using PCA space.
# FindClusters          Find clusters within the neighbor graph.
# RunUMAP               Cartoonishly squash 50 dimensions into 2 dimensions
#                         for visualization purposes.
#

# dimensionality reduction (get PCA),
# also we need "clusters" to remove doublets.
# This will be run again after doublet removal.
obj[["Data"]] <- NULL
obj <- obj %>%
        NormalizeData() %>% # why do we need to do this again? This is "cell-level" normalization, after all.
        FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>% # Redo since some variance lost to filtering
        ScaleData(features = row.names(.)) %>%                               # ditto
        RunPCA( reduction.name = "pca.filtered") %>%
        FindNeighbors(reduction = "pca.filtered", dims = 1:30) %>%
        FindClusters(resolution = 0.4, reduction = "pca.filtered", cluster.name = "clusters.unintegrated") %>%
        RunUMAP(reduction = "pca.filtered",
                dims = 1:30,
                reduction.name = "umap.unintegrated")

# OR equivalently, without pipes %>%
# obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000)
# obj <- ScaleData(obj, features = row.names(obj))
# obj <- RunPCA(obj)
# obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30)
# obj <- FindClusters(obj, resolution = 0.4, cluster.name = "clusters.unintegrated")
# obj <- RunUMAP(obj,
#                reduction = "pca",
#                dims = 1:30,
#                reduction.name = "umap.unintegrated")

# Now we rejoin to run scDblFinder
obj=JoinLayers(obj)

### scDblFinder works on sce (not Seurat) object, so we need to make a new object and convert.
obj.sce=as.SingleCellExperiment(obj)

### Run scDblFinder on the sce object to identify which are singlets and which are doublets
obj.default <- scDblFinder(
                   obj.sce,
                   cluster = obj.sce@colData@listData[["clusters.unintegrated"]],
                   samples = obj$sample
               )

### We assign the singlet/doublet status to the seurat object meta data.
obj@meta.data$doublet <- obj.default@colData@listData[["scDblFinder.class"]]

### Subset to only the singlets moving forward (Only overwrite if you are SURE you won't need old data)
obj <- subset(obj, subset = doublet == "singlet")

message("Find clusters after removing doublets")

# Now we reprocess the data (because it's different now)
obj <- obj %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>% # Refind variable features since some variableness may be lost with doublet removal
  ScaleData(features = row.names(.)) %>%                               # Recale features for the same reason
  RunPCA(reduction.name = "pca.doublets_removed") %>%
  FindNeighbors(reduction = "pca.doublets_removed", dims = 1:30) %>%
  FindClusters(resolution = 0.4, cluster.name = "clusters.doublets_removed") %>%
  RunUMAP(reduction = "pca.doublets_removed",
          dims = 1:30,
          reduction.name = "umap.doublets_removed")

### The data should be ready for "integration".
### Which is basically fancy batch correction.
### We need to split samples for integration
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$sample)

obj <- IntegrateLayers(object = obj,
                       method = CCAIntegration,
                       orig.reduction = "pca.doublets_removed",
                       new.reduction = "cca",
                       verbose = TRUE)

obj <- IntegrateLayers(object = obj, 
                       method = HarmonyIntegration,
                       orig.reduction = "pca.doublets_removed",
                       new.reduction = "harmony",
                       verbose = TRUE)

obj <- JoinLayers(obj)

# Now find clusters on integrated data (two ways)
obj <- obj %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%
  FindClusters(resolution = 0.4, cluster.name = "clusters.harmony") %>%
  RunUMAP(reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony") %>%
  FindNeighbors(reduction = "cca", dims = 1:30) %>%
  FindClusters(resolution = 0.4, cluster.name = "clusters.cca") %>%
  RunUMAP(reduction = "cca",     dims = 1:30, reduction.name = "umap.cca")

saveRDS(obj, "DoubletsRemoved_filtered_integrated.rds")
