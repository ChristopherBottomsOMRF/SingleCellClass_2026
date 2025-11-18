#!/bin/env -S Rscript --vanilla
#SBATCH --cpus-per-task 4
#SBATCH --mem-per-cpu 20G

library(Seurat)
library(dplyr) # Provide pipe operator %>%
library(ggplot2)

### read in the seurat object with all 4 integrations
obj=readRDS("DoubletsRemoved_filtered_integrated.rds")

### Up to this point, we have only created the integrations based on PCs.
### We have no projected these new integrations into a 2d plane or viewing (UMAP)
#obj=RunUMAP(obj,reduction = "cca", dims = 1:30, reduction.name = "umap.cca")
#obj=RunUMAP(obj,reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")

pdf("UMAPS.pdf")
DimPlot(obj, reduction="umap.unintegrated", group.by = "group")
DimPlot(obj, reduction="umap.cca",          group.by = "group")
DimPlot(obj, reduction="umap.harmony",      group.by = "group")
DimPlot(obj, reduction="umap.unintegrated", group.by = "sample")
DimPlot(obj, reduction="umap.cca",          group.by = "sample")
DimPlot(obj, reduction="umap.harmony",      group.by = "sample")
dev.off()

#saveRDS(obj,"/hpc-prj/cbds/SingleCellClass/Objects/DoubletsRemoved_filtered_integrated4ways_withUMAPS.rds")













