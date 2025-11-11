#!/bin/env -S Rscript --vanilla
#SBATCH --cpus-per-task 20
#SBATCH --mem-per-cpu 4G

library(Seurat)
library(ggplot2)
### Read in the unfiltered, merged object.
### You will need to change the file path to wherever 
### your object is located.
### You may need to "mount" or "map" the HPC to your computer
unfiltered_merged=readRDS("Merged_Samples.rds")

pdf("plots.pdf")
### change the yintercepts values to help visualize appropriate cut points
print(
    VlnPlot(unfiltered_merged, features="nFeature_RNA", group.by = "group", pt.size = 0) +
    NoLegend() +
    geom_hline(yintercept = 500, col = "red") +
    geom_hline(yintercept = 9000, col = "blue")
)

print(
    VlnPlot(unfiltered_merged,features="nCount_RNA",group.by="group",pt.size=0) +
    NoLegend() +
    geom_hline(yintercept = 1000, col = "red") +
    geom_hline(yintercept = 90000, col = "blue")
)

print(
    VlnPlot(unfiltered_merged,features = "percent.mt", group.by = "group", pt.size = 0) +
    NoLegend() +
    geom_hline(yintercept = 5, col = "red")
)

print(
    VlnPlot(unfiltered_merged,features = "percent.ribo", group.by = "group", pt.size = 0) +
    NoLegend() +
    geom_hline(yintercept = 10, col = "red")
)

### nFeatures and nCount should be highly correlated
print(
    FeatureScatter(unfiltered_merged,
                   feature1 = "nFeature_RNA",
                   feature2 = "nCount_RNA",
                   pt.size = 0)
) # We expect high correlation

### show data with threshold set to get a sense of how much we are losing
print(FeatureScatter(unfiltered_merged,
                     feature1 = "percent.mt",
                     feature2 = "percent.ribo",
                     pt.size = 0) +
      NoLegend() +
      geom_hline(yintercept = 4, col = "red") +
      geom_vline(xintercept = 15,col = "blue")
)

dev.off()

### We want to keep that upper left hand corner: Low mitochondrial reads, high ribosomal reads
### We can do this with the "subset" function.

### Usage: new_obj=subset(old_obj, subset=meta_data_name *relation* value)
### example: keloid_subjects=subset(unfiltered_merged, subset=group=="keloid")
filtered_merged_unintegrated <- subset(unfiltered_merged,
                                       subset = nFeature_RNA > 200 &
                                                nFeature_RNA < 9000 &
                                                percent.mt < 15 &
                                                percent.ribo > 4)

### save the filtered seurat object to file for later use
saveRDS(filtered_merged_unintegrated, "filtered_merged_unintegrated.rds")
