#!/bin/env Rscript
#SBATCH --cpus-per-task 5
#SBATCH --mem-per-cpu 4G

library(Seurat)

# Read in CSV file (expecting columns: Fastq_ID,group,replicate)
id_groups <- read.csv("samplesheet.csv")

# Create sample name by a combination of group and replicate
id_groups$sample=paste0(id_groups$group,"_",id_groups$replicate)

# Create empty list, which will contain all of our sample objects
all_samples <- list()

# Process all of the samples individually
for (samp_name in as.character(id_groups$sample)) {

    # Get the group that corresponds to this sample
    group <- id_groups$group[which(id_groups$sample == samp_name)]

    # Piece together the expected location of this sample's count and bc matrices
    matrix_dir <- paste0(samp_name,
                         "/outs/filtered_feature_bc_matrix/")

    # Read in the count data for this sample
    sample_counts <- Read10X(matrix_dir)

    # Create a Seurat object for this sample
    all_samples[[samp_name]] <-
        CreateSeuratObject(
            counts = sample_counts,
            project = samp_name, # metadata "orig.ident" field
            min.cells = 1,     # features found in fewer cells are excluded
            min.features = 200 # cells with fewer features (genes) are excluded
        )

    # Store sample and group names in the sample object itself
    all_samples[[samp_name]]$sample <- samp_name
    all_samples[[samp_name]]$group <- group

    # Determine mitochondrial reads per cell
    # make sure pattern matches mitochondrial gene names for the genome annotation (e.g. Humans may be "^MT-")
    all_samples[[samp_name]][["percent.mt"]] <-
        PercentageFeatureSet(all_samples[[samp_name]], pattern = "^MT-")

    # Determine ribosomal reads per cell
    # make sure pattern matches ribosome names for genome annotation
    all_samples[[samp_name]][["percent.ribo"]] <-
        PercentageFeatureSet(all_samples[[samp_name]], pattern = "^RP[LS]")

    # Save this sample's data
    saveRDS(all_samples[[samp_name]], file = paste0(samp_name, "_prefiltered.rds"))
}

# Start off with the first item as "merged". Others will be merged into it.
merged <- all_samples[[1]] # critical: use double square brackets

# Get all samples after the first
remaining_samples <- unlist(all_samples[2:length(all_samples)])

merged <- merge(merged, y = remaining_samples, add.cell.ids = as.character(id_groups$sample))

merged <- JoinLayers(merged)
saveRDS(merged,"Merged_Samples.rds")
