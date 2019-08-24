library(mTEC.10x.pipeline)

library(Seurat)
library(dplyr)
library(scran)
library(DDRTree)

load("/home/kwells4/mTEC_dev/mtec_snakemake/aireTrace/analysis_outs/seurat_aireTrace_empty.rda")
data_dir <- "/home/kwells4/mTEC_dev/data/"
load(paste0(data_dir, "TFs.rda"))
load(paste0(data_dir, "gene_to_ensembl.rda"))

stage_color_df <- data.frame("Cortico_medullary" = "#CC6600",
                             "Ccl21a_high" = "#009933",
                             "Early_Aire" = "#0066CC",
                             "Aire_positive" = "#660099",
                             "Late_Aire" = "#FF0000",
                             "Tuft" = "#990000",
                             "unknown" = "#666666")

stage_color <- t(stage_color_df)[ , 1]

stage_levels <- c("Cortico_medullary", "Ccl21a_high", "Early_Aire",
                  "Aire_positive", "Late_Aire", "Tuft", "unknown")


pdf("/home/kwells4/mTEC_dev/mtec_snakemake/reproducability_test.pdf")

# Add mitochondiral percent to the meta data
mtec <- add_perc_mito(mtec)

# Plot quality plots
qc_plot(mtec)

# remove low quality cells, normalize and scale data, find variable genes, and perform PCA 
mtec <- process_cells(mtec)

# Plot a PCA and determine dimensions to use by looking at heatmaps and elbow plots
PC_plots(mtec, jackstraw = FALSE, test_pcs = 1:20)

PC_plots(mtec, jackstraw = TRUE, test_pcs = 1:20)

# Decide PCs based on output of jackstraw plot
mtec <- group_cells(mtec, dims_use = 1:12)

plot <- tSNE_PCA(mtec, "cluster")

# Name clusters (these were determined originally based on gene expression)
stage_list <- c("0" = "Ccl21a_high", "1" = "Ccl21a_high", "2" = "Ccl21a_high",
                "3" = "Late_Aire", "4" = "Cortico_medullary",
                "5" = "Early_Aire", "6" = "Aire_positive",
                "7" = "Tuft", "8" = "unknown")

mtec <- set_stage(mtec, stage_list)

mtec@meta.data$stage <- factor(mtec@meta.data$stage,
                               levels = stage_levels)

plot_2 <- tSNE_PCA(mtec, "stage", color = stage_color)

dev.off()