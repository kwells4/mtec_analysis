# This script is used to analyze the isotype control taken
# at the first timepoint.

library(Seurat)
library(mTEC.10x.pipeline)
library(dplyr)
library(scran)
library(DDRTree)

######## if on local, change this to current path to directory! ###########
#setwd("/Users/kristen/Documents/sshfs_qian/isoControlBeg")

# Grab name from snakemake object
seurat_mtec <- snakemake@input[[1]]
qc_pdf <- snakemake@output[[1]]
seurat_output <- snakemake@output[[2]]
image_pdf <- snakemake@params[[1]]
data_dir <- snakemake@params[[2]]
monocle_output <- snakemake@params[[3]]

# Load file

######### If on local, comment out ###########
load(seurat_mtec)

load(paste0(data_dir, "TFs.rda"))
load(paste0(data_dir, "gene_to_ensembl.rda"))
pdf(qc_pdf)

####### If on local, uncomment
#load("analysis_outs/seurat_isoControlBeg_empty.rda")
#pdf("analysis_outs/isoControlBeg_images.pdf")

##############################################################
# Set stage colors this is specific to the mtec dev project! #
##############################################################

# Must be the same order as stage_levels
stage_color_df <- data.frame("Cortico_medullary" = "#CC6600",
                             "Ccl21a_high" = "#009933",
                             "Early_Aire" = "#0066CC",
                             "Aire_positive" = "#660099",
                             "Late_Aire" = "#FF0000",
                             "Tuft" = "#990000",
                             "other" = "#666666")

stage_color <- t(stage_color_df)[ , 1]

stage_levels <- c("Cortico_medullary", "Ccl21a_high", "Early_Aire",
                  "Aire_positive", "Late_Aire", "Tuft", "other")

######################
# Initial processing #
######################
mtec <- add_perc_mito(mtec)
qc_plot(mtec)

mtec <- process_cells(mtec)

PC_plots(mtec, jackstraw = TRUE, test_pcs = 1:20)

# Determine dims_use from the jackstraw output
mtec <- group_cells(mtec, dims_use = 1:12)

# Change this based on your data!!!
mtec <- change_clus_ids(mtec, new_ids = c(3, 1, 2, 4, 0, 5, 6))

# Change this based on your data!!!
stage_list <- c("0" = "Ccl21a_high", "1" = "Ccl21a_high", "2" = "Early_Aire",
                "3" = "Aire_positive" , "4" = "Late_Aire", "5" = "Tuft",
                "6" = "Late_Aire")  

mtec <- set_stage(mtec, stage_list)

mtec@meta.data$stage <- factor(mtec@meta.data$stage,
                               levels = stage_levels) 

mtec <- Seurat::StashIdent(mtec, save.name = "seurat_cluster")

mtec <- Seurat::SetAllIdent(mtec, id = "stage")

mtec@ident <- factor(mtec@ident, levels = stage_levels)

# Find markers for each cluster
mtec.markers <- FindAllMarkers(object = mtec, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
top10 <- mtec.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = mtec, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)

# Find genes that are DE between each cluster (pairwise comparison)
print("find markers")
mtec <- significant_markers(mtec)

print("run cyclone")
mtec <- run_cyclone(mtec, gene_to_ensembl)

# Run monocle for pseudotime
pseudotime_mtec <- run_monocle(mtec)
save(pseudotime_mtec, file = monocle_output)

tSNE_PCA(mtec, "seurat_cluster", PCA = TRUE)
tSNE_PCA(mtec, "cluster", PCA = TRUE, color = stage_color)
tSNE_PCA(mtec, "nUMI")
tSNE_PCA(mtec, "nGene")

dev.off()

pdf(image_pdf)

tSNE_PCA(mtec, "seurat_cluster")
tSNE_PCA(mtec, "cluster", color = stage_color)
tSNE_PCA(mtec, "Aire")
tSNE_PCA(mtec, "Mki67")
tSNE_PCA(mtec, "Ascl1")
tSNE_PCA(mtec, "Hmgb2")
tSNE_PCA(mtec, "Dclk1")
tSNE_PCA(mtec, "Trpm5")
tSNE_PCA(mtec, "Tspan8")
tSNE_PCA(mtec, "Krt5")
tSNE_PCA(mtec, "Krt10")
tSNE_PCA(mtec, "Il25")
tSNE_PCA(mtec, "Il33")
tSNE_PCA(mtec, "Il13")
tSNE_PCA(mtec, "Top2a")
tSNE_PCA(mtec, "Hes6")
tSNE_PCA(mtec, "Hes1")
tSNE_PCA(mtec, "Fezf2")
tSNE_PCA(mtec, "cycle_phase", color = c("black", "red", "purple"))

PCA_loadings(mtec)
PCA_loadings(mtec, PC_val = "PC2")
TFs_all <- c("H2afz", "Top2a", "Hmgb1", "Hmgn1", "H2afx", as.character(TFs))
plot_heatmap(mtec, subset_list = TFs_all, cell_color = stage_color)

plot_heatmap(mtec, subset_list = TFs_all,
	color_list = c("Cdx1", "Utf1", "Tcf7", "Spib", "Cdk4", "Ptma",
		"H2afx", "Hmgb1"),
	color_list2 = c("Aire", "Irf7", "Cited2", "Spib", "Hes1", "Pax1",
                "Relb", "Lmo4", "Pou2f3"),
	cell_color = stage_color)

dev.off()

print(sessionInfo())

save(mtec, file = seurat_output)
