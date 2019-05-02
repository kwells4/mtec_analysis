library(mTEC.10x.pipeline)
library(slingshot)
library(dplyr)
library(RColorBrewer)

input_object <- snakemake@input[[1]]
output_pdf <- snakemake@params[[1]]
data_dir <- snakemake@params[[2]]
output_object <- snakemake@output[[1]]
slingshot_obj <- snakemake@output[[2]]

# Load in data
mtec_wt <- get(load(input_object))

load(paste0(data_dir, "TFs.rda"))
load(paste0(data_dir, "gene_to_ensembl.rda"))

stage_list <- c("0" = "Aire_positive", "1" = "Late_Aire", "2" = "Early_Aire",
                  "3" = "Aire_positive" , "4" = "Ccl21a_high", "5" = "Cortico_medullary",
                  "6" = "Late_Aire","7" = "Tuft", "8" = "unknown",
                  "9" = "unknown")


mtec_wt <- set_stage(mtec_wt, stage_list)

stage_levels <- c("Cortico_medullary", "Ccl21a_high", "Early_Aire",
	              "Aire_positive", "Late_Aire", "Tuft", "unknown")

mtec_wt@meta.data$stage <- factor(mtec_wt@meta.data$stage,
                               levels = stage_levels)

# Colors for plotting
stage_color_df <- data.frame("Cortico_medullary" = "#CC6600", "Ccl21a_high" = "#009933",
                            "Early_Aire" = "#0066CC", "Aire_positive" = "#660099",
                            "Late_Aire" = "#FF0000", "Tuft" = "#990000",
                            "unknown" = "#FDFBFB")

stage_color <- t(stage_color_df)[ , 1]

pdf(output_pdf)
mtec_markers <- Seurat::FindAllMarkers(object = mtec_wt, only.pos = TRUE,
	                                   min.pct = 0.25, thresh.use = 0.25)

top10 <- mtec_markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
Seurat::DoHeatmap(object = mtec_wt, genes.use = top10$gene, slim.col.label = TRUE,
	      remove.key = TRUE)

# Find genes that are DE between each cluster (pairwise comparison)
print("find markers")
mtec_wt <- significant_markers(mtec_wt)

print("run cyclone")
mtec_wt <- run_cyclone(mtec_wt, gene_to_ensembl)

# Run slingshot
print("run slingshot")
# Remove unknown cells because we know these are not in the same lineage as mTECs
mtec_wt_slingshot <- mtec_wt
mtec_wt_slingshot@assay$DE <- NULL
mtec_wt_slingshot <- Seurat::SetAllIdent(mtec_wt_slingshot, "stage")
mtec_wt_slingshot <- Seurat::SubsetData(mtec_wt_slingshot, ident.remove = "unknown")

# Pull out dataframes for slingshot to use UMAP embeddings
cluster_labels <- mtec_wt_slingshot@meta.data$stage
slingshot_umap <- mtec_wt_slingshot@dr$umap@cell.embeddings

sce.umap <- slingshot(slingshot_umap, clusterLabels = cluster_labels,
	start.clus = "Early_Aire", extend = "n")

umap.lineages <- getLineages(slingshot_umap, clusterLabels = cluster_labels,
  start.clus = "Early_Aire")

# Pull out pseudotime information
pseudotime <- data.frame(slingPseudotime(sce.umap))
pseudotime_cols1 <- pseudotime$curve1
pseudotime_cols2 <- pseudotime$curve2
pseudotime_cols3 <- pseudotime$curve3

# Plot pseudotime
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
# plot(slingshot_umap, col = stage_color, pch=16, asp = 1)
# lines(sce.umap, lwd=2, type = "curves")

# plot(slingshot_umap, col = stage_color, pch=16, asp = 1)
# lines(umap.lineages, lwd=2)

plot(slingshot_umap, col = colors[cut(pseudotime_cols1,breaks=100)], pch=16, asp = 1)
lines(umap.lineages, lwd=2)

plot(slingshot_umap, col = colors[cut(pseudotime_cols1,breaks=100)], pch=16, asp = 1)
lines(sce.umap, lwd=2, type = "curves")

plot(slingshot_umap, col = colors[cut(pseudotime_cols2,breaks=100)], pch=16, asp = 1)
lines(sce.umap, lwd=2, type = "curves")

plot(slingshot_umap, col = colors[cut(pseudotime_cols3,breaks=100)], pch=16, asp = 1)
lines(sce.umap, lwd=2, type = "curves")



tSNE_PCA(mtec_wt, "Aire")
tSNE_PCA(mtec_wt, "Mki67")
tSNE_PCA(mtec_wt, "Ascl1")
tSNE_PCA(mtec_wt, "Hmgb2")
tSNE_PCA(mtec_wt, "Dclk1")
tSNE_PCA(mtec_wt, "Trpm5")
tSNE_PCA(mtec_wt, "Tspan8")
tSNE_PCA(mtec_wt, "Krt5")
tSNE_PCA(mtec_wt, "Krt10")
tSNE_PCA(mtec_wt, "Il25")
tSNE_PCA(mtec_wt, "Il33")
tSNE_PCA(mtec_wt, "Il13")
tSNE_PCA(mtec_wt, "Top2a")
tSNE_PCA(mtec_wt, "Hes6")
tSNE_PCA(mtec_wt, "Hes1")
tSNE_PCA(mtec_wt, "Fezf2")
tSNE_PCA(mtec_wt, "Ccl21a")
tSNE_PCA(mtec_wt, "stage", color = stage_color)
tSNE_PCA(mtec_wt, "cycle_phase", color = c("black", "red", "purple"))

dev.off()
save(sce.umap, file = slingshot_obj)
save(mtec_wt, file = output_object)
