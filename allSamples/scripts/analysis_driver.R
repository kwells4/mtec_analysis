library(mTEC.10x.pipeline)
library(slingshot)
library(RColorBrewer)

set_stage <- function(mtec, current_cluster){
  stage_name <- sort(table(mtec@meta.data$at_stage[mtec@meta.data$res.0.6 ==
    current_cluster]), decreasing = TRUE)[1]
  stage <- names(stage_name)
  names(stage) <- current_cluster
  return(stage)
}

input_object <- snakemake@input[[1]]
mapping_object <- snakemake@input[[2]]
output_pdf <- snakemake@params[[1]]
data_dir <- snakemake@params[[2]]
output_object <- snakemake@output[[1]]
slingshot_obj <- snakemake@output[[2]]

pdf(output_pdf)

# Set colors for the plots
stages_colors <- data.frame("Cortico_medullary" = "#CC6600", "Ccl21a_high" = "#009933",
                            "Early_Aire" = "#0066CC", "Aire_positive" = "#660099",
                            "Late_Aire" = "#FF0000", "Tuft" = "#990000",
                            "unknown" = "#666666")

# Make the colors into a readable form for plotting
stage_color <- t(stages_colors)[, 1]

# Load in data
mtecCombined <- get(load(input_object))

aireTrace <- get(load(mapping_object))

load(paste0(data_dir, "TFs.rda"))
load(paste0(data_dir, "gene_to_ensembl.rda"))

aireTrace@meta.data$stage <- factor(aireTrace@meta.data$stage, levels = c(levels(aireTrace@meta.data$stage), "unknown"))
aireTrace@meta.data$stage[is.na(aireTrace@meta.data$stage)] <- "unknown"

aire_trace_stage_orig <- data.frame(at_stage = aireTrace@meta.data$stage,
	                                row.names = rownames(aireTrace@meta.data))

rownames(aire_trace_stage_orig) <- paste0("aireTrace_",
	                                      rownames(aire_trace_stage_orig))

new_meta_data <- mtecCombined@meta.data

new_meta_data <- merge(new_meta_data, aire_trace_stage_orig,
	                   by = "row.names", all.x = TRUE)

rownames(new_meta_data) <- new_meta_data$Row.names
new_meta_data$Row.names <- NULL

if (!identical(rownames(mtecCombined@meta.data), rownames(new_meta_data))) {
  print("must reorder genes")
  new_meta_data <- new_meta_data[match(rownames(mtecCombined@meta.data),
                               rownames(new_meta_data)), , drop = FALSE]
}

mtecCombined@meta.data <- new_meta_data


# Name the clusters based on the aire trace data. This blindly names based on
# the most common stage associated with each cluster, but this should be visually
# checked to ensure accuracy
stage_names <- lapply(unique(mtecCombined@meta.data$res.0.6), function(x)
	set_stage(mtecCombined, x))
stage_names <- unlist(stage_names)

mtecCombined@meta.data$stage <- stage_names[mtecCombined@meta.data$res.0.6]

mtecCombined@meta.data$stage <- factor(mtecCombined@meta.data$stage,
                                     levels = c("Cortico_medullary", "Ccl21a_high",
                                                "Early_Aire", "Aire_positive",
                                                "Late_Aire", "Tuft",
                                                "unknown"))

print("run_cyclone")
mtecCombined <- run_cyclone(mtecCombined, gene_to_ensembl)

print("run slingshot")
# Remove unknown cells
if (!identical(rownames(mtecCombined@meta.data), names(mtecCombined@ident))) {
  print("must reorder genes")
  mtecCombined@meta.data <- mtecCombined@meta.data[match(rownames(mtecCombined@ident),
                               rownames(mtecCombined@meta.data)), , drop = FALSE]
}

mtecCombined_full <- mtecCombined

mtecCombined <- Seurat::SetAllIdent(mtecCombined, id = "stage")
mtecCombined <- Seurat::SubsetData(mtecCombined, ident.remove = "unknown")

# Pull out dataframes for slingshot to use UMAP embeddings
cluster_labels <- mtecCombined@meta.data$stage

slingshot_umap <- mtecCombined@dr$umap@cell.embeddings

sce.umap <- slingshot(slingshot_umap, clusterLabels = cluster_labels,
  start.clus = "Early_Aire")



# Pull out pseudotime information
pseudotime <- data.frame(slingPseudotime(sce.umap))
pseudotime_cols1 <- pseudotime$curve1
pseudotime_cols2 <- pseudotime$curve2
pseudotime_cols3 <- pseudotime$curve3

# Plot pseudotime
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plot(slingshot_umap, col = stage_color, pch=16, asp = 1)
lines(sce.umap, lwd=2, type = "curves")

plot(slingshot_umap, col = colors[cut(pseudotime_cols1,breaks=100)], pch=16, asp = 1)
lines(sce.umap, lwd=2, type = "curves")

plot(slingshot_umap, col = colors[cut(pseudotime_cols2,breaks=100)], pch=16, asp = 1)
lines(sce.umap, lwd=2, type = "curves")

plot(slingshot_umap, col = colors[cut(pseudotime_cols3,breaks=100)], pch=16, asp = 1)
lines(sce.umap, lwd=2, type = "curves")

tSNE_PCA(mtecCombined, "Aire")
tSNE_PCA(mtecCombined, "Mki67")
tSNE_PCA(mtecCombined, "Ascl1")
tSNE_PCA(mtecCombined, "Hmgb2")
tSNE_PCA(mtecCombined, "Dclk1")
tSNE_PCA(mtecCombined, "Trpm5")
tSNE_PCA(mtecCombined, "Tspan8")
tSNE_PCA(mtecCombined, "Krt5")
tSNE_PCA(mtecCombined, "Krt10")
tSNE_PCA(mtecCombined, "Il25")
tSNE_PCA(mtecCombined, "Il33")
tSNE_PCA(mtecCombined, "Il13")
tSNE_PCA(mtecCombined, "Top2a")
tSNE_PCA(mtecCombined, "Hes6")
tSNE_PCA(mtecCombined, "Hes1")
tSNE_PCA(mtecCombined, "Fezf2")
tSNE_PCA(mtecCombined, "Ccl21a")
tSNE_PCA(mtecCombined, "stage", color = stage_color)
tSNE_PCA(mtecCombined, "cycle_phase", color = c("black", "red", "purple"))

dev.off()

save(sce.umap, file = slingshot_obj)
save(mtecCombined_full, file = output_object)
