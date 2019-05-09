library(mTEC.10x.pipeline)

working_dir <- "/home/kwells4/mTEC_dev/mtec_snakemake/"

out_file <- paste0(working_dir, "not_yet_included/seurat_info.csv")

load(paste0(working_dir, "controls/analysis_outs/seurat_controls_merged.rda"))

clusters <- data.frame(clusters = mtec_wt@meta.data$stage,
	                   row.names = rownames(mtec_wt@meta.data))
umap_data <- data.frame(mtec_wt@dr$umap@cell.embeddings)

seurat_data <- merge(clusters, umap_data, by = "row.names")

name_1 <- regmatches(seurat_data$Row.names, regexpr(".+_", seurat_data$Row.names))
name_2 <- regmatches(seurat_data$Row.names, regexpr("_.+", seurat_data$Row.names))
name_2 <- sub("_", ":", name_2)
rownames(seurat_data) <- paste0(name_1, "count", name_2, "x")

seurat_data$Row.names <- NULL

write.csv(seurat_data, out_file)