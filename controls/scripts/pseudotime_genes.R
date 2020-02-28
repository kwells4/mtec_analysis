library(mTEC.10x.pipeline)
library(slingshot)
library(gam)
library(pheatmap)
library(viridis)

load("/home/kwells4/mTEC_dev/mtec_snakemake/controls/analysis_outs/controls_merged_slingshot.rda")
load("/home/kwells4/mTEC_dev/mtec_snakemake/controls/analysis_outs/seurat_controls_merged.rda")

stage_color_df <- data.frame("Cortico_medullary" = "#CC6600", "Ccl21a_high" = "#009933",
                            "Early_Aire" = "#0066CC", "Aire_positive" = "#660099",
                            "Late_Aire" = "#FF0000", "Tuft" = "#990000",
                            "unknown" = "#FDFBFB")

stage_color <- t(stage_color_df)[ , 1]

mtec_wt_slingshot@assay$DE <- NULL
mtec_wt_slingshot <- Seurat::SetAllIdent(mtec_wt_slingshot, "stage")
mtec_wt_slingshot <- Seurat::SubsetData(mtec_wt_slingshot, ident.remove = "unknown")

pseudotime <- data.frame(slingPseudotime(sce.umap))

pdf("/home/kwells4/mTEC_dev/mtec_snakemake/controls/analysis_outs/pseudotime_plots.pdf")

# Aire branch
pseudotime_1 <- pseudotime[!is.na(pseudotime$curve1), ]
t_1 <- pseudotime_1$curve1

mtec_1 <- Seurat::SubsetData(mtec_wt_slingshot, cells.use = rownames(pseudotime_1))

Y_1 <- mtec_1@data

var3000 <- names(sort(apply(Y_1,1,var),decreasing = TRUE))[1:3000]

Y_1 <- Y_1[var3000,]

clusters_1 <- data.frame(stage = mtec_1@ident)



gam.pval <- apply(Y_1,1,function(z){
    d <- data.frame(z=z, t=t_1)
    suppressWarnings({
      tmp <- suppressWarnings(gam(z ~ lo(t_1), data=d))
    })
    p <- summary(tmp)[3][[1]][2,3]
    p
})

topgenes <- names(sort(gam.pval, decreasing = FALSE))[1:100]

heatdata <- Y_1[topgenes, order(t_1, na.last = NA)]
heatclus <- clusters_1[order(t_1, na.last = NA), , drop = FALSE]
heatmap_colors <- list(stage = stage_color)
pdf("/home/kwells4/mTEC_dev/mtec_snakemake/controls/analysis_outs/pseudotime_aire_branch.pdf")

p_1 <- pheatmap(scale(heatdata, center = TRUE),
	            color = viridis(10),
	            annotation_colors = heatmap_colors,
	            annotation_col = heatclus,
	            cluster_cols = FALSE,
	            cluster_rows = TRUE,
	            show_rownames = TRUE,
	            show_colnames = FALSE)

p_1

plot_list <- lapply(topgenes, function(x) plot_sling_pseudotime(
  seurat_object = mtec_wt_slingshot, sling_object = sce.umap, y_val = x,
  col_by = "stage", pseudotime_curve = "curve1", color = stage_color,
  range = c(0, 16)))

plot_list

dev.off()


# Ccl21a branch
pseudotime_2 <- pseudotime[!is.na(pseudotime$curve2), ]
t_2 <- pseudotime_2$curve2

mtec_2 <- Seurat::SubsetData(mtec_wt_slingshot, cells.use = rownames(pseudotime_2))

Y_2 <- mtec_2@data

var3000 <- names(sort(apply(Y_2,1,var),decreasing = TRUE))[1:3000]

Y_2 <- Y_2[var3000,]

clusters_2 <- data.frame(stage = mtec_2@ident)



gam.pval <- apply(Y_2,1,function(z){
    d <- data.frame(z=z, t=t_2)
    suppressWarnings({
      tmp <- suppressWarnings(gam(z ~ lo(t_2), data=d))
    })
    p <- summary(tmp)[3][[1]][2,3]
    p
})

topgenes <- names(sort(gam.pval, decreasing = FALSE))[1:100]

heatdata <- Y_2[topgenes, order(t_2, na.last = NA)]
heatclus <- clusters_2[order(t_2, na.last = NA), , drop = FALSE]
heatmap_colors <- list(stage = stage_color)
pdf("/home/kwells4/mTEC_dev/mtec_snakemake/controls/analysis_outs/pseudotime_ccl21a_branch.pdf")

p_2 <- pheatmap(scale(heatdata, center = TRUE),
	            color = viridis(10),
	            annotation_colors = heatmap_colors,
	            annotation_col = heatclus,
	            cluster_cols = FALSE,
	            cluster_rows = TRUE,
	            show_rownames = TRUE,
	            show_colnames = FALSE)

p_2

plot_list <- lapply(topgenes, function(x) plot_sling_pseudotime(
  seurat_object = mtec_wt_slingshot, sling_object = sce.umap, y_val = x,
  col_by = "stage", pseudotime_curve = "curve2", color = stage_color,
  range = c(0, 16)))

plot_list

dev.off()


#Tuft branch
pseudotime_3 <- pseudotime[!is.na(pseudotime$curve3), ]
t_3 <- pseudotime_3$curve3

mtec_3 <- Seurat::SubsetData(mtec_wt_slingshot, cells.use = rownames(pseudotime_3))

Y_3 <- mtec_3@data

var3000 <- names(sort(apply(Y_3,1,var),decreasing = TRUE))[1:3000]

Y_3 <- Y_3[var3000,]

clusters_3 <- data.frame(stage = mtec_3@ident)



gam.pval <- apply(Y_3,1,function(z){
    d <- data.frame(z=z, t=t_3)
    suppressWarnings({
      tmp <- suppressWarnings(gam(z ~ lo(t_3), data=d))
    })
    p <- summary(tmp)[3][[1]][2,3]
    p
})

topgenes <- names(sort(gam.pval, decreasing = FALSE))[1:100]

heatdata <- Y_3[topgenes, order(t_3, na.last = NA)]
heatclus <- clusters_3[order(t_3, na.last = NA), , drop = FALSE]
heatmap_colors <- list(stage = stage_color)
pdf("/home/kwells4/mTEC_dev/mtec_snakemake/controls/analysis_outs/pseudotime_tuft_branch.pdf")

p_3 <- pheatmap(scale(heatdata, center = TRUE),
	            color = viridis(10),
	            annotation_colors = heatmap_colors,
	            annotation_col = heatclus,
	            cluster_cols = FALSE,
	            cluster_rows = TRUE,
	            show_rownames = TRUE,
	            show_colnames = FALSE)

p_3

plot_list <- lapply(topgenes, function(x) plot_sling_pseudotime(
  seurat_object = mtec_wt_slingshot, sling_object = sce.umap, y_val = x,
  col_by = "stage", pseudotime_curve = "curve3", color = stage_color,
  range = c(0, 16)))

plot_list

dev.off()
pdf("/home/kwells4/mTEC_dev/mtec_snakemake/controls/analysis_outs/pseudotime_interest_genes.pdf")
Ccl21a_genes <- c("Ltb", "Ccl21a", "Ascl1", "Hmgb2", "Stmn1")

aire_genes <- c("Tnfrsf11a", "Traf3", "Traf6", "Nfkb1", "Creb1", "Nfatc1", "Aire", "Fezf2",
	"Ascl1", "Hmgb2", "Ascl1")

plot_list <- lapply(Ccl21a_genes, function(x) plot_sling_pseudotime(
  seurat_object = mtec_wt_slingshot, sling_object = sce.umap, y_val = x,
  col_by = "stage", pseudotime_curve = "curve2", color = stage_color,
  range = c(0, 16)))
plot_list

plot_list <- lapply(aire_genes, function(x) plot_sling_pseudotime(
  seurat_object = mtec_wt_slingshot, sling_object = sce.umap, y_val = x,
  col_by = "stage", pseudotime_curve = "curve1", color = stage_color,
  range = c(0, 16)))
plot_list

dev.off()