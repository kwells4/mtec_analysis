library(mTEC.10x.pipeline)
source("/home/kwells4/mTEC_dev/mtec_snakemake/scripts/figure_funcs.R")

ggplot2::theme_set(ggplot2::theme_classic(base_size = 18))

#load("/home/kwells4/mTEC_dev/mtec_snakemake/controls/analysis_outs/seurat_controls_merged.rda")
mtec_wt <- readRDS("/home/kwells4/mTEC_dev/mtec_snakemake/controls/analysis_outs/seurat_controls_merged.rda")
percent_ident <- function(seurat_object, data_set, meta_data_col, ident){
  cells_use <- rownames(seurat_object@meta.data)[
    seurat_object@meta.data[[meta_data_col]] == data_set]
  new_seurat <- Seurat::SubsetData(seurat_object, cells.use = cells_use)
  ident_cells <- table(new_seurat@meta.data[[ident]])
  print(ident_cells)
  ident_percent <- (ident_cells[["TRUE"]])/nrow(new_seurat@meta.data) * 100
  #names(ident_percent) <- data_set
  return(ident_percent)
}

stage_color_df <- data.frame("Cortico_medullary" = "#CC6600", "Ccl21a_high" = "#009933",
                            "Early_Aire" = "#0066CC", "Aire_positive" = "#660099",
                            "Late_Aire" = "#FF0000", "Tuft" = "#990000",
                            "unknown" = "#FDFBFB")

stage_color <- t(stage_color_df)[ , 1]


pdf("/home/kwells4/mTEC_dev/mtec_snakemake/controls/analysis_outs/revision_plots.pdf")

#Subset out "unknown" cells

mtec_wt <- Seurat::SetAllIdent(mtec_wt, "stage")
mtec_wt@assay$DE <- NULL
mtec_wt <- Seurat::SubsetData(mtec_wt, ident.remove = "unknown")

# In WT do Mki67 percent in Aire low (0) and Aire high (>0)


plotDimRed(mtec_wt, col_by = "Mki67")
plotDimRed(mtec_wt, col_by = "Pdpn")
plotDimRed(mtec_wt, col_by = "Hmmr")
plotDimRed(mtec_wt, col_by = "stage", color = stage_color, show_legend = FALSE)


mtec_wt@meta.data$Ccl21a <- mtec_wt@data["Ccl21a", ]
mtec_wt@meta.data$Aire <- mtec_wt@data["Aire", ]
mtec_wt@meta.data$Mki67 <- mtec_wt@data["Mki67", ]
mtec_wt@meta.data$Hmgb2 <- mtec_wt@data["Hmgb2", ]
mtec_wt@meta.data$Stmn1 <- mtec_wt@data["Stmn1", ]
mtec_wt@meta.data$Top2a <- mtec_wt@data["Top2a", ]

mtec_wt@meta.data$db_pos <- mtec_wt@meta.data$Ccl21a > 4 & 
  mtec_wt@meta.data$Aire > 0.5


mtec_wt@meta.data$Mki67_pos <- mtec_wt@meta.data$Mki67 > 0.25

plotDimRed(mtec_wt, col_by = "db_pos", color = c("#000000", "#FF007F"))





plot_trio_subset <- function(seurat_object, subset_on = NULL, ...) {
	if (!is.null(subset_on)) {
		seurat_object <- Seurat::SubsetData(seurat_object, ident.use = subset_on)
	}
	trio_plots(seurat_object, plot_name = subset_on, ...)
}

geneset <- c("Ccl21a", "Ccl21b.1", "Ccl21a.1")
geneset_2 <- c("Aire", "Ccl21a", "Mki67")



trio_plots(mtec_wt, plot_name = "Ccl21 genes", geneset = geneset,
  cell_cycle = FALSE, plot_violin = TRUE, plot_jitter = FALSE, jitter_and_violin = FALSE,
  color = stage_color, sep_by = "stage")

trio_plots(mtec_wt, plot_name = "Ccl21 genes", geneset = geneset_2,
  cell_cycle = FALSE, plot_violin = TRUE, plot_jitter = FALSE, jitter_and_violin = FALSE,
  color = stage_color, sep_by = "stage")


trio_plots(mtec_wt, plot_name = "Cycling cells", geneset = c("Mki67"),
  cell_cycle = FALSE, plot_violin = TRUE, plot_jitter = FALSE, jitter_and_violin = FALSE,
  color = stage_color, sep_by = "stage")

plotting_data <- mtec_wt@meta.data

plotDimRed(mtec_wt, col_by = "cycle_phase", color = c("black", "red", "purple"))

plotDimRed(mtec_wt, col_by = "cycle_phase", color = c("black", "red", "purple"),
  show_legend = FALSE)

aire_ccl21a <- ggplot2::ggplot(plotting_data, ggplot2::aes(x = Aire, y = Ccl21a,
  color = stage)) +
  ggplot2::geom_point() +
  ggplot2::scale_color_manual(values = stage_color)

ki67_aire <- ggplot2::ggplot(plotting_data, ggplot2::aes(x = Mki67, y = Aire,
  color = stage)) +
  ggplot2::geom_point() +
  ggplot2::scale_color_manual(values = stage_color)

ki67_ccl21a <- ggplot2::ggplot(plotting_data, ggplot2::aes(x = Mki67, y = Ccl21a,
  color = stage)) +
  ggplot2::geom_point() +
  ggplot2::scale_color_manual(values = stage_color)


aire_ccl21a
ki67_aire
ki67_ccl21a

ggplot2::ggplot(plotting_data, ggplot2::aes(x = Hmgb2, y = Ccl21a,
  color = stage)) +
  ggplot2::geom_point() +
  ggplot2::scale_color_manual(values = stage_color)

ggplot2::ggplot(plotting_data, ggplot2::aes(x = Stmn1, y = Ccl21a,
  color = stage)) +
  ggplot2::geom_point() +
  ggplot2::scale_color_manual(values = stage_color)

ggplot2::ggplot(plotting_data, ggplot2::aes(x = Top2a, y = Ccl21a,
  color = stage)) +
  ggplot2::geom_point() +
  ggplot2::scale_color_manual(values = stage_color)


mtec_wt_TAC <- Seurat::SetAllIdent(mtec_wt, "stage")
mtec_wt_TAC <- Seurat::SubsetData(mtec_wt_TAC, ident.use = "Early_Aire")

plotting_data_TAC <- mtec_wt_TAC@meta.data
ggplot2::ggplot(plotting_data_TAC, ggplot2::aes(x = Aire, y = Ccl21a,
  color = exp)) +
  ggplot2::geom_point() +
  ggplot2::scale_color_manual(values = c("#FC8D62", "#8DA0CB"))

dev.off()

