library(mTEC.10x.pipeline)
source("/home/kwells4/mTEC_dev/mtec_snakemake/scripts/figure_funcs.R")

load("/home/kwells4/mTEC_dev/mtec_snakemake/allSamples/analysis_outs/seurat_allSamples_combined.rda")

timepoints <- c("isoControlBeg", "isoControlEnd", "timepoint1",
	"timepoint2", "timepoint3", "timepoint5")

percent_ident <- function(seurat_object, data_set, meta_data_col, ident){
  cells_use <- rownames(seurat_object@meta.data)[
    seurat_object@meta.data[[meta_data_col]] == data_set]
  new_seurat <- Seurat::SubsetData(seurat_object, cells.use = cells_use)
  ident_cells <- table(new_seurat@meta.data[[ident]])
  print(ident_cells)
  if("TRUE" %in% names(ident_cells)){
    ident_percent <- (ident_cells[["TRUE"]])/nrow(new_seurat@meta.data) * 100
  }
  else{
    ident_percent <- 0
  }
  #names(ident_percent) <- data_set
  return(ident_percent)
}



timecourse_color <- c("#FC8D62", "#8DA0CB", "#E78AC3", "#65C42D",
  "#FFB62F", "#E5C494")

pdf("/home/kwells4/mTEC_dev/mtec_snakemake/allSamples/analysis_outs/revision_plots.pdf")

#Subset out "unknown" cells

mtecCombined_full <- Seurat::SetAllIdent(mtecCombined_full, "stage")
mtecCombined <- Seurat::SubsetData(mtecCombined_full, ident.remove = "unknown")

# In WT do Mki67 percent in Aire low (0) and Aire high (>0)

lapply(timepoints, function(x) full_umap(mtecCombined,
  data_set = x, col_by = "Mki67",
  show_legend = TRUE))

plotDimRed(mtecCombined, col_by = "Mki67")

mtecCombined@meta.data$Ccl21a <- mtecCombined@data["Ccl21a", ]
mtecCombined@meta.data$Aire <- mtecCombined@data["Aire", ]
mtecCombined@meta.data$Mki67 <- mtecCombined@data["Mki67", ]
mtecCombined@meta.data$Hmgb2 <- mtecCombined@data["Hmgb2", ]
mtecCombined@meta.data$Stmn1 <- mtecCombined@data["Stmn1", ]
mtecCombined@meta.data$Top2a <- mtecCombined@data["Top2a", ]

mtecCombined@meta.data$db_pos <- mtecCombined@meta.data$Ccl21a > 4 & 
  mtecCombined@meta.data$Aire > 0.5

mtecCombined@meta.data$Aire_exp <- mtecCombined@meta.data$Aire > 0 &
  mtecCombined@meta.data$Ccl21a <= 4
mtecCombined@meta.data$Ccl21a_exp <- mtecCombined@meta.data$Ccl21a > 4 &
  mtecCombined@meta.data$Aire <= 0
mtecCombined@meta.data$Aire_and_Ccl21a <- mtecCombined@meta.data$Aire > 0 &
  mtecCombined@meta.data$Ccl21a > 4
mtecCombined@meta.data$db_neg <- mtecCombined@meta.data$Ccl21a <= 4 &
  mtecCombined@meta.data$Aire <= 0

mtecCombined@meta.data$Mki67_pos <- mtecCombined@meta.data$Mki67 > 0.25

plotDimRed(mtecCombined, col_by = "db_pos", color = c("#000000", "#FF007F"))

lapply(timepoints, function(x) full_umap(mtecCombined,
  data_set = x, col_by = "db_pos",
  show_legend = TRUE, color = c("#000000", "#FF007F")))

percent_db_pos <- sapply(timepoints, USE.NAMES = TRUE,
  function(x) percent_ident(mtecCombined,
  data_set = x, meta_data_col = "exp", ident = "db_pos"))
percent_db_pos <- data.frame(percent_db_pos)
percent_db_pos$experiment <- rownames(percent_db_pos)

ggplot2::ggplot(percent_db_pos, ggplot2::aes(x = experiment,
                                              y = percent_db_pos)) +
  ggplot2::geom_bar(ggplot2::aes(fill = experiment),
    stat = "identity") + 
  ggplot2::scale_fill_manual(values = timecourse_color)

percent_ki67 <- sapply(timepoints, USE.NAMES = TRUE,
  function(x) percent_ident(mtecCombined,
  data_set = x, meta_data_col = "exp", ident = "Mki67_pos"))
percent_ki67 <- data.frame(percent_ki67)
percent_ki67$experiment <- rownames(percent_ki67)

ggplot2::ggplot(percent_ki67, ggplot2::aes(x = experiment,
                                              y = percent_ki67)) +
  ggplot2::geom_bar(ggplot2::aes(fill = experiment),
    stat = "identity") + 
  ggplot2::scale_fill_manual(values = timecourse_color)

stage_color_df <- data.frame("Cortico_medullary" = "#CC6600", "Ccl21a_high" = "#009933",
                            "Early_Aire" = "#0066CC", "Aire_positive" = "#660099",
                            "Late_Aire" = "#FF0000", "Tuft" = "#990000",
                            "unknown" = "#FDFBFB")

stage_color <- t(stage_color_df)[ , 1]

plot_trio_subset <- function(seurat_object, subset_on = NULL, ...) {
	if (!is.null(subset_on)) {
		seurat_object <- Seurat::SubsetData(seurat_object, ident.use = subset_on)
	}
	trio_plots(seurat_object, plot_name = subset_on, ...)
}

geneset <- c("Ccl21a", "Ccl21b.1", "Ccl21a.1")
geneset_2 <- c("Aire", "Ccl21a", "Mki67")


mtecCombined <- Seurat::SetAllIdent(mtecCombined, "exp")

trio_plots(mtecCombined, plot_name = "Ccl21 genes", geneset = geneset,
  cell_cycle = FALSE, plot_violin = TRUE, plot_jitter = FALSE, jitter_and_violin = FALSE,
  color = stage_color, sep_by = "stage")
lapply(timepoints, function(x) plot_trio_subset(mtecCombined, subset_on = x,
  geneset = geneset, cell_cycle = FALSE, plot_violin = TRUE, plot_jitter = FALSE,
  jitter_and_violin = FALSE, color = stage_color, sep_by = "stage"))

lapply(timepoints, function(x) plot_trio_subset(mtecCombined, subset_on = x,
  geneset = geneset_2, cell_cycle = FALSE, plot_violin = TRUE, plot_jitter = FALSE,
  jitter_and_violin = FALSE, color = stage_color, sep_by = "stage"))

trio_plots(mtecCombined, plot_name = "Cycling cells", geneset = c("Mki67"),
  cell_cycle = FALSE, plot_violin = TRUE, plot_jitter = FALSE, jitter_and_violin = FALSE,
  color = stage_color, sep_by = "stage")

plotting_data <- mtecCombined@meta.data

plotDimRed(mtecCombined, col_by = "cycle_phase", color = c("black", "red", "purple"))

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

mtecCombined_TAC <- Seurat::SubsetData(mtecCombined, ident.remove = "aireTrace")
mtecCombined_TAC <- Seurat::SetAllIdent(mtecCombined_TAC, "stage")
mtecCombined_TAC <- Seurat::SubsetData(mtecCombined_TAC, ident.use = "Early_Aire")
percent_Aire <- sapply(timepoints, USE.NAMES = TRUE,
  function(x) percent_ident(mtecCombined_TAC,
  data_set = x, meta_data_col = "exp", ident = "Aire_exp"))
percent_Ccl21 <- sapply(timepoints, USE.NAMES = TRUE,
  function(x) percent_ident(mtecCombined_TAC,
  data_set = x, meta_data_col = "exp", ident = "Ccl21a_exp"))
percent_dp_pos <- sapply(timepoints, USE.NAMES = TRUE,
  function(x) percent_ident(mtecCombined_TAC,
  data_set = x, meta_data_col = "exp", ident = "Aire_and_Ccl21a"))
percent_db_neg <- sapply(timepoints, USE.NAMES = TRUE,
  function(x) percent_ident(mtecCombined_TAC,
  data_set = x, meta_data_col = "exp", ident = "db_neg"))


percent_all <- data.frame(percent_Aire = percent_Aire,
                          percent_Ccl21a = percent_Ccl21,
                          percent_db_pos = percent_dp_pos,
                          percent_db_neg = percent_db_neg)
percent_all$experiment <- rownames(percent_all)\

# Start here
percent_all_m <- tidyr::gather(data = percent_all, key = classification, percent,
  value = percent_Aire, percent_Ccl21a, percent_db_pos, percent_db_neg,
  factor_key = TRUE)


ggplot2::ggplot(percent_Aire, ggplot2::aes(x = experiment,
                                              y = percent_Aire)) +
  ggplot2::geom_bar(ggplot2::aes(fill = experiment),
    stat = "identity") + 
  ggplot2::scale_fill_manual(values = timecourse_color)

plotting_data_TAC <- mtecCombined_TAC@meta.data

ggplot2::ggplot(plotting_data_TAC, ggplot2::aes(x = Aire, y = Ccl21a,
  color = exp)) +
  ggplot2::geom_point() +
  ggplot2::scale_color_manual(values = timecourse_color)
dev.off()
