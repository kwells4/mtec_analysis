library(mTEC.10x.pipeline)
source("/home/kwells4/mTEC_dev/mtec_snakemake/scripts/figure_funcs.R")

load("/home/kwells4/mTEC_dev/mtec_snakemake/allSamples/analysis_outs/seurat_allSamples_combined.rda")


timepoints <- c("isoControlBeg", "isoControlEnd", "timepoint1",
	"timepoint2", "timepoint3", "timepoint5")

pdf("/home/kwells4/mTEC_dev/mtec_snakemake/allSamples/analysis_outs/alice_umap.pdf")

lapply(timepoints, function(x) full_umap(mtecCombined,
  data_set = x, col_by = "5330417C22Rik",
  show_legend = TRUE))

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

geneset <- c("5330417C22Rik", "Aire", "Trpm5")

mtecCombined <- Seurat::SetAllIdent(mtecCombined, "exp")

lapply(timepoints, function(x) plot_trio_subset(mtecCombined, subset_on = x,
	geneset = geneset, cell_cycle = FALSE, plot_violin = TRUE, plot_jitter = FALSE,
	jitter_and_violin = FALSE, color = stage_color, sep_by = "stage"))

dev.off()
