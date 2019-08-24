library(mTEC.10x.pipeline)
source("/home/kwells4/mTEC_dev/mtec_snakemake/scripts/figure_funcs.R")

load("/home/kwells4/mTEC_dev/mtec_snakemake/allSamples/analysis_outs/seurat_allSamples_combined.rda")


timepoints <- c("isoControlBeg", "isoControlEnd", "timepoint1",
	"timepoint2", "timepoint3", "timepoint5")

pdf("/home/kwells4/mTEC_dev/mtec_snakemake/allSamples/analysis_outs/mark_umap.pdf")

lapply(timepoints, function(x) full_umap(mtecCombined,
  data_set = x, col_by = "Aire",
  show_legend = TRUE))

lapply(timepoints, function(x) full_umap(mtecCombined,
  data_set = x, col_by = "Fezf2",
  show_legend = TRUE))

dev.off()