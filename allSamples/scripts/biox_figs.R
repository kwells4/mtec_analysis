library(mTEC.10x.pipeline)
library(dplyr)
library(slingshot)
#library(svglite)
library(ggplot2)
library(gplots)
library(reshape)
source("/home/kwells4/mTEC_dev/mtec_snakemake/scripts/figure_funcs.R")

allSamples <- "/home/kwells4/mTEC_dev/mtec_snakemake/allSamples/analysis_outs/seurat_allSamples_combined.rda"
data_directory <- "/home/kwells4/mTEC_dev/data/"
save_dir <- "/home/kwells4/mTEC_dev/mtec_snakemake/allSamples/analysis_outs/"

stage_color_df <- data.frame("cTEC" = "#CC6600", "Ccl21a_high" = "#009933",
                            "TAC_TEC" = "#0066CC", "Aire_positive" = "#660099",
                            "Late_Aire" = "#FF0000", "Tuft" = "#990000",
                            "unknown" = "#FDFBFB")

stage_color <- t(stage_color_df)[ , 1]

mtecCombined <- get(load(allSamples))

load(paste0(data_directory, "gene_lists.rda"))

ggplot2::theme_set(ggplot2::theme_classic(base_size = 18))

stage_levels <- c("cTEC", "Ccl21a_high", "TAC_TEC",
                  "Aire_positive", "Late_Aire", "Tuft", "other")

stage_levels_wt <- c("cTEC", "Ccl21a_high", "TAC_TEC",
                  "Aire_positive", "Late_Aire", "Tuft", "unknown")

new_exp_names <- c(aireTrace = "Trace exp",
                  isoControlBeg = "Ctl wk 2",
                  isoControlEnd = "Ctl wk 10",
                  timepoint1 = "wk 2",
                  timepoint2 = "wk 4",
                  timepoint3 = "wk 6",
                  timepoint5 = "wk 10")

mtecCombined_with_un <- mtecCombined
mtecCombined@assay$DE <- NULL
# mtecCombined <- Seurat::SetAllIdent(mtecCombined, id = "stage")
# mtecCombined <- Seurat::SubsetData(mtecCombined, ident.remove = "unknown")

# Rename everything
mtecCombined <- Seurat::SetAllIdent(mtecCombined, id = "stage")
mtecCombined <- Seurat::SubsetData(mtecCombined, ident.remove = "unknown")
mtecCombined@meta.data$pub_exp <- new_exp_names[mtecCombined@meta.data$exp]
mtecCombined@meta.data$pub_exp <- factor(mtecCombined@meta.data$pub_exp,
                                         levels = unname(new_exp_names))
mtecCombined@meta.data$stage <- as.character(mtecCombined@meta.data$stage)
mtecCombined@meta.data$stage[mtecCombined@meta.data$stage ==
  "Cortico_medullary"] <- "cTEC"
mtecCombined@meta.data$stage[mtecCombined@meta.data$stage ==
  "Early_Aire"] <- "TAC_TEC"
mtecCombined@meta.data$stage <- factor(mtecCombined@meta.data$stage,
  levels = stage_levels)
mtecCombined@meta.data$stage_exp <- paste0(mtecCombined@meta.data$stage,
  "_", mtecCombined@meta.data$exp)

mtecCombined@meta.data$at_stage <- as.character(mtecCombined@meta.data$at_stage)
mtecCombined@meta.data$at_stage[mtecCombined@meta.data$at_stage ==
  "Cortico_medullary"] <- "cTEC"
mtecCombined@meta.data$at_stage[mtecCombined@meta.data$at_stage ==
  "Early_Aire"] <- "TAC_TEC"
mtecCombined@meta.data$at_stage <- factor(mtecCombined@meta.data$at_stage,
  levels = stage_levels)

timecourse_color <- c("#FC8D62", "#8DA0CB", "#E78AC3", "#65C42D",
  "#FFB62F", "#E5C494")
#timecourse_color <- RColorBrewer::brewer.pal(8, "Set2")
#timecourse_color <- c(timecourse_color[2:7])

data_sets <- unique(mtecCombined@meta.data$pub_exp)
data_sets <- data_sets[data_sets != new_exp_names['aireTrace']]

stage_color_df_3 <- data.frame("cTEC" = "#CC6600",
                               "Ccl21aHigh" = "#009933",
                               "EarlymTEC" = "#0066CC",
                               "AirePositive" = "#660099",
                               "LateAire" = "#FF0000",
                               "Tuft" = "#990000",
                               "unknown" = "#D3D3D3")

stage_color3 <- t(stage_color_df_3)[ , 1]



# Figure 4c Barplots of recovery
mtecCombSub <- mtecCombined
mtecCombSub@assay$ablation_DE <- NULL



# Figure 4f Timecourse in pseudotime
mtec_no_at <- mtecCombSub
mtec_no_at <- Seurat::SetAllIdent(mtec_no_at, id = "exp")
mtec_no_at <- Seurat::SubsetData(mtec_no_at, ident.remove = "aireTrace")

average_gene_list <- c("Aire", "Fezf2", "Gapdh", "Emc7", "Tnfrsf11a")

mtec_aire_positive <- Seurat::SubsetData(mtecCombSub, ident.use = "Aire_positive",
  subset.raw = TRUE)
cells_use <- rownames(mtec_aire_positive@meta.data)[mtec_aire_positive@meta.data$exp !=
                                                 "aireTrace"]

no_at_mtec_aire <- Seurat::SubsetData(mtec_aire_positive, cells.use = cells_use)


# Figure 6c violin plot of average expression of gene set across ablation

mtecCombined_all <- no_at_mtec_aire
# This is an okay place for a for loop (recursion) 
# http://adv-r.had.co.nz/Functionals.html
for (gene_set in names(gene_lists)) {
  mtecCombined_all <- plot_gene_set(mtecCombined_all,
    gene_lists[[gene_set]], gene_set, make_plot = FALSE)
}


trio_plots_median(mtecCombined_all, geneset = c("Aire", "aire_genes", "tra_fantom"),
  cell_cycle = FALSE, plot_violin = TRUE, jitter_and_violin = FALSE,
  plot_jitter = FALSE, color = timecourse_color, sep_by = "pub_exp",
  save_plot = paste0(save_dir, "/biox_violin_set1.pdf"))

trio_plots_median(mtecCombined_all, geneset = c("Fezf2", "fezf2_genes", "all_other_genes"),
  cell_cycle = FALSE, plot_violin = TRUE, jitter_and_violin = FALSE,
  plot_jitter = FALSE, color = timecourse_color, sep_by = "pub_exp",
  save_plot = paste0(save_dir, "/biox_violin_set2.pdf"))
