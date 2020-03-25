library(mTEC.10x.pipeline)
library(dplyr)
library(slingshot)
#library(svglite)
library(ggplot2)
library(gplots)
library(reshape)
source("/home/kwells4/mTEC_dev/mtec_snakemake/scripts/figure_funcs.R")



##################################################################################

# Get files from Snakemake
aireTrace <- "aireTrace/analysis_outs/seurat_aireTrace.rda"
controls <- "controls/analysis_outs/seurat_controls_merged.rda"
allSamples <- "allSamples/analysis_outs/seurat_allSamples_combined.rda"
controls_slingshot <- "controls/analysis_outs/controls_merged_slingshot.rda"
allSamples_slingshot <- "allSamples/analysis_outs/allSamples_combined_slingshot.rda"
early_aire_mtec <- "allSamples/analysis_outs/seurat_allSamples_Early_Aire_combined.rda"
save_file <- "figure_output/complete_figs.txt"
data_directory <- "/home/kwells4/mTEC_dev/data/"
save_dir <- "figure_output"


# Colors for plotting
stage_color_df <- data.frame("cTEC" = "#CC6600", "Ccl21a_high" = "#009933",
                            "TAC_TEC" = "#0066CC", "Aire_positive" = "#660099",
                            "Late_Aire" = "#FF0000", "Tuft" = "#990000",
                            "unknown" = "#FDFBFB")

stage_color <- t(stage_color_df)[ , 1]

# Load in data
mtec <- get(load(aireTrace))

mtec_wt <- readRDS(controls)

mtecCombined <- get(load(allSamples))

wt_slingshot <- get(load(controls_slingshot))

all_slingshot <- get(load(allSamples_slingshot))

progenitor_mtec <- get(load(early_aire_mtec))

fig_list <- list()

TFs <- get(load(paste0(data_directory, "TFs.rda")))

TFs_all <- c("H2afz", "Top2a", "Hmgb1", "Hmgn1", "H2afx", as.character(TFs))

bootstrap <- FALSE

stage_levels <- c("cTEC", "Ccl21a_high", "TAC_TEC",
                  "Aire_positive", "Late_Aire", "Tuft", "other")

stage_levels_wt <- c("cTEC", "Ccl21a_high", "TAC_TEC",
                  "Aire_positive", "Late_Aire", "Tuft", "unknown")

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 18))

mtec_wt@meta.data$stage <- as.character(mtec_wt@meta.data$stage)


mtec_wt@meta.data$stage[mtec_wt@meta.data$stage ==
  "Cortico_medullary"] <- "cTEC"
mtec_wt@meta.data$stage[mtec_wt@meta.data$stage ==
  "Early_Aire"] <- "TAC_TEC"

mtec_wt@meta.data$stage <- factor(mtec_wt@meta.data$stage, levels = stage_levels_wt)

############
# Figure 1 #
############
print("Figure 1")
# Figure 1a
# Graphical outline

# Figure 1b UMAP of controls
# Remove unknown cells
mtec_wt_plot <- mtec_wt
mtec_wt_plot@assay$DE <- NULL

mtec_wt_plot <- Seurat::SetAllIdent(mtec_wt_plot, id = "stage")
mtec_wt_plot <- Seurat::SubsetData(mtec_wt_plot, ident.remove = "unknown")

mtec_wt_plot@ident <- factor(mtec_wt_plot@ident, levels = stage_levels)

# Make a umap
tSNE_PCA(mtec_wt_plot, "stage", color = stage_color, show_legend = FALSE,
  save_plot = paste0(save_dir, "/figure_1bI.pdf"))

tSNE_PCA(mtec_wt_plot, "stage", color = stage_color, show_legend = TRUE,
  save_plot = paste0(save_dir, "/figure_1bII.pdf"))

# Figure 1c violin plots of markers
trio_plots_median(mtec_wt_plot, geneset = c("Ackr4", "Ccl21a", "Trpm5"),
  cell_cycle = FALSE, plot_violin = TRUE, jitter_and_violin = FALSE,
  plot_jitter = FALSE, color = stage_color, sep_by = "cluster",
  save_plot = paste0(save_dir, "/figure_1cI.pdf"))
trio_plots_median(mtec_wt_plot, geneset = c("Aire", "Trpm5", "Dclk1"),
  cell_cycle = FALSE, plot_violin = TRUE, jitter_and_violin = FALSE,
  plot_jitter = FALSE, color = stage_color, sep_by = "cluster",
  save_plot = paste0(save_dir, "/figure_1cII.pdf"))
trio_plots_median(mtec_wt_plot, geneset = c("Prss16", "Krt5", "Dclk1"),
  cell_cycle = FALSE, plot_violin = TRUE, jitter_and_violin = FALSE,
  plot_jitter = FALSE, color = stage_color, sep_by = "cluster",
  save_plot = paste0(save_dir, "/figure_1cIII.pdf"))
trio_plots_median(mtec_wt_plot, geneset = c("Fezf2", "Dclk1", "Ccl21a"),
  cell_cycle = FALSE, plot_violin = TRUE, jitter_and_violin = FALSE,
  plot_jitter = FALSE, color = stage_color, sep_by = "cluster",
  save_plot = paste0(save_dir, "/figure_1cIV.pdf"))

# Figure 1d RNA velocity

fig_list <- c(fig_list, "figure_1")
############
# Figure 2 #
############
print("Figure 2")
load(paste0(data_directory, "gene_lists.rda"))



# Figure 2b Dot plot of genes in TAC-TECs
comparison_list <- names(mtec_wt@assay$DE)
clusters <- lapply(comparison_list, get_slots)
clusters <- unique(unlist(clusters))
all_clusters <- sapply(clusters, cluster_gene_list,
                       cluster_list = comparison_list,
                       seurat_object = mtec_wt,
                       USE.NAMES = TRUE)

all_clusters_but_ea <- all_clusters
all_clusters_but_ea$Early_Aire <- NULL
ea_genes <- all_clusters$Early_Aire


genes_no_ea <- unique(unlist(all_clusters_but_ea))
genes_ea_overlap <- intersect(ea_genes, genes_no_ea)

genes_ea_unique <- setdiff(ea_genes, genes_ea_overlap)

print(genes_ea_unique)

TA_5 <- read.table("/home/kwells4/mTEC_dev/data/neural_transit_amplifying_cluster5.txt")
TA_8 <- read.table("/home/kwells4/mTEC_dev/data/neural_transit_amplifying_cluster8.txt")

names(TA_5) <- "gene"
names(TA_8) <- "gene"
TA_5_gene <- gsub("__.*", "", TA_5$gene)
TA_8_gene <- gsub("__.*", "", TA_8$gene)

TA_genes_all <- unique(c(TA_5_gene, TA_8_gene))

text_color <- ifelse(genes_ea_unique %in% TA_genes_all, "red", "black")

pdf(paste0(save_dir, "/figure_2bI.pdf"))
dot_plot <- Seurat::DotPlot(mtec_wt_plot, genes.plot = genes_ea_unique,
                            cols.use = c("blue", "red"), x.lab.rot = T,
                            plot.legend = F, dot.scale = 8, do.return = T)
dev.off()


pdf(paste0(save_dir, "/figure_2bII.pdf"), width = 12, height = 8)

dot_plot <- dot_plot + theme(axis.text.x = element_text(colour = text_color))

dot_plot

dev.off()

pdf(paste0(save_dir, "/figure_2bIII.pdf"))
dot_plot <- Seurat::DotPlot(mtec_wt_plot, genes.plot = genes_ea_unique,
                            cols.use = c("blue", "red"), x.lab.rot = T,
                            plot.legend = T, dot.scale = 8, do.return = T)
dev.off()

# Figure 2c Jitter plots of Cycling
trio_plots(mtec_wt_plot, geneset = c("Hmgb2", "Tubb5", "Stmn1"),
  cell_cycle = TRUE, jitter_and_violin = FALSE, plot_jitter = TRUE,
  color = stage_color, sep_by = "cluster",
  save_plot = paste0(save_dir, "/figure_2cI.pdf"), group_color = FALSE)
trio_plots(mtec_wt_plot, geneset = c("Aire", "Ccl21a", "Fezf2"),
  cell_cycle = TRUE, jitter_and_violin = FALSE, plot_jitter = TRUE,
  color = stage_color, sep_by = "cluster",
  save_plot = paste0(save_dir, "/figure_2cII.pdf"), group_color = FALSE)

# Figure 2d Flow plot of Ki67

# Figure 2e pseudotime analysis
# UMAP of cell stage with slingshot overlay
##########
# This is my crummy workaround until I get my new package running. If plotDimRed
# is returning a ggplot object, I'm golden
dims <- seq_len(2)
umap_coords <- data.frame(mtec_wt_plot@dr$umap@cell.embeddings)
umap_coords$stage <- mtec_wt_plot@meta.data$stage
base_plot <- ggplot2::ggplot(umap_coords,
  ggplot2::aes(x = UMAP1, y = UMAP2, color = stage)) +
  ggplot2::geom_point() + 
  ggplot2::scale_color_manual(values = stage_color) +
  ggplot2::theme(legend.position = "none")

##############
# Keep this
c <- slingCurves(wt_slingshot)[[1]]
curve1_coord <- data.frame(c$s[c$ord, dims])
curve1_coord$stage <- "line"

# This line cuts off the long tail... Probably a better option is
# to just remove unknown before running slingshot.
curve1_coord <- curve1_coord[curve1_coord$UMAP2 > -5, ]
base_plot <- base_plot + ggplot2::geom_path(data = curve1_coord,
  ggplot2::aes(UMAP1, UMAP2), color = "black", size = 1)

ggplot2::ggsave(paste0(save_dir, "/figure_2e.pdf"), plot = base_plot)


gene_names <- c("Aire", "Fezf2")
plot_sets <- c("tra_fantom", "aire_genes", "fezf2_genes")

# Figure 2f
# Pseudotime of genes
plot_names <- list(protein_coding = paste0(save_dir, "/figure_2fI.pdf"),
  tra_fantom = paste0(save_dir, "/figure_2fII.pdf"),
  aire_genes = paste0(save_dir, "/figure_2fIV.pdf"),
  Aire = paste0(save_dir, "/figure_2fIII.pdf"),
  fezf2_genes = paste0(save_dir, "/figure_2fVI.pdf"),
  Fezf2 = paste0(save_dir, "/figure_2fV.pdf"))


for (gene_set in names(gene_lists)) {
  mtec_wt_plot <- plot_gene_set(mtec_wt_plot,
    gene_lists[[gene_set]], gene_set, make_plot = FALSE)
}



# mtec_wt_plot@meta.data$stage <- factor(mtec_wt_plot@meta.data$stage,
#   levels = stage_levels)

# Plot each of the genes and gene sets in pseudotime, end at 16 because the
# is where the "unknown" cells are
plot_list <- lapply(names(plot_names), function(x) plot_sling_pseudotime(
  seurat_object = mtec_wt_plot, sling_object = wt_slingshot, y_val = x,
  col_by = "stage", pseudotime_curve = "curve1", color = stage_color,
  range = c(0, 16), save_plot = plot_names[[x]]))

fig_list <- c(fig_list, "figure_2")


############
# Figure 3 #
############
print("Figure 3")

mtec@meta.data$stage <- as.character(mtec@meta.data$stage)

mtec@meta.data$stage[mtec@meta.data$stage ==
  "Cortico_medullary"] <- "cTEC"
mtec@meta.data$stage[mtec@meta.data$stage ==
  "Early_Aire"] <- "TAC_TEC"
# Make plots from aire trace
mtec_no_un <- mtec
mtec_no_un@assay$DE <- NULL
########################################################################
# Fix this section by reruning aireTrace analysis driver
mtec_no_un@meta.data$stage[is.na(mtec_no_un@meta.data$stage)] <- "other"



mtec_no_un@meta.data$stage <- factor(mtec_no_un@meta.data$stage,
  levels = stage_levels)
mtec_no_un <- Seurat::SetAllIdent(mtec_no_un, id = "res.0.6")
mtec_no_un <- Seurat::SetAllIdent(mtec_no_un, id = "stage")
mtec_no_un@ident <- factor(mtec_no_un@ident, levels = stage_levels)
#########################################################################
mtec_no_un <- Seurat::SubsetData(mtec_no_un, ident.remove = "other")



# Figure 3b UMAPs of Aire trace
tSNE_PCA(mtec_no_un, "stage", color = stage_color, show_legend = FALSE,
  save_plot = paste0(save_dir, "/figure_3b.pdf"))
tSNE_PCA(mtec_no_un, "stage", color = stage_color,
  save_plot = paste0(save_dir, "/figure_3bI.pdf"))

# Figure 3c Violin plots of marker genes
trio_plots_median(mtec_no_un, geneset = c("Aire", "Ackr4", "Trpm5"),
  cell_cycle = FALSE, plot_violin = TRUE, jitter_and_violin = FALSE,
  plot_jitter = FALSE, color = stage_color, sep_by = "cluster",
  save_plot = paste0(save_dir, "/figure_3cI.pdf"))
trio_plots_median(mtec_no_un, geneset = c("Ccl21a", "Trpm5", "Dclk1"),
  cell_cycle = FALSE, plot_violin = TRUE, jitter_and_violin = FALSE,
  plot_jitter = FALSE, color = stage_color, sep_by = "cluster",
  save_plot = paste0(save_dir, "/figure_3cII.pdf"))
trio_plots_median(mtec_no_un, geneset = c("GFP", "Prss16","Dclk1"),
  cell_cycle = FALSE, plot_violin = TRUE, jitter_and_violin = FALSE,
  plot_jitter = FALSE, color = stage_color, sep_by = "cluster",
  save_plot = paste0(save_dir, "/figure_3cIII.pdf"))
trio_plots_median(mtec_no_un, geneset = c("Krt5", "Dclk1", "Ccl21a"),
  cell_cycle = FALSE, plot_violin = TRUE, jitter_and_violin = FALSE,
  plot_jitter = FALSE, color = stage_color, sep_by = "cluster",
  save_plot = paste0(save_dir, "/figure_3cIV.pdf"))

# Figure 3d UMAPs of Ccl21a, Aire, and GFP
tSNE_PCA(mtec_no_un, "Ccl21a", save_plot = paste0(save_dir, "/figure_3dI.pdf"))
tSNE_PCA(mtec_no_un, "Aire", save_plot = paste0(save_dir, "/figure_3dII.pdf"))
tSNE_PCA(mtec_no_un, "GFP", save_plot = paste0(save_dir, "/figure_3dIII.pdf"))

fig_list <- c(fig_list, "figure_3")
############
# Figure 4 #
############

print("Figure 4")

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

# Figure 4b umaps of Aire trace and all samples
tSNE_PCA(mtecCombined, "stage", save_plot = paste0(save_dir, "/figure_4bI.pdf"),
  color = stage_color, show_legend = FALSE)

# UMAP highlighting aire_trace cells colored by aire_trace labels
full_umap(mtecCombined, "aireTrace", col_by = "at_stage", color = stage_color,
  save_plot = paste0(save_dir, "/figure_4bII.pdf"), show_legend = FALSE)

tSNE_PCA(mtecCombined, "stage", save_plot = paste0(save_dir, "/figure_4bIII.pdf"),
  color = stage_color, show_legend = TRUE)

# Figure 4c Barplots of recovery
mtecCombSub <- mtecCombined
mtecCombSub@assay$ablation_DE <- NULL

stage_list_all <- lapply(data_sets, function(x) populations_dfs_new(mtecCombSub,
                         x, subsample = TRUE, subsample_by = "pub_exp"))
stage_df_all <- do.call("rbind", stage_list_all)

stage_df_all$sample <- factor(stage_df_all$sample, levels = unname(new_exp_names))

population_plots(stage_df_all, color = stage_color,
  save_plot = paste0(save_dir, "/figure_4c.pdf"))

# Figure 4d Flow of week 4

# Figure 4e Umap with slingshot overlay
##########
# This is my crummy workaround until I get my new package running. If plotDimRed
# is returning a ggplot object, I'm golden
dims <- seq_len(2)
umap_coords <- data.frame(mtecCombined@dr$umap@cell.embeddings)
umap_coords$stage <- mtecCombined@meta.data$stage
base_plot <- ggplot2::ggplot(umap_coords,
  ggplot2::aes(x = UMAP1, y = UMAP2, color = stage)) +
  ggplot2::geom_point() + 
  ggplot2::scale_color_manual(values = stage_color) +
  ggplot2::theme(legend.position = "none")

##############
# Keep this
c <- slingCurves(all_slingshot)[[3]]
curve1_coord <- data.frame(c$s[c$ord, dims])
curve1_coord$stage <- "line"

# This line cuts off the long tail... Probably a better option is
# to just remove unknown before running slingshot.
curve1_coord <- curve1_coord[curve1_coord$UMAP2 > -5, ]
base_plot <- base_plot + ggplot2::geom_path(data = curve1_coord,
  ggplot2::aes(UMAP1, UMAP2), color = "black", size = 1)

ggplot2::ggsave(paste0(save_dir, "/figure_3d.pdf"), plot = base_plot)


# Figure 4f Timecourse in pseudotime
mtec_no_at <- mtecCombSub
mtec_no_at <- Seurat::SetAllIdent(mtec_no_at, id = "exp")
mtec_no_at <- Seurat::SubsetData(mtec_no_at, ident.remove = "aireTrace")
plot_sling_pseudotime(
  seurat_object = mtec_no_at, sling_object = all_slingshot, y_val = "pub_exp",
  col_by = "pub_exp", pseudotime_curve = "curve3", color = timecourse_color,
  plot_type = "density", width = 7, height = 7,
  save_plot = paste0(save_dir, "/figure_3e.pdf"))

fig_list <- c(fig_list, "figure_4")

# ############
# # Figure 5 #
# ############
print("Figure 5")
reanalysis_colors <- c("#603E95", "#009DA1", "#FAC22B", "#D7255D")

projenitor_mtec <- mtec_no_at

projenitor_mtec <- Seurat::SetAllIdent(projenitor_mtec, id = "stage")
projenitor_mtec <- Seurat::SubsetData(projenitor_mtec, ident.use = "TAC_TEC")

progenitor_mtec@meta.data$pub_exp <- new_exp_names[progenitor_mtec@meta.data$exp]
progenitor_mtec@meta.data$pub_exp <- factor(progenitor_mtec@meta.data$pub_exp,
                                         levels = unname(new_exp_names))
cells_use <- rownames(progenitor_mtec@meta.data)[progenitor_mtec@meta.data$exp !=
                                                 "aireTrace"]

no_at_mtec <- Seurat::SubsetData(progenitor_mtec, cells.use = cells_use)

# Figure 5a Highlight just the early aire cells 
highlight_one_group(mtecCombined, meta_data_col = "stage", group = "TAC_TEC",
  color_df = stage_color, show_legend = FALSE,
  save_plot = paste0(save_dir, "/figure_5a.pdf"))

highlight_one_group(mtecCombined, meta_data_col = "stage", group = "TAC_TEC",
  color_df = stage_color, show_legend = TRUE,
  save_plot = paste0(save_dir, "/figure_5aI.pdf"))

# Figure 5b % of cells in G2
cell_cycle <- mtecCombined@meta.data[rownames(mtecCombined@meta.data) %in% 
  rownames(no_at_mtec@meta.data), ]
if (!identical(rownames(cell_cycle), rownames(no_at_mtec@meta.data))) {
  print("must reorder cells")
  cell_cycle <- cell_cycle[match(rownames(no_at_mtec@meta.data),
                                     rownames(cell_cycle)), , drop = FALSE]
}
no_at_mtec@meta.data$cycle_phase <- cell_cycle$cycle_phase
percent_cycling <- sapply(data_sets, USE.NAMES = TRUE,
  function(x) percent_cycling_cells(projenitor_mtec,
  data_set = x, meta_data_col = "pub_exp"))
percent_cycling <- data.frame(percent_cycling)
percent_cycling$experiment <- rownames(percent_cycling)
percent_cycling$percent_cycling <- percent_cycling$percent_cycling * 100
percent_cycling$experiment <- factor(percent_cycling$experiment,
                                     levels = unname(new_exp_names))
pdf(paste0(save_dir, "/figure_5b.pdf"))

ggplot2::ggplot(percent_cycling, ggplot2::aes(x = experiment,
                                              y = percent_cycling)) +
  ggplot2::geom_bar(ggplot2::aes(fill = experiment),
    stat = "identity") + 
  ggplot2::scale_fill_manual(values = timecourse_color)

dev.off()

stat_comparisons <- list(c("Ctl wk 10", "Ctl wk 2"),
                         c("Ctl wk 10", "wk 2"),
                         c("Ctl wk 10", "wk 4"),
                         c("Ctl wk 10", "wk 6"),
                         c("Ctl wk 10", "wk 10"))

# Figure 5c Violin plots of genes of interest
trio_plots_median(projenitor_mtec, geneset = c("Aire", "Ccl21a", "Fezf2"),
  cell_cycle = FALSE, plot_violin = TRUE, jitter_and_violin = FALSE,
  plot_jitter = FALSE, color = timecourse_color, sep_by = "pub_exp",
  save_plot = paste0(save_dir, "/figure_5cII.pdf"))


# Figure 5d UMAP of reanalysis of early aire cells
tSNE_PCA(no_at_mtec, "cluster", color = reanalysis_colors,
  save_plot = paste0(save_dir, "/figure_5dI.pdf"))

tSNE_PCA(no_at_mtec, "cluster", color = reanalysis_colors, show_legend = FALSE,
  save_plot = paste0(save_dir, "/figure_5d.pdf"))

# Figure 5e Violin plots of genes of interest
trio_plots_median(no_at_mtec, geneset = c("Aire", "Ccl21a", "Fezf2"), cell_cycle = FALSE,
  plot_violin = TRUE, jitter_and_violin = FALSE, plot_jitter = FALSE,
  sep_by = "cluster", color = reanalysis_colors,
  save_plot = paste0(save_dir, "/figure_5eII.pdf"))

# Figure 5f Percent of cycling cells in each new cluster
clusters <- unique(no_at_mtec@meta.data$res.0.6)
percent_cycling <- sapply(clusters, USE.NAMES = TRUE,
  function(x) percent_cycling_cells(no_at_mtec,
  data_set = x, meta_data_col = "res.0.6"))
percent_cycling <- data.frame(percent_cycling)
percent_cycling$cluster <- sub("\\d\\.", "", rownames(percent_cycling))
pdf(paste0(save_dir, "/figure_5f.pdf"))

ggplot2::ggplot(percent_cycling, ggplot2::aes(x = cluster,
                                              y = percent_cycling)) +
  ggplot2::geom_bar(ggplot2::aes(fill = cluster),
    stat = "identity") + 
  ggplot2::scale_fill_manual(values = reanalysis_colors)

dev.off()

# Figure 5g Population plot of new clusters over time.
stage_list_all <- lapply(data_sets, function(x) populations_dfs_new(no_at_mtec,
                         x, subsample = TRUE, subsample_by = "pub_exp",
                         meta_data_col = "res.0.6"))
stage_df_all <- do.call("rbind", stage_list_all)
stage_df_all$sample <- factor(stage_df_all$sample, levels = unname(new_exp_names))


population_plots(stage_df_all, color = reanalysis_colors,
  save_plot = paste0(save_dir, "/figure_5g.pdf"))



fig_list <- c(fig_list, "figure_5")

############
# Figure 6 #
############
print("Figure 6")

lowest_UMI_exp <- "timepoint2"

downsample_UMI <- TRUE

if (downsample_UMI) {
  
  lowest_UMI <- get_umi(mtecCombSub, subset_seurat = TRUE, subset_by = "stage_exp",
    subset_val = paste0("Aire_positive_", lowest_UMI_exp))
} else {
  lowest_UMI <- NULL
}

# Figure 6a Start with full UMAP only coloring the Aire positive cluster
highlight_one_group(mtecCombined, meta_data_col = "stage", group = "Aire_positive",
  color_df = stage_color, show_legend = FALSE,
  save_plot = paste0(save_dir, "/figure_6a.pdf"))

highlight_one_group(mtecCombined, meta_data_col = "stage", group = "Aire_positive",
  color_df = stage_color, show_legend = TRUE,
  save_plot = paste0(save_dir, "/figure_6aI.pdf"))

# Figure 6b Violin plots of expression across ablation
average_gene_list <- c("Aire", "Fezf2", "Gapdh", "Emc7", "Tnfrsf11a")

mtec_aire_positive <- Seurat::SubsetData(mtecCombSub, ident.use = "Aire_positive",
  subset.raw = TRUE)
cells_use <- rownames(mtec_aire_positive@meta.data)[mtec_aire_positive@meta.data$exp !=
                                                 "aireTrace"]

no_at_mtec_aire <- Seurat::SubsetData(mtec_aire_positive, cells.use = cells_use)

trio_plots_median(no_at_mtec_aire, geneset = c("Aire", "Fezf2", "Gapdh"),
  cell_cycle = FALSE, jitter_and_violin = FALSE, plot_jitter = FALSE,
  plot_violin = TRUE, color = timecourse_color, sep_by = "pub_exp",
  save_plot = paste0(save_dir, "/figure_6bI.pdf"))

# Figure 6c violin plot of average expression of gene set across ablation

mtecCombined_all <- no_at_mtec_aire
# This is an okay place for a for loop (recursion) 
# http://adv-r.had.co.nz/Functionals.html
for (gene_set in names(gene_lists)) {
  mtecCombined_all <- plot_gene_set(mtecCombined_all,
    gene_lists[[gene_set]], gene_set, make_plot = FALSE)
}


gene_sets <- c("all_other_genes", "tra_fantom", "aire_genes")
trio_plots_median(mtecCombined_all, geneset = gene_sets,
  cell_cycle = FALSE, plot_violin = TRUE, jitter_and_violin = FALSE,
  plot_jitter = FALSE, color = timecourse_color, sep_by = "pub_exp",
  save_plot = paste0(save_dir, "/figure_6c.pdf"))

# Figure 6d percent of genes across ablation 
old_data_sets <- unique(mtecCombined@meta.data$exp)
old_data_sets <- old_data_sets[old_data_sets != "aireTrace"]

percents_counts_all <- sapply(old_data_sets, function(x) percents_and_counts(mtecCombSub,
  gene_lists = gene_lists, batch_name = paste0("Aire_positive_", x),
  downsample_UMI = downsample_UMI, one_batch = TRUE, batch = "stage_exp",
  lowest_UMI = lowest_UMI, count = c("genes", "percent", "umi")))

percents <- sapply(names(percents_counts_all), function(x) 
  get_perc_count(percents_counts_all, x, data_type = "percent"), USE.NAMES = TRUE)

counts <- lapply(names(percents_counts_all), function(x)
  get_perc_count(percents_counts_all, x, data_type = "counts"))

counts_df <- do.call(rbind, counts)
counts_df$exp <- sub("Aire_positive_", "", counts_df$exp)
counts_df$pub_exp <- new_exp_names[counts_df$exp]
counts_df$pub_exp <- factor(counts_df$pub_exp,
                            levels = unname(new_exp_names))
counts_df_m <- reshape2::melt(counts_df, variable.name = "gene_list",
  value.name = "gene_count")

# I Changed this on 082219 to add co-expression data
to_plot <- c("tra_fantom", "all_other_genes", "aire_genes", "fezf2_genes",
  "co_expr_A", "co_expr_C", "co_expr_D", "co_expr_E")
short_list <- c("tra_fantom", "aire_genes", "fezf2_genes")

counts_df_plot <- counts_df_m[counts_df_m$gene_list %in% to_plot, ]
counts_df_short <- counts_df_m[counts_df_m$gene_list %in% short_list, ]

colnames(percents) <- new_exp_names[colnames(percents)]
percents_m <- reshape2::melt(percents)
names(percents_m) <- c("gene_list", "pub_exp", "percent_of_genes")

# I suspect this line is the problem...

percents_m$pub_exp <- factor(percents_m$pub_exp,
                             levels = unname(new_exp_names))
percents_plot <- percents_m[percents_m$gene_list %in% to_plot, ]

pdf(paste0(save_dir, "/figure_6d.pdf"))
ggplot2::ggplot(percents_plot, ggplot2::aes(x = pub_exp, y = percent_of_genes,
                                         group = gene_list, color = gene_list)) +
  ggplot2::geom_line(size = 2) +
  ggplot2::ylim(0,1) +
  ggplot2::scale_color_brewer(palette = "Dark2")
  #ggplot2::theme_classic()

dev.off()


fig_list <- c(fig_list, "figure_6")

############
# Figure 7 #
############
print("Figure 7")

fezf2_plot_names <- c(isoControlBeg = paste0(save_dir, "/figure_7aI.pdf"),
                     isoControlEnd = paste0(save_dir, "/figure_7aII.pdf"),
                     timepoint1 = paste0(save_dir, "/figure_7aIII.pdf"),
                     timepoint2 = paste0(save_dir, "/figure_7aIV.pdf"),
                     timepoint3 = paste0(save_dir, "/figure_7aV.pdf"),
                     timepoint5 = paste0(save_dir, "/figure_7aVI.pdf"))

lapply(names(fezf2_plot_names), function(x) full_umap(mtecCombined,
  data_set = x, col_by = "Fezf2", show_legend = TRUE,
  save_plot = fezf2_plot_names[[x]]))

heatmap_aire_mtec <- Seurat::SetAllIdent(no_at_mtec_aire, id = "exp")
mtec.markers <- Seurat::FindAllMarkers(object = mtec_wt_plot, only.pos = TRUE,
  min.pct = 0.25, thresh.use = 0.25)
top30 <- mtec.markers %>% group_by(cluster) %>% top_n(30, avg_logFC)
top30 <- as.data.frame(top30)

aire_positive_genes <- top30[top30$cluster == "Aire_positive", ]$gene

aire_positive_genes <- c(aire_positive_genes, "Fezf2")


pdf(paste0(save_dir, "/figure_7b.pdf"))

print(plot_heatmap_new(heatmap_aire_mtec, cell_color = timecourse_color,
  subset_list = aire_positive_genes))

dev.off()

fig_list <- c(fig_list, "figure_7")

############################################################################

########################
# Supplemental Figures #
########################

#########################
# Supplemental Figure 1 #
#########################
print("supplemental_figure_1")

# Supplemental Figure 1a sorting strategy

# Supplemental Figure 1b umap colored by original cluster
tSNE_PCA(mtec_wt, "res.0.6",show_legend = FALSE,
  save_plot = paste0(save_dir, "/figure_s1bI.pdf"))

tSNE_PCA(mtec_wt, "res.0.6", show_legend = TRUE,
  save_plot = paste0(save_dir, "/figure_s1bII.pdf"))

# Supplemental Figure 1c umaps by gene
tSNE_PCA(mtec_wt, "Ackr4", save_plot = paste0(save_dir, "/figure_s1cI.pdf"))
tSNE_PCA(mtec_wt, "Prss16", save_plot = paste0(save_dir, "/figure_s1cII.pdf"))
tSNE_PCA(mtec_wt, "Ccl21a", save_plot = paste0(save_dir, "/figure_s1cIII.pdf"))
tSNE_PCA(mtec_wt, "Krt5", save_plot = paste0(save_dir, "/figure_s1cIV.pdf"))
tSNE_PCA(mtec_wt, "Hmgb2", save_plot = paste0(save_dir, "/figure_s1cV.pdf"))
tSNE_PCA(mtec_wt, "Stmn1", save_plot = paste0(save_dir, "/figure_s1cVI.pdf"))
tSNE_PCA(mtec_wt, "Fezf2", save_plot = paste0(save_dir, "/figure_s1cVII.pdf"))
tSNE_PCA(mtec_wt, "Aire", save_plot = paste0(save_dir, "/figure_s1cVIII.pdf"))
tSNE_PCA(mtec_wt, "Tnfrsf11a", save_plot = paste0(save_dir, "/figure_s1cIX.pdf"))
tSNE_PCA(mtec_wt, "Krt10", save_plot = paste0(save_dir, "/figure_s1cX.pdf"))
tSNE_PCA(mtec_wt, "Cldn4", save_plot = paste0(save_dir, "/figure_s1cXI.pdf"))
tSNE_PCA(mtec_wt, "Trpm5", save_plot = paste0(save_dir, "/figure_s1cXII.pdf"))
tSNE_PCA(mtec_wt, "Dclk1", save_plot = paste0(save_dir, "/figure_s1cXIII.pdf"))

# Supplemental Figure 1d heatmap of DE genes
# Heatmap of all TFs with interesting TFs highlighted
pdf(paste0(save_dir, "/figure_s1d.pdf"))

heatmap_mtec <- mtec_wt_plot

heatmap_mtec@assay$DE <- mtec_wt@assay$DE

plot_heatmap(heatmap_mtec, subset_list = TFs_all,
  color_list = c("Cdx1", "Utf1", "Tcf7", "Spib", "Cdk4", "Ptma",
    "H2afx", "Hmgb1"),
  color_list2 = c("Aire", "Irf7", "Cited2", "Spib", "Hes1", "Pax1",
                "Relb", "Lmo4", "Pou2f3"),
  cell_color = stage_color)
dev.off()

fig_list <- c(fig_list, "figure_s1")

#########################
# Supplemental Figure 2 #
#########################
print("supplemental_figure_2")


# Supplemental Figure 2a Heatmap of all early Aire genes
comparison_list <- names(mtec_wt@assay$DE)
clusters <- lapply(comparison_list, get_slots)
clusters <- unique(unlist(clusters))
all_clusters <- sapply(clusters, cluster_gene_list,
                       cluster_list = comparison_list,
                       seurat_object = mtec_wt,
                       USE.NAMES = TRUE)

all_clusters_but_ea <- all_clusters
all_clusters_but_ea$Early_Aire <- NULL
ea_genes <- all_clusters$Early_Aire

pdf(paste0(save_dir, "/Figure_s2a.pdf"), width = 8, height = 15)

plot_heatmap(heatmap_mtec, subset_list = ea_genes,
  cell_color = stage_color)

dev.off()

# Supplemental Figure 2b Heatmap of unique Early Aire genes labeled by TA overlap
genes_no_ea <- unique(unlist(all_clusters_but_ea))
genes_ea_overlap <- intersect(ea_genes, genes_no_ea)

genes_ea_unique <- setdiff(ea_genes, genes_ea_overlap)

TA_5 <- read.table("/home/kwells4/mTEC_dev/data/neural_transit_amplifying_cluster5.txt")
TA_8 <- read.table("/home/kwells4/mTEC_dev/data/neural_transit_amplifying_cluster8.txt")

names(TA_5) <- "gene"
names(TA_8) <- "gene"
TA_5_gene <- gsub("__.*", "", TA_5$gene)
TA_8_gene <- gsub("__.*", "", TA_8$gene)

TA_genes_all <- unique(c(TA_5_gene, TA_8_gene))


pdf(paste0(save_dir, "/figure_s2b.pdf"))

plot_heatmap(heatmap_mtec, subset_list = genes_ea_unique,
  color_list = TA_genes_all,
  cell_color = stage_color)


dev.off()


# Supplemental Figure 2c violin plots of stem cell genes
trio_plots_median(mtec_wt_plot, geneset = c("Pou5f1", "Sox2", "Krt15"),
  cell_cycle = FALSE, plot_violin = TRUE, plot_jitter = FALSE,
  color = stage_color, sep_by = "stage",
  save_plot = paste0(save_dir, "/figure_s2c.pdf"))


# Supplemental Figure 2d bar plots of cycling counts
mtec_wt_plot@meta.data$Ccl21a <- mtec_wt_plot@data["Ccl21a", ]
mtec_wt_plot@meta.data$Aire <- mtec_wt_plot@data["Aire", ]
mtec_wt_plot@meta.data$db_pos <- mtec_wt_plot@meta.data$Ccl21a > 4 & 
  mtec_wt_plot@meta.data$Aire > 0.5

mtec_wt_plot@meta.data$Aire_exp <- mtec_wt_plot@meta.data$Aire > 0.5 &
  mtec_wt_plot@meta.data$Ccl21a <= 4
mtec_wt_plot@meta.data$Ccl21a_exp <- mtec_wt_plot@meta.data$Ccl21a > 4 &
  mtec_wt_plot@meta.data$Aire <= 0.5
mtec_wt_plot@meta.data$Aire_and_Ccl21a <- mtec_wt_plot@meta.data$Aire > 0.5 &
  mtec_wt_plot@meta.data$Ccl21a > 4
mtec_wt_plot@meta.data$db_neg <- mtec_wt_plot@meta.data$Ccl21a <= 4 &
  mtec_wt_plot@meta.data$Aire <= 0.5

mtec_tac <- mtec_wt_plot

mtec_tac <- Seurat::SetAllIdent(mtec_tac, id = "stage")
mtec_tac <- Seurat::SubsetData(mtec_tac, ident.use = "TAC_TEC")

cell_cycle_phase <- c("G1", "G2M", "S")

count_Aire <- sapply(cell_cycle_phase, USE.NAMES = TRUE,
  function(x) percent_ident(mtec_tac,
  data_set = x, meta_data_col = "cycle_phase", ident = "Aire_exp", count = TRUE))
count_Ccl21 <- sapply(cell_cycle_phase, USE.NAMES = TRUE,
  function(x) percent_ident(mtec_tac,
  data_set = x, meta_data_col = "cycle_phase", ident = "Ccl21a_exp", count = TRUE))
count_dp_pos <- sapply(cell_cycle_phase, USE.NAMES = TRUE,
  function(x) percent_ident(mtec_tac,
  data_set = x, meta_data_col = "cycle_phase", ident = "Aire_and_Ccl21a", count = TRUE))
count_db_neg <- sapply(cell_cycle_phase, USE.NAMES = TRUE,
  function(x) percent_ident(mtec_tac,
  data_set = x, meta_data_col = "cycle_phase", ident = "db_neg", count = TRUE))

blue_palette <- c(count_Aire = "#111E6C",
                  count_Ccl21a = "#008081",
                  count_db_pos = "#0080FF",
                  count_db_neg = "#4C516D")

count_all <- data.frame(count_Aire = count_Aire,
                        count_Ccl21a = count_Ccl21,
                        count_db_pos = count_dp_pos,
                        count_db_neg = count_db_neg)

count_all$cycle_phase <- rownames(count_all)


count_all_m <- tidyr::gather(data = count_all, key = classification,
  value = percent, count_Aire:count_db_neg,
  factor_key = TRUE)
pdf(paste0(save_dir, "/figure_s2d.pdf"))
ggplot2::ggplot(count_all_m, ggplot2::aes(x = cycle_phase, y = percent,
  fill = classification)) +
  ggplot2::geom_bar(position = "dodge", stat = "identity") +
  ggplot2::scale_fill_manual(values = blue_palette)
dev.off()


fig_list <- c(fig_list, "figure_s2")

######################### 
# Supplemental Figure 3 #
#########################
print("supplemental_figure_3")

mtec_wt <- Seurat::SetAllIdent(mtec_wt, id = "stage")
# Supplemental Figure 1b umap colored by original cluster
tSNE_PCA(mtec, "res.0.6",show_legend = FALSE,
  save_plot = paste0(save_dir, "/figure_s3aIII.pdf"))

tSNE_PCA(mtec, "res.0.6", show_legend = TRUE,
  save_plot = paste0(save_dir, "/figure_s3aIV.pdf"))


# Supplental Figure 3a Dotplot of Aire trace
stage_levels <- c("cTEC", "Ccl21a_high", "TAC_TEC",
                  "Aire_positive", "Late_Aire", "Tuft")

mtec_no_un@ident <- factor(mtec_no_un@ident,
   levels = stage_levels)
pdf(paste0(save_dir, "/figure_s3aI.pdf"))
dot_plot <- Seurat::DotPlot(mtec_no_un, genes.plot = genes_ea_unique,
                            cols.use = c("blue", "red"), x.lab.rot = T,
                            plot.legend = F, dot.scale = 8, do.return = T)


dev.off()

pdf(paste0(save_dir, "/figure_s3aII.pdf"))
dot_plot <- Seurat::DotPlot(mtec_no_un, genes.plot = genes_ea_unique,
                            cols.use = c("blue", "red"), x.lab.rot = T,
                            plot.legend = T, dot.scale = 8, do.return = T)
dev.off()

pdf(paste0(save_dir, "/figure_s3aIII.pdf"), width = 12, height = 8)

dot_plot <- dot_plot + theme(axis.text.x = element_text(colour = text_color))

dot_plot

dev.off()

# Supplemental Figure 3b Jitter plots of cycling with markers from AT
trio_plots(mtec_no_un, geneset = c("Hmgb2", "Tubb5", "Stmn1"),
  cell_cycle = TRUE, jitter_and_violin = FALSE, plot_jitter = TRUE,
  color = stage_color, sep_by = "stage",
  save_plot = paste0(save_dir, "/figure_s3bI.pdf"), group_color = FALSE)
trio_plots(mtec_no_un, geneset = c("Aire", "Ccl21a", "Fezf2"),
  cell_cycle = TRUE, jitter_and_violin = FALSE, plot_jitter = TRUE,
  color = stage_color, sep_by = "stage",
  save_plot = paste0(save_dir, "/figure_s3bII.pdf"), group_color = FALSE)

# Supplemental Figure 3c Correlation of WT with Aire Trace
pdf(paste0(save_dir, "/figure_s3c.pdf"))
master_plot(mtec, "aire_trace", mtec_wt, "wt", stage_color_df)
dev.off()

fig_list <- c(fig_list, "supplemental_figure_3")
#########################
# Supplemental Figure 4 #
#########################
print("supplemental_figure_4")

# Supplemental Figure 4a UMAP of original cluster
tSNE_PCA(mtecCombined_with_un, "res.0.6", save_plot = paste0(save_dir, "/figure_s4a.pdf"))

# Supplemental Figure 4b Umap of stage recovery
plot_names <- list(isoControlBeg = paste0(save_dir, "/figure_s4bI.pdf"),
  isoControlEnd = paste0(save_dir, "/figure_s4bII.pdf"),
  timepoint1 = paste0(save_dir, "/figure_s4bIII.pdf"),
  timepoint2 = paste0(save_dir, "/figure_s4bIV.pdf"),
  timepoint3 = paste0(save_dir, "/figure_s4bV.pdf"),
  timepoint5 = paste0(save_dir, "/figure_s4bVI.pdf"))

lapply(names(plot_names), function(x) full_umap(mtecCombined,
  data_set = x, col_by = "stage", color = stage_color,
  save_plot = plot_names[[x]], show_legend = FALSE))

# Supplemental Figure 4c UMAP of marker genes
tSNE_PCA(mtecCombined_with_un, "Ackr4", save_plot = paste0(save_dir, "/figure_s4cI.pdf"))
tSNE_PCA(mtecCombined_with_un, "Prss16", save_plot = paste0(save_dir, "/figure_s4cII.pdf"))
tSNE_PCA(mtecCombined_with_un, "Ccl21a", save_plot = paste0(save_dir, "/figure_s4cIII.pdf"))
tSNE_PCA(mtecCombined_with_un, "Krt5", save_plot = paste0(save_dir, "/figure_s4cIV.pdf"))
tSNE_PCA(mtecCombined_with_un, "Hmgb2", save_plot = paste0(save_dir, "/figure_s4cV.pdf"))
tSNE_PCA(mtecCombined_with_un, "Stmn1", save_plot = paste0(save_dir, "/figure_s4cVI.pdf"))
tSNE_PCA(mtecCombined_with_un, "Fezf2", save_plot = paste0(save_dir, "/figure_s4cVII.pdf"))
tSNE_PCA(mtecCombined_with_un, "Aire", save_plot = paste0(save_dir, "/figure_s4cVIII.pdf"))
tSNE_PCA(mtecCombined_with_un, "Tnfrsf11a", save_plot = paste0(save_dir, "/figure_s4cIX.pdf"))
tSNE_PCA(mtecCombined_with_un, "Krt10", save_plot = paste0(save_dir, "/figure_s4cX.pdf"))
tSNE_PCA(mtecCombined_with_un, "Trpm5", save_plot = paste0(save_dir, "/figure_s4cXI.pdf"))
tSNE_PCA(mtecCombined_with_un, "Dclk1", save_plot = paste0(save_dir, "/figure_s4cXII.pdf"))
tSNE_PCA(mtecCombined_with_un, "Cldn4", save_plot = paste0(save_dir, "/figure_s4cXIII.pdf"))
tSNE_PCA(mtecCombined_with_un, "GFP", save_plot = paste0(save_dir, "/figure_s4cXIV.pdf"))

fig_list <- c(fig_list, "supplemental_figure_4")
#########################
# Supplemental Figure 5 #
#########################
print("supplemental_figure_5")

# Supplemental Figure 5a Violin markers for all as in S2
trio_plots_median(mtecCombined, geneset = c("Ackr4", "Ccl21a", "Aire"),
  cell_cycle = FALSE, plot_violin = TRUE, jitter_and_violin = FALSE,
  plot_jitter = FALSE, color = stage_color, sep_by = "stage",
  save_plot = paste0(save_dir, "/figure_s5aI.pdf"))
trio_plots_median(mtecCombined, geneset = c("Krt10", "Trpm5", "GFP"),
  cell_cycle = FALSE, plot_violin = TRUE, jitter_and_violin = FALSE,
  plot_jitter = FALSE, color = stage_color, sep_by = "stage",
  save_plot = paste0(save_dir, "/figure_s5aII.pdf"))

# Supplemental Figure 5b Correlation of TP5 with WT
mtecCombExp <- Seurat::SetAllIdent(mtecCombSub, id = "exp")
mtec_tp5 <- Seurat::SubsetData(mtecCombExp, ident.use = "timepoint5")
mtec_end <- Seurat::SubsetData(mtecCombExp, ident.use = "isoControlEnd")
pdf(paste0(save_dir, "/figure_s5b.pdf"))
master_plot(mtec_tp5, "wk_10", mtec_end, "ctl_wk_10", stage_color_df)
dev.off()

# Supplemental Figure 5c expression of other Ccl21 genes (Ccl21c is not in any gene list)
trio_plots_median(mtecCombined, geneset = c("Ccl21a", "Ccl21a.1", "Ccl21b.1"),
  cell_cycle = FALSE, plot_violin = TRUE, jitter_and_violin = FALSE,
  plot_jitter = FALSE, color = stage_color, sep_by = "stage",
  save_plot = paste0(save_dir, "/figure_s5c.pdf"))



fig_list <- c(fig_list, "supplemental_figure_5")
#########################
# Supplemental Figure 6 #
#########################
print("supplemental_figure_6")
# Supplemental Figure 6a violin plots of extra genes of interest
trio_plots_median(projenitor_mtec, geneset = c("Hmgb2", "Tubb5", "Stmn1"),
  cell_cycle = FALSE, plot_violin = TRUE, jitter_and_violin = FALSE,
  plot_jitter = FALSE, color = timecourse_color, sep_by = "pub_exp",
  save_plot = paste0(save_dir, "/figure_s6a.pdf"))

# Supplemental Figure 6b violing plots of extra genes of interest, extra analysis
trio_plots_median(no_at_mtec, geneset = c("Hmgb2", "Tubb5", "Stmn1"), cell_cycle = FALSE,
  plot_violin = TRUE, jitter_and_violin = FALSE, plot_jitter = FALSE,
  sep_by = "cluster", color = reanalysis_colors,
  save_plot = paste0(save_dir, "/figure_s6b.pdf"))

# Supplemental Figure 6c bar plots of cycling mTECS
mtecCombExp@meta.data$Ccl21a <- mtecCombExp@data["Ccl21a", ]
mtecCombExp@meta.data$Aire <- mtecCombExp@data["Aire", ]
mtecCombExp@meta.data$db_pos <- mtecCombExp@meta.data$Ccl21a > 4 & 
  mtecCombExp@meta.data$Aire > 0.5

mtecCombExp@meta.data$Aire_exp <- mtecCombExp@meta.data$Aire > 0.5 &
  mtecCombExp@meta.data$Ccl21a <= 4
mtecCombExp@meta.data$Ccl21a_exp <- mtecCombExp@meta.data$Ccl21a > 4 &
  mtecCombExp@meta.data$Aire <= 0.5
mtecCombExp@meta.data$Aire_and_Ccl21a <- mtecCombExp@meta.data$Aire > 0.5 &
  mtecCombExp@meta.data$Ccl21a > 4
mtecCombExp@meta.data$db_neg <- mtecCombExp@meta.data$Ccl21a <= 4 &
  mtecCombExp@meta.data$Aire <= 0.5

mtecComb_cycle <- Seurat::SubsetData(mtecCombExp, ident.remove = "aireTrace")
# cells_keep <- rownames(mtecComb_cycle@meta.data[
#   mtecComb_cycle@meta.data$cycle_phase == "G2M", ])
# mtecComb_cycle <- Seurat::SubsetData(mtecComb_cycle, cells.use = cells_keep)

mtecComb_cycle <- Seurat::SetAllIdent(mtecComb_cycle, id = "stage")
mtecComb_cycle <- Seurat::SubsetData(mtecComb_cycle, ident.use = "TAC_TEC")

blue_palette <- c(percent_Aire = "#111E6C",
                  percent_Ccl21a = "#008081",
                  percent_db_pos = "#0080FF",
                  percent_db_neg = "#4C516D")

timepoints <- c("isoControlBeg", "isoControlEnd", "timepoint1",
                "timepoint2", "timepoint3", "timepoint5")

percent_Aire <- sapply(timepoints, USE.NAMES = TRUE,
  function(x) percent_ident(mtecComb_cycle,
  data_set = x, meta_data_col = "exp", ident = "Aire_exp"))
percent_Ccl21 <- sapply(timepoints, USE.NAMES = TRUE,
  function(x) percent_ident(mtecComb_cycle,
  data_set = x, meta_data_col = "exp", ident = "Ccl21a_exp"))
percent_dp_pos <- sapply(timepoints, USE.NAMES = TRUE,
  function(x) percent_ident(mtecComb_cycle,
  data_set = x, meta_data_col = "exp", ident = "Aire_and_Ccl21a"))
percent_db_neg <- sapply(timepoints, USE.NAMES = TRUE,
  function(x) percent_ident(mtecComb_cycle,
  data_set = x, meta_data_col = "exp", ident = "db_neg"))



percent_all <- data.frame(percent_Aire = percent_Aire,
                          percent_Ccl21a = percent_Ccl21,
                          percent_db_pos = percent_dp_pos,
                          percent_db_neg = percent_db_neg)
percent_all$experiment <- rownames(percent_all)


percent_all_m <- tidyr::gather(data = percent_all, key = classification,
  value = percent, percent_Aire:percent_db_neg,
  factor_key = TRUE)

pdf(paste0(save_dir, "/figure_s6c.pdf"))
ggplot2::ggplot(percent_all_m, ggplot2::aes(x = experiment, y = percent,
  fill = classification)) +
  ggplot2::geom_bar(position = "stack", stat = "identity") +
  ggplot2::scale_fill_manual(values = blue_palette)
dev.off()

# Supplemental Figure 6d bar plots of cycling counts
count_Aire <- sapply(timepoints, USE.NAMES = TRUE,
  function(x) percent_ident(mtecComb_cycle,
  data_set = x, meta_data_col = "exp", ident = "Aire_exp", count = TRUE))
count_Ccl21 <- sapply(timepoints, USE.NAMES = TRUE,
  function(x) percent_ident(mtecComb_cycle,
  data_set = x, meta_data_col = "exp", ident = "Ccl21a_exp", count = TRUE))
count_dp_pos <- sapply(timepoints, USE.NAMES = TRUE,
  function(x) percent_ident(mtecComb_cycle,
  data_set = x, meta_data_col = "exp", ident = "Aire_and_Ccl21a", count = TRUE))
count_db_neg <- sapply(timepoints, USE.NAMES = TRUE,
  function(x) percent_ident(mtecComb_cycle,
  data_set = x, meta_data_col = "exp", ident = "db_neg", count = TRUE))

blue_palette <- c(count_Aire = "#111E6C",
                  count_Ccl21a = "#008081",
                  count_db_pos = "#0080FF",
                  count_db_neg = "#4C516D")

count_all <- data.frame(count_Aire = count_Aire,
                        count_Ccl21a = count_Ccl21,
                        count_db_pos = count_dp_pos,
                        count_db_neg = count_db_neg)

count_all$experiment <- rownames(count_all)


count_all_m <- tidyr::gather(data = count_all, key = classification,
  value = percent, count_Aire:count_db_neg,
  factor_key = TRUE)
pdf(paste0(save_dir, "/figure_s6d.pdf"))
ggplot2::ggplot(count_all_m, ggplot2::aes(x = experiment, y = percent,
  fill = classification)) +
  ggplot2::geom_bar(position = "dodge", stat = "identity") +
  ggplot2::scale_fill_manual(values = blue_palette)
dev.off()

# Supplemental Figure 6e Flow of Ki67 cells in ablation

fig_list <- c(fig_list, "supplemental_figure_6")

##########################################################################################

#########################
# Supplemental Figure 7 #
#########################
print("supplemental_figure_7")

# Supplementary Figure 7a UMAPs of TRAs in recovery
limit_list <- list(tra_fantom = c(0, 0.100),
                   aire_genes = c(0, 0.100),
                   fezf2_genes = c(0, 0.300))

plot_names_fig5 <- list(isoControlBeg = paste0(save_dir, "/figure_s7aI.pdf"),
  isoControlEnd = paste0(save_dir, "/figure_s7aII.pdf"),
  timepoint1 = paste0(save_dir, "/figure_s7aIII.pdf"),
  timepoint2 = paste0(save_dir, "/figure_s7aIV.pdf"),
  timepoint3 = paste0(save_dir, "/figure_s7aV.pdf"),
  timepoint5 = paste0(save_dir, "/figure_s7aVI.pdf"))

names(plot_names_fig5) <- new_exp_names[names(plot_names_fig5)]

lapply(names(plot_names_fig5), function(x) plot_gene_set(mtecCombined,
                                            gene_set = gene_lists[["tra_fantom"]],
                                            plot_name = "tra_fantom",
                                            one_dataset = FALSE,
                                            data_set = x,
                                            meta_data_col = "pub_exp",
                                            limits = limit_list[["tra_fantom"]],
                                            save_plot = plot_names_fig5[[x]]))



# Supplementary Figure 7b cumulative number of genes in WT not downsampled
wt_aire <- Seurat::SetAllIdent(mtec_aire_positive, id = "exp")
wt_aire <- Seurat::SubsetData(wt_aire,
  ident.use = c("isoControlBeg", "isoControlEnd"),
  subset.raw = TRUE)
wt_matrix <- as.matrix(wt_aire@raw.data)

list_names <- c("all_other_genes", "tra_fantom", "aire_genes", "fezf2_genes")

cumFreqAll <- sapply(list_names, function(x) cumFreqFunc(wt_matrix,
  gene_lists[[x]], x))

cumFreq_df <- do.call("rbind", cumFreqAll)

cumFreqPlot <- ggplot2::ggplot(cumFreq_df, ggplot2::aes(x = ID, y = list_percent,
                                         group = gene_list,
                                         color = gene_list)) +
  ggplot2::geom_line() +
  ggplot2::scale_color_brewer(palette = "Set2")

ggplot2::ggsave(paste0(save_dir, "/figure_s7b.pdf"), plot = cumFreqPlot)

to_plot_sup <- c("tra_fantom", "all_other_genes", "aire_genes", "fezf2_genes")

counts_df_plot_sup <- counts_df_m[counts_df_m$gene_list %in% to_plot_sup, ]

# Supplemental Figure 7c box plots for gene counts after downsampling
full_plot <- ggplot2::ggplot(counts_df_plot_sup, ggplot2::aes(x = gene_list,
                                                          y = gene_count,
                                                          fill = pub_exp)) +
  ggplot2::geom_boxplot() +
  ggplot2::scale_fill_manual(values = timecourse_color) +
  #ggplot2::theme_classic() +
  ggpubr::stat_compare_means(method = "anova", size = 2, label.y = 6150)

zoom_plot <- ggplot2::ggplot(counts_df_short, ggplot2::aes(x = gene_list,
                                                           y = gene_count,
                                                           fill = pub_exp)) +
  ggplot2::geom_boxplot(show.legend = FALSE) +
  ggplot2::scale_fill_manual(values = timecourse_color) +
  #ggplot2::theme_classic() +
  ggplot2::theme(panel.background = ggplot2::element_blank(),
                 axis.title.x = ggplot2::element_blank(),
                 axis.title.y = ggplot2::element_blank(),
                 panel.border = ggplot2::element_rect(color = "black",
                                                      fill = NA,
                                                      size = 1)) +
  ggpubr::stat_compare_means(method = "anova", size = 2, label.y = 300)

zoom_plot_g <- ggplot2::ggplotGrob(zoom_plot)

all_plots <- full_plot + ggplot2::annotation_custom(grob = zoom_plot_g,
                                                    xmin = 1.5,
                                                    xmax = Inf,
                                                    ymin = 1000,
                                                    ymax = Inf) +
  ggplot2::annotation_custom(grob = grid::rectGrob(gp = grid::gpar(fill = NA)),
                             xmin = 1.5,
                             xmax = Inf,
                             ymin = -Inf,
                             ymax = 500)
pdf(paste0(save_dir, "/figure_s7c.pdf"))
all_plots

dev.off()

# Supplemental Figure 7d histograms of genes after downsampling
if (bootstrap) {
  source(bootstrap_script)

}

fig_list <- c(fig_list, "supplemental_figure_7")


#########################
# Supplemental Figure 8 #
#########################
print("supplemental_figure_8")
heatmap_names <- list(isoControlBeg = paste0(save_dir, "/figure_s8I.pdf"),
  isoControlEnd = paste0(save_dir, "/figure_s8II.pdf"),
  timepoint1 = paste0(save_dir, "/figure_s8III.pdf"),
  timepoint2 = paste0(save_dir, "/figure_s8IV.pdf"),
  timepoint3 = paste0(save_dir, "/figure_s8V.pdf"),
  timepoint5 = paste0(save_dir, "/figure_s8VI.pdf"))

# Supplemental Figure 8 heatmaps of marker genes for all samples
mtec_wt_plot <- mtec_wt
mtec_wt_plot@assay$DE <- NULL
mtec_wt_plot <- Seurat::SubsetData(mtec_wt_plot, ident.remove = "unknown")
mtec.markers <- Seurat::FindAllMarkers(object = mtec_wt_plot, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
top10 <- mtec.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
print(head(top10))

lapply(names(plot_names), function(x) plot_marker_heatmap(mtecCombSub, subset_val = x,
  gene_df <- top10, save_plot = heatmap_names[[x]]))

fig_list <- c(fig_list, "supplemental_figure_8")

#########################
# Supplemental Figure 9 #
#########################


print("supplemental_figure_9")

# Supplementary Figure 9a nGene and nUMI in all cells before correction
at_timecourse <- c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#65C42D",
  "#FFB62F", "#E5C494")
#at_timecourse <- RColorBrewer::brewer.pal(7, "Set2")


mtec_meta_data <- mtecCombined@meta.data

plot_all_cell_count <- ggplot2::ggplot(mtec_meta_data,
  ggplot2::aes(x = pub_exp, y = nGene, group = pub_exp,
               fill = pub_exp)) + 
  ggplot2::geom_violin(scale = "width") +
  ggplot2::scale_fill_manual(values = at_timecourse) + 
  ggplot2::ylab("nGene")

ggplot2::ggsave(paste0(save_dir, "/figure_s9aI.pdf"), plot = plot_all_cell_count)

plot_all_umi_count <- ggplot2::ggplot(mtec_meta_data,
  ggplot2::aes(x = pub_exp, y = nUMI, group = pub_exp,
               fill = pub_exp)) + 
  ggplot2::geom_violin(scale = "width") +
  ggplot2::scale_fill_manual(values = at_timecourse) + 
  ggplot2::ylab("nUMI")

ggplot2::ggsave(paste0(save_dir, "/figure_s9aII.pdf"), plot = plot_all_umi_count)


# Supplemetary Figure 9b nGene and nUMI just Aire Positive before correction

aire_positive_meta <- mtec_meta_data[mtec_meta_data$stage == "Aire_positive", ]

plot_all_cell_count_ap <- ggplot2::ggplot(aire_positive_meta,
  ggplot2::aes(x = pub_exp, y = nGene, group = pub_exp,
               fill = pub_exp)) + 
  ggplot2::geom_violin(scale = "width") +
  ggplot2::scale_fill_manual(values = at_timecourse) + 
  ggplot2::ylab("nGene") +
  ggplot2::xlab("experiment")

ggplot2::ggsave(paste0(save_dir, "/figure_s9bI.pdf"), plot = plot_all_cell_count_ap)

plot_all_umi_count_ap <- ggplot2::ggplot(aire_positive_meta,
  ggplot2::aes(x = pub_exp, y = nUMI, group = pub_exp,
               fill = pub_exp)) + 
  ggplot2::geom_violin(scale = "width") +
  ggplot2::scale_fill_manual(values = at_timecourse, name = "experiment") + 
  ggplot2::ylab("nUMI") +
  ggplot2::xlab("experiment")

ggplot2::ggsave(paste0(save_dir, "/figure_s9bII.pdf"), plot = plot_all_umi_count_ap)

# Supplementary Figure 9c nGene and nUMI just aire positive after correction
umis <- lapply(names(percents_counts_all), function(x)
  get_perc_count(percents_counts_all, x, data_type = "umi"))

umis_df <- do.call(rbind, umis)
umis_df$exp <- sub("Aire_positive_", "", umis_df$exp)
umis_df$pub_exp <- new_exp_names[umis_df$exp]
umis_df$pub_exp <- factor(umis_df$pub_exp,
                          levels = unname(new_exp_names))
umis_df_m <- reshape2::melt(umis_df, variable.name = "gene_list",
  value.name = "gene_count")


umis_df_all <- umis_df_m[umis_df_m$gene_list == "protein_coding", ]
counts_df_all <- counts_df_m[counts_df_m$gene_list == "protein_coding", ]

plot_corrected_gene <- ggplot2::ggplot(counts_df_all,
  ggplot2::aes(x = pub_exp, y = gene_count, 
               group = pub_exp, fill = pub_exp)) +
  ggplot2::geom_violin(scale = "width") +
  ggplot2::scale_fill_manual(values = timecourse_color, name = "experiment") + 
  ggplot2::ylab("nGene") +
  ggplot2::xlab("experiment")

ggplot2::ggsave(paste0(save_dir, "/figure_s9cI.pdf"), plot = plot_corrected_gene)


plot_corrected_umi <- ggplot2::ggplot(umis_df_all,
  ggplot2::aes(x = pub_exp, y = gene_count, 
               group = pub_exp, fill = pub_exp)) +
  ggplot2::geom_violin(scale = "width") +
  ggplot2::scale_fill_manual(values = timecourse_color) + 
  ggplot2::ylab("nUMI") +
  ggplot2::xlab("experiment")

ggplot2::ggsave(paste0(save_dir, "/figure_s9cII.pdf"), plot = plot_corrected_umi)

# Supplementary figure 9d Dropouts of different house keeping genes before and after
experiments <- as.character(unique(no_at_mtec_aire@meta.data$exp))

# Housekeeping genes to test
gene_list <- c("Chmp2a", "Emc7", "Psmb2", "Psmb4", "Vcp", "Gapdh")

# Before downsample
dropout_list <- sapply(experiments, function(x) get_cell_matrix(no_at_mtec_aire,
                                                              gene_list = gene_list,
                                                              subset = x))

colnames(dropout_list) <- new_exp_names[colnames(dropout_list)]

# After downsample
dropout_downsample <- sapply(experiments, function(x) 
                           get_cell_matrix(no_at_mtec_aire, gene_list = gene_list,
                                         subset = x, downsample_matrix = TRUE,
                                         lowest_UMI = lowest_UMI))

colnames(dropout_downsample) <- new_exp_names[colnames(dropout_downsample)]

# Change dfs into form for ggplot2
dropout_list_m <- reshape2::melt(dropout_list)

dropout_downsample_m <- reshape2::melt(dropout_downsample)

# Make plots
dropout_list_m$Var2 <- factor(dropout_list_m$Var2,
                              levels = unname(new_exp_names))

dropout_plot <- ggplot2::ggplot(dropout_list_m,
  ggplot2::aes(x = Var1, y = value, fill = Var2)) +
  ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge()) +
  ggplot2::scale_fill_manual(values = timecourse_color, name = "experiment") +
  ggplot2::xlab("housekeeping gene") +
  ggplot2::ylab("dropout percent")

ggplot2::ggsave(paste0(save_dir, "/figure_s9dI.pdf"), plot = dropout_plot)

dropout_downsample_m$Var2 <- factor(dropout_downsample_m$Var2,
                              levels = unname(new_exp_names))

dropout_downsample_plot <- ggplot2::ggplot(dropout_downsample_m,
  ggplot2::aes(x = Var1, y = value, fill = Var2)) +
  ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge()) +
  ggplot2::scale_fill_manual(values = timecourse_color, name = "experiment") +
  ggplot2::xlab("housekeeping gene") +
  ggplot2::ylab("dropout percent")

ggplot2::ggsave(paste0(save_dir, "/figure_s9dII.pdf"), plot = dropout_downsample_plot)


fig_list <- c(fig_list, "supplemental_figure_6")

write.table(fig_list, save_file)
