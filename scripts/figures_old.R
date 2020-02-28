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
aireTrace <- snakemake@input[[1]]
controls <- snakemake@input[[2]]
allSamples <- snakemake@input[[3]]
controls_slingshot <- snakemake@input[[4]]
allSamples_slingshot <- snakemake@input[[5]]
early_aire_mtec <- snakemake@input[[6]]
save_file <- snakemake@output[[1]]
data_directory <- snakemake@params[[1]]
save_dir <- snakemake@params[[2]]

# Colors for plotting
stage_color_df <- data.frame("cTEC" = "#CC6600", "Ccl21a_high" = "#009933",
                            "Early_mTEC" = "#0066CC", "Aire_positive" = "#660099",
                            "Late_Aire" = "#FF0000", "Tuft" = "#990000",
                            "unknown" = "#FDFBFB")

stage_color <- t(stage_color_df)[ , 1]

# Load in data
mtec <- get(load(aireTrace))

mtec_wt <- get(load(controls))

mtecCombined <- get(load(allSamples))

wt_slingshot <- get(load(controls_slingshot))

all_slingshot <- get(load(allSamples_slingshot))

progenitor_mtec <- get(load(early_aire_mtec))

fig_list <- list()

TFs <- get(load(paste0(data_directory, "TFs.rda")))

TFs_all <- c("H2afz", "Top2a", "Hmgb1", "Hmgn1", "H2afx", as.character(TFs))

bootstrap <- FALSE

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 18))

############
# Figure 1 #
############
print("Figure 1")

mtec@meta.data$stage <- as.character(mtec@meta.data$stage)

mtec@meta.data$stage[mtec@meta.data$stage ==
  "Cortico_medullary"] <- "cTEC"
mtec@meta.data$stage[mtec@meta.data$stage ==
  "Early_Aire"] <- "Early_mTEC"
# Make plots from aire trace
mtec_no_un <- mtec
mtec_no_un@assay$DE <- NULL
########################################################################
# Fix this section by reruning aireTrace analysis driver
mtec_no_un@meta.data$stage[is.na(mtec_no_un@meta.data$stage)] <- "other"


stage_levels <- c("cTEC", "Ccl21a_high", "Early_mTEC",
                  "Aire_positive", "Late_Aire", "Tuft", "other")
mtec_no_un@meta.data$stage <- factor(mtec_no_un@meta.data$stage,
  levels = stage_levels)
mtec_no_un <- Seurat::SetAllIdent(mtec_no_un, id = "res.0.6")
mtec_no_un <- Seurat::SetAllIdent(mtec_no_un, id = "stage")
mtec_no_un@ident <- factor(mtec_no_un@ident, levels = stage_levels)
#########################################################################
mtec_no_un <- Seurat::SubsetData(mtec_no_un, ident.remove = "other")



# Figure 1b
tSNE_PCA(mtec_no_un, "stage", color = stage_color, show_legend = FALSE,
  save_plot = paste0(save_dir, "/figure_1b.pdf"))
tSNE_PCA(mtec_no_un, "stage", color = stage_color,
  save_plot = paste0(save_dir, "/figure_1bI.pdf"))

# Figure 1c
# Violin plots of marker genes, pick a few more here
trio_plots_median(mtec_no_un, geneset = c("Aire", "Ccl21a", "Trpm5"),
  cell_cycle = FALSE, plot_violin = TRUE, jitter_and_violin = FALSE,
  plot_jitter = FALSE, color = stage_color, sep_by = "cluster",
  save_plot = paste0(save_dir, "/figure_1cI.pdf"))
trio_plots_median(mtec_no_un, geneset = c("GFP", "Ackr4", "Krt10"),
  cell_cycle = FALSE, plot_violin = TRUE, jitter_and_violin = FALSE,
  plot_jitter = FALSE, color = stage_color, sep_by = "cluster",
  save_plot = paste0(save_dir, "/figure_1cII.pdf"))

# Figure 1d
# UMAPs of Ccl21a, Aire, and GFP
tSNE_PCA(mtec_no_un, "Ccl21a", save_plot = paste0(save_dir, "/figure_1dI.pdf"))
tSNE_PCA(mtec_no_un, "Aire", save_plot = paste0(save_dir, "/figure_1dII.pdf"))
tSNE_PCA(mtec_no_un, "GFP", save_plot = paste0(save_dir, "/figure_1dIII.pdf"))

fig_list <- c(fig_list, "figure_1")
############
# Figure 2 #
############

print("Figure 2")
load(paste0(data_directory, "gene_lists.rda"))


# Figure 2a
# Dot plot of most interesting markers and stem cell markers
markers_to_plot_full <- c("Ackr4", "Psmb11", "Krt5", "Ccl21a", "Aire", "Fezf2",
  "Cldn4", "Spink5", "Trpm5", "Dclk1", "Hmgb2", "Hmgn2", "Hmgb1", "H2afx",
  "Stmn1", "Tubb5")

pdf(paste0(save_dir, "/figure_2aI.pdf"))
dot_plot <- Seurat::DotPlot(mtec_no_un, genes.plot = rev(markers_to_plot_full),
                            cols.use = c("blue", "red"), x.lab.rot = T,
                            plot.legend = F, dot.scale = 8, do.return = T)


dev.off()


pdf(paste0(save_dir, "/figure_2aII.pdf"))
dot_plot <- Seurat::DotPlot(mtec_no_un, genes.plot = rev(markers_to_plot_full),
                            cols.use = c("blue", "red"), x.lab.rot = T,
                            plot.legend = T, dot.scale = 8, do.return = T)
dev.off()


# Figure 2b
# Jitter plots of chromatin modifiers overlayed with cell cycle state
trio_plots(mtec_no_un, geneset = c("Hmgb2", "Tubb5", "Stmn1"),
  cell_cycle = TRUE, jitter_and_violin = FALSE, plot_jitter = TRUE,
  color = stage_color, sep_by = "cluster",
  save_plot = paste0(save_dir, "/figure_2bI.pdf"), group_color = FALSE)
trio_plots(mtec_no_un, geneset = c("Aire", "Ccl21a", "Fezf2"),
  cell_cycle = TRUE, jitter_and_violin = FALSE, plot_jitter = TRUE,
  color = stage_color, sep_by = "cluster",
  save_plot = paste0(save_dir, "/figure_2bII.pdf"), group_color = FALSE)

# Figure 2c
# Flow of Mki76/Aire

# Figure 2d
# SC velocity

# Figure 2e
mtec_wt@meta.data$stage <- as.character(mtec_wt@meta.data$stage)


mtec_wt@meta.data$stage[mtec_wt@meta.data$stage ==
  "Cortico_medullary"] <- "cTEC"
mtec_wt@meta.data$stage[mtec_wt@meta.data$stage ==
  "Early_Aire"] <- "Early_mTEC"


# Remove unknown cells
mtec_wt_plot <- mtec_wt
mtec_wt_plot@assay$DE <- NULL

mtec_wt_plot <- Seurat::SetAllIdent(mtec_wt_plot, id = "stage")
mtec_wt_plot <- Seurat::SubsetData(mtec_wt_plot, ident.remove = "unknown")
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



mtec_wt_plot@meta.data$stage <- factor(mtec_wt_plot@meta.data$stage,
  levels = stage_levels)

# Plot each of the genes and gene sets in pseudotime, end at 16 because the
# is where the "unknown" cells are
plot_list <- lapply(names(plot_names), function(x) plot_sling_pseudotime(
  seurat_object = mtec_wt_plot, sling_object = wt_slingshot, y_val = x,
  col_by = "stage", pseudotime_curve = "curve1", color = stage_color,
  range = c(0, 16), save_plot = plot_names[[x]]))

fig_list <- c(fig_list, "figure_2")

# ############
# # Figure 3 #
# ############

print("Figure 3")

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


mtecCombined <- Seurat::SetAllIdent(mtecCombined, id = "stage")
mtecCombined <- Seurat::SubsetData(mtecCombined, ident.remove = "unknown")
mtecCombined@meta.data$pub_exp <- new_exp_names[mtecCombined@meta.data$exp]
mtecCombined@meta.data$pub_exp <- factor(mtecCombined@meta.data$pub_exp,
                                         levels = unname(new_exp_names))
mtecCombined@meta.data$stage <- as.character(mtecCombined@meta.data$stage)
mtecCombined@meta.data$stage[mtecCombined@meta.data$stage ==
  "Cortico_medullary"] <- "cTEC"
mtecCombined@meta.data$stage[mtecCombined@meta.data$stage ==
  "Early_Aire"] <- "Early_mTEC"
mtecCombined@meta.data$stage <- factor(mtecCombined@meta.data$stage,
  levels = stage_levels)
mtecCombined@meta.data$stage_exp <- paste0(mtecCombined@meta.data$stage,
  "_", mtecCombined@meta.data$exp)

mtecCombined@meta.data$at_stage <- as.character(mtecCombined@meta.data$at_stage)
mtecCombined@meta.data$at_stage[mtecCombined@meta.data$at_stage ==
  "Cortico_medullary"] <- "cTEC"
mtecCombined@meta.data$at_stage[mtecCombined@meta.data$at_stage ==
  "Early_Aire"] <- "Early_mTEC"
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

# Figure 2B
# UMAP of all cells Either put the key on both or remove the key from both
tSNE_PCA(mtecCombined, "stage", save_plot = paste0(save_dir, "/figure_3bI.pdf"),
  color = stage_color, show_legend = FALSE)

# UMAP highlighting aire_trace cells colored by aire_trace labels
full_umap(mtecCombined, "aireTrace", col_by = "at_stage", color = stage_color,
  save_plot = paste0(save_dir, "/figure_3bII.pdf"), show_legend = FALSE)

tSNE_PCA(mtecCombined, "stage", save_plot = paste0(save_dir, "/figure_3bIII.pdf"),
  color = stage_color, show_legend = TRUE)

# Figure 3C
# Barplots of recovery
mtecCombSub <- mtecCombined
mtecCombSub@assay$ablation_DE <- NULL

stage_list_all <- lapply(data_sets, function(x) populations_dfs_new(mtecCombSub,
                         x, subsample = TRUE, subsample_by = "pub_exp"))
stage_df_all <- do.call("rbind", stage_list_all)

stage_df_all$sample <- factor(stage_df_all$sample, levels = unname(new_exp_names))

population_plots(stage_df_all, color = stage_color,
  save_plot = paste0(save_dir, "/figure_3c.pdf"))

# Figure 3D
# Umap with slingshot overlay
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


# Figure 3E
# Timecourse in pseudotime
mtec_no_at <- mtecCombSub
mtec_no_at <- Seurat::SetAllIdent(mtec_no_at, id = "exp")
mtec_no_at <- Seurat::SubsetData(mtec_no_at, ident.remove = "aireTrace")
plot_sling_pseudotime(
  seurat_object = mtec_no_at, sling_object = all_slingshot, y_val = "pub_exp",
  col_by = "pub_exp", pseudotime_curve = "curve3", color = timecourse_color,
  plot_type = "density", width = 7, height = 7,
  save_plot = paste0(save_dir, "/figure_3e.pdf"))

fig_list <- c(fig_list, "figure_3")
# ###############################################################################

# ############
# # Figure 4 #
# ############
print("Figure 4")


lowest_UMI_exp <- "timepoint2"

downsample_UMI <- TRUE

if (downsample_UMI) {
  
  lowest_UMI <- get_umi(mtecCombSub, subset_seurat = TRUE, subset_by = "stage_exp",
    subset_val = paste0("Aire_positive_", lowest_UMI_exp))
} else {
  lowest_UMI <- NULL
}

# Figure 4a
# Start with full UMAP only coloring the Aire positive cluster
highlight_one_group(mtecCombined, meta_data_col = "stage", group = "Aire_positive",
  color_df = stage_color, show_legend = FALSE,
  save_plot = paste0(save_dir, "/figure_4a.pdf"))

highlight_one_group(mtecCombined, meta_data_col = "stage", group = "Aire_positive",
  color_df = stage_color, show_legend = TRUE,
  save_plot = paste0(save_dir, "/figure_4aI.pdf"))

# Figure 4b
# Average expression of set of genes over ablation
average_gene_list <- c("Aire", "Fezf2", "Gapdh", "Emc7", "Tnfrsf11a")

mtec_aire_positive <- Seurat::SubsetData(mtecCombSub, ident.use = "Aire_positive",
  subset.raw = TRUE)
cells_use <- rownames(mtec_aire_positive@meta.data)[mtec_aire_positive@meta.data$exp !=
                                                 "aireTrace"]

no_at_mtec_aire <- Seurat::SubsetData(mtec_aire_positive, cells.use = cells_use)

# Make the lines thicker here
plot_avg_exp_genes(no_at_mtec_aire, average_gene_list,
                   save_plot = paste0(save_dir, "/figure_4b.pdf"),
                   avg_expr_id = "pub_exp")

trio_plots_median(no_at_mtec_aire, geneset = c("Aire", "Fezf2", "Gapdh"),
  cell_cycle = FALSE, jitter_and_violin = FALSE, plot_jitter = FALSE,
  plot_violin = TRUE, color = timecourse_color, sep_by = "pub_exp",
  save_plot = paste0(save_dir, "/figure_4bII.pdf"))

trio_plots_median(no_at_mtec_aire, geneset = c("Aire", "Fezf2", "Gapdh"),
  cell_cycle = FALSE, jitter_and_violin = FALSE, plot_jitter = FALSE,
  plot_violin = TRUE, color = timecourse_color, sep_by = "pub_exp",
  save_plot = paste0(save_dir, "/figure_4bIII.pdf"))

mtecCombined_all <- no_at_mtec_aire
# This is an okay place for a for loop (recursion) 
# http://adv-r.had.co.nz/Functionals.html
for (gene_set in names(gene_lists)) {
  mtecCombined_all <- plot_gene_set(mtecCombined_all,
    gene_lists[[gene_set]], gene_set, make_plot = FALSE)
}


# Figure 4c
# Violin plot of average expression of gene sets in all cells separated
# by experiment
gene_sets <- c("all_other_genes", "tra_fantom", "aire_genes")
trio_plots_median(mtecCombined_all, geneset = gene_sets,
  cell_cycle = FALSE, plot_violin = TRUE, jitter_and_violin = FALSE,
  plot_jitter = FALSE, color = timecourse_color, sep_by = "pub_exp",
  save_plot = paste0(save_dir, "/figure_4c.pdf"))


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
#percents_m$pub_exp <- new_exp_names[percents_m$exp]
percents_m$pub_exp <- factor(percents_m$pub_exp,
                             levels = unname(new_exp_names))
percents_plot <- percents_m[percents_m$gene_list %in% to_plot, ]

# Figure 4d
# Percent of gene lists
pdf(paste0(save_dir, "/figure_4d.pdf"))
ggplot2::ggplot(percents_plot, ggplot2::aes(x = pub_exp, y = percent_of_genes,
                                         group = gene_list, color = gene_list)) +
  ggplot2::geom_line(size = 2) +
  ggplot2::ylim(0,1) +
  ggplot2::scale_color_brewer(palette = "Dark2")
  #ggplot2::theme_classic()

dev.off()

# Figure 4e
# Number of genes per cell
full_plot <- ggplot2::ggplot(counts_df_plot, ggplot2::aes(x = gene_list,
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
pdf(paste0(save_dir, "/figure_4e.pdf"))
all_plots

dev.off()

fig_list <- c(fig_list, "figure_4")

######################
# Figure 4 extension #
######################

aire_plot_names <- c(isoControlBeg = paste0(save_dir, "/figure_e4Ia.pdf"),
                     isoControlEnd = paste0(save_dir, "/figure_e4If.pdf"),
                     timepoint1 = paste0(save_dir, "/figure_e4Ib.pdf"),
                     timepoint2 = paste0(save_dir, "/figure_e4Ic.pdf"),
                     timepoint3 = paste0(save_dir, "/figure_e4Id.pdf"),
                     timepoint5 = paste0(save_dir, "/figure_e4Ie.pdf"))

lapply(names(aire_plot_names), function(x) full_umap(mtecCombined,
  data_set = x, col_by = "Aire", show_legend = TRUE,
  save_plot = aire_plot_names[[x]]))

fezf2_plot_names <- c(isoControlBeg = paste0(save_dir, "/figure_e4IIa.pdf"),
                      isoControlEnd = paste0(save_dir, "/figure_e4IIf.pdf"),
                      timepoint1 = paste0(save_dir, "/figure_e4IIb.pdf"),
                      timepoint2 = paste0(save_dir, "/figure_e4IIc.pdf"),
                      timepoint3 = paste0(save_dir, "/figure_e4IId.pdf"),
                      timepoint5 = paste0(save_dir, "/figure_e4IIe.pdf"))

lapply(names(fezf2_plot_names), function(x) full_umap(mtecCombined,
  data_set = x, col_by = "Fezf2", show_legend = TRUE,
  save_plot = fezf2_plot_names[[x]]))


# ############
# # Figure 5 #
# ############

# This still needs some work

print("Figure 5")
reanalysis_colors <- c("#603E95", "#009DA1", "#FAC22B", "#D7255D")

projenitor_mtec <- mtec_no_at

projenitor_mtec <- Seurat::SetAllIdent(projenitor_mtec, id = "stage")
projenitor_mtec <- Seurat::SubsetData(projenitor_mtec, ident.use = "Early_mTEC")

progenitor_mtec@meta.data$pub_exp <- new_exp_names[progenitor_mtec@meta.data$exp]
progenitor_mtec@meta.data$pub_exp <- factor(progenitor_mtec@meta.data$pub_exp,
                                         levels = unname(new_exp_names))
cells_use <- rownames(progenitor_mtec@meta.data)[progenitor_mtec@meta.data$exp !=
                                                 "aireTrace"]

no_at_mtec <- Seurat::SubsetData(progenitor_mtec, cells.use = cells_use)

# Figure 5a
# Highlight just the early aire cells 
highlight_one_group(mtecCombined, meta_data_col = "stage", group = "Early_mTEC",
  color_df = stage_color, show_legend = FALSE,
  save_plot = paste0(save_dir, "/figure_5a.pdf"))

highlight_one_group(mtecCombined, meta_data_col = "stage", group = "Early_mTEC",
  color_df = stage_color, show_legend = TRUE,
  save_plot = paste0(save_dir, "/figure_5aI.pdf"))

# Figure 5b
# % of cells in G2
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

# Figure 5c
# Violin plots of genes of interest
trio_plots_median(projenitor_mtec, geneset = c("Hmgb2", "Tubb5", "Stmn1"),
  cell_cycle = FALSE, plot_violin = TRUE, jitter_and_violin = FALSE,
  plot_jitter = FALSE, color = timecourse_color, sep_by = "pub_exp",
  save_plot = paste0(save_dir, "/figure_5cI.pdf"))

trio_plots_median(projenitor_mtec, geneset = c("Aire", "Ccl21a", "Fezf2"),
  cell_cycle = FALSE, plot_violin = TRUE, jitter_and_violin = FALSE,
  plot_jitter = FALSE, color = timecourse_color, sep_by = "pub_exp",
  save_plot = paste0(save_dir, "/figure_5cII.pdf"))


# Figure 5d
# UMAP of reanalysis of early aire cells
tSNE_PCA(no_at_mtec, "cluster", color = reanalysis_colors,
  save_plot = paste0(save_dir, "/figure_5dI.pdf"))

tSNE_PCA(no_at_mtec, "cluster", color = reanalysis_colors, show_legend = FALSE,
  save_plot = paste0(save_dir, "/figure_5d.pdf"))

# Figure 5e
# Violin plots of genes of interest
trio_plots_median(no_at_mtec, geneset = c("Hmgb2", "Tubb5", "Stmn1"), cell_cycle = FALSE,
  plot_violin = TRUE, jitter_and_violin = FALSE, plot_jitter = FALSE,
  sep_by = "cluster", color = reanalysis_colors,
  save_plot = paste0(save_dir, "/figure_5eI.pdf"))

trio_plots_median(no_at_mtec, geneset = c("Aire", "Ccl21a", "Fezf2"), cell_cycle = FALSE,
  plot_violin = TRUE, jitter_and_violin = FALSE, plot_jitter = FALSE,
  sep_by = "cluster", color = reanalysis_colors,
  save_plot = paste0(save_dir, "/figure_5eII.pdf"))

# Figure 5f
# Percent of cycling cells in each new cluster
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

# Figure 5g
# Population plot of new clusters over time.
stage_list_all <- lapply(data_sets, function(x) populations_dfs_new(no_at_mtec,
                         x, subsample = TRUE, subsample_by = "pub_exp",
                         meta_data_col = "res.0.6"))
stage_df_all <- do.call("rbind", stage_list_all)
stage_df_all$sample <- factor(stage_df_all$sample, levels = unname(new_exp_names))


population_plots(stage_df_all, color = reanalysis_colors,
  save_plot = paste0(save_dir, "/figure_5g.pdf"))



fig_list <- c(fig_list, "figure_5")



############################################################################

########################
# Supplemental Figures #
########################

#########################
# Supplemental Figure 1 #
#########################
print("supplemental_figure_1")

# S1a
# Heatmap of all TFs with interesting TFs highlighted
pdf(paste0(save_dir, "/figure_s1a.pdf"))

heatmap_mtec <- mtec_no_un

heatmap_mtec@assay$DE <- mtec@assay$DE

plot_heatmap(heatmap_mtec, subset_list = TFs_all,
  color_list = c("Cdx1", "Utf1", "Tcf7", "Spib", "Cdk4", "Ptma",
    "H2afx", "Hmgb1"),
  color_list2 = c("Aire", "Irf7", "Cited2", "Spib", "Hes1", "Pax1",
                "Relb", "Lmo4", "Pou2f3"),
  cell_color = stage_color)
dev.off()

# S1b
# Marker genes on UMAP
tSNE_PCA(mtec, "Ackr4", save_plot = paste0(save_dir, "/figure_s1bI.pdf"))
tSNE_PCA(mtec, "Psmb11", save_plot = paste0(save_dir, "/figure_s1bII.pdf"))
tSNE_PCA(mtec, "Ccl21a", save_plot = paste0(save_dir, "/figure_s1bIII.pdf"))
tSNE_PCA(mtec, "Ascl1", save_plot = paste0(save_dir, "/figure_s1bIV.pdf"))
tSNE_PCA(mtec, "Hmgb2", save_plot = paste0(save_dir, "/figure_s1bV.pdf"))
tSNE_PCA(mtec, "Stmn1", save_plot = paste0(save_dir, "/figure_s1bVI.pdf"))
tSNE_PCA(mtec, "Fezf2", save_plot = paste0(save_dir, "/figure_s1bVII.pdf"))
tSNE_PCA(mtec, "Aire", save_plot = paste0(save_dir, "/figure_s1bVIII.pdf"))
#tSNE_PCA(mtec, "Tnfrsf11a", save_plot = paste0(save_dir, "/figure_s1bIX.pdf"))
tSNE_PCA(mtec, "GFP", save_plot = paste0(save_dir, "/figure_s1bIX.pdf"))
tSNE_PCA(mtec, "Spink5", save_plot = paste0(save_dir, "/figure_s1bX.pdf"))
tSNE_PCA(mtec, "Cldn4", save_plot = paste0(save_dir, "/figure_s1bXI.pdf"))
tSNE_PCA(mtec, "Trpm5", save_plot = paste0(save_dir, "/figure_s1bXII.pdf"))
tSNE_PCA(mtec, "Dclk1", save_plot = paste0(save_dir, "/figure_s1bXIII.pdf"))
tSNE_PCA(mtec, "res.0.6", save_plot = paste0(save_dir, "/figure_s1bXIV.pdf"))

fig_list <- c(fig_list, "supplemental_figure_1")

#S1d Heatmap of previous markers
# pdf(paste0(save_dir, "/figure_s1d.pdf"))

# previous_marker_list <- c("Bmp4", "Kitl", "Il7", "Atg5", "Fgfr2", "Ackr4",
#   "Pax9", "Pax1", "Dll4", "Cxcl12", "Ly75", "Psmb11", "Ccl25", "Map3k14",
#   "Prss16", "Wnt4", "Cbx4", "Cd83", "Foxn1", "Hoxa3", "Ltbr", "Ctsl", "Enpep",
#   "Trp63", "Kremen1", "Rela", "Tbata", "Egfr", "Sirt1", "H2-Aa", "Cd74", "H2-Eb1",
#   "Relb", "Traf6", "Krt17", "Krt14", "Cd40", "Ccl19", "Ascl1", "Skint1", "Ctss",
#   "Aire", "Tnfrsf11a", "Fgf21", "Fezf2", "Spib", "Ehf", "Ctsz", "Ctsh", "Hdac3",
#   "Krt5", "Cldn4", "Hey1", "Cldn3", "Plet1", "Six1", "Epcam", "Chuk", "Bmi1",
#   "Six4", "Eya1", "Sox4", "H2-ab1", "Sbsn", "Trpm5", "Gng13", "Ly6a", "Spink5",
#   "Hpgds", "Avil", "Pex1", "Ccl21a", "Ivl", "Il25", "Dclk1", "Pdpn", "Krt10",
#   "Cd80", "Ptsg1", "Il17rb")

# plot_heatmap(mtec_no_un, subset_list = previous_marker_list,
#   cell_color = stage_color)

# dev.off()

# Figure s1c
comparison_list <- names(mtec@assay$DE)
clusters <- lapply(comparison_list, get_slots)
clusters <- unique(unlist(clusters))
all_clusters <- sapply(clusters, cluster_gene_list,
                       cluster_list = comparison_list, USE.NAMES = TRUE)

all_clusters_but_ea <- all_clusters
all_clusters_but_ea$Early_Aire <- NULL
ea_genes <- all_clusters$Early_Aire


genes_no_ea <- unique(unlist(all_clusters_but_ea))
genes_ea_overlap <- intersect(ea_genes, genes_no_ea)

genes_ea_unique <- setdiff(ea_genes, genes_ea_overlap)

print(genes_ea_unique)

cell_cycle_genes <- c("Mif", "Birc5", "Cks2", "Cdk1", "Cdk4", "Cdc20",
                      "Ranbp1", "Npm1", "Ccnd2", "Ube2c", "Erh", "Cks1b",
                      "Top2a", "Tubb2b", "H2afx", "Ran", "Tuba1b",
                      "Stmn1", "Tubb5", "Ube2s")

chromosome_organization <- c("Ptma", "Hmgb2", "Cdk1", "Cdc20", "Npm1",
                             "H2afx", "Set", "Ran", "H2afz", "Top2a")

regulation_of_proliferation <- c("Mif", "Hmgb2", "Cdk1", "Eif5a", "Cdk4",
                                 "Cdc20", "Npm1", "Ccnd2", "Vim")

regulation_of_differentiation <- c("Mif", "Tubb2b", "Hmgb2", "Hes6", "Eif5a",
                                   "Cdc20", "Ranbp1", "Vim")

gene_hits <- unique(c(cell_cycle_genes, chromosome_organization,
                      regulation_of_proliferation,
                      regulation_of_differentiation))

no_hit <- setdiff(genes_ea_unique, gene_hits)

# Make a heatmap
ea_gene_df <- as.data.frame(genes_ea_unique)
rownames(ea_gene_df) <- genes_ea_unique
ea_gene_df$cell_cycle <- 0
ea_gene_df$cell_cycle[rownames(ea_gene_df) %in% cell_cycle_genes] <- 1
ea_gene_df$chromosome_organization <- 0
ea_gene_df$chromosome_organization[rownames(ea_gene_df) %in% chromosome_organization] <- 1
ea_gene_df$regulation_of_proliferation <- 0
ea_gene_df$regulation_of_proliferation[rownames(ea_gene_df) %in% regulation_of_proliferation] <- 1
ea_gene_df$regulation_of_differentiation <- 0
ea_gene_df$regulation_of_differentiation[rownames(ea_gene_df) %in% regulation_of_differentiation] <- 1
ea_gene_df$genes_ea_unique <- NULL

pdf(paste0(save_dir, "/figure_s1c.pdf"))

ea_gene_df$gene <- rownames(ea_gene_df)
ea_gene_df_m <- melt(ea_gene_df)
ea_gene_df_m[ea_gene_df_m$value == 0, ]$value <- "gene not in GO term"
ea_gene_df_m[ea_gene_df_m$value == 1, ]$value <- "gene in GO term"
ea_gene_df_m$value <- factor(ea_gene_df_m$value)


ggplot(ea_gene_df_m, aes(x = variable, y = gene, fill = value)) +
  geom_tile(color = "white") +
  coord_equal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("black", "white")) +
  xlab("GO term")


dev.off()


# Figure s1e
# Jitter plots of chromatin modifiers overlayed with cell cycle state
trio_plots_median(mtec_no_un, geneset = c("Pou5f1", "Sox2", "Krt15"),
  cell_cycle = FALSE, plot_violin = TRUE, plot_jitter = FALSE,
  color = stage_color, sep_by = "stage",
  save_plot = paste0(save_dir, "/figure_s1e.pdf"))


# Figure_s1_extra
# Heatmap of all early Aire genes
pdf(paste0(save_dir, "/Figure_s1_extra.pdf"), width = 8, height = 15)

plot_heatmap(heatmap_mtec, subset_list = ea_genes,
  cell_color = stage_color)

dev.off()

######################### 
# Supplemental Figure 2 #
#########################
print("supplemental_figure_2")

mtec_wt <- Seurat::SetAllIdent(mtec_wt, id = "stage")

# S2a
# Violin plots of marker genes for WT
trio_plots_median(mtec_wt_plot, geneset = c("Ackr4", "Ccl21a", "Aire"),
  cell_cycle = FALSE, plot_violin = TRUE, jitter_and_violin = FALSE,
  plot_jitter = FALSE, color = stage_color, sep_by = "stage",
  save_plot = paste0(save_dir, "/figure_s2aI.pdf"))
trio_plots_median(mtec_wt_plot, geneset = c("Krt10", "Trpm5", "Ascl1"),
  cell_cycle = FALSE, plot_violin = TRUE, jitter_and_violin = FALSE,
  plot_jitter = FALSE, color = stage_color, sep_by = "stage",
  save_plot = paste0(save_dir, "/figure_s2aII.pdf"))

# S2b
# Dotplot of WT
# markers_to_plot_full <- c("Krt5", "Ccl21a", "Ascl1", "Hes1", "Hmgb2", "Hmgn2",
#   "Hmgb1", "H2afx", "Stmn1", "Tubb5", "Mki67", "Ptma", "Aire", "Utf1", "Fezf2", 
#   "Krt10", "Nupr1", "Cebpb", "Trpm5", "Pou2f3", "Dclk1")
stage_levels <- c("cTEC", "Ccl21a_high", "Early_mTEC",
                  "Aire_positive", "Late_Aire", "Tuft")

mtec_wt_plot@ident <- factor(mtec_wt_plot@ident,
   levels = stage_levels)
pdf(paste0(save_dir, "/figure_s2bI.pdf"))
dot_plot <- Seurat::DotPlot(mtec_wt_plot, genes.plot = rev(markers_to_plot_full),
                            cols.use = c("blue", "red"), x.lab.rot = T,
                            plot.legend = F, dot.scale = 8, do.return = T)


dev.off()

pdf(paste0(save_dir, "/figure_s2bII.pdf"))
dot_plot <- Seurat::DotPlot(mtec_wt_plot, genes.plot = rev(markers_to_plot_full),
                            cols.use = c("blue", "red"), x.lab.rot = T,
                            plot.legend = T, dot.scale = 8, do.return = T)
dev.off()

# S2c
# Jitter plots of cycling with markers from AT
trio_plots(mtec_wt_plot, geneset = c("Hmgb2", "Tubb5", "Stmn1"),
  cell_cycle = TRUE, jitter_and_violin = FALSE, plot_jitter = TRUE,
  color = stage_color, sep_by = "stage",
  save_plot = paste0(save_dir, "/figure_s2cI.pdf"), group_color = FALSE)
trio_plots(mtec_wt_plot, geneset = c("Aire", "Ccl21a", "Fezf2"),
  cell_cycle = TRUE, jitter_and_violin = FALSE, plot_jitter = TRUE,
  color = stage_color, sep_by = "stage",
  save_plot = paste0(save_dir, "/figure_s2cII.pdf"), group_color = FALSE)

# S2d
# Correlation of WT with Aire Trace
pdf(paste0(save_dir, "/figure_s2d.pdf"))
master_plot(mtec, "aire_trace", mtec_wt, "wt", stage_color_df)
dev.off()

# S2f
# number of early mtecs

# S2g
# Percent Ki67 positive

fig_list <- c(fig_list, "supplemental_figure_2")
#########################
# Supplemental Figure 3 #
#########################
print("supplemental_figure_3")
# S3a
# Umap of stage recovery
plot_names <- list(isoControlBeg = paste0(save_dir, "/figure_s3aI.pdf"),
  isoControlEnd = paste0(save_dir, "/figure_s3aII.pdf"),
  timepoint1 = paste0(save_dir, "/figure_s3aIII.pdf"),
  timepoint2 = paste0(save_dir, "/figure_s3aIV.pdf"),
  timepoint3 = paste0(save_dir, "/figure_s3aV.pdf"),
  timepoint5 = paste0(save_dir, "/figure_s3aVI.pdf"))

lapply(names(plot_names), function(x) full_umap(mtecCombined,
  data_set = x, col_by = "stage", color = stage_color,
  save_plot = plot_names[[x]], show_legend = FALSE))

# S3b
# Violin markers for all as in S2
trio_plots_median(mtecCombined, geneset = c("Ackr4", "Ccl21a", "Aire"),
  cell_cycle = FALSE, plot_violin = TRUE, jitter_and_violin = FALSE,
  plot_jitter = FALSE, color = stage_color, sep_by = "stage",
  save_plot = paste0(save_dir, "/figure_s3bI.pdf"))
trio_plots_median(mtecCombined, geneset = c("Krt10", "Trpm5", "GFP"),
  cell_cycle = FALSE, plot_violin = TRUE, jitter_and_violin = FALSE,
  plot_jitter = FALSE, color = stage_color, sep_by = "stage",
  save_plot = paste0(save_dir, "/figure_s3bII.pdf"))

# S3c Correlation of TP5 with WT
mtecCombExp <- Seurat::SetAllIdent(mtecCombSub, id = "exp")
mtec_tp5 <- Seurat::SubsetData(mtecCombExp, ident.use = "timepoint5")
mtec_end <- Seurat::SubsetData(mtecCombExp, ident.use = "isoControlEnd")
pdf(paste0(save_dir, "/figure_s3c.pdf"))
master_plot(mtec_tp5, "wk_10", mtec_end, "ctl_wk_10", stage_color_df)
dev.off()

# S3d expression of other Ccl21 genes (Ccl21c is not in any gene list)
trio_plots_median(mtecCombined, geneset = c("Ccl21a", "Ccl21a.1", "Ccl21b.1"),
  cell_cycle = FALSE, plot_violin = TRUE, jitter_and_violin = FALSE,
  plot_jitter = FALSE, color = stage_color, sep_by = "pub_exp",
  save_plot = paste0(save_dir, "/figure_s3d.pdf"))



fig_list <- c(fig_list, "supplemental_figure_3")
#########################
# Supplemental Figure 4 #
#########################
print("supplemental_figure_4")

# S4a
# Marker genes on UMAP
tSNE_PCA(mtecCombined_with_un, "Ackr4", save_plot = paste0(save_dir, "/figure_s4aI.pdf"))
tSNE_PCA(mtecCombined_with_un, "Psmb11", save_plot = paste0(save_dir, "/figure_s4aII.pdf"))
tSNE_PCA(mtecCombined_with_un, "Ccl21a", save_plot = paste0(save_dir, "/figure_s4aIII.pdf"))
tSNE_PCA(mtecCombined_with_un, "Krt5", save_plot = paste0(save_dir, "/figure_s4aIV.pdf"))
tSNE_PCA(mtecCombined_with_un, "Krt8", save_plot = paste0(save_dir, "/figure_s4aV.pdf"))
tSNE_PCA(mtecCombined_with_un, "Ascl1", save_plot = paste0(save_dir, "/figure_s4aVI.pdf"))
tSNE_PCA(mtecCombined_with_un, "Fezf2", save_plot = paste0(save_dir, "/figure_s4aVII.pdf"))
tSNE_PCA(mtecCombined_with_un, "Aire", save_plot = paste0(save_dir, "/figure_s4aVIII.pdf"))
tSNE_PCA(mtecCombined_with_un, "Tnfrsf11a", save_plot = paste0(save_dir, "/figure_s4aIX.pdf"))
tSNE_PCA(mtecCombined_with_un, "Spink5", save_plot = paste0(save_dir, "/figure_s4aX.pdf"))
tSNE_PCA(mtecCombined_with_un, "Trpm5", save_plot = paste0(save_dir, "/figure_s4aXI.pdf"))
tSNE_PCA(mtecCombined_with_un, "Dclk1", save_plot = paste0(save_dir, "/figure_s4aXII.pdf"))
tSNE_PCA(mtecCombined_with_un, "Pou2f3", save_plot = paste0(save_dir, "/figure_s4aXIII.pdf"))
tSNE_PCA(mtecCombined_with_un, "GFP", save_plot = paste0(save_dir, "/figure_s4aXIV.pdf"))
tSNE_PCA(mtecCombined_with_un, "res.0.6", save_plot = paste0(save_dir, "/figure_s4XV.pdf"))

plot_marker_heatmap <- function(mtec, subset_val, gene_df, subset_by = "exp", save_plot = NULL){
  if (!(is.null(save_plot))){
    extension <- substr(save_plot, nchar(save_plot)-2, nchar(save_plot))
    if (extension == "pdf"){
      pdf(save_plot)
    } else if (extension == "png") {
      png(save_plot)
    } else {
      print("save plot must be .png or .pdf")
    }
  }
  print(save_plot)
  print(subset_val)
  mtec <- Seurat::SetAllIdent(mtec, id = subset_by)
  mtec_sub <- Seurat::SubsetData(mtec, ident.use = subset_val)
  mtec_sub <- Seurat::SetAllIdent(mtec_sub, id = "stage")
  print(Seurat::DoHeatmap(object = mtec_sub, genes.use = gene_df$gene, slim.col.label = TRUE,
    remove.key = TRUE, group.label.rot = TRUE))
  #heatmap <- heatmap + theme(axis.text.x = element_text(angle = 45))
  if (!(is.null(save_plot))){
    dev.off()
  }
}

heatmap_names <- list(isoControlBeg = paste0(save_dir, "/figure_s4bI.pdf"),
  isoControlEnd = paste0(save_dir, "/figure_s4bII.pdf"),
  timepoint1 = paste0(save_dir, "/figure_s4bIII.pdf"),
  timepoint2 = paste0(save_dir, "/figure_s4bIV.pdf"),
  timepoint3 = paste0(save_dir, "/figure_s4bV.pdf"),
  timepoint5 = paste0(save_dir, "/figure_s4bVI.pdf"))

# Heatmaps here
mtec_wt_plot <- mtec_wt
mtec_wt_plot@assay$DE <- NULL
mtec_wt_plot <- Seurat::SubsetData(mtec_wt_plot, ident.remove = "unknown")
mtec.markers <- Seurat::FindAllMarkers(object = mtec_wt_plot, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
top10 <- mtec.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
print(head(top10))

lapply(names(plot_names), function(x) plot_marker_heatmap(mtecCombSub, subset_val = x,
  gene_df <- top10, save_plot = heatmap_names[[x]]))

fig_list <- c(fig_list, "supplemental_figure_4")
#########################
# Supplemental Figure 5 #
#########################
print("supplemental_figure_5")

# S5a 
# TRA recovery UMAPs
# TRAs in recovery
limit_list <- list(tra_fantom = c(0, 0.100),
                   aire_genes = c(0, 0.100),
                   fezf2_genes = c(0, 0.300))

plot_names_fig5 <- list(isoControlBeg = paste0(save_dir, "/figure_s5aI.pdf"),
  isoControlEnd = paste0(save_dir, "/figure_s5aII.pdf"),
  timepoint1 = paste0(save_dir, "/figure_s5aIII.pdf"),
  timepoint2 = paste0(save_dir, "/figure_s5aIV.pdf"),
  timepoint3 = paste0(save_dir, "/figure_s5aV.pdf"),
  timepoint5 = paste0(save_dir, "/figure_s5aVI.pdf"))

names(plot_names_fig5) <- new_exp_names[names(plot_names_fig5)]

lapply(names(plot_names_fig5), function(x) plot_gene_set(mtecCombined,
                                            gene_set = gene_lists[["tra_fantom"]],
                                            plot_name = "tra_fantom",
                                            one_dataset = FALSE,
                                            data_set = x,
                                            meta_data_col = "pub_exp",
                                            limits = limit_list[["tra_fantom"]],
                                            save_plot = plot_names_fig5[[x]]))



# S5b
# Number of protein coding genes seen in WT not downsampled
# Double check this is correct
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

ggplot2::ggsave(paste0(save_dir, "/figure_s5b.pdf"), plot = cumFreqPlot)

# S5c
# Bootstrap downsample plots
if (bootstrap) {
  source(bootstrap_script)

}

fig_list <- c(fig_list, "supplemental_figure_5")
#########################
# Supplemental Figure 6 #
#########################


print("supplemental_figure_6")

# S6a
# nGene and nUMI in all cells before correction
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

ggplot2::ggsave(paste0(save_dir, "/figure_s6aI.pdf"), plot = plot_all_cell_count)

plot_all_umi_count <- ggplot2::ggplot(mtec_meta_data,
  ggplot2::aes(x = pub_exp, y = nUMI, group = pub_exp,
               fill = pub_exp)) + 
  ggplot2::geom_violin(scale = "width") +
  ggplot2::scale_fill_manual(values = at_timecourse) + 
  ggplot2::ylab("nUMI")

ggplot2::ggsave(paste0(save_dir, "/figure_s6aII.pdf"), plot = plot_all_umi_count)


# S6b
# nGene and nUMI just Aire Positive before correction

aire_positive_meta <- mtec_meta_data[mtec_meta_data$stage == "Aire_positive", ]

plot_all_cell_count_ap <- ggplot2::ggplot(aire_positive_meta,
  ggplot2::aes(x = pub_exp, y = nGene, group = pub_exp,
               fill = pub_exp)) + 
  ggplot2::geom_violin(scale = "width") +
  ggplot2::scale_fill_manual(values = at_timecourse) + 
  ggplot2::ylab("nGene") +
  ggplot2::xlab("experiment")

ggplot2::ggsave(paste0(save_dir, "/figure_s6bI.pdf"), plot = plot_all_cell_count_ap)

plot_all_umi_count_ap <- ggplot2::ggplot(aire_positive_meta,
  ggplot2::aes(x = pub_exp, y = nUMI, group = pub_exp,
               fill = pub_exp)) + 
  ggplot2::geom_violin(scale = "width") +
  ggplot2::scale_fill_manual(values = at_timecourse, name = "experiment") + 
  ggplot2::ylab("nUMI") +
  ggplot2::xlab("experiment")

ggplot2::ggsave(paste0(save_dir, "/figure_s6bII.pdf"), plot = plot_all_umi_count_ap)

# S6c
# nGene and nUMI just aire positive after correction
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

ggplot2::ggsave(paste0(save_dir, "/figure_s6cI.pdf"), plot = plot_corrected_gene)


plot_corrected_umi <- ggplot2::ggplot(umis_df_all,
  ggplot2::aes(x = pub_exp, y = gene_count, 
               group = pub_exp, fill = pub_exp)) +
  ggplot2::geom_violin(scale = "width") +
  ggplot2::scale_fill_manual(values = timecourse_color) + 
  ggplot2::ylab("nUMI") +
  ggplot2::xlab("experiment")

ggplot2::ggsave(paste0(save_dir, "/figure_s6cII.pdf"), plot = plot_corrected_umi)

# S6d
# Dropouts of different house keeping genes before and after
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

ggplot2::ggsave(paste0(save_dir, "/figure_s6dI.pdf"), plot = dropout_plot)

dropout_downsample_m$Var2 <- factor(dropout_downsample_m$Var2,
                              levels = unname(new_exp_names))

dropout_downsample_plot <- ggplot2::ggplot(dropout_downsample_m,
  ggplot2::aes(x = Var1, y = value, fill = Var2)) +
  ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge()) +
  ggplot2::scale_fill_manual(values = timecourse_color, name = "experiment") +
  ggplot2::xlab("housekeeping gene") +
  ggplot2::ylab("dropout percent")

ggplot2::ggsave(paste0(save_dir, "/figure_s6dII.pdf"), plot = dropout_downsample_plot)

# Heatmaps here
mtec_wt_genes <- mtec_wt
mtec_wt_genes@assay$DE <- NULL
mtec_wt_genes <- Seurat::SetAllIdent(mtec_wt_genes, id = "stage")
mtec_wt_genes <- Seurat::SubsetData(mtec_wt_genes, ident.remove = "unknown")
mtec.markers <- Seurat::FindAllMarkers(object = mtec_wt_genes, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
top30 <- mtec.markers %>% group_by(cluster) %>% top_n(30, avg_logFC)
top30 <- as.data.frame(top30)

aire_positive_genes <- top30[top30$cluster == "Aire_positive", ]$gene

aire_positive_genes <- c(aire_positive_genes, "Fezf2")

aire_positive_matrix <- as.matrix(no_at_mtec_aire@data)
heatmap_aire_mtec <- Seurat::SetAllIdent(no_at_mtec_aire, id = "exp")

plot_heatmap_new <- function(mtec, cell_color = NULL, subset_list = NULL,
  color_list = NULL, color_list2 = NULL, order_cells = TRUE,
  seed = 0){
  mtec_data <- mtec@data
  
  
  # Subset the list if desired (ie by a list of specific genes)
  if (!is.null(subset_list)){
    mtec_data <- mtec_data[rownames(mtec_data) %in% subset_list, ]
  }
  mtec_data <- as.matrix(mtec_data)
  
  # Center values to plot on heatmap
  mtec_data_heatmap <- t(scale(t(mtec_data), scale = FALSE))
  cluster <- as.data.frame(mtec@ident)
  names(cluster) <- "cluster_val"
  
  # Order cells by cluster
  if (order_cells){
    cluster <- cluster[order(cluster$cluster_val), , drop=FALSE]
    mtec_data_heatmap <- mtec_data_heatmap[, match(rownames(cluster),
                                                   colnames(mtec_data_heatmap))]
  }

  colors <- as.numeric(cluster$cluster_val)
  if (!is.null(cell_color)) {
    col1 <- cell_color
  } else {
    col1 <- RColorBrewer::brewer.pal(length(levels(cluster$cluster_val)), "Set1")
  }
  
  cols <- rep("black", nrow(mtec_data_heatmap))
  
  # Color some text red if desired
  if (!is.null(color_list)){
    cols[row.names(mtec_data_heatmap) %in% color_list] <- "red"
  }
  if (!is.null(color_list2)){
    cols[row.names(mtec_data_heatmap) %in% color_list2] <- "blue"
  }
  
  sep_list <- lapply(1:length(unique(colors)), function(x) grep(x, colors)[1])
  sep_list <- unlist(sep_list)

  # Seed for reporducibility
  set.seed(seed)
  gplots::heatmap.2(mtec_data_heatmap,
                    density.info  = "none",
                    labCol        = FALSE,
                    Colv          = !order_cells,
                    colRow        = cols,
                    ColSideColors = col1[colors],
                    colsep        = sep_list,
                    sepcolor      = "white",
                    trace         = "none",
                    col           = grDevices::colorRampPalette(c("blue", "yellow")),
                    dendrogram    = "row")


}

pdf(paste0(save_dir, "/figure_s6e.pdf"))

print(plot_heatmap_new(heatmap_aire_mtec, cell_color = timecourse_color,
  subset_list = aire_positive_genes))

dev.off()

fig_list <- c(fig_list, "supplemental_figure_6")

write.table(fig_list, save_file)
