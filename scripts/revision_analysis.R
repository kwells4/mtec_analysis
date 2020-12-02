library(mTEC.10x.pipeline)
library(dplyr)
library(slingshot)
#library(svglite)
library(ggplot2)
library(gplots)
library(reshape)
source("/home/kwells4/mTEC_dev/mtec_snakemake/scripts/figure_funcs.R")
library(Seurat)
library(pheatmap)
library(ggforce)

process_cells_var <- function (mtec, filter_by = c("nGene", "percent_mito"),
	low_thresholds = c(200, -Inf), high_thresholds = c(7500, 0.1),
	scale_factor = 10000, filter_cells = TRUE, normalize_data = TRUE,
	find_var_genes = TRUE, scale_data = TRUE, PCA = TRUE,
	var_genes = NULL){
    if (filter_cells) {
        mtec <- Seurat::FilterCells(object = mtec, subset.names = filter_by, 
            low.thresholds = low_thresholds, high.thresholds = high_thresholds)
    }
    if (normalize_data) {
        mtec <- Seurat::NormalizeData(object = mtec, normalization.method = "LogNormalize", 
            scale.factor = scale_factor)
    }
    if (find_var_genes) {
        mtec <- Seurat::FindVariableGenes(object = mtec, mean.function = ExpMean, 
            dispersion.function = LogVMR, x.low.cutoff = 0.0125, 
            x.high.cutoff = 3, y.cutoff = 0.5)
    }
    if (scale_data) {
        mtec <- Seurat::ScaleData(object = mtec)
    }
    if (PCA) {
    	if(is.null(var_genes)){
    		var_genes = mtec@var.genes
    	}
        mtec <- Seurat::RunPCA(object = mtec, pc.genes = var_genes, 
            do.print = TRUE, pcs.print = 1:5, genes.print = 5, 
            seed.use = 42)
        mtec <- Seurat::ProjectPCA(object = mtec, do.print = FALSE)
    }
    return(mtec)
}

group_cells_new <- function (mtec, dims_use = 1:10, random_seed = 0, resolution = 0.6, 
    cluster = TRUE, tSNE = TRUE, UMAP = TRUE) 
{
    if (cluster) {
        mtec <- Seurat::FindClusters(object = mtec, reduction.type = "pca", 
            dims.use = dims_use, resolution = resolution, print.output = 0, 
            save.SNN = TRUE, random.seed = random_seed, force.recalc = TRUE)
    }
    if (tSNE) {
        mtec <- Seurat::RunTSNE(object = mtec, dims.use = dims_use, 
            do.fast = TRUE, seed.use = random_seed)
    }
    if (UMAP) {
        mtec <- Seurat::RunUMAP(object = mtec, dims.use = dims_use, 
            reduction.use = "pca", sed.usee = random_seed)
    }
    return(mtec)
}

trio_plots_new <- function (seurat_object, geneset, cell_cycle = FALSE, plot_jitter = TRUE, 
    plot_violin = FALSE, jitter_and_violin = FALSE, color = NULL, 
    sep_by = "cluster", save_plot = NULL, nrow = NULL, ncol = NULL, 
    group_color = TRUE, plot_name = NULL) 
{
    gene_list_stage <- c()
    if (plot_jitter) {
        if (group_color) {
            for (gene in geneset) {
                gene_stage <- jitter_plot(seurat_object, gene, 
                  sep_by, color = color)
                gene_list_stage[[gene]] <- gene_stage
            }
            full_plot <- gridExtra::grid.arrange(grobs = gene_list_stage, 
                nrow = length(geneset), top = grid::textGrob(plot_name, 
                  gp = grid::gpar(fontsize = 20, font = 3)))
        }
        if (cell_cycle) {
            gene_list_cycle <- c()
            for (gene in geneset) {
                gene_cycle <- jitter_plot(seurat_object, gene, 
                  "stage", "cycle_phase", color = c("black", 
                    "#FF8C00", "#4169E1"))
                gene_list_cycle[[gene]] <- gene_cycle
            }
            full_plot <- gridExtra::grid.arrange(grobs = gene_list_cycle, 
                nrow = length(geneset))
        }
    }
    if (plot_violin || jitter_and_violin) {
        for (gene in geneset) {
            gene_stage <- violin_plot(seurat_object, gene, sep_by, 
                color = color, plot_jitter = jitter_and_violin)
            gene_list_stage[[gene]] <- gene_stage
        }
        if (is.null(nrow)) {
            nrow <- length(geneset)
        }
        if (is.null(ncol)) {
            ncol <- 1
        }
        full_plot <- gridExtra::grid.arrange(grobs = gene_list_stage, 
            nrow = nrow, ncol = ncol, top = grid::textGrob(plot_name, 
                gp = grid::gpar(fontsize = 20, font = 3)))
    }
    if (!(is.null(save_plot))) {
        ggplot2::ggsave(save_plot, plot = full_plot)
    }
    return(full_plot)
}


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
revision_dir <- "revision_plot/"

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

load(paste0(data_directory, "gene_lists.rda"))


stage_levels <- c("cTEC", "Ccl21a_high", "TAC_TEC",
                  "Aire_positive", "Late_Aire", "Tuft", "other")

stage_levels_wt <- c("cTEC", "Ccl21a_high", "TAC_TEC",
                  "Aire_positive", "Late_Aire", "Tuft", "unknown")

# rename cell types in objects

# WT
mtec_wt@meta.data$stage <- as.character(mtec_wt@meta.data$stage)

mtec_wt@meta.data$stage[mtec_wt@meta.data$stage ==
  "Cortico_medullary"] <- "cTEC"
mtec_wt@meta.data$stage[mtec_wt@meta.data$stage ==
  "Early_Aire"] <- "TAC_TEC"

mtec_wt@meta.data$stage <- factor(mtec_wt@meta.data$stage, levels = stage_levels_wt)

mtec_wt_plot <- mtec_wt
mtec_wt_plot@assay$DE <- NULL

mtec_wt_plot <- Seurat::SetAllIdent(mtec_wt_plot, id = "stage")
mtec_wt_plot <- Seurat::SubsetData(mtec_wt_plot, ident.remove = "unknown")

mtec_wt_plot@ident <- factor(mtec_wt_plot@ident, levels = stage_levels)


# Aire Trace 

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

# Combined cells 

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

# First investigate if TSAs are in the pseudotime analysis
var_genes <- mtec_wt@var.genes

tras <- gene_lists$tra_fantom

# Find common genes
overlap <- intersect(var_genes, tras)

# 348/3298 genes overlap


# dim.scores <- GetDimReduction(object = object, reduction.type = reduction.type, 
#         slot = slot.use)
# PrintDim(object = object, dims.print = pcs.print, genes.print = genes.print, 
#             reduction.type = reduction.name)

# Repeat analysis removing them
new_var_genes <- var_genes[!(var_genes %in% tras)]
mtec_wt_repeat <- mtec_wt
mtec_wt_repeat@assay$DE <- NULL

mtec_wt_repeat <- process_cells_var(mtec_wt_repeat, var_genes = new_var_genes)
mtec_wt_repeat@dr$umap <- NULL
mtec_wt_repeat <- group_cells_new(mtec_wt_repeat, dims_use = 1:13)

pdf(paste0(revision_dir, "new_wt.pdf"))
plotDimRed(mtec_wt_repeat, "stage", color = stage_color)
dev.off()

mtec_wt_repeat <- Seurat::SetAllIdent(mtec_wt_repeat, "stage")
mtec_wt_repeat <- Seurat::SubsetData(mtec_wt_repeat, ident.remove = "unknown")


# Pull out dataframes for slingshot to use UMAP embeddings
cluster_labels <- mtec_wt_repeat@meta.data$stage
slingshot_umap <- mtec_wt_repeat@dr$umap@cell.embeddings

sce.umap <- slingshot(slingshot_umap, clusterLabels = cluster_labels,
	start.clus = "TAC_TEC", extend = "n")

plot_names <- list(protein_coding = paste0(revision_dir, "/protein_coding.pdf"),
  tra_fantom = paste0(revision_dir, "/tra_fantom.pdf"),
  aire_genes = paste0(revision_dir, "/aire_genes.pdf"),
  Aire = paste0(revision_dir, "/aire.pdf"),
  fezf2_genes = paste0(revision_dir, "/fezf2_genes.pdf"),
  Fezf2 = paste0(revision_dir, "/fezf2.pdf"))


for (gene_set in names(gene_lists)) {
  mtec_wt_repeat <- plot_gene_set(mtec_wt_repeat,
    gene_lists[[gene_set]], gene_set, make_plot = FALSE)
}



# mtec_wt_plot@meta.data$stage <- factor(mtec_wt_plot@meta.data$stage,
#   levels = stage_levels)

# Plot each of the genes and gene sets in pseudotime, end at 16 because the
# is where the "unknown" cells are
plot_list <- lapply(names(plot_names), function(x) plot_sling_pseudotime(
  seurat_object = mtec_wt_repeat, sling_object = sce.umap, y_val = x,
  col_by = "stage", pseudotime_curve = "curve1", color = stage_color,
  range = c(0, 16), save_plot = plot_names[[x]]))

# Determine the # of TRAs expressed per cell
tras_plot <- tras[tras %in% rownames(mtec_wt_plot@data)]
tra_df <- data.frame(FetchData(object = mtec_wt_plot, vars.all = tras_plot))

number_tras_cell <- rowSums(tra_df > 0)

percent_tras_cell <- number_tras_cell/length(tras)

mtec_wt_plot@meta.data$tra_percent <- percent_tras_cell

plot_sling_pseudotime( seurat_object = mtec_wt_plot, sling_object = wt_slingshot,
	y_val = "tra_percent", col_by = "stage", pseudotime_curve = "curve1",
	color = stage_color, range = c(0, 16),
	save_plot = paste0(revision_dir, "/tra_percent.pdf"))

total_tra_expression <- rowSums(tra_df)
avg_tra_expression_expressed <- total_tra_expression/number_tras_cell

mtec_wt_plot@meta.data$avg_tra_expression_expressed <- avg_tra_expression_expressed

plot_sling_pseudotime( seurat_object = mtec_wt_plot, sling_object = wt_slingshot,
	y_val = "avg_tra_expression_expressed", col_by = "stage", pseudotime_curve = "curve1",
	color = stage_color, range = c(0, 16),
	save_plot = paste0(revision_dir, "/avg_tra_expression_expressed.pdf"))

# Make heatmap of other markers suggested by reviewers.
heatmap_list <- c("Tnfrsf11a", #RANK
			   "Bmp4",
			   "Kitl",
			   "Il7",
			   "Fgfr2",
			   "Ackr4",
			   "Pax1",
			   "Cxcl12",
			   "Ly75",
			   "Psmb11",
			   "Ccl25",
			   "Prss16",
			   "Foxn1",
			   "Ltbr",
			   "Ctsl",
			   "Enpep",
			   "Trp63",
			   "Rela",
			   "Egfr",
			   "Sirt1",
			   "H2-Aa",
			   "Cd74",
			   "H2-Eb1",
			   "Relb",
			   "Traf6",
			   "Krt17",
			   "Krt14",
			   "Cd40",
			   "Ccl19",
			   "Ascl1",
			   "Skint1",
			   "Aire",
			   "Fezf2",
			   "Spib",
			   "Krt5",
			   "Cldn4",
			   "Hey1",
			   "Cldn3",
			   "Epcam",
			   "Ccl21a", # All above genes from the Maehr paper
			   "Itgb4",
			   "Itga6",
			   "Sox4",
			   "Pdpn",
			   "Cd80", # From reviewers
         "Dclk1",
         "Trpm5",
         "Pou2f3",
         "Alox5",
         "Il25",
         "Ivl",
         "Lor",
         "Krt10",
         "Krt1",
         "Cd86",
         #"Cd24", # not in the data set
         "Ska1",
         "Dsc2" # From Mark
			   )

# Find average expression of the genes
stage_levels_new <- stage_levels[1:6]
mtec_wt_plot@ident <- factor(mtec_wt_plot@ident, levels = stage_levels_new)

# Find average expression
gene_expression <- AverageExpression(mtec_wt_plot, genes.use = heatmap_list)

# Scale
gene_expression <- t(scale(t(gene_expression)))

gene_expression <- MinMax(gene_expression, min = -2.5, max = 2.5)

column_color <- data.frame(stage = names(stage_color))
rownames(column_color) <- column_color$stage

column_color_code <- list(stage = stage_color[1:6])

# I am taking a color palette from ArchR
blueYellow <- c("#352A86", "#343DAE", "#0262E0", "#1389D2", "#2DB7A3",
	"#A5BE6A", "#F8BA43", "#F6DA23", "#F8FA0D") 
palOut <- colorRampPalette(blueYellow)(256)

# Plot heatmap
pdf(paste0(revision_dir, "control_curated_genes2.pdf"),
	height = 12, width = 8)

pheatmap(gene_expression, annotation_colors = column_color_code,
	annotation_col = column_color, color = palOut, cluster_cols = FALSE,
  cellwidth = 12, cellheight = 8)

dev.off()

# Further analysis of the TAC-TEC population

# First plot zsGreen expression here.
reanalysis_colors <- c("0" = "#603E95", "1" = "#009DA1", "2" = "#FAC22B")

progenitor_mtec@meta.data$pub_exp <- new_exp_names[progenitor_mtec@meta.data$exp]
progenitor_mtec@meta.data$pub_exp <- factor(progenitor_mtec@meta.data$pub_exp,
                                         levels = unname(new_exp_names))

cells_use <- rownames(progenitor_mtec@meta.data)[progenitor_mtec@meta.data$exp ==
                                                 "aireTrace"]

at_progenitor <- Seurat::SubsetData(progenitor_mtec, cells.use = cells_use)

pdf(paste0(revision_dir, "GFP_in_progenitor.pdf"), height = 4, width = 8)
plot(trio_plots_median(at_progenitor, geneset = c("GFP"),
  cell_cycle = FALSE, plot_violin = TRUE, jitter_and_violin = FALSE,
  plot_jitter = FALSE, color = reanalysis_colors, sep_by = "cluster"))


dev.off()

# Make a bimodal plot of Aire and zsGreen in the progenitor
aire_df <- data.frame(FetchData(object = at_progenitor,
  vars.all = c("Aire", "GFP", "res.0.6")))

pdf(paste0(revision_dir, "bipotent_GFP_in_progenitor.pdf"))
ggplot(aire_df, aes(x = Aire, y = GFP, color = res.0.6)) +
  geom_point(size = 2) +
  scale_color_manual(values = reanalysis_colors) +
  labs(x = "Aire expression", y = "zsGreen expression", color = "cluster")

dev.off()

# Then repeat differential gene expression and make a heatmap
progenitor_DEG <- FindAllMarkers(progenitor_mtec, logfc.threshold = 0.5)

progenitor_DEG <- progenitor_DEG[progenitor_DEG$p_val_adj < 0.05, ]



# Find average expression
prog_gene_expression <- AverageExpression(progenitor_mtec,
  genes.use = unique(progenitor_DEG$gene))

# Scale
prog_gene_expression <- t(scale(t(prog_gene_expression)))

prog_gene_expression <- MinMax(prog_gene_expression, min = -2.5, max = 2.5)

column_color_prog <- data.frame(cluster = c("0", "1", "2"))
rownames(column_color_prog) <- column_color_prog$cluster

column_color_code_prog <- list(cluster = reanalysis_colors)

# I am taking a color palette from ArchR
blueYellow <- c("#352A86", "#343DAE", "#0262E0", "#1389D2", "#2DB7A3",
  "#A5BE6A", "#F8BA43", "#F6DA23", "#F8FA0D") 
palOut <- colorRampPalette(blueYellow)(256)

heatmap_1 <- pheatmap(prog_gene_expression, annotation_colors = column_color_code_prog,
  annotation_col = column_color_prog, color = palOut, cluster_cols = FALSE,
  cellwidth = 45, cellheight = 8)


# Plot heatmap
pdf(paste0(revision_dir, "prog_DE_genes.pdf"), height = 16.5, width = 5)

heatmap_1
dev.off()


# Find average expression
all_gene_expression <- AverageExpression(mtec_wt_plot,
  genes.use = unique(progenitor_DEG$gene))

# Scale
all_gene_expression <- t(scale(t(all_gene_expression)))

all_gene_expression <- MinMax(all_gene_expression, min = -2.5, max = 2.5)

# Find order from the previous heatmap
rowclust <- hclust(dist(prog_gene_expression))
all_gene_expression <- all_gene_expression[rowclust$order,]

# I am taking a color palette from ArchR
blueYellow <- c("#352A86", "#343DAE", "#0262E0", "#1389D2", "#2DB7A3",
  "#A5BE6A", "#F8BA43", "#F6DA23", "#F8FA0D") 
palOut <- colorRampPalette(blueYellow)(256)

# Plot heatmap
heatmap_2 <-  pheatmap(all_gene_expression, annotation_colors = column_color_code,
  annotation_col = column_color, color = palOut, cluster_cols = FALSE,
  cluster_rows = FALSE, cellwidth = 45, cellheight = 8)

pdf(paste0(revision_dir, "prog_DE_genes_all_cells.pdf"), height = 17, width = 7)

heatmap_2

dev.off()
library(gridExtra)
library(grid)
pdf(paste0(revision_dir, "prog_DE_genes_both_plots.pdf"), height = 16, width = 7)
do.call(grid.arrange, list(heatmap_1[[4]], heatmap_2[[4]]))

dev.off()

# Aire trace
# Bipotent plots of Aire trace
at_df <- data.frame(FetchData(object = mtec_no_un,
  vars.all = c("Aire", "GFP", "Ccl21a", "Ackr4", "Trpm5", "Krt10", "Fezf2", "stage")))
stage_color_short <- stage_color[1:6]
pdf(paste0(revision_dir, "bipotent_at_plots.pdf"))
ggplot(at_df, aes(x = Aire, y = GFP, color = stage)) +
  geom_point(size = 2) +
  scale_color_manual(values = stage_color_short) +
  #scale_fill_manual(values = stage_color_short) +
  labs(x = "Aire expression", y = "zsGreen expression", color = "population") 
  #ggforce::geom_mark_ellipse(expand = 0, aes(fill = stage))

ggplot(at_df, aes(x = Ccl21a, y = GFP, color = stage)) +
  geom_point(size = 2) +
  scale_color_manual(values = stage_color_short) +
  #scale_fill_manual(values = stage_color_short) +
  labs(x = "Ccl21a expression", y = "zsGreen expression", color = "population") 
  #geom_mark_ellipse(expand = 0, aes(fill = stage))


ggplot(at_df, aes(x = Ackr4, y = GFP, color = stage)) +
  geom_point(size = 2) +
  scale_color_manual(values = stage_color_short) +
  #scale_fill_manual(values = stage_color_short) +
  labs(x = "Ackr4 expression", y = "zsGreen expression", color = "population") 
  #geom_mark_ellipse(expand = 0, aes(fill = stage))


ggplot(at_df, aes(x = Trpm5, y = GFP, color = stage)) +
  geom_point(size = 2) +
  scale_color_manual(values = stage_color_short) +
  #scale_fill_manual(values = stage_color_short) +
  labs(x = "Trpm5 expression", y = "zsGreen expression", color = "population") 
  #geom_mark_ellipse(expand = 0, aes(fill = stage))


ggplot(at_df, aes(x = Krt10, y = GFP, color = stage)) +
  geom_point(size = 2) +
  scale_color_manual(values = stage_color_short) +
  scale_fill_manual(values = stage_color_short) +
  labs(x = "Krt10 expression", y = "zsGreen expression", color = "population") 
  #geom_mark_ellipse(expand = 0, aes(fill = stage))

ggplot(at_df, aes(x = Fezf2, y = GFP, color = stage)) +
  geom_point(size = 2) +
  scale_color_manual(values = stage_color_short) +
  scale_fill_manual(values = stage_color_short) +
  labs(x = "Fezf2 expression", y = "zsGreen expression", color = "population") 
  #geom_mark_ellipse(expand = 0, aes(fill = stage))

dev.off()


# Table of percent of cells in each category
ccl21_cutoff <- 4
Ackr4 cutoff <- 0
Aire_cutoff <- 0
zsGreen_cutoff <- 0
Trpm5_cutoff <- 0
Fezf2_cutoff <- 0
Krt10_cutoff <- 1
cutoff_list <- c(Aire = 0,
                 Ackr4 = 0,
                 Trpm5 = 0,
                 Fezf2 = 0,
                 Krt10 = 2,
                 Ccl21a = 4)

all_percents <- lapply(names(stage_color[1:6]), function(identity){
  at_df_sub <- at_df[at_df$stage == identity, ]
  percents <- lapply(names(cutoff_list), function(x){
    double_positive <- nrow(at_df_sub[at_df_sub[[x]] > cutoff_list[x] &
      at_df_sub$GFP > zsGreen_cutoff, ])/nrow(at_df_sub) * 100
    single_pos1 <- nrow(at_df_sub[at_df_sub[[x]] <= cutoff_list[x] &
      at_df_sub$GFP > zsGreen_cutoff, ])/nrow(at_df_sub) *100
    single_pos2 <- nrow(at_df_sub[at_df_sub[[x]] > cutoff_list[x] &
      at_df_sub$GFP <= zsGreen_cutoff, ])/nrow(at_df_sub) *100
    double_neg <- nrow(at_df_sub[at_df_sub[[x]] <= cutoff_list[x] &
      at_df_sub$GFP <= zsGreen_cutoff, ])/nrow(at_df_sub) * 100
    percent_list <- c(double_positive, single_pos1, single_pos2, double_neg)
    name1 <- paste0(x, "+_zsgreen+")
    name2 <- paste0(x, "-_zsgreen+")
    name3 <- paste0(x, "+_zsgreen-")
    name4 <- paste0(x, "-_zsgreen-")
    names(percent_list) <- c(name1, name2, name3, name4)
    return(percent_list)

    })
  percent_df <- data.frame(unlist(percents))
  colnames(percent_df) <- identity
  return(percent_df)
})

all_df <- do.call(cbind, all_percents)

write.csv(all_df, paste0(revision_dir, "zsgreen_percents.csv"))


# Remake plots for cell cell_cycle
trio_plots_new(mtec_wt_plot, geneset = c("Hmgb2", "Tubb5", "Stmn1"),
  cell_cycle = TRUE, jitter_and_violin = FALSE, plot_jitter = TRUE,
  color = stage_color, sep_by = "cluster",
  save_plot = paste0(save_dir, "/figure_2cI.pdf"), group_color = FALSE)
trio_plots_new(mtec_wt_plot, geneset = c("Aire", "Ccl21a", "Fezf2"),
  cell_cycle = TRUE, jitter_and_violin = FALSE, plot_jitter = TRUE,
  color = stage_color, sep_by = "cluster",
  save_plot = paste0(save_dir, "/figure_2cII.pdf"), group_color = FALSE)


# Supplemental Figure 3b Jitter plots of cycling with markers from AT
trio_plots_new(mtec_no_un, geneset = c("Hmgb2", "Tubb5", "Stmn1"),
  cell_cycle = TRUE, jitter_and_violin = FALSE, plot_jitter = TRUE,
  color = stage_color, sep_by = "stage",
  save_plot = paste0(save_dir, "/figure_s3bI.pdf"), group_color = FALSE)
trio_plots_new(mtec_no_un, geneset = c("Aire", "Ccl21a", "Fezf2"),
  cell_cycle = TRUE, jitter_and_violin = FALSE, plot_jitter = TRUE,
  color = stage_color, sep_by = "stage",
  save_plot = paste0(save_dir, "/figure_s3bII.pdf"), group_color = FALSE)
