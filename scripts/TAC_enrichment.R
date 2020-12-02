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
controls <- "controls/analysis_outs/seurat_controls_merged.rda"
data_directory <- "/home/kwells4/mTEC_dev/data/"
save_dir <- "figure_output"


# Colors for plotting
stage_color_df <- data.frame("cTEC" = "#CC6600", "Ccl21a_high" = "#009933",
                            "TAC_TEC" = "#0066CC", "Aire_positive" = "#660099",
                            "Late_Aire" = "#FF0000", "Tuft" = "#990000",
                            "unknown" = "#FDFBFB")

stage_color <- t(stage_color_df)[ , 1]



mtec_wt <- readRDS(controls)

bootstrap <- FALSE

stage_levels <- c("cTEC", "Ccl21a_high", "TAC_TEC",
                  "Aire_positive", "Late_Aire", "Tuft", "other")

stage_levels_wt <- c("cTEC", "Ccl21a_high", "TAC_TEC",
                  "Aire_positive", "Late_Aire", "Tuft", "unknown")


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



TA_5 <- read.table("/home/kwells4/mTEC_dev/data/neural_transit_amplifying_cluster5.txt")
TA_8 <- read.table("/home/kwells4/mTEC_dev/data/neural_transit_amplifying_cluster8.txt")

names(TA_5) <- "gene"
names(TA_8) <- "gene"
TA_5_gene <- gsub("__.*", "", TA_5$gene)
TA_8_gene <- gsub("__.*", "", TA_8$gene)

TA_genes_all <- unique(c(TA_5_gene, TA_8_gene))

get_unique_genes <- function(cluster){
  all_clusters_but_one <- all_clusters
  all_clusters_but_one[[cluster]] <- NULL
  gene_set <- all_clusters[[cluster]]


  genes_all_but_one <- unique(unlist(all_clusters_but_one))
  genes_overlap <- intersect(gene_set, genes_all_but_one)

  genes_unique <- setdiff(gene_set, genes_overlap)
  print(cluster)
  print(length(genes_unique))
  print(table(genes_unique %in% TA_genes_all))
  print(genes_unique[genes_unique %in% TA_genes_all])
  return(genes_unique)
}

unique_clusters <- sapply(clusters, function(x) 
  get_unique_genes(x))

text_color <- ifelse(genes_ea_unique %in% TA_genes_all, "red", "black")