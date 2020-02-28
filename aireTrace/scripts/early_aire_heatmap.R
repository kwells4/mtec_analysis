library(mTEC.10x.pipeline)
library(gplots)
library(ggplot2)
library(reshape)
load("/home/kwells4/mTEC_dev/mtec_snakemake/aireTrace/analysis_outs/seurat_aireTrace.rda")

get_slots <- function(comparison){
  return(strsplit(comparison, "v(?=[A-Z])", perl = TRUE))
}

get_genes <- function(comparison, cluster){
  split_comparison <- strsplit(comparison, "v(?=[A-Z])", perl = TRUE)
  if(split_comparison[[1]][1] == cluster){
    gene_table <- mtec@assay$DE[[comparison]]
    gene_table <- gene_table[gene_table$avg_logFC > 0, ]
    gene_list <- rownames(gene_table)
  } else if(split_comparison[[1]][2] == cluster) {
    gene_table <- mtec@assay$DE[[comparison]]
    gene_table <- gene_table[gene_table$avg_logFC < 0, ] 
    gene_list <- rownames(gene_table)
  } else {
    gene_list <- NULL
  }
  
  return(gene_list)
}

cluster_gene_list <- function(cluster, cluster_list){
  full_list <- lapply(cluster_list, get_genes, cluster = cluster)
  short_list <- unique(unlist(full_list))
  return(short_list)
}

stage_color_df <- data.frame("cTEC" = "#CC6600", "Ccl21a_high" = "#009933",
                            "Early_mTEC" = "#0066CC", "Aire_positive" = "#660099",
                            "Late_Aire" = "#FF0000", "Tuft" = "#990000",
                            "unknown" = "#FDFBFB")

stage_color <- t(stage_color_df)[ , 1]

comparison_list <- names(mtec@assay$DE)
clusters <- lapply(comparison_list, get_slots)
clusters <- unique(unlist(clusters))
all_clusters <- sapply(clusters, cluster_gene_list,
                       cluster_list = comparison_list, USE.NAMES = TRUE)

ea_genes <- all_clusters$Early_Aire



pdf("/home/kwells4/mTEC_dev/mtec_snakemake/figure_output/Figure_s1c.pdf")

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


heatmap_mtec <- mtec_no_un

heatmap_mtec@assay$DE <- mtec@assay$DE


plot_heatmap(heatmap_mtec, subset_list = TFs_all,
  color_list = c("Cdx1", "Utf1", "Tcf7", "Spib", "Cdk4", "Ptma",
    "H2afx", "Hmgb1"),
  color_list2 = c("Aire", "Irf7", "Cited2", "Spib", "Hes1", "Pax1",
                "Relb", "Lmo4", "Pou2f3"),
  cell_color = stage_color)


dev.off()


