library(mTEC.10x.pipeline)
library(gplots)
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

pdf("/home/kwells4/mTEC_dev/mtec_snakemake/aireTrace/analysis_outs/gene_list.pdf")
heatmap.2(as.matrix(ea_gene_df), Rowv = FALSE, Colv=FALSE, col=c("white", "black"),
  trace = "none", density.info = "none", dendrogram = "none", cexRow = 0.5)
dev.off()

# Make venn diagram
#area1 <- length(genes_ea_unique)
#area2 <- length(cell_cycle_genes)
#area3 <- length(chromosome_organization)
#area4 <- length(regulation_of_proliferation)
#area5 <- length(regulation_of_differentiation)

