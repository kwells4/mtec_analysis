


##########
# Add in all the steps to get to the mtec_combined from all figures

mtecCombined_all

mtec_aire_positive <- Seurat::SubsetData(mtecCombSub, ident.use = "Aire_positive",
  subset.raw = TRUE)

week_4 <- Seurat::SetAllIdent(mtecCombined_all, id = "exp")

week_4 <- Seurat::SubsetData(week_4, ident.use = "timepoint2")

counts <- week_4@data

counts <- counts[rownames(counts) %in% gene_lists$tra_fantom,]

counts <- counts[rowSums(data.frame(counts)) > 0, ]

cell_count <- apply(counts, 1, function(x) sum(x>0))

expr_genes <- names(cell_count[cell_count > 0])

mtecCombined_all <- plot_gene_set(mtecCombined_all,
    expr_genes, "expressed_tras", make_plot = FALSE)

gene_sets <- c("expressed_tras", "tra_fantom", "aire_genes")

trio_plots_median(mtecCombined_all, geneset = gene_sets,
  cell_cycle = FALSE, plot_violin = TRUE, jitter_and_violin = FALSE,
  plot_jitter = FALSE, color = timecourse_color, sep_by = "pub_exp",
  save_plot = "allSamples/analysis_outs/expressed_tras.pdf")