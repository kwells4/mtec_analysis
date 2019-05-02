library(dropseqr)
library(mTEC.10x.pipeline)

run_bootstrap <- function(seurat_obj, compare_type = "exp",
	                      compare_against = "isoControlBeg",
	                      test_samples = c("timepoint1", "timepoint2"),
	                      subset_mtec = FALSE, subset_by = "stage",
	                      subset_val = "Aire_positive",
	                      downsample_UMI = FALSE){
  seurat_obj$subset_by <- paste0(seurat_obj@meta.data[[compare_type]],
  	"_", seurat_obj@meta.data[[subset_by]])
  compare_against_full <- paste0(compare_against, "_", subset_val)
  seurat_obj <- Seurat::SetAllIdent(seurat_obj, id = "subset_by")
  compare_obj <- Seurat::SubsetData(seurat_obj, ident.use = compare_against_full,
  	subset.raw = TRUE)
  if (downsample_UMI)
}

if (run_bootstrap){
  if (downsample_UMI) {
    raw_isoBeg <- subset_matrix(seurat_obj = mtecCombined,
                                start_matrix = mtecCombined@raw.data,
                                subset_by = "stage_exp",
                                subset_val = paste0(mtec_subset, "_isoControlBeg"))

    data_umi_isoBeg <- median(colSums(as.matrix(raw_isoBeg)))

    factor <- lowest_UMI/data_umi_isoBeg

    set.seed(0)
    isoBeg_cell_matrix <- DropletUtils::downsampleMatrix(raw_isoBeg, prop = factor)
   } else {
    isoBeg_cell_matrix <- subset_matrix(seurat_obj = mtecCombined,
                                        start_matrix = mtecCombined@raw.data,
                                        subset_by = "stage_exp",
                                        subset_val = paste0(mtec_subset, "_isoControlBeg"))
   } 
  test_samples <- c("timepoint1", "timepoint2", "timepoint3", "timepoint5", "isoControlEnd")
  boostrap_list <- c("tra_fantom", "all_other_genes", "aire_genes", "fezf2_genes")
  pdf(paste0(working_dir,
  	"combinedAnalysis/images/Test_data_qual/bootstrap_downsample_tra2.pdf"))
  
  cell_matrix <- isoBeg_cell_matrix
  for (sample in test_samples){
    percents <- NULL
    print(sample)
    ncells <- sum(mtecCombined@meta.data$stage_exp ==
                  paste0("AirePositive_", sample))

    for (gene_set in boostrap_list){
      percents <- NULL
      gene_list <- gene_lists[[gene_set]]
      data_percents <- all_percent_list[[gene_set]]
      data_percents <- data.frame(data_percents)
      data_percents <- as.data.frame(t(data_percents))
      print(data_percents)
      subsample_list <- cell_sampler(cell_matrix, ncells,
                                     percent_list_function,
                                     n = n, replace = FALSE,
                                     funArgs = list(gene_list = gene_list))

      z <- (data_percents[[sample]] - mean(subsample_list)) / sd(subsample_list)
      pval <- 2*pnorm(-abs(z))
      pval <- scales::scientific(pval, digits = 3)
      print(pval)
      grob <- grobTree(textGrob(pval, x = 0.85, y = 0.95, hjust = 0,
                                gp = gpar(col = "black", fontsize = 12)))
      subsample_percents <- data.frame(subsample_list)
     
      p <- ggplot2::ggplot(subsample_percents, ggplot2::aes(x = subsample_list)) +
                ggplot2::geom_histogram() +
                ggplot2::theme_classic() +
                ggplot2::geom_vline(data = data_percents,
                                    ggplot2::aes_string(xintercept = sample),
                                    color = "red") +
                ggplot2::ggtitle(paste0(sample, ", n = ", ncells, " ",
                                        gene_set)) +
                ggplot2::annotation_custom(grob)
      print(p)
    }

  }
  dev.off()
}
