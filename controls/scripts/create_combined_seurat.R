# This script is used to analyze the controls. This only combines two samples
library(Seurat)
library(mTEC.10x.pipeline)

#############
# Functions #
#############

list_seurat_obj <- function(file_path, name){
  print(file_path)
  seurat_obj <- get(load(file_path))
  seurat_obj <- add_perc_mito(seurat_obj)
  seurat_obj <- process_cells(seurat_obj)
  seurat_obj@meta.data$exp <- name
  return(seurat_obj)
}

# Grab names from snakemake object
data_sets <- snakemake@input[["data_list"]]
pdf_file <- snakemake@params[[1]]
data_names <- snakemake@params[[2]]
output_file <- snakemake@output[[1]]
print(data_sets)
names(data_sets) <- data_names

seurat_obj_list <- lapply(names(data_sets), function(x)
  list_seurat_obj(data_sets[[x]], x))

mtec_wt <- Seurat::MergeSeurat(object1 = seurat_obj_list[[1]],
	                           object2 = seurat_obj_list[[2]],
	                           add.cell.id1 = names(data_sets[1]),
	                           add.cell.id2 = names(data_sets[2]))

pdf(pdf_file)

mtec_wt <- add_perc_mito(mtec_wt)

qc_plot(mtec_wt)

mtec_wt <- process_cells(mtec_wt)

PC_plots(mtec_wt, jackstraw = TRUE, test_pcs = 1:20)

# Determine dims use from jackstraw output
mtec_wt <- group_cells(mtec_wt, dims_use = 1:13)

tSNE_PCA(mtec_wt, "exp")
tSNE_PCA(mtec_wt, "cluster")
tSNE_PCA(mtec_wt, "cluster", PCA = TRUE, tSNE = TRUE, UMAP = FALSE)

dev.off()

save(mtec_wt, file = output_file)
