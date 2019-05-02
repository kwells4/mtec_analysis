# This script is used to analyze all data files together
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

top_genes <- function(seurat_obj){
  top_genes <- head(rownames(seurat_obj@hvg.info), 1000)
  return(top_genes)
}

get_genes <- function(gene_list, seurat_obj){
  gene_use <- intersect(gene_list, rownames(seurat_obj@scale.data))
  return(gene_use)
}

# Grab names from snakemake object
data_sets <- snakemake@input[["data_list"]]
pdf_file <- snakemake@params[[1]]
data_names <- snakemake@params[[2]]
output_file <- snakemake@output[[1]]
print(data_sets)
names(data_sets) <- data_names

print(data_names)
print(data_sets)

print(output_file)

pdf(pdf_file)

# Make a list of seurat objects that have been normalized and have a list of hvgs
seurat_obj_list <- sapply(names(data_sets), function(x) list_seurat_obj(
  data_sets[[x]], x))

# Pull out top 1000 hvgs from all seurat objects
gene_list <- lapply(seurat_obj_list, function(x) top_genes(x))

# Combine all hvgs into one list
gene_list <- unique(unlist(gene_list))

# Ensure all genes used for downstream analysis are in all samples
genes_use <- lapply(seurat_obj_list, function(x) get_genes(gene_list, x))

# Make a list intersecting these genes
genes_use <- Reduce(intersect, genes_use)


print("multi CCA")
# Run multi CCA to determine CCA space
mtecCombined <- RunMultiCCA(object.list = seurat_obj_list, add.cell.ids = data_names,
                            genes.use = genes_use, num.ccs = 30)

# Plot the CCA coordinates
DimPlot(object = mtecCombined, reduction.use = "cca", group.by = "exp", 
        pt.size = 0.5, do.return = TRUE)

# Plot violin plot of each experiment along CC1
VlnPlot(object = mtecCombined, features.plot = "CC1", group.by = "exp", 
        do.return = TRUE)

# Print genes associated with each dimension
PrintDim(object = mtecCombined, reduction.type = "cca", dims.print = 1:2, 
         genes.print = 10)

# Test significance for each dim
MetageneBicorPlot(mtecCombined, grouping.var = "exp", dims.eval = 1:30, 
                  display.progress = FALSE)

print("align subspace")
# Align based on CCA coordinates
mtecCombined <- AlignSubspace(mtecCombined, reduction.type = "cca",
                              grouping.var = "exp", 
                              dims.align = 1:20)

# Determine tSNE coordinates based on CCA space (rather than PCA space)
mtecCombined <- RunTSNE(mtecCombined, reduction.use = "cca.aligned", dims.use = 1:20, 
                        do.fast = T, seed.use = 0)

# Determine UMAP coordinates based on CCA space (rather than PCA space)
mtecCombined <- RunUMAP(mtecCombined, reduction.use = "cca.aligned", dims.use = 1:20,
                        seed.use = 0)

# Determine clusters based on CCA space (rather than PCA space)
mtecCombined <- FindClusters(mtecCombined, reduction.type = "cca.aligned", 
                             resolution = 0.6, dims.use = 1:20,
                             random.seed = 0)

# Plot the clusters
tSNE_PCA(mtecCombined, "cluster")

dev.off()

save(mtecCombined, file = output_file)