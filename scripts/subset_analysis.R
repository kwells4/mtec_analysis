# This script is used to analyze one population of cells in all samples combined
library(mTEC.10x.pipeline)
library(dplyr)
library(Seurat)

#############
# Functions #
#############

list_seurat_obj <- function(seurat_obj, name, project_name){
  seurat_obj <- Seurat::SetAllIdent(seurat_obj, id = "exp")
  seurat_obj <- Seurat::SubsetData(seurat_obj, ident.use = name, subset.raw = TRUE)
  new_seurat_obj <- Seurat::CreateSeuratObject(raw.data = seurat_obj@raw.data,
                                               min.cells = 3,
                                               min.genes = 200,
                                               project = project_name)
  new_seurat_obj <- add_perc_mito(new_seurat_obj)

  # The default for PCA is 20 PCs. This is impossible with so few cells.
  new_seurat_obj <- process_cells(new_seurat_obj, PCA = FALSE)
  # Run PCA with fewer PCs computed
  new_seurat_obj <- Seurat::RunPCA(object = new_seurat_obj,
    pc.genes = new_seurat_obj@var.genes, do.print = TRUE, pcs.print = 1:5,
    genes.print = 5, seed.use = 42, pcs.compute = 5)
  new_seurat_obj@meta.data$exp <- name
  return(new_seurat_obj)
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
input_file <- snakemake@input[[1]]
subset_name <- snakemake@params[[1]]
pdf_file <- snakemake@params[[2]]
output_file <- snakemake@output[[1]]

mtecCombined <- get(load(input_file))

# Subset the Seurat object to be only the cell population for futher analysis
mtecSub <- mtecCombined
mtecSub@assay$ablation_DE <- NULL
mtec_subset <- Seurat::SubsetData(mtecSub, ident.use = subset_name, subset.raw = TRUE)

data_sets <- unique(mtec_subset@meta.data$exp)


pdf(pdf_file)

# Make a list of seurat objects that have been normalized and have a list of hvgs
seurat_obj_list <- sapply(data_sets, function(x) list_seurat_obj(
  mtec_subset, x, subset_name))

# Pull out top 1000 hvgs from all seurat objects
gene_list <- lapply(seurat_obj_list, function(x) top_genes(x))

# Combine all hvgs into one list
gene_list <- unique(unlist(gene_list))

# Ensure all genes used for downstream analysis are in all samples
genes_use <- lapply(seurat_obj_list, function(x) get_genes(gene_list, x))

# Make a list intersecting these genes
genes_use <- Reduce(intersect, genes_use)

print("multi CCA")
# Run multi CCA to determine CCA space Rcude num.ccs to 5 because there are
# so few cells per sample
mtecSubComb <- RunMultiCCA(object.list = seurat_obj_list,
                            genes.use = genes_use, num.ccs = 10)

# Plot the CCA coordinates
DimPlot(object = mtecSubComb, reduction.use = "cca", group.by = "exp", 
        pt.size = 0.5, do.return = TRUE)

# Plot violin plot of each experiment along CC1
VlnPlot(object = mtecSubComb, features.plot = "CC1", group.by = "exp", 
        do.return = TRUE)

# Print genes associated with each dimension
PrintDim(object = mtecSubComb, reduction.type = "cca", dims.print = 1:2, 
         genes.print = 10)

# Test significance for each dim
MetageneBicorPlot(mtecSubComb, grouping.var = "exp", dims.eval = 1:10, 
                  display.progress = FALSE)

# Align based on CCA coordinates
print("align subspace")
mtecSubComb <- AlignSubspace(mtecSubComb, reduction.type = "cca",
                              grouping.var = "exp", 
                              dims.align = 1:10)

# Determine tSNE coordinates based on CCA space (rather than PCA space)
mtecSubComb <- RunTSNE(mtecSubComb, reduction.use = "cca.aligned", dims.use = 1:7, 
                        do.fast = T, seed.use = 0)

# Determine UMAP coordinates based on CCA space (rather than PCA space)
mtecSubComb <- RunUMAP(mtecSubComb, reduction.use = "cca.aligned", dims.use = 1:7,
                        seed.use = 0)

# Determine clusters based on CCA space (rather than PCA space)
mtecSubComb <- FindClusters(mtecSubComb, reduction.type = "cca.aligned", 
                             resolution = 0.6, dims.use = 1:7,
                             random.seed = 0)

save(mtecSubComb, file = output_file)