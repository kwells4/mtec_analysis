library(Seurat)

sample_dir <- snakemake@input[[1]]
project_name <- snakemake@params[[1]]
data_path <- snakemake@params[[2]]
output_obj <- snakemake@output[[1]]

sequence_dir <- paste0(sample_dir, data_path)
print(sequence_dir)

mtec_data <- Read10X(data.dir = sequence_dir)
mtec <- CreateSeuratObject(raw.data = mtec_data, min.cells = 3, min.genes = 200,
                           project = project_name)

save(mtec, file = output_obj)
