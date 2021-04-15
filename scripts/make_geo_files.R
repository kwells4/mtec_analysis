library(mTEC.10x.pipeline)

file_paths <- "/home/kwells4/mTEC_dev/mtec_snakemake/"

save_path <- "/home/kwells4/mTEC_dev/geo_files/"

file_names <- c(aireTrace = "aireTrace",
	            isoCtl_wk2 = "isoControlBeg",
	            isoCtl_wk10 = "isoControlEnd",
	            RANKL_wk2 = "timepoint1",
	            RANKL_wk4 = "timepoint2",
	            RANKL_wk6 = "timepoint3",
	            RANKL_wk10 = "timepoint5",
	            allSamples = "allSamples",
	            combinedControl = "controls")



save_matrix <- function(seurat_name, save_name){
  object_path <- paste0(file_paths, seurat_name, "/analysis_outs/seurat_",
  	seurat_name, ".rda")
  print(object_path)
  seurat_object <- get(load(object_path))
  data_matrix <- as.data.frame(as.matrix(seurat_object@data))
  write.table(data_matrix, file = paste0(save_path, save_name, ".csv"),
  	sep = ",", row.names = TRUE, col.names = TRUE)
}


lapply(names(file_names), function(x) save_matrix(file_names[x], x))