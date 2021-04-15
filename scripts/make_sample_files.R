library(Seurat)

aireTrace <- "aireTrace/analysis_outs/seurat_aireTrace.rda"
controls <- "controls/analysis_outs/seurat_controls_merged.rda"
allSamples <- "allSamples/analysis_outs/seurat_allSamples_combined.rda"


# Load in data
mtec <- get(load(aireTrace))

mtec_wt <- readRDS(controls)

mtecCombined <- get(load(allSamples))

mtec_meta <- mtec@meta.data

mtec_wt_meta <- mtec_wt@meta.data

mtecCombined_meta <- mtecCombined@meta.data

write.csv(mtec_meta, "aireTrace/analysis_outs/aire_trace_meta.csv")
write.csv(mtec_wt_meta, "controls/analysis_outs/controls_meta.csv")
write.csv(mtecCombined_meta, "allSamples/analysis_outs/allSamples_meta.csv")
write.table(mtec_meta, "aireTrace/analysis_outs/aire_trace_meta.txt", sep = "\t")
write.table(mtec_wt_meta, "controls/analysis_outs/controls_meta.txt", sep = "\t")
write.table(mtecCombined_meta, "allSamples/analysis_outs/allSamples_meta.txt", sep = "\t")
