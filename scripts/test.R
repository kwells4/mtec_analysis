no_at_mtec_aire <- Seurat::SetAllIdent(no_at_mtec_aire, "exp")

Seurat::FindMarkers(no_at_mtec_aire, ident.1 = "isoControlBeg",
	ident.2 = "isoControlEnd", genes.use = c("Gapdh", "Aire", "Fezf2"),
	logfc.threshold = 0, min.pct = 0)

Seurat::FindMarkers(progenitor_mtec, ident.1 = "isoControlBeg",
	ident.2 = "timepoint3",
	genes.use = c("Aire", "Ccl21a", "Fezf2", "Hmgb2", "Tubb5", "Stmn1"),
	logfc.threshold = 0, min.pct = 0)


