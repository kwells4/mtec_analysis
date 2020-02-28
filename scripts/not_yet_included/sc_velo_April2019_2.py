import scvelo as scv
import matplotlib
matplotlib.use('agg')
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


scv.settings.set_figure_params('scvelo')

working_dir = "/home/kwells4/mTEC_dev/mtec_snakemake/"

seurat_cells = working_dir + "not_yet_included/seurat_info_new.csv"
output = working_dir + "/figure_output/figure_2d_new2.pdf"
input_loom = working_dir + "not_yet_included/wt_velocyto.loom"
seurat_df = pd.read_csv(seurat_cells, index_col = 0)
supplement_output = working_dir +"/figure_output/supplemental_scvelo_oct_2019.pdf"

seurat_cell_list = list(seurat_df.index)
adata = scv.read(input_loom, sparse = True, cache = True)

adata.var_names_make_unique()
adata = adata[adata.obs.index.isin(seurat_cell_list)]
scv.utils.show_proportions(adata)
scv.utils.cleanup(adata, clean='all')



scv.pp.filter_and_normalize(adata, min_counts=20, min_counts_u=10, n_top_genes=3000)

scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

scv.tl.velocity(adata)

scv.tl.velocity_graph(adata)

new_seurat = seurat_df.reindex(adata.obs.index)

# adata.obs = pd.merge(adata.obs, seurat_df, right_index = True, left_index = True, how = "left")
adata.obs["clusters"] = new_seurat["clusters"]
umap_coord = list(zip(new_seurat["UMAP1"], new_seurat["UMAP2"]))
umap_coord = np.asarray(umap_coord)

# umap_coord = list(zip(adata.obs["UMAP1"], adata.obs["UMAP2"]))
# umap_coord = np.asarray(umap_coord)
adata.obsm["X_umap"] = umap_coord

scv.tl.velocity_embedding(adata, basis='umap')


# change #666666 to #ffffff as a cheat
with PdfPages(output) as pdf:
	adata.uns["clusters_colors"] = ("#660099", "#009933", "#CC6600", "#0066CC", "#FF0000", "#990000", "#FFFFFF")
	plt.figure(None,(8,8))
	scv.pl.velocity_embedding_stream(adata, legend_loc = "none", arrow_color = "black")
	pdf.savefig()
	plt.close()

with PdfPages(supplement_output) as pdf:
	plt.figure(None,(8,8))
	scv.pl.velocity_embedding(adata, basis='umap', dpi=200)
	pdf.savefig()
	plt.close()
	