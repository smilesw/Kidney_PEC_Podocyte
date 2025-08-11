import warnings
warnings.filterwarnings('ignore')

import dynamo as dyn

dyn.get_all_dependencies_version()

adata=dyn.read_loom('/mnt/d/home/smilesw/project/mouse_scRNA-seq_diabetic_kidney_cellbender/scVelo/velocyto_merged.loom')

dyn.configuration.set_figure_params('dynamo', background='white')

import pandas as pd
sample_obs = pd.read_csv("/mnt/d/home/smilesw/project/mouse_scRNA-seq_diabetic_kidney_cellbender/scVelo/barcode_dbm.csv")
umap_cord = pd.read_csv("/mnt/d/home/smilesw/project/mouse_scRNA-seq_diabetic_kidney_cellbender/scVelo/umap_dbm.csv")
cell_clusters = pd.read_csv("/mnt/d/home/smilesw/project/mouse_scRNA-seq_diabetic_kidney_cellbender/scVelo/cluster_dbm.csv")

import numpy as np
#Transfer seurat info to vlm
dict_group_to_id = {'dbm_glom_1':"dbm_glom_1", 'dbm_glom_2':"dbm_glom_2", 'dbm_nonglom_1':"dbm_nonglom_1", 'dbm_nonglom_2':"dbm_nonglom_2",
                   'dbdb_glom_1':"dbdb_glom_1", 'dbdb_glom_2':"dbdb_glom_2", 'dbdb_nonglom_1':"dbdb_nonglom_1", 'dbdb_nonglom_2':"dbdb_nonglom_2"}

def get_seurat_cell_id(cell):
    lib, barcode = cell.split(':')
    return '%s_%s_1' % (dict_group_to_id[lib], barcode.replace('x', ''))

adata.obs.index = np.array([get_seurat_cell_id(cell) for cell in adata.obs.index])

bc_intersect = set(adata.obs.index) & set(sample_obs['barcode'])
keep_cell = [cell in bc_intersect for cell in adata.obs.index]
adata = adata[keep_cell]

sample_one_index = pd.DataFrame(adata.obs.index)
sample_one = sample_one_index.rename(columns = {0:'Cell ID'})

umap_cord = umap_cord.rename(columns = {'Unnamed: 0':'Cell ID'})
umap_ordered = sample_one.merge(umap_cord, on = "Cell ID")

umap_ordered = umap_ordered.iloc[:,1:]
adata.obs["_X"] = np.array(umap_ordered["UMAP_1"])
adata.obs["_Y"] = np.array(umap_ordered["UMAP_2"])

cell_clusters = cell_clusters.rename(columns = {'Unnamed: 0':'Cell ID'})

clusters_ordered = sample_one.merge(cell_clusters, on = "Cell ID")

adata.obs["Clusters"] = np.array(clusters_ordered["cluster"])

dyn.pp.recipe_monocle(adata)

dyn.tl.dynamics(adata, model='stochastic', cores=3)

dyn.tl.reduceDimension(adata)

dyn.pl.umap(adata, color='Clusters')

adata.obsm["X_umap"]=np.array(umap_ordered)

dyn.pl.umap(adata, color='Clusters')

dyn.tl.gene_wise_confidence(adata, group='Clusters', lineage_dict={'PEC1': ['Podocyte1']})

dyn.tl.cell_velocities(adata, method='pearson', other_kernels_dict={'transform': 'sqrt'})

dyn.tl.cell_wise_confidence(adata)

dyn.tl.confident_cell_velocities(adata, group='Clusters', lineage_dict={'PEC1': ['Podocyte1']},)

dyn.pl.cell_wise_vectors(adata, color=['Clusters'], basis='umap', show_legend='on data', quiver_length=6, quiver_size=6, pointsize=0.1, show_arrowed_spines=False)

dyn.pl.streamline_plot(adata, color=['Clusters'], basis='umap', show_legend='', figsize=(8,6), linewidth=2, show_arrowed_spines=True, save_show_or_return='return')
plt.savefig("streamline_dbm.png")
plt.show()

dyn.vf.VectorField(adata, basis='umap', M=1000, pot_curl_div=True)

dyn.pl.plot_energy(adata, basis='umap')

dyn.pl.topography(adata, basis='umap', background='white', color=['ntr', 'Clusters'], streamline_color='black', show_legend='', frontier=True)

dyn.pl.umap(adata,  color='umap_ddhodge_potential', frontier=True, figsize=(4,3))

import matplotlib.pyplot as plt
fig1, f1_axes = plt.subplots(ncols=2, nrows=2, constrained_layout=True, figsize=(10, 8))
f1_axes
f1_axes[0, 0] = dyn.pl.cell_wise_vectors(adata, color='speed_pca', pointsize=0.5, alpha = 0.7, ax=f1_axes[0, 0], quiver_length=6, quiver_size=6, save_show_or_return='return')
f1_axes[0, 1] = dyn.pl.grid_vectors(adata, color='divergence_pca', ax=f1_axes[0, 1], quiver_length=12, quiver_size=12, save_show_or_return='return')
f1_axes[1, 0] = dyn.pl.streamline_plot(adata, color='acceleration_pca', ax=f1_axes[1, 0], save_show_or_return='return')
f1_axes[1, 1] = dyn.pl.streamline_plot(adata, color='curvature_pca', ax=f1_axes[1, 1], save_show_or_return='return')
plt.savefig("figure.png")
plt.show()
