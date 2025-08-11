import scvelo as scv
scv.settings.verbosity = 3 
scv.settings.presenter_view = True  
scv.set_figure_params('scvelo')  
import loompy

files = ["/home/smilesw/DKD/dbm_glom_1/dbm_glom_1/velocyto/dbm_glom_1.loom","/home/smilesw/DKD/dbm_glom_2/dbm_glom_2/velocyto/dbm_glom_2.loom",
         "/home/smilesw/DKD/dbm_nonglom_1/dbm_nonglom_1/velocyto/dbm_nonglom_1.loom","/home/smilesw/DKD/dbm_nonglom_2/dbm_nonglom_2/velocyto/dbm_nonglom_2.loom",
         "/home/smilesw/DKD/dbdb_glom_1/dbdb_glom_1/velocyto/dbdb_glom_1.loom","/home/smilesw/DKD/dbdb_glom_2/dbdb_glom_2/velocyto/dbdb_glom_2.loom",
         "/home/smilesw/DKD/dbdb_nonglom_1/dbdb_nonglom_1/velocyto/dbdb_nonglom_1.loom","/home/smilesw/DKD/dbdb_nonglom_2/dbdb_nonglom_2/velocyto/dbdb_nonglom_2.loom"]

loompy.combine(files, output_file= "/home/smilesw/DKD/velocyto_merged.loom")

adata = scv.read('/mnt/d/dbdb_looom/velocyto_merged.loom', cache=True)

scv.pl.proportions(adata)

import pandas as pd
sample_obs = pd.read_csv("/mnt/z/project/mouse_scRNA-seq_diabetic_kidney_cellbender/scVelo/barcode_dbm.csv")
umap_cord = pd.read_csv("/mnt/z/project/mouse_scRNA-seq_diabetic_kidney_cellbender/scVelo/umap_dbm.csv")
cell_clusters = pd.read_csv("/mnt/z/project/mouse_scRNA-seq_diabetic_kidney_cellbender/scVelo/cluster_dbm.csv")

import numpy as np
#Transfer seurat info to vlm
dict_group_to_id = {'dbm_glom_1':"dbm_glom_1", 'dbm_glom_2':"dbm_glom_2", 'dbm_nonglom_1':"dbm_nonglom_1", 'dbm_nonglom_2':"dbm_nonglom_2",
                   'dbdb_glom_1':"dbdb_glom_1", 'dbdb_glom_2':"dbdb_glom_2", 'dbdb_nonglom_2':"dbdb_nonglom_2"}

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

import scanpy 

scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)

scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

scv.tl.velocity(adata)

scv.tl.velocity_graph(adata)

adata.obsm['X_umap'] = np.array(umap_ordered)
scv.pl.velocity_embedding_stream(adata, basis='umap', color='Clusters', palette = ('#00B0F6','#CDB38B','#9932CC','#FF8C00'), dpi=600)

#Dynamical Modeling
scv.tl.recover_dynamics(adata)

scv.tl.velocity(adata, mode='dynamical')
scv.tl.velocity_graph(adata)

scv.pl.velocity_embedding_stream(adata, basis='umap', color='Clusters', arrow_size=2.5, palette = ('#00B0F6','#228B22','#CDB38B','#9932CC','#FF6347','#FF8C00'), dpi=600, legend_loc='none', save="a.png")

df = adata.var
df = df[(df['fit_likelihood'] > .1) & df['velocity_genes'] == True]

kwargs = dict(xscale='log', fontsize=16)
with scv.GridSpec(ncols=3) as pl:
    pl.hist(df['fit_alpha'], xlabel='transcription rate', **kwargs)
    pl.hist(df['fit_beta'] * df['fit_scaling'], xlabel='splicing rate', xticks=[.1, .4, 1], **kwargs)
    pl.hist(df['fit_gamma'], xlabel='degradation rate', xticks=[.1, .4, 1], **kwargs)

scv.get_df(adata, 'fit*', dropna=True).head()

scv.tl.latent_time(adata)
scv.pl.scatter(adata, color='latent_time', color_map='gnuplot', size=80, save="latent_dbm.png")

scv.tl.terminal_states(adata)

scv.pl.scatter(adata, color=['root_cells', 'end_points'])

var_names = ['Srgap1', 'Thsd7a', 'Robo2', 'Dpp4','Dach1','Magi2']
scv.pl.scatter(adata, var_names, frameon=False,  color='Clusters', palette = ('#00B0F6','#228B22','#CDB38B','#9932CC','#FF6347','#FF8C00'))
scv.pl.scatter(adata, x='latent_time', y=var_names,  color='Clusters', frameon=False)
