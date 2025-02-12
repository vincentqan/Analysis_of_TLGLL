###Import dependent python packages
import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import anndata
import magic
import scprep
from matplotlib import rcParams
import scanpy as sc
import decoupler as dc
import scirpy as ir
import infercnvpy as cnv
from matplotlib import pyplot as plt, cm as mpl_cm
from cycler import cycler
## RNA velocity
import scvelo as scv
import loompy as lp

scv.logging.print_version()
sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
sc.settings.set_figure_params(dpi_save=300)

#######doublets detection##########
###### use P1_Dx as an example######
adata = sc.read('./P1_Dx.h5ad')
######remove doublet cells detected by scrublet######
counts_matrix =  adata.to_df().astype(pd.SparseDtype(int, fill_value=0)).sparse.to_coo().tocsc()
print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))

scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.076)
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2,
                                                          min_cells=3,
                                                          min_gene_variability_pctl=85,
                                                          n_prin_comps=30)
scrub.plot_histogram();
threshold_set=0.25

scrub.call_doublets(threshold=threshold_set)
scrub.plot_histogram();

adata.obs["doublet_score"] = doublet_scores
adata.obs["doublet_detected"] = ((doublet_scores > threshold_set) + 0)

###rewrite 'doublets' into meta.data
adata.obs.to_csv("P1_Dx_doublet_annot.csv")
##### subsample data randomly #####
sc.pp.subsample(adata_use, n_obs=100)

###########################################################
#######Integrated-object construct########
adata = sc.read('./adata_all.h5ad')
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
adata.raw = adata
adata = adata[:, adata.var.highly_variable]
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata,n_comps=50)
sc.external.pp.bbknn(adata, batch_key='sample')
sc.tl.umap(adata)
sc.tl.leiden(adata)

###################################
#######Figure1######
##### volocity analysis ########
#run velocyto run10x first
loompy.combine(files, output_filename, key="Accession")
scv.pp.filter_and_normalize(adata_ve, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata_ve)
scv.tl.velocity(adata_ve)
scv.tl.velocity_graph(adata_ve)
scv.pl.velocity_embedding_stream(adata_ve, basis='umap', color="celltypes",frameon=True,legend_loc= 'right margin')
##### integrate celltypes sankey plot #####
sankey.sankey(
    df['CellType_source'], df['integrate_celltypes'],
    leftLabels=leftLabels[::-1],
    rightLabels=rightLabels[::-1],
    aspect=20, colorDict=colorDict,
    fontsize=12, figure_name="RNA_celltypes_to_integrate_celltypes",
    closePlot=False,
)

#############################
#######Figure2##############
##### imputation and infercnv analysis #####
emt_data = adata.to_df('normalized_counts')
emt_data.head()
#Setting the MAGIC operator parameters
magic_op = magic.MAGIC(n_jobs=6)
#run MAGIC
emt_magic = magic_op.fit_transform(emt_data,genes="all_genes")
adata_raw_imputataion.layers["counts_global_magic_count"] = np.array(emt_magic)

# run infercnv analysis
cnv.io.genomic_position_from_gtf(
    '../data/reference/refdata-gex-GRCh38-2020-A/genes/genes.gtf',
    adata,
    gtf_gene_id = 'gene_id',
    adata_gene_id = 'gene_ensemble'
)
cnv.tl.infercnv(
    adata,
    reference_key="cell_type",
    reference_cat=reference_celltypes,
    window_size=100,
    step=1,
    n_jobs=8,
)
cnv.tl.cnv_score(adata, groupby='cells')
sc.pl.violin(adata, ['cnv_score'], groupby='cell_types',rotation=90 )

#############################
#######Figure3##############
##### gene signatures scoring of scRNA-seq #####
exhausted_genes = ['PDCD1','TIGIT','LAG3','HAVCR2','CTLA4']
cytotoxic_genes = ['NKG7','CCL4','CST7','PRF1','GZMA','GZMB','IFNG','CCL3','FGFBP2']
naive_genes = ['CCR7','TCF7','LEF1','SELL']
proliferation_genes = ['MKI67','TYMS','TOP2A']

sc.tl.score_genes(adata, exhausted_genes, ctrl_size=50, gene_pool=None, n_bins=25, score_name='exhausted_genes_score', random_state=0, copy=False, use_raw=True)
sc.tl.score_genes(adata, cytotoxic_genes, ctrl_size=50, gene_pool=None, n_bins=25, score_name='cytotoxic_genes_score', random_state=0, copy=False, use_raw=True)
sc.tl.score_genes(adata, naive_genes, ctrl_size=50, gene_pool=None, n_bins=25, score_name='naive_genes_score', random_state=0, copy=False, use_raw=True)
sc.tl.score_genes(adata, proliferation_genes, ctrl_size=50, gene_pool=None, n_bins=25, score_name='proliferation_genes_score', random_state=0, copy=False, use_raw=True)

#############################
#######Figure5##############
##### Calculate the activity of pathways #####
dc.run_wsum(mat=adata, net=model, source='source', target='target', weight='weight',
            #times=10,
              verbose=True,use_raw=True)

#############################
#######Sfigure5##############
##### corrlationship plot #####
def fun(cluster):
    frame = pd.DataFrame(adata_plot.obsm['X_pca'],
                         index=adata_plot.obs.index).T.astype(float)[
        adata_plot.obs[adata_plot.obs.Cell_types == cluster].cells].sum(axis=1)
    return (frame)

cor_data = pd.Series(list(adata_plot.obs.Cell_types.unique()),
                     index=list(adata_plot.obs.Cell_types.unique())).apply(fun).T.corr(method='pearson',
                                                                                       min_periods=1)
g = sb.clustermap(cor_data, cmap="vlag", yticklabels=True, xticklabels=True,
                  cbar_pos=(0.02, 0.8, 0.03, 0.1), figsize=(5, 5),
                  vmin=-1, vmax=1, center=0,
                  )
g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize=10,)
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize=10, rotation=0 )

#############################
#######Sfigure6##############
##### TCR analysis of scRNA-seq #####
ir.pp.merge_with_ir(adata, adata_tcr)
ir.tl.chain_qc(adata)
adata = adata[adata.obs["chain_pairing"] != "multichain", :].copy()
adata = adata[~adata.obs["chain_pairing"].isin(["orphan VDJ", "orphan VJ"]), :].copy()
ir.pp.ir_dist(
    adata,
    metric="alignment",
    sequence="aa",
    cutoff=15,
)
ir.tl.define_clonotype_clusters(
    adata, sequence="aa", metric="alignment", receptor_arms="all", dual_ir="any"
)
ir.tl.clonotype_network(adata, min_cells=3, sequence="aa", metric="alignment")