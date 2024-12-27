import scanpy as sc
import numpy as np
import pandas as pd
import random

def qc_metrics(adata, path_RB=''):
    mito_genes = adata.var_names.str.startswith('MT-')
    adata.var['mito'] = mito_genes
    mito_genes = np.array(adata.var.index)[mito_genes]
    
    with open(path_RB,'r') as file:
        lines = file.readlines()
    RB_genes = [x.rstrip('\n') for x in lines]
    data_genes = list(adata.var.index)
    RB_genes_in_data = set(data_genes).intersection(RB_genes)
    RB_genes_in_data = list(RB_genes_in_data)
    
    adata.var['ribo'] = False
    adata.var.loc[RB_genes_in_data, 'ribo'] = True
    
    sc.pp.calculate_qc_metrics(adata, qc_vars=('mito', 'ribo'), inplace=True, layer='raw_counts')
    
    return adata

def preprocessing(adata, n_pcs = 45, n_hvgs = 4500):
    # Normalize
    adata.layers['raw_counts'] = adata.X
    adata.layers['median'] = adata.layers['raw_counts'].copy()
    sc.pp.normalize_total(adata, layer='median')
    adata.layers['log'] = adata.layers['median'].copy()
    sc.pp.log1p(adata, layer='log')
    
    # HVG
    sc.pp.highly_variable_genes(adata, layer='raw_counts', flavor='seurat_v3', max_mean=np.inf, n_top_genes=n_hvgs)
    
    # PCA
    adata.X = adata.layers['log']
    sc.tl.pca(adata, n_comps=n_pcs, use_highly_variable=True)
    adata.obsm['X_pca'] = adata.obsm['X_pca'].copy()
    
    # UMAP
    sc.pp.neighbors(adata, method='umap', n_neighbors = 30, use_rep='X_pca')
    sc.tl.umap(adata, min_dist = 0.1, random_state=5)
    
    # Phenograph
    k = 30
    communities, _, _ = sc.external.tl.phenograph(
        pd.DataFrame(adata.obsm['X_pca']), k=k,
        nn_method = 'brute', njobs = -1,
    )
    adata.obs['PhenoGraph_clusters'] = pd.Categorical(communities)

    # randomly shuffle indexes
    barcodes = list(adata.obs_names)
    random.shuffle(barcodes)
    adata = adata[barcodes, :]
    return adata