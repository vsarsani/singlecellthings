import os
import argparse
import scanpy as sc
import pandas as pd
import numpy as np
import scanpy.external as sce
import anndata as ad

def main(args):
    sc.settings.n_jobs = 64
    
    # Read the input h5ad file
    adata = sc.read_h5ad(args.input_file)
    
    # Add the new category 'Mathys' to the 'Cohort' column
    adata.obs['Cohort'] = adata.obs['Cohort'].cat.add_categories('Mathys')
    adata.obs['Cohort'].fillna('Mathys', inplace=True)
    
    # Remove MT, ribosomal, and hemoglobin genes
    mito_genes = adata.var_names.str.startswith('MT-')
    ribo_genes = adata.var_names.str.startswith(("RPS", "RPL"))
    hb_genes = adata.var_names.str.contains("^HB[^(P)]")
    remove = np.logical_or.reduce((mito_genes, ribo_genes, hb_genes))
    keep = np.invert(remove)
    adata = adata[:, keep]
    
    # Normalize, log transform, and find highly variable genes
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=args.n_top_genes, batch_key=args.batch_key)
    sc.pp.scale(adata)
    
    # Perform PCA and neighbors
    sc.tl.pca(adata, svd_solver='arpack', n_comps=args.n_comps)
    sc.pp.neighbors(adata, n_neighbors=args.n_neighbors, n_pcs=args.n_comps)
    
    # Harmony integration
    sce.pp.harmony_integrate(adata, ['Cohort', 'Sample'], verbose=1, max_iter_harmony=50)
    adata.obsm['X_pca'] = adata.obsm['X_pca_harmony']
    sc.pp.neighbors(adata, n_neighbors=args.n_neighbors, n_pcs=args.n_comps)
    
    # UMAP and Leiden clustering
    sc.tl.umap(adata)
    
    resolutions = [float(r) for r in args.resolutions.split(',')]
    for res in resolutions:
        sc.tl.leiden(adata, resolution=res, key_added=f'leiden_{res:.2f}')
    
    # Create a trimmed down h5ad file
    adata_new = sc.read_h5ad(args.input_file)
    adata_w = ad.AnnData(adata_new.X)
    adata_w.obs_names = adata_new.obs_names
    adata_w.var_names = adata_new.var_names
    
    columns_to_copy = [
        'Sample', 'nCount_RNA', 'nFeature_RNA', 'mitoRatio',
        'old.predicted.Subclass', 'old.predicted.Supertype', 'Cohort'
    ]
    columns_to_copy.extend([f'leiden_{res:.2f}' for res in resolutions])
    
    for column in columns_to_copy:
        adata_w.obs[column] = adata.obs[column].copy()
    
    adata_w.obsm['X_umap'] = adata.obsm['X_umap']
    
    # Write the output h5ad file
    adata_w.write_h5ad(args.output_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process single-cell RNA-seq data.')
    parser.add_argument('--input_file', type=str, required=True, help='Input h5ad file')
    parser.add_argument('--batch_key', type=str, required=True, help='Batch key for highly variable genes')
    parser.add_argument('--n_top_genes', type=int, required=True, help='Number of top genes')
    parser.add_argument('--n_neighbors', type=int, required=True, help='Number of neighbors for UMAP and Leiden')
    parser.add_argument('--n_comps', type=int, required=True, help='Number of principal components')
    parser.add_argument('--resolutions', type=str, required=True, help='Comma-separated Leiden clustering resolutions')
    parser.add_argument('--output_file', type=str, required=True, help='Output h5ad file')
    
    args = parser.parse_args()
    main(args)
