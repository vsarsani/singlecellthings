import os
import logging
import argparse
import scanpy as sc
import pandas as pd
import numpy as np
import scrublet as scr
import celltypist
from celltypist import models
import scanpy.external as sce
from sklearn.metrics.pairwise import cosine_distances

# Set up logging
logging.basicConfig(filename='sample_preprocessing.log', level=logging.INFO)

# Preprocessing
def preprocess(sample_file, reference_file, sampleprefix):
    sample = sc.read_h5ad(sample_file, chunk_size=1000)
    reference = sc.read_h5ad(reference_file, chunk_size=10000)

    # Extract metadata
    metadata = sample.obs.copy()

    # Define function to calculate IQR bounds
    def calculate_iqr_bounds(series):
        q1 = series.quantile(0.25)
        q3 = series.quantile(0.75)
        iqr = q3 - q1
        lower_bound = q1 - 3 * iqr
        upper_bound = q3 + 3 * iqr
        return lower_bound, upper_bound

    # Calculate IQR bounds
    ncount_lower, ncount_upper = calculate_iqr_bounds(metadata['nCount_RNA'])
    nfeature_lower, nfeature_upper = calculate_iqr_bounds(metadata['nFeature_RNA'])
    mito_lower, mito_upper = calculate_iqr_bounds(metadata['mitoRatio'])

    # Filter cells based on IQR bounds
    filtered_cells = metadata[
        (metadata['nCount_RNA'] >= ncount_lower) & (metadata['nCount_RNA'] <= ncount_upper) &
        (metadata['nFeature_RNA'] >= nfeature_lower) & (metadata['nFeature_RNA'] <= nfeature_upper) &
        (metadata['mitoRatio'] >= mito_lower) & (metadata['mitoRatio'] <= mito_upper)
    ].index

    # Subset the dataset
    sample_filtered = sample[filtered_cells, :].copy()

    logging.info(f"Filtered dataset has {sample_filtered.shape[0]} cells.")

    # Calculate doublet rate
    num_cells = sample_filtered.shape[0]
    doublet_rate = (0.066 + 0.000757 * num_cells) / 100

    # Run Scrublet
    #scrub = scr.Scrublet(sample_filtered.X, expected_doublet_rate=doublet_rate)
    #doublet_scores, predicted_doublets = scrub.scrub_doublets()

    # Add Scrublet results to sample
    #sample_filtered.obs['doublet_scores'] = doublet_scores
    #sample_filtered.obs['predicted_doublets'] = predicted_doublets

    # Handle None values in 'predicted_doublets'
    #if 'predicted_doublets' in sample_filtered.obs:
        #if sample_filtered.obs['predicted_doublets'].isnull().all():
           # logging.warning("Column 'predicted_doublets' is all None, proceeding without filtering.")
        #else:
            #sample_filtered = sample_filtered[~sample_filtered.obs['predicted_doublets']].copy()

    # Filter out predicted doublets
    #sample_filtered = sample_filtered[~sample_filtered.obs['predicted_doublets']].copy()

    # Save filtered dataset
    #sample_filtered.write_h5ad(f"{sampleprefix}_filtered.h5ad")

    logging.info(f"Filtered dataset saved with {sample_filtered.shape[0]} cells.")
    sc.pp.normalize_total(sample_filtered, target_sum=1e4)
    sc.pp.log1p(sample_filtered)
    sc.pp.highly_variable_genes(sample_filtered, min_mean=0.0125, max_mean=3, min_disp=0.25)
    sample_filtered.raw = sample_filtered
    sc.pp.scale(sample_filtered, max_value=10)
    sc.tl.pca(sample_filtered, svd_solver='arpack')
    sc.pp.neighbors(sample_filtered, n_neighbors=10, n_pcs=50)
    sc.tl.umap(sample_filtered)
    sample_filtered = sample_filtered[:, reference.var_names]
    return sample_filtered

# Ingestion
def ingest(sample, reference_file):
    reference = sc.read_h5ad(reference_file, chunk_size=10000)
    sc.tl.ingest(sample, reference, obs=["Class_8", "Subclass"])
    return sample

# Integration and label transfer
def integrate_and_transfer(sample, reference_file):
    logging.info("Reading reference data...")
    reference = sc.read_h5ad(reference_file, chunk_size=10000)
    # Concatenate and integrate
    logging.info("Concatenating and integrating...")
    adata_concat = reference.concatenate(sample, batch_categories=["ref", "new"])    
    sc.tl.pca(adata_concat)
    sc.pp.neighbors(adata_concat, n_pcs=30, use_rep='X_pca')
    sce.pp.harmony_integrate(adata_concat, 'batch')

    # Calculate distances
    logging.info("Calculating distances...")
    distances_harmony = 1 - cosine_distances(
        adata_concat[adata_concat.obs.batch == "ref"].obsm["X_pca_harmony"],
        adata_concat[adata_concat.obs.batch == "new"].obsm["X_pca_harmony"],
    )

    def label_transfer(dist, labels, index):
        lab = pd.get_dummies(labels)
        class_prob = lab.to_numpy().T @ dist
        norm = np.linalg.norm(class_prob, 2, axis=0)
        class_prob = class_prob / norm
        class_prob = (class_prob.T - class_prob.min(1)) / class_prob.ptp(1)
        cp_df = pd.DataFrame(class_prob, columns=lab.columns)
        cp_df.index = index
        m = cp_df.idxmax(axis=1)
        return m

    # Label transfer
    logging.info("Performing label transfer...")
    class_def = label_transfer(distances_harmony, reference.obs.Class_8, sample.obs.index)
    sample.obs['predicted_class'] = class_def
    subclass_def = label_transfer(distances_harmony, reference.obs.Subclass, sample.obs.index)
    sample.obs['predicted_subclass'] = subclass_def
    
    logging.info("Integration and transfer completed.")
    
    return sample

# Run Celltypist
def run_celltypist(sample, pickl_ref1, pickl_ref2):
    celltypist1 = models.Model.load(pickl_ref1)
    celltypist2 = models.Model.load(pickl_ref2)
    result1 = celltypist.annotate(sample, model=celltypist1, majority_voting = True)
    result2 = celltypist.annotate(sample, model=celltypist2, majority_voting = True)
    return result1, result2

# Main function
def main(sample_file, reference_file, pickl_ref1, pickl_ref2, sampleprefix):
    sample_filtered = preprocess(sample_file, reference_file, sampleprefix)
    sample_added_annots = ingest(sample_filtered, reference_file)
    sample_added_annots_integ = integrate_and_transfer(sample_added_annots, reference_file)
    celltypist1, celltypist2 = run_celltypist(sample_added_annots_integ, pickl_ref1, pickl_ref2)
    sample_added_annots_integ.obs['Class_celltypist_label'] = celltypist1.predicted_labels['majority_voting']
    sample_added_annots_integ.obs['Subclass_celltypist_label'] = celltypist2.predicted_labels['majority_voting']

    # Save the annotated AnnData object
    sample_added_annots_integ.write_h5ad(f"{sampleprefix}_annotated.h5ad")

    # Save the obs DataFrame to a CSV file
    sample_added_annots_integ.obs.to_csv(f"{sampleprefix}_annotated.csv")

    # Define the list of columns to plot
    colors = ['Class_8', 'predicted_class', 'Class_celltypist_label', 'Subclass', 'predicted_subclass', 'Subclass_celltypist_label']

    # Plot UMAP with specified colors and layout, and save to PDF
    sc.pl.umap(sample_added_annots_integ, color=colors, wspace=0.5, ncols=3, save=f"{sampleprefix}_plot.pdf")

    logging.info(f"Saved annotated data to {sampleprefix}_annotated.h5ad, {sampleprefix}_annotated.csv, and UMAP plot to {sampleprefix}_plot.pdf.")

if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Single-cell RNA-seq preprocessing, integration, and annotation")
    parser.add_argument("--sample_file", type=str, help="Path to the sample h5ad file")
    parser.add_argument("--reference_file", type=str, help="Path to the reference h5ad file")
    parser.add_argument("--pickl_ref1", type=str, help="Path to the first Celltypist model pickle file")
    parser.add_argument("--pickl_ref2", type=str, help="Path to the second Celltypist model pickle file")
    parser.add_argument("--sampleprefix", type=str, help="Prefix for the output files")

    args = parser.parse_args()

    main(args.sample_file, args.reference_file, args.pickl_ref1, args.pickl_ref2, args.sampleprefix)
