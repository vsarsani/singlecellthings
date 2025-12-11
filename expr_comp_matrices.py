#!/usr/bin/env python3
import os
import sys
import scanpy as sc
import pandas as pd
import decoupler as dc
from sccoda.util import cell_composition_data as dat

sc.settings.n_jobs = 8

# Expect: expr_comp_matrices.py <workdir> <h5ad_path> <cell_type_cluster>
if len(sys.argv) != 4:
    raise SystemExit(
        f"Usage: {sys.argv[0]} <workdir> <h5ad_path> <cell_type_cluster>\n"
        f"Got: {sys.argv[1:]}"
    )

workdir = sys.argv[1]
h5ad_path = sys.argv[2]
cluster_name = sys.argv[3]

# Use the arguments instead of hard-coded local paths
input_files = [h5ad_path]
cluster_names = [cluster_name]

# Process each file
for file, cluster in zip(input_files, cluster_names):
    print(f"Processing {file}...")

    # Load the h5ad file
    adata = sc.read_h5ad(file)

    data_scanpy_1 = dat.from_scanpy(
        adata,
        cell_type_identifier=cluster,
        sample_identifier="Sample"
    )
    # print(data_scanpy_1)

    combined_df_corrected = pd.concat(
        [
            data_scanpy_1.obs,
            pd.DataFrame(
                data_scanpy_1.X,
                index=data_scanpy_1.obs.index,
                columns=data_scanpy_1.var.index,
            ),
        ],
        axis=1,
    )

    if "leiden" in cluster:
        new_column_names = ['cluster' + str(col) for col in data_scanpy_1.var.index]
    else:
        new_column_names = list(data_scanpy_1.var.index)

    # Update the relevant columns in combined_df_corrected with the new names
    combined_df_corrected.columns = list(data_scanpy_1.obs.columns) + new_column_names

    # Pseudobulk aggregation
    pdata = dc.pp.pseudobulk(
        adata,
        sample_col="Sample",  # Use "Sample" column for aggregation
        groups_col="Sample",
        skip_checks=True,
        mode="sum",
        empty=True
    )
    dc.pp.filter_samples(pdata, min_cells=10, min_counts=100)

    # Normalize total counts and log-transform
    sc.pp.normalize_total(pdata, target_sum=1e6)
    sc.pp.log1p(pdata)

    # Filter genes with mean log-transformed expression >= 3
    gene_filter = pdata.X.mean(axis=0) >= 3
    pdata = pdata[:, gene_filter]

    # Scale the data to a maximum value of 10
    sc.pp.scale(pdata.copy(), max_value=10)

    # Save the processed data as a CSV file
    data = pd.DataFrame(pdata.X, index=pdata.obs_names, columns=pdata.var_names)

    output_file1 = os.path.join(
        f"{workdir}/processed_matrices/{os.path.splitext(os.path.basename(file))[0]}_expression_matrix_ds.csv"
    )
    output_file2 = os.path.join(
        f"{workdir}/processed_matrices/{os.path.splitext(os.path.basename(file))[0]}_composition_matrix_ds.csv"
    )

    # Align samples between expression and composition
    common_samples = data.index.intersection(combined_df_corrected.index)
    data_aligned = data.loc[common_samples]
    combined_df_corrected_aligned = combined_df_corrected.loc[common_samples]

    data_aligned.to_csv(output_file1)
    combined_df_corrected_aligned.to_csv(output_file2)

    print(f"Saved {output_file1}")
    print(f"Saved {output_file2}")

print("Processing completed.")
