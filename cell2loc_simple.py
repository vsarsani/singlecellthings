import numpy as np
from sklearn.neighbors import KernelDensity
from sklearn.metrics import pairwise_distances_argmin_min
from scipy.spatial.distance import cdist
import anndata
from sklearn.model_selection import GridSearchCV

class CorticalThicknessPredictor:
    def __init__(self, adata, bandwidth=0.1):
        
        self.adata = adata
        self.bandwidth = bandwidth  # bandwidth for Gaussian kernel
        
        # Extract the reference and new data for harmony integrated space and other variables
        self.X_ref = adata[adata.obs['batch'] == 'ref'].obsm['X_pca_harmony']
        self.X_new = adata[adata.obs['batch'] == 'new'].obsm['X_pca_harmony']
        
        # Extract cortical thickness ('i') and cell-type clusters ('seurat_clusters') from reference
        self.i_ref = adata[adata.obs['batch'] == 'ref'].obs['i_smooth']
        self.clusters_ref = adata[adata.obs['batch'] == 'ref'].obs['subtype']

        # Extract new cells' anndata object
        self.new_cells_adata = adata[adata.obs['batch'] == 'new', :].copy()

        print("Initialized CorticalThicknessPredictor with reference and new data.")
    
    def predict_i_and_clusters(self):
        """
        Predict the cortical thickness 'i' and clusters for the new data cells
        using the nearest neighbors from the reference data.
        """
        print("Calculating nearest neighbors and predicting 'i' and clusters for new cells...")
        
        # Find nearest neighbors in reference data for each new cell
        nearest_ref_idx, _ = pairwise_distances_argmin_min(self.X_new, self.X_ref)
        
        # Predict 'i' and clusters for new cells based on nearest neighbors
        self.predicted_i_new = self.i_ref.iloc[nearest_ref_idx].values
        self.predicted_clusters_new = self.clusters_ref.iloc[nearest_ref_idx].values
        
        print(f"Predicted 'i' for {len(self.predicted_i_new)} new cells.")
        
        # Add 'predicted_i_nn' to the new cells' anndata object (nearest neighbor)
        self.new_cells_adata.obs['predicted_i_nn'] = self.predicted_i_new

    def estimate_cell_proportions(self):
        """
        Estimate the cell-type proportions for density-based bins of i in the reference dataset.
        """
        print("Estimating cell-type proportions based on density bins of 'i' in the reference data...")
        
        # We will store cell-type proportions for different density "bins" of i
        self.cell_proportions_by_density = {}
        
        # Sort reference 'i' values and create density bins
        sorted_i_ref = np.sort(self.i_ref)
        density_bins = np.percentile(sorted_i_ref, np.linspace(0, 100, 20 + 1))  # Quartiles
        
        # Estimate proportions for each density bin
        for i in range(len(density_bins) - 1):
            bin_start = density_bins[i]
            bin_end = density_bins[i + 1]
            in_bin = (self.i_ref >= bin_start) & (self.i_ref < bin_end)
            
            if in_bin.sum() == 0:
                continue  # Skip if no cells fall in this range
            
            cell_composition = {}
            for cluster in set(self.clusters_ref):
                # Proportion of this cluster in the bin
                cell_composition[cluster] = np.mean(self.clusters_ref[in_bin] == cluster)
            
            # Store the proportions for this density bin
            self.cell_proportions_by_density[(bin_start, bin_end)] = cell_composition
        
        print("Estimated cell-type proportions by density bins of 'i'.")

    def gaussian_kernel(self, dist, bandwidth):
        """
        Gaussian kernel function for smoothing.
        """
        return np.exp(-0.5 * (dist / bandwidth) ** 2)
    
    def gaussian_kernel_smoothing(self, X_train, y_train, X_test):
        """
        Perform Gaussian kernel smoothing to predict 'i' values based on cell-type proportions.
        """
        # Compute the pairwise distances between training and test data
        distances = cdist(X_test, X_train, metric='euclidean')
        
        # Apply Gaussian kernel to distances
        weights = self.gaussian_kernel(distances, self.bandwidth)
        
        # Normalize the weights
        weights /= weights.sum(axis=1, keepdims=True)
        
        # Predict 'i' as weighted average of training 'i' values
        predictions = np.dot(weights, y_train)
        return predictions

    def tune_bandwidth(self, X_train, y_train):
        """
        Tune the bandwidth for the KDE to match the density of the reference 'i' values.
        """
        print("Tuning bandwidth for Gaussian kernel smoothing...")

        # Use GridSearchCV to find the best bandwidth for kernel smoothing
        params = {'bandwidth': np.logspace(-2, 0, 10)}
        grid_search = GridSearchCV(KernelDensity(kernel='gaussian'), param_grid=params)
        
        # Fit the KDE model on the training data and find the best bandwidth
        grid_search.fit(X_train)
        best_bandwidth = grid_search.best_params_['bandwidth']
        
        print(f"Best bandwidth found: {best_bandwidth}")
        return best_bandwidth

    def train_model_for_i(self):
        """
        Prepare the training data for Gaussian kernel smoothing.
        """
        print("Preparing cell-type proportions for Gaussian kernel smoothing...")
        
        # Create feature matrix (proportions by density bins for reference cells)
        self.X_train = np.array([
            [self.cell_proportions_by_density[bin_range][cluster] for bin_range in self.cell_proportions_by_density]
            for cluster in self.clusters_ref
        ])
        
        self.y_train = self.i_ref.values  # Reference 'i' values
        
        print(f"Prepared training data for {len(self.X_train)} reference cells.")
        
        # Tune bandwidth to match the density of reference i values
        self.bandwidth = self.tune_bandwidth(self.X_train, self.y_train)

    def predict_i_using_model(self):
        """
        Use Gaussian kernel smoothing to predict 'i' for all new cells based on cell-type proportions.
        """
        print(f"Predicting 'i' for {len(self.X_new)} new cells using Gaussian kernel smoothing with bandwidth {self.bandwidth}...")
        
        # Create test feature matrix (proportions by density bins for new cells)
        X_test = np.array([
            [self.cell_proportions_by_density[bin_range][cluster] for bin_range in self.cell_proportions_by_density]
            for cluster in self.predicted_clusters_new
        ])
        
        # Predict 'i' for the new cells using Gaussian kernel smoothing
        self.predicted_i_model = self.gaussian_kernel_smoothing(self.X_train, self.y_train, X_test)
        
        # Add 'predicted_i_ctype' to the new cells' anndata object (cell-type proportion model)
        self.new_cells_adata.obs['predicted_i_ctype'] = self.predicted_i_model
        
        print(f"Predicted 'i' for {len(self.predicted_i_model)} new cells using Gaussian kernel smoothing.")

    def predict_final(self, alpha=0.9):
        """
        Use the default alpha value of 0.5 to make final predictions for 'i' in new cells.
        """
        print(f"Making final predictions for new cells using alpha = {alpha}...")
        
        # Combine the predictions using the specified alpha
        self.combine_predictions(alpha)
        
        # Add the 'predicted_i_combined' column to the anndata object for new cells
        self.new_cells_adata.obs['predicted_i_combined'] = self.predicted_i_final
        
        print(f"Final predictions completed for {len(self.predicted_i_final)} new cells.")
        
        return self.new_cells_adata

    


    def run_all(self, alpha=0.9):
        """
        Run the entire prediction workflow from nearest neighbor predictions to final output.
        
        Parameters:
        - n_neighbors: Number of nearest neighbors to consider for predictions.
        - alpha: Weight for combining predictions from nearest neighbors and model predictions.
        - n_bins: Number of density bins to create for estimating cell-type proportions.
        
        Returns:
        - new_cells_with_predictions: Updated AnnData object with predicted 'i' values for new cells.
        """
        print("Starting the full prediction workflow...")

        # Step 1: Predict nearest neighbor 'i' and clusters for new cells
        self.predict_i_and_clusters()

        # Step 2: Estimate cell-type proportions for the reference dataset by density bins
        self.estimate_cell_proportions()

        # Step 3: Train the model based on the reference dataset
        self.train_model_for_i()

        # Step 4: Predict 'i' using the trained model for new cells
        self.predict_i_using_model()

        # Step 5: Make the final prediction for 'i'
        return self.predict_final(alpha=alpha)


