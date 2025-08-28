# Required imports
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import fdrcorrection
import matplotlib.pyplot as plt
import seaborn as sns
from cliffs_delta import cliffs_delta


def calculate_eigenspecies(expr_data, species_FRC, group):
    """
    Calculate eigenspecies for a specific group
    
    Parameters:
    -----------
    expr_data : DataFrame
        Expression data for the group
    species_FRC : DataFrame
        Species clustering information
    group : str
        Group name
    
    Returns:
    --------
    list
        List of eigenspecies results
    """
    eigenspecies_results = []
    
    for cluster in species_FRC['cluster'].unique():
        cluster_speciess = species_FRC[species_FRC['cluster'] == cluster]['species'].tolist()
        cluster_speciess = [species for species in cluster_speciess if species in expr_data.index and sum(expr_data.loc[species] != 0) > 0.1 * len(expr_data.columns)]
        if len(cluster_speciess) > 0:
            cluster_expr = expr_data.loc[cluster_speciess]
            pca = PCA(n_components=1)
            eigenvectors = pca.fit_transform(cluster_expr.T)
            eigenspecies = eigenvectors[:, 0]
            for sample, sample_id in enumerate(expr_data.columns):
                eigenspecies_results.append({
                    'sample': sample_id,
                    'group': group,
                    'cluster': cluster,
                    'eigenspecies': eigenspecies[sample]
                })
    
    return eigenspecies_results


def calculate_eigenspecies_together(expr_data, species_FRC, meta_df, g1_samples, g2_samples, g1, g2):
    """
    Calculate eigenspecies for all samples together
    
    Parameters:
    -----------
    expr_data : DataFrame
        Expression data for all samples
    species_FRC : DataFrame
        Species clustering information
    meta_df : DataFrame
        Metadata information
    g1_samples : list
        List of group 1 sample IDs
    g2_samples : list
        List of group 2 sample IDs
    g1 : str
        Group 1 name
    g2 : str
        Group 2 name
    
    Returns:
    --------
    list
        List of eigenspecies results
    """
    eigenspecies_results = []
    
    # Merge all samples
    all_samples = g1_samples + g2_samples

    # Ensure expr_data only contains the samples we need
    expr_data = expr_data.loc[:, all_samples]
    
    for cluster in species_FRC['cluster'].unique():
        cluster_speciess = species_FRC[species_FRC['cluster'] == cluster]['species'].tolist()
        # filter species based on expr
        cluster_speciess = [species for species in cluster_speciess if species in expr_data.index and sum(expr_data.loc[species] != 0) > 0.1 * len(expr_data.columns)]
        
        if len(cluster_speciess) > 0:
            cluster_expr = expr_data.loc[cluster_speciess]
            pca = PCA(n_components=1)
            eigenvectors = pca.fit_transform(cluster_expr.T)
            eigenspecies = eigenvectors[:, 0]
            
            for sample, sample_id in enumerate(expr_data.columns):
                # Determine group based on sample ID
                group = g1 if sample_id in g1_samples else g2
                
                eigenspecies_results.append({
                    'sample': sample_id,
                    'group': group,
                    'cluster': cluster,
                    'eigenspecies': eigenspecies[sample]
                })
    
    return eigenspecies_results




def eigenspecies_correlation_network(eigenspecies_df, group_value, prefix):
    """
    Create and save eigenspecies network for a specific group
    
    Parameters:
    -----------
    eigenspecies_df : DataFrame
        DataFrame containing eigenspecies data
    group_value : str
        The group value to filter by (e.g., 'Health', 'Disease')
    prefix : str
        Prefix for output files
        
    Returns:
    --------
    network : DataFrame
        Correlation network
    sample_cluster_matrix : DataFrame
        Matrix of eigenspecies values
    """
    
    # 1. Filter data for the specific group
    filtered_eigenspecies = eigenspecies_df[eigenspecies_df['group'] == group_value]
    
    # 2. Prepare matrices
    eigenspecies_matrix = []
    cluster_list = []
    
    # 3. Get unique samples from the filtered data
    unique_samples = filtered_eigenspecies['sample'].unique()
    
    # 4. Process each cluster
    for cluster in filtered_eigenspecies['cluster'].unique():
        # Get values for this cluster
        cluster_values = filtered_eigenspecies[filtered_eigenspecies['cluster'] == cluster]['eigenspecies'].tolist()
        
        # Only include clusters with complete data and variation
        if len(cluster_values) == len(unique_samples) and len(set(cluster_values)) > 1:
            # Calculate standard deviation to ensure it's not zero
            std_dev = np.std(cluster_values)
            if std_dev > 0:  # Only add clusters with non-zero standard deviation
                eigenspecies_matrix.append(cluster_values)
                cluster_list.append(cluster)
            else:
                print(f"Warning: Cluster {cluster} has zero standard deviation (all values identical). Skipping.")
    
    # 5. Create sample-cluster matrix
    sample_cluster_matrix = pd.DataFrame(
        eigenspecies_matrix, 
        columns=unique_samples, 
        index=cluster_list
    )
    
    # 6. Calculate correlation network and remove rows/columns with all NA values
    network = sample_cluster_matrix.T.corr(method='pearson')
    na_rows = network.isna().all(axis=1)
    na_cols = network.isna().all(axis=0)
    network = network.loc[~na_rows, ~na_cols]
    
    return network, sample_cluster_matrix
    
def get_preserv_matrix(network1, network2):
    """
    Calculate the Preserv(1,2) matrix between two networks
    
    Parameters:
    -----------
    network1 : DataFrame
        Eigenspecies network for the first group
    network2 : DataFrame
        Eigenspecies network for the second group
    
    Returns:
    --------
    DataFrame
        Preserv(1,2) matrix
    """
    # Get common modules between the two networks
    common_clusters = list(set(network1.index) & set(network2.index))
    
    # Create Preserv(1,2) matrix
    preserv_matrix = pd.DataFrame(index=common_clusters, columns=common_clusters)
    
    for i, cluster1 in enumerate(common_clusters):
        for j, cluster2 in enumerate(common_clusters):
            #preserv_matrix.iloc[i, j] = 1 - max(abs(network1.loc[cluster1, cluster2] - network2.loc[cluster1, cluster2]),
                                              #abs(network1.loc[cluster1, cluster2] - network2.loc[cluster2, cluster1]))
            preserv_matrix.iloc[i, j] = 1 - (max(network1.loc[cluster1, cluster2], network2.loc[cluster1, cluster2]) - min(network1.loc[cluster1, cluster2], network2.loc[cluster1, cluster2]))
    
    return preserv_matrix


    
def compare_eigenspecies_networks(g1_network, g2_network):
    """
    Compare eigenspecies differences between common modules in g1_network and g2_network
    
    Parameters:
    -----------
    g1_network : DataFrame
        Eigenspecies network for the first group (e.g., healthy)
    g2_network : DataFrame
        Eigenspecies network for the second group (e.g., disease)
    
    Returns:
    --------
    DataFrame
        Comparison results containing module names, mean eigenspecies for both groups,
        p-values, FDR-corrected p-values, and effect sizes
    """
    # Find common modules between the two networks
    common_clusters = list(set(g1_network.index) & set(g2_network.index))
    
    # Initialize results list
    results = []
    
    # Compare each common module
    for cluster in common_clusters:
        g1_eigenspecies = g1_network.loc[cluster]
        g2_eigenspecies = g2_network.loc[cluster]
        
        # Calculate Mann-Whitney U test p-value
        u, p_value = mannwhitneyu(g1_eigenspecies, g2_eigenspecies)
        effect_size, eff_p_value = cliffs_delta(g1_eigenspecies, g2_eigenspecies)
        
        # Record results
        results.append({
            'cluster': cluster,
            'Health Mean': g1_eigenspecies.mean(),
            'Disease Mean': g2_eigenspecies.mean(),
            'p-value': p_value,
            'effect_size': effect_size
        })
    
    # Convert results to DataFrame
    results_df = pd.DataFrame(results)
    
    # Apply FDR correction to p-values
    _, fdr_p_values = fdrcorrection(results_df['p-value'].values)
    results_df['FDR p-value'] = fdr_p_values
    
    return results_df