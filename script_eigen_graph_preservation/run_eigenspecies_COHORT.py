# This script computes eigenspecies (by phenotype) and compare eigenspecies networks between health and disease groups.
# This script also computes preservation.
import pandas as pd
import numpy as np
import seaborn as sns
from sklearn.decomposition import PCA
from matplotlib import pyplot as plt
import os
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import fdrcorrection
from cliffs_delta import cliffs_delta

import argparse

module_df = pd.read_csv('../result/GCN_fix_tree/leaves_cluster.tsv', sep='\t')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Compare eigensp networks')
    parser.add_argument('--prefix', type=str, required=True, help='Output file prefix')
    parser.add_argument('--pheno', type=str, required=True, help='Phenotype')
    parser.add_argument('--idir', type=str, required=True, help='input dir')
    parser.add_argument('--odir', type=str, required=True, help='output dir')
    parser.add_argument('--cluster', type=str, required=True, help='cluster leaves')
    args = parser.parse_args()
    prefix = args.prefix
    pheno = args.pheno
    idir = args.idir
    odir = args.odir
    module_df = pd.read_csv(args.cluster, sep='\t')
    expr_df = pd.read_csv(os.path.join(idir,pheno,prefix,'abd.tsv'), sep='\t', index_col=0)
    meta_df = pd.read_csv(os.path.join(idir,pheno,prefix,'metadata.tsv'), sep='\t')


# extract species abundance profile 
expr_df.index = expr_df.index.str.split('|').str[-1]


# get sample ids of different phenotypes
health_samples = meta_df[meta_df['disease'] == 'Health']['sample_id'].tolist()
disease_samples = meta_df[meta_df['disease'] != 'Health']['sample_id'].tolist()

health_expr = expr_df.loc[:, health_samples]
disease_expr = expr_df.loc[:, disease_samples]

# compute eigensp
eigensp_results = []
for expr_df, group in [(health_expr, 'Health'), (disease_expr, 'Disease')]:
    for module in module_df['cluster'].unique():
        module_genes = module_df[module_df['cluster'] == module]['species'].tolist()
        module_genes = [gene for gene in module_genes if gene in expr_df.index and sum(expr_df.loc[gene] != 0) > 0.2 * len(expr_df.columns)]
        if len(module_genes) > 0:
            module_expr = expr_df.loc[module_genes]
            pca = PCA(n_components=1)
            eigenvectors = pca.fit_transform(module_expr.T)
            eigensp = eigenvectors[:, 0]
            for sample, sample_id in enumerate(expr_df.columns):
                eigensp_results.append({
                    'sample': sample_id,
                    'group': group,
                    'module': module,
                    'eigensp': eigensp[sample]
                })

# summary to DataFrame
eigensp_df = pd.DataFrame(eigensp_results)
eigensp_out = os.path.join(odir, pheno, "{}.eigensp.csv".format(prefix))

eigensp_df.to_csv(eigensp_out, sep='\t')


# plot eigen sp network heatmap for health group
health_eigensp = eigensp_df[eigensp_df['group'] == 'Health']
health_eigensp_matrix = []
health_module_list = []
for module in health_eigensp['module'].unique():
    module_eigensp = health_eigensp[health_eigensp['module'] == module]['eigensp'].tolist()
    if len(module_eigensp) == len(health_eigensp['sample'].unique()) and len(set(module_eigensp)) > 1:
        health_eigensp_matrix.append(module_eigensp)
        health_module_list.append(module)
        #print(module_eigensp)

health_sample_cluster_matrix = pd.DataFrame(health_eigensp_matrix, columns=health_eigensp['sample'].unique(), index=health_module_list)
health_network = health_sample_cluster_matrix.T.corr(method='pearson')
health_network = health_network.dropna(axis=0, how='any')
health_network = health_network.dropna(axis=1, how='any')


health_network.to_csv(os.path.join(odir, pheno, "{}.health.eigensp_cor.tsv".format(prefix)),sep='\t')

plt.figure(figsize=(12, 10))
sns.heatmap(health_network, annot=True, cmap='YlOrRd')
plt.title('Eigensp Matrix')
plt.savefig(os.path.join(odir, pheno, "{}.health.eigensp_cor.png".format(prefix)))

#sns.heatmap(health_network, annot=True, cmap='YlOrRd')


# plot eigen sp network heatmap for disease group
disease_eigensp = eigensp_df[eigensp_df['group'] != 'Health']
disease_eigensp_matrix = []
disease_module_list = []
for module in disease_eigensp['module'].unique():
    module_eigensp = disease_eigensp[disease_eigensp['module'] == module]['eigensp'].tolist()
    if len(module_eigensp) == len(disease_eigensp['sample'].unique()) and len(set(module_eigensp)) > 1:
        disease_eigensp_matrix.append(module_eigensp)
        disease_module_list.append(module)

disease_sample_cluster_matrix = pd.DataFrame(disease_eigensp_matrix, columns=disease_eigensp['sample'].unique(), index=disease_module_list)
disease_network = disease_sample_cluster_matrix.T.corr(method='pearson')
disease_network = disease_network.dropna(axis=0, how='any')
disease_network = disease_network.dropna(axis=1, how='any')
#disease_network.to_csv('disease_eigensp_network.csv')
#sns.heatmap(disease_network, annot=True, cmap='YlOrRd')



disease_network.to_csv(os.path.join(odir, pheno, "{}.disease.eigensp_cor.tsv".format(prefix)),sep='\t')

plt.figure(figsize=(12, 10))
sns.heatmap(disease_network, annot=True, cmap='YlOrRd')
plt.title('Eigensp Matrix')
plt.savefig(os.path.join(odir, pheno, "{}.disease.eigensp_cor.png".format(prefix)))

sns.heatmap(disease_network, annot=True, cmap='YlOrRd')

# compute preservation
def preserv_matrix(network1, network2):
    """
    compute preserv matrix for two eigensp networks Preserv(1,2)
    
    parameters:
    network1 (pd.DataFrame): eigensp network in control
    network2 (pd.DataFrame): eigensp network in disease
    
    return:
    pd.DataFrame: Preserv(1,2) matrix
    """
    # common modules
    common_modules = list(set(network1.index) & set(network2.index))
    
    # initialize preserv matrix
    preserv_matrix = pd.DataFrame(index=common_modules, columns=common_modules)
    
    for i, module1 in enumerate(common_modules):
        for j, module2 in enumerate(common_modules):
            #preserv_matrix.iloc[i, j] = 1 - max(abs(network1.loc[module1, module2] - network2.loc[module1, module2]),
                                              #abs(network1.loc[module1, module2] - network2.loc[module2, module1]))
            preserv_matrix.iloc[i, j] = 1 - (max(network1.loc[module1, module2], network2.loc[module1, module2]) - min(network1.loc[module1, module2], network2.loc[module1, module2]))
    
    return preserv_matrix

# compute preserv matrix Preserv(1,2)
preserv_matrix = preserv_matrix(health_network, disease_network)


preserv_matrix
preserv_matrix = preserv_matrix.astype(float)
#sns.clustermap(preserv_matrix, annot=True, cmap='YlOrRd')


preserv_matrix.to_csv(os.path.join(odir, pheno, "{}.preserv_matrix.tsv").format(prefix),sep='\t')

plt.figure(figsize=(12, 10))
sns.clustermap(preserv_matrix, annot=True, cmap='YlOrRd')
plt.title('Eigensp Matrix')
plt.savefig(os.path.join(odir, pheno, "{}.preserv_matrix.png".format(prefix)))


# differential test for health network and disease_network
def compare_eigensp_networks(health_network, disease_network):
    """
    compare health_network and disease_network for common module of the eigensp difference
    
    parameters:
    health_network (pd.DataFrame)
    disease_network (pd.DataFrame)
    return:
    pd.DataFrame: module, eigensp, p value, p adj, effect size 
    """
    # common modules
    common_modules = list(set(health_network.index) & set(disease_network.index))
    
    results = []
    
    # compare eigensp difference
    for module in common_modules:
        health_eigensp = health_network.loc[module]
        disease_eigensp = disease_network.loc[module]
        
        # Mann-Whitney U test
        u, p_value = mannwhitneyu(health_eigensp, disease_eigensp)
        
        # effect size
        #effect_size = calculate_effect_size(u, len(health_eigensp), len(disease_eigensp))
        #health_mean = health_eigensp.mean()
        #disease_mean = disease_eigensp.mean()
        #health_std = health_eigensp.std()
        #effect_size = (health_mean - disease_mean) / health_std

        effect_size, eff_p_value = cliffs_delta(health_eigensp, disease_eigensp)
        
        # result
        results.append({
            'Module': module,
            'Health Mean': health_eigensp.mean(),
            'Disease Mean': disease_eigensp.mean(),
            'p-value': p_value,
            'effect_size': effect_size
        })
    
    # to DataFrame
    results_df = pd.DataFrame(results)
    
    # FDR adj
    _, fdr_p_values = fdrcorrection(results_df['p-value'].values)
    results_df['FDR p-value'] = fdr_p_values
    
    return results_df


results = compare_eigensp_networks(health_sample_cluster_matrix,disease_sample_cluster_matrix)
results

results.to_csv(os.path.join(odir, pheno, "{}.compare_eigensp_networks.tsv".format(prefix)), sep='\t')