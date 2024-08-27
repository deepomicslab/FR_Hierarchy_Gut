import pandas as pd
import numpy as np
import seaborn as sns
from sklearn.decomposition import PCA
from matplotlib import pyplot as plt

from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import fdrcorrection
from cliffs_delta import cliffs_delta

import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Compare eigengene networks')
    parser.add_argument('--prefix', type=str, required=True, help='input file prefix')
    parser.add_argument('--pheno', type=str, required=True, help='Phenotype')
    parser.add_argument('--oprefix', type=str, required=True, help='output prefix')
    args = parser.parse_args()
    prefix = args.prefix
    pheno = args.pheno
    expr_df = pd.read_csv(f'../data/{pheno}/{prefix}/abd.tsv', sep='\t', index_col=0)
    meta_df = pd.read_csv(f'../data/{pheno}/{prefix}/metadata.tsv', sep='\t')
    prefix = args.oprefix

# read data
module_df = pd.read_csv('../result/GCN_fix_tree/leaves_cluster.tsv', sep='\t')

# extract gene
expr_df.index = expr_df.index.str.split('|').str[-1]

health_samples = meta_df[meta_df['disease'] == 'Health']['sample_id'].tolist()
disease_samples = meta_df[meta_df['disease'] != 'Health']['sample_id'].tolist()


health_expr = expr_df.loc[:, health_samples]
disease_expr = expr_df.loc[:, disease_samples]

# compute eigengene
eigengene_results = []
for expr_df, group in [(health_expr, 'Health'), (disease_expr, 'Disease')]:
    for module in module_df['cluster'].unique():
        module_genes = module_df[module_df['cluster'] == module]['species'].tolist()
        module_genes = [gene for gene in module_genes if gene in expr_df.index and sum(expr_df.loc[gene] != 0) > 0.2 * len(expr_df.columns)]
        if len(module_genes) > 0:
            module_expr = expr_df.loc[module_genes]
            pca = PCA(n_components=1)
            eigenvectors = pca.fit_transform(module_expr.T)
            eigengene = eigenvectors[:, 0]
            for sample, sample_id in enumerate(expr_df.columns):
                eigengene_results.append({
                    'sample': sample_id,
                    'group': group,
                    'module': module,
                    'eigengene': eigengene[sample]
                })

# summary to DataFrame
eigengene_df = pd.DataFrame(eigengene_results)
eigengene_out = "{}.eigengene.csv".format(prefix)

eigengene_df.to_csv(eigengene_out, sep='\t')



health_eigengene = eigengene_df[eigengene_df['group'] == 'Health']
health_eigengene_matrix = []
health_module_list = []
for module in health_eigengene['module'].unique():
    module_eigengene = health_eigengene[health_eigengene['module'] == module]['eigengene'].tolist()
    if len(module_eigengene) == len(health_eigengene['sample'].unique()) and len(set(module_eigengene)) > 1:
        health_eigengene_matrix.append(module_eigengene)
        health_module_list.append(module)
        #print(module_eigengene)

health_sample_cluster_matrix = pd.DataFrame(health_eigengene_matrix, columns=health_eigengene['sample'].unique(), index=health_module_list)
health_network = health_sample_cluster_matrix.T.corr(method='pearson')
health_network = health_network.dropna(axis=0, how='any')
health_network = health_network.dropna(axis=1, how='any')


health_network.to_csv("{}.health.eigengene_cor.tsv".format(prefix),sep='\t')

plt.figure(figsize=(12, 10))
sns.heatmap(health_network, annot=True, cmap='YlOrRd')
plt.title('Eigengene Matrix')
plt.savefig("{}.health.eigengene_cor.png".format(prefix))

#sns.heatmap(health_network, annot=True, cmap='YlOrRd')


# compute eigengene network for disease group
disease_eigengene = eigengene_df[eigengene_df['group'] != 'Health']
disease_eigengene_matrix = []
disease_module_list = []
for module in disease_eigengene['module'].unique():
    module_eigengene = disease_eigengene[disease_eigengene['module'] == module]['eigengene'].tolist()
    if len(module_eigengene) == len(disease_eigengene['sample'].unique()) and len(set(module_eigengene)) > 1:
        disease_eigengene_matrix.append(module_eigengene)
        disease_module_list.append(module)

disease_sample_cluster_matrix = pd.DataFrame(disease_eigengene_matrix, columns=disease_eigengene['sample'].unique(), index=disease_module_list)
disease_network = disease_sample_cluster_matrix.T.corr(method='pearson')
disease_network = disease_network.dropna(axis=0, how='any')
disease_network = disease_network.dropna(axis=1, how='any')
#disease_network.to_csv('disease_eigengene_network.csv')
#sns.heatmap(disease_network, annot=True, cmap='YlOrRd')



disease_network.to_csv("{}.disease.eigengene_cor.tsv".format(prefix),sep='\t')

plt.figure(figsize=(12, 10))
sns.heatmap(disease_network, annot=True, cmap='YlOrRd')
plt.title('Eigengene Matrix')
plt.savefig("{}.disease.eigengene_cor.png".format(prefix))

sns.heatmap(disease_network, annot=True, cmap='YlOrRd')


def preserv_matrix(network1, network2):
    """
    compute preserv matrix for two eigengene networks Preserv(1,2)
    
    parameters:
    network1 (pd.DataFrame): eigengene network in control
    network2 (pd.DataFrame): eigengene network in disease
    
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


preserv_matrix.to_csv("{}.preserv_matrix.tsv".format(prefix),sep='\t')

plt.figure(figsize=(12, 10))
sns.clustermap(preserv_matrix, annot=True, cmap='YlOrRd')
plt.title('Eigengene Matrix')
plt.savefig("{}.preserv_matrix.png".format(prefix))



def compare_eigengene_networks(health_network, disease_network):
    """
    compare health_network and disease_network for common module of the eigengene difference
    
    parameters:
    health_network (pd.DataFrame)
    disease_network (pd.DataFrame)
    return:
    pd.DataFrame: module, eigengene, p value, p adj, effect size 
    """
    # common modules
    common_modules = list(set(health_network.index) & set(disease_network.index))
    
    results = []
    
    # compare eigengene difference
    for module in common_modules:
        health_eigengene = health_network.loc[module]
        disease_eigengene = disease_network.loc[module]
        
        # Mann-Whitney U test
        u, p_value = mannwhitneyu(health_eigengene, disease_eigengene)
        
        # effect size
        #effect_size = calculate_effect_size(u, len(health_eigengene), len(disease_eigengene))
        #health_mean = health_eigengene.mean()
        #disease_mean = disease_eigengene.mean()
        #health_std = health_eigengene.std()
        #effect_size = (health_mean - disease_mean) / health_std

        effect_size, eff_p_value = cliffs_delta(health_eigengene, disease_eigengene)
        
        # result
        results.append({
            'Module': module,
            'Health Mean': health_eigengene.mean(),
            'Disease Mean': disease_eigengene.mean(),
            'p-value': p_value,
            'effect_size': effect_size
        })
    
    # to DataFrame
    results_df = pd.DataFrame(results)
    
    # FDR adj
    _, fdr_p_values = fdrcorrection(results_df['p-value'].values)
    results_df['FDR p-value'] = fdr_p_values
    
    return results_df


results = compare_eigengene_networks(health_sample_cluster_matrix,disease_sample_cluster_matrix)
results

results.to_csv("{}.compare_eigengene_networks.tsv".format(prefix), sep='\t')