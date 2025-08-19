# GCN checking and functions according to GCN input
import pandas as pd
from distance import *
import numpy as np
import copy
from scipy.spatial.distance import pdist, squareform

# GCN format: pandas dataframe, row_index=functional gene, col_index=Species, sep default= '\t', transfer=False
def input_GCN(GCN_path, sep='\t', transfer=False):
    ref_GCN= pd.read_csv(GCN_path, sep=sep, index_col=0, header=0)
    if transfer:
        ref_GCN = ref_GCN.T
    if check(ref_GCN):
        return ref_GCN
    else:
        exit(1)

def check(ref_GCN):
    # some checking may be needed here
    return True

# compute distance by different functions (default is jaccard distance)
def sp_d(ref_GCN, distance_method='JSD'):
    sp_list = list(ref_GCN.columns)
    if distance_method == 'SORENSEN':
        distance_compressed = pdist(ref_GCN.values.T, sorensen_distance)
    elif distance_method == 'COR':
        distance_compressed = pdist(ref_GCN.values.T, corr_distance)
    else: 
        distance_compressed = pdist(ref_GCN.values.T, Jaccard_distance)
    distance_matrix = squareform(distance_compressed)
    distance_df = pd.DataFrame(distance_matrix, columns = sp_list, index = sp_list)
    return distance_df

# compute functional redundancy index
def fr_df(profile, distance_df):
    sp_list = list(distance_df.columns)
    distance_matrix = distance_df.values
    n_sp = len(sp_list)
    fr_matrix = np.zeros(shape=(n_sp, n_sp))
    max_fr = -np.Inf
    min_fr = np.Inf
    for i in range(n_sp):
        for j in range(i+1, n_sp):
            sp1 = sp_list[i]
            sp2 = sp_list[j]
            x = np.array(profile[sp1])
            y = np.array(profile[sp2])
            FR = np.dot(x, y)*(1-distance_matrix[i][j])
            fr_matrix[i][j] = FR
            # update max and min
            if FR > max_fr:
                max_fr = FR
            if FR < min_fr and FR > 0:
                min_fr = FR
    #print(fr_matrix)
    #print(max_fr, min_fr)            
    log_max = np.log10(max_fr)
    log_min = np.log10(min_fr)
    fr_matrix = fr_normalize(fr_matrix,  log_max, log_min)
    fr_matrix += fr_matrix.T
    row, col = np.diag_indices_from(fr_matrix)
    fr_matrix[row, col] = 1
    fr_df = pd.DataFrame(fr_matrix, columns=sp_list, index=sp_list)
    return fr_df, log_max, log_min
    
# normalize by log10 and rescale to (0,1) divide sum
def fr_normalize(fr_matrix, log_max, log_min):
    range_interval = log_max - log_min
    fr_matrix = np.array(list(map(np.log10, fr_matrix)))
    fr_matrix -= log_min
    fr_matrix /= range_interval
    fr_matrix[fr_matrix == -np.inf] = 0
    return fr_matrix

# update GCN with cluster result
def cluster_GCN(ref_GCN, cluster_sp):
    for cluster, cluster_sp_list in cluster_sp.items():
        ref_GCN['cluster_'+str(cluster)] = np.median(ref_GCN[cluster_sp_list].values, axis=1)
        ref_GCN.drop(columns=cluster_sp_list, inplace=True)
    return ref_GCN


def reciprocal(fr_df):
    value_df = fr_df.replace(0, 1000)
    value_matrix = 1./(value_df.values)
    return value_matrix

# used when hierarchical structure from SEAT
# split GCN for each cluster
def level_GCN(cluster_sp_dict, ori_GCN):
    level_result = copy.deepcopy(ori_GCN)
    for name, cluster_sp in cluster_sp_dict.items():
        cluster_GCN = level_result[cluster_sp]
        sum = np.median(cluster_GCN.values, axis=1)
        level_result.drop(columns=cluster_sp, inplace = True)
        level_result[name] = sum
    return level_result
 