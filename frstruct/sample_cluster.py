# function to cluster the samples according to profile
import numpy as np
import scipy.spatial.distance as scipyDistance
from sklearn_extra.cluster import KMedoids
from sklearn.metrics import calinski_harabasz_score as CH
import copy
from pyseat.SEAT import SEAT
from pyseat.HierachicalEmbedding import HierachicalEmbedding
import pandas as pd

def sample_cluster(profile, params):
    #print("clustering")
    #print(profile.shape)
    if params['sample cluster'] == 'seat':
        return seat(profile)
    else:
        if params['max_cluster'] != None:
            return enterotype(profile, params['max_cluster'])
        else:
            return enterotype(profile)

def enterotype(profile, max_cluster=20):
    sample_list = list(profile.index)
    n_sample = len(sample_list)
    profile = np.array(profile.values)
    # compute distance matrix
    JSD_matrix = np.zeros(shape=(n_sample, n_sample))
    for i in range(n_sample):
        for j in range(i+1, n_sample):
            sample_vector1 = profile[i]
            sample_vector2 = profile[j]
            JSD = scipyDistance.jensenshannon(sample_vector1, sample_vector2)
            JSD_matrix[i][j] = JSD
    JSD_matrix += JSD_matrix.T
    # print(JSD_matrix)
    max_score = 0
    best_label = None
    n_cluster = 0
    if max_cluster < 2:
        print("There must be more than 2 cluster")
        max_cluster = 2
    max_cluster = max(min(max_cluster, int(profile.shape[0]/2)), 2)
    for i in range(2, max_cluster+1):
        round_score, round_label = enterotype_cluster(profile, JSD_matrix, i)
        if round_score > max_score:
            max_score = round_score
            best_label = copy.deepcopy(round_label)
            n_cluster = i
        # print(i, round_score)
    return n_cluster, max_score, best_label


def enterotype_cluster(profile, distance, n_cluster):
    # PAM algorithm
    km_model = KMedoids(n_clusters = n_cluster, random_state = 0, metric = 'precomputed', method = 'pam', init =  'k-medoids++').fit(distance)  
    cluster_label = km_model.labels_
    # compute CH score
    score = CH(profile, cluster_label)
    return score, cluster_label
    

def seat(profile, fn):
    seat = SEAT(affinity="gaussian_kernel",
                sparsification="knn_neighbors",
                objective="SE",
                strategy="bottom_up",
                max_k = 5)
    seat.fit_predict(profile.values)
    cluster_label = seat.labels_
    n_cluster = len(set(cluster_label))
    y = pd.DataFrame({
    'Global': [1]*len(seat.labels_),
    'Subpopulation': seat.labels_,
    'Club': seat.clubs})
    y.index = profile.index
    HE = HierachicalEmbedding(device='cpu', n_epochs=200, init='spectral', random_state=0,min_dist=0.1,n_components=2)
    embed = HE.fit_transform(seat.aff_m, y, thetas=[1, 1, 1])
    HE.viz_fit(fn=fn)
    return n_cluster, cluster_label
    
def seat_plot(profile, params, fn):
    seat = SEAT(affinity=params['aff'],
            sparsification=params['sf'],
            objective=params['ob'],
            strategy=params['st'],
            max_k = params['maxk'])
    seat.fit_predict(profile.values)
    cluster_label = seat.labels_
    # compute CH score
    score = CH(profile, cluster_label)
    n_cluster = len(set(cluster_label))
    y = pd.DataFrame({
    'Global': [1]*len(seat.labels_),
    'Subpopulation': seat.labels_,
    'Club': seat.clubs})
    y.index = profile.index
    HE = HierachicalEmbedding(device='cpu', n_epochs=200, init='spectral', random_state=0,min_dist=0.1,n_components=2)
    embed = HE.fit_transform(seat.aff_m, y, thetas=[1, 1, 1])
    print('SEAT hierachical embedding: ', embed.shape)
    HE.viz_fit(fn=fn)
    return n_cluster, score, cluster_label
 