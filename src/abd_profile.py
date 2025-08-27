# functions to process profile data
import pandas as pd
import numpy as np
import sklearn
import copy

# input format : col_index=sp, row_index=sname, sep default='\t', transfer=False

def input_profile(profile_path, sep='\t', transfer=False):
    profile= pd.read_csv(profile_path, sep=sep, index_col=0, header=0)
    if transfer:
        profile = profile.T
    return profile

# check species not in GCN
def check(profile, ref_GCN, name_dict = None):
    # check if all species is in ref GCN
    profile = rename_s_level(profile)
    if name_dict:
        profile = profile.rename(columns=name_dict)
    diff = set(profile.columns).difference(set(ref_GCN.columns))
    if len(diff) > 0:
        print("Species not in GCN2008: {}".format(diff))
        print('!!!species are deleted')
        profile.drop(columns=list(diff), inplace=True)
    profile = clean(profile)
    return profile

# get species level profile
def rename_s_level(profile):
    # choose species level
    species = {}
    for tax in list(profile.columns):
        if tax.split('|')[-1][0] == "s":
            species[tax] = tax.split('|')[-1]
    profile = profile[list(species.keys())].rename(columns=species)
    profile = profile.loc[(profile > 0).any(axis=1), :]
    return profile
    

# normalize delete 0 columns
def clean(profile):
    # delete 0 columns and 0 rows
    profile = delete_zero_col(profile)
    profile = delete_zero_row(profile)

    # normalize
    profile = row_normalize(profile)
    return profile

def delete_zero_col(profile):
    profile = profile.loc[:, (profile > 0).any(axis=0)]
    return profile

def delete_zero_row(profile):
    profile = profile.loc[(profile > 0).any(axis=1), :]
    return profile

# used after sample cluster
def split_to_clusters(profile, cluster_labels):
    cluster_profiles = {}
    cluster_samples = {}
    sample_list = list(profile.index)
    for label in set(cluster_labels):
        cluster_samples[label] = []
        
    for i, label in enumerate(cluster_labels):
        cluster_samples[label].append(sample_list[i])

    for label, v in cluster_samples.items():
        cluster_profiles[label] = delete_zero_col(profile.loc[v])
    return cluster_profiles

# used when hierarchical structure from SEAT
# split profile to each clsuter
def level_profile(cluster_sp_dict, ori_profile):
    level_result = copy.deepcopy(ori_profile)
    for name, cluster_sp in cluster_sp_dict.items():
        cluster_abd = level_result[cluster_sp]
        sum = np.sum(cluster_abd, axis=1)
        cluster_sp.sort()
        level_result.drop(columns = cluster_sp, inplace= True)
        level_result[name] = sum
    return level_result

# total abundance computation for each cluster
def level_profile_stree(cluster_sp_dict, ori_profile):
    level_result = copy.deepcopy(ori_profile)
    for name, cluster_sp in cluster_sp_dict.items():
        cluster_common = list(set(cluster_sp).intersection(set(ori_profile.columns)))
        cluster_abd = level_result[cluster_common]
        sum = np.sum(cluster_abd, axis=1)
        #cluster_sp.sort()
        level_result.drop(columns = cluster_common, inplace= True)
        level_result[name] = sum
    return level_result

# convert taxonomy profile to KO profile
def ko_profile(tax_profile, ref_GCN):
    # tax_profile: columns = sp, rows = sample
    # ref_GCN: columns = ko, rows = sp
    short_GCN = ref_GCN.loc[tax_profile.columns, :]
    value = np.dot(tax_profile.values, short_GCN.values)
    ko_profile = pd.DataFrame(value, columns = short_GCN.columns, index=tax_profile.index)
    ko_profile = row_normalize(ko_profile)
    ko_profile = delete_zero_col(ko_profile)
    ko_profile = delete_zero_row(ko_profile)
    return ko_profile

def row_normalize(df):
    new_df = df.div(df.sum(axis=1), axis=0)
    return new_df

def col_normalize(df):
    new_df = df.div(df.sum(axis=0), axis=1)
    return new_df