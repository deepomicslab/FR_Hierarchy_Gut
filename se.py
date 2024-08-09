import pandas as pd
import numpy as np
import copy

# input is node:parent, node:leaves, node:direct childrens, node: childrens 

def node_se(edge_df, sp_list, parent_list):
    degree = edge_df.sum()
    #print(degree)
    v_g = degree.sum()
    v_parent = degree[parent_list].sum()
    v_node = degree[sp_list].sum()
    v_outside = v_node - edge_df.loc[sp_list, sp_list].sum().sum()
    if v_parent*v_node == 0:
        return 0, v_outside
    value = - v_outside/v_g*np.log2(v_node/v_parent)
    return value, v_outside

def all_gamma(direct_children_dict, out_v_dict, node_leaves, m2):
    gamma_dict = {} 
    s_dict = {}
    log_s_dict = {}
    for node in node_leaves.keys():
        value = 0
        direct_children = direct_children_dict[node]
        for c in direct_children:
            if c in out_v_dict:
                value += out_v_dict[c]
        gamma_dict[node] = value
        s_dict[node] = value - out_v_dict[node]
        if s_dict[node] <= 0:
            log_s_dict[node] = 0
        else:
            log_s_dict[node] = -s_dict[node]/m2*np.log2(s_dict[node]/m2)
    return copy.deepcopy(log_s_dict)
        
def se_sum(se_dict, children_list, node):
    value = 0
    for c in children_list:
        if c in se_dict.keys():
            value += se_dict[c]
    if node in children_list:
        return value
    value += se_dict[node]
    return value

def adj_sum(log_s_dict, children_list, node):
    value = 0
    for c in children_list:
        if c in log_s_dict.keys():
            value += log_s_dict[c]
    if node in children_list:
        return value
    value += log_s_dict[node]
    return value

def all_node_se(edge_df, parent_dict, node_leaves):
    se_dict = {'root': 0}
    outside_dict = {'root': 0}
    valid_leaves = set(edge_df.index)
    for l in valid_leaves:
        sp_list = [l]
        parent = parent_dict[l]
        if parent not in parent_dict.keys():
            parent_list = list(valid_leaves)
        else:
            parent_list = list(set(node_leaves[parent]).intersection(valid_leaves))
        value, out_v = node_se(edge_df, sp_list, parent_list)
        se_dict[l] = value
        outside_dict[l] = out_v

    for node, sp_list in node_leaves.items():
        sp_list = list(set(sp_list).intersection(valid_leaves))
        if node not in parent_dict.keys():
            continue
        parent = parent_dict[node]
        if parent not in parent_dict:
            parent_list = list(valid_leaves)
        else:
            parent_list = list(set(node_leaves[parent]).intersection(valid_leaves))
        value, out_v = node_se(edge_df, sp_list, parent_list)
        se_dict[node] = value
        outside_dict[node] = out_v
    return copy.deepcopy(se_dict), copy.deepcopy(outside_dict)

def subtree_se_adj(edge_df, parent_dict, node_leaves, child_dict, direct_children_dict, param):
    v_g = edge_df.sum().sum()
    #v_g = 1
    se_dict, outside_dict = all_node_se(edge_df, parent_dict, node_leaves)
    #print(len(outside_dict), len(outside_dict))
    log_s_dict = all_gamma(direct_children_dict, outside_dict, node_leaves, v_g)
    # print(len(se_dict), len(log_s_dict))
    result = {}
    for node in node_leaves.keys():
        children_list = child_dict[node]
        value = se_sum(se_dict, children_list, node) + param*adj_sum(log_s_dict, children_list, node)
        result[node] = value
    return copy.deepcopy(result)

def subtree_se(edge_df, parent_dict, node_leaves, child_dict):
    se_dict, outside_dict = all_node_se(edge_df, parent_dict, node_leaves)
    result = {}
    for node in se_dict.keys():
        children_list = child_dict[node]
        value = se_sum(se_dict, children_list, node)
        result[node] = value
    return copy.deepcopy(result)