import abd_profile
import detect_keystone
import copy
import FR
import GCN
import sample_cluster
import tree_util
import warnings
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import json
import itol_util
import scipy.stats

# avoid log warning
warnings.filterwarnings('ignore')

def level_process(profile, ref_GCN):
    level_result = {"keystone_node": None, 'fr': None, 'scores': None}
    distance = GCN.sp_d(ref_GCN)
    fr_df, log_max, log_min = FR.fr_df(profile, distance)
    cluster_idx, scores = detect_keystone.pr(fr_df.values)
    keystone_cluster = list((fr_df.columns))[cluster_idx]
    level_result["keystone_node"] = keystone_cluster
    level_result['fr'] = fr_df
    level_result['scores'] = scores
    return level_result

def cluster_process(cluster_profile, distance_df, ref_GCN, params):
    result = {'newick': None, "level_result": {}, "reflection": {}, 'seat': None}
    # level = {"keystone_node": None, 'fr': None, }
    # print('profile_size :{}'.format(cluster_profile.shape))
    fr_d, log_max, log_min = FR.fr_df(cluster_profile, distance_df)
    # print('fr_d size {}'.format(fr_d.shape))
    # test normal
    fr_list = []
    for i in range(fr_d.shape[0]):
        for j in range(i+1, fr_d.shape[0]):
            if fr_d.iloc[i, j] > 0:
                fr_list.append(fr_d.iloc[i, j])

    tree = tree_util.make_tree(fr_d)
    # result['seat'] = copy.deepcopy(tree)
    result['seat'] = copy.deepcopy(tree.Z_)
    newick_tree = tree.newick
    
    name_dict, reverse_dict = tree_util.name_reflection(fr_d)
    result['reflection']['name_dict'] = copy.deepcopy(name_dict)
    result['reflection']['reverse_dict'] = copy.deepcopy(reverse_dict)
    json_tree = tree_util.parse(newick_tree)
    largest = {'largest': 0}
    leaf_list, l = tree_util.recu_compute(json_tree, 0, largest)
    largest_level = largest['largest']
    if "max_depth" in list(params.keys()):
        top_list = []
        tree_util.depth_limit_newick(json_tree, params['max_depth'], top_list)
        newick_tree = tree_util.limit_newick_last(top_list)
        json_tree = tree_util.parse(newick_tree)
        nlayer = params['max_depth']
    else:
        nlayer = largest_level
    result['newick'] = newick_tree  

    leaf_list, l = tree_util.recu_compute(json_tree, 0, largest)
    layer_leaves_dict = tree_util.make_layer_dict(nlayer)
    
    tree_util.recu_layer(json_tree, layer_leaves_dict)
    tree_util.to_layer_leaves(layer_leaves_dict, nlayer)
    renamed_GCN = ref_GCN[list(fr_d.columns)].rename(columns=name_dict)
    renamed_profile = cluster_profile.rename(columns=name_dict)
    
    for i in range(1, nlayer):
        if len(list(layer_leaves_dict[i].keys())) > 2:
            level_profile = abd_profile.level_profile(layer_leaves_dict[i], renamed_profile)
            level_GCN = GCN.level_GCN(layer_leaves_dict[i], renamed_GCN)
            result['level_result'][i] = copy.deepcopy(level_process(level_profile, level_GCN))
        else:
            result['level_result'][i] = "parent of level {}".format(i+1)
    
    # compute leaf layer
    cluster_idx, scores = detect_keystone.pr(fr_d.values)
    keystone_cluster = name_dict[list((fr_d.columns))[cluster_idx]]
    result['level_result'][0] = {}
    result['level_result'][0]["keystone_node"] = keystone_cluster
    result['level_result'][0]['fr'] = copy.deepcopy(fr_d)
    result['level_result'][0]['scores'] = scores
    result['leaves_dict'] = copy.deepcopy(layer_leaves_dict)

    return result

def post_process(cluster_profiles, final_result, outdir):
    # draw the tree 
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    leaf_list = []
    for k, v in final_result.items():
        cluster_dir = os.path.join(outdir, 'cluster_{}'.format(k))
        if not os.path.exists(cluster_dir):
            os.mkdir(cluster_dir)
        cluster_profiles[k].to_csv(os.path.join(cluster_dir, 'cluster_profile.tsv'), sep='\t')

        newick = v['newick']
        reverse_dict = v['reflection']['reverse_dict']
        json_tree = tree_util.parse(v['newick'])
        renamed_tree = copy.deepcopy(json_tree)
        newick = tree_util.call_newick(renamed_tree)
        leaf_list = [x.replace('_', '-') for x in reverse_dict.values()]
        newick = newick.replace('_', '-')
        
        for node_name, real_name in reverse_dict.items():
            pos = 0
            while pos < len(newick):
                pos = newick.find(node_name, pos)
                if pos >=0:
                    next_char = newick[pos + len(node_name)]
                    if not next_char.isdigit() and newick[pos-1] in [',', ')', '(']:
                        # replace
                        newick = newick[:pos] + real_name.replace('_', '-') + newick[(pos + len(node_name)):]
                        leaf_list.append(real_name.replace('_', '-'))
                        
                        break
                    else:
                        pos += 1
                else:
                    break
        
        # add parent information
        
        parent_dict = {}
        adding = {}
        tree_util.parents(json_tree, parent_dict) 
        for n, p in parent_dict.items():
            if n in list(reverse_dict.keys()):
                adding[reverse_dict[n]] = p
        for n, p in adding.items():
            parent_dict[n] = p
        for n, p in parent_dict.items():
            if p in list(reverse_dict.keys()):
                parent_dict[n] = reverse_dict[p]
     
        opath1 = os.path.join(cluster_dir, 'tree.newick')
        with open(opath1, 'w') as fp:
            fp.write(newick)

        opath2 = os.path.join(cluster_dir, "keystone_node.tsv")
        node_list = []
        keystone_df = pd.DataFrame(columns=['node', 'layer', 'PR_score', 'is_keystone', 'leaves', 'parent'])
        for l, value in v['level_result'].items():
            if type(value) != str:
                level_dir = os.path.join(cluster_dir, 'layer_{}'.format(l))
                if not os.path.exists(level_dir):
                    os.mkdir(level_dir)
                value['fr'].to_csv(os.path.join(level_dir, 'fr.tsv'), sep = '\t')
                level_pr_path= os.path.join(level_dir, 'PR_scores.tsv')
                level_pr_df = pd.DataFrame(columns=['Node', 'PR'])
                
                level_node_list = list(value['fr'].columns)
                
                scores = value['scores']
                for nid, node in enumerate(level_node_list):
                    if node in list(reverse_dict.keys()):
                        new_node_name = reverse_dict[node]
                    else:
                        new_node_name = node
                    level_pr_df.loc[new_node_name, 'PR'] = scores[nid]
                    level_pr_df.loc[new_node_name, 'Node'] = new_node_name

                    if not (node.startswith('s__') or node.startswith('K')):
                        leaves_list = v['leaves_dict'][l][node]
                        new_leaves_list = [reverse_dict[leaf] for leaf in leaves_list]
                        leaves_s = ','.join(new_leaves_list)
                    else:
                        leaves_s = new_node_name
                    
                    if node == value['keystone_node'] or (l==0 and node == reverse_dict[value['keystone_node']]):
                        is_keystone = True
                    else:
                        is_keystone = False
                    parent = parent_dict[node]
                    keystone_df.loc[keystone_df.shape[0]] = [new_node_name, l, scores[nid], is_keystone, leaves_s, parent]
                level_pr_df = level_pr_df.sort_values(by='PR', ascending=False)
                level_pr_df.to_csv(level_pr_path, sep='\t', index=None)
                if (value["keystone_node"]) in list(reverse_dict.keys()):
                    value['keystone_node'] = reverse_dict[value["keystone_node"]]
                node_list.append(value['keystone_node'])
        keystone_df = keystone_df.sort_values(by=['layer', 'PR_score'], ascending=[True, False])
        keystone_df.to_csv(opath2, sep='\t', index=None)
        opath3 = os.path.join(cluster_dir, 'itol_conf.txt')
        with open(opath3, 'w') as fp:
            # print(node_list)
            node_list[-1] = node_list[-1].replace('_', '-')
            fp.write(itol_util.branch_style(node_list, node_list[-1], leaf_list))
        
        seat_Z = v['seat']
        fr = v['level_result'][0]['fr']

        plt.figure(figsize = (18, 15))
        sns.clustermap(fr,
                #row_linkage=seat.Z_,
                row_linkage=seat_Z,
                #col_linkage=seat.Z_,
                col_linkage=seat_Z,
                cmap='YlGnBu',
                cbar_pos = (1, 0.5, 0.05, 0.3),
                yticklabels=False,
                xticklabels=False )
        plt.savefig(os.path.join(cluster_dir, "heatmap.jpg"))
        plt.close()
        
def main(ori_GCN, ori_profile, outdir, params, d_df, cluster_profiles=None):
    if not cluster_profiles:
        # step1: sample cluster
        n_cluster, score, cluster_labels = sample_cluster.sample_cluster(ori_profile, params)
        # step2: split samples according to cluster
        cluster_profiles = abd_profile.split_to_clusters(ori_profile, cluster_labels)
    # step3: compute each cluster result
    final_result = {}
    for k, cluster_profile in cluster_profiles.items():
        if cluster_profile.shape[0] < 2:
            continue
        final_result[k] = copy.deepcopy(cluster_process(cluster_profile, d_df, ori_GCN, params))
    post_process(cluster_profiles, final_result, outdir)
    return final_result, cluster_profiles