import re
import copy
from pyseat.SEAT import SEAT

# Use seat to make tree
def make_tree(fr_df):
    seat = SEAT(affinity="gaussian_kernel",
            sparsification="knn_neighbors",
            objective="SE",
            strategy="bottom_up")
    seat.fit(fr_df)
    return seat

# Seat default node name reflect to species
def name_reflection(fr_d):
    name_dict = {}
    reverse_dict = {}
    for i, old_name in enumerate(list(fr_d.columns)):
        new_name = "n{}".format(i)
        name_dict[old_name] = new_name
        reverse_dict[new_name] = old_name
    return name_dict, reverse_dict

# covert newick to json tree structure
def parse(newick):
    tokens = re.finditer(r"([^:;,()\s]*)(?:\s*:\s*([\d.]+)\s*)?([,);])|(\S)", newick+";")

    def recurse(nextid = 0, parentid = -1): # one node
        thisid = nextid;
        children = []

        name, length, delim, ch = next(tokens).groups(0)
        if ch == "(":
            while ch in "(,":
                node, ch, nextid = recurse(nextid+1, thisid)
                children.append(node)
            name, length, delim, ch = next(tokens).groups(0)
        node = {"id": thisid, "name": name, "length": float(length) if length else None, 
                "parentid": parentid, "children": children, "taken": False}
        return copy.deepcopy(node), delim, nextid

    return recurse()[0]

# find parent of each node/leaf
def parents(node, parent_dict):
    for child in node['children']:
        node_name = child['name']
        parent_dict[node_name] = node['name']
        parents(child, parent_dict)


# compute the leaves and level of each node, and clade deepth
def recu_compute(node, level=0, largest = {'largest': 0}):
    node['leaves'] = []
    node['level'] = level
    if len(node['children']) == 0:
        node['clade_depth'] = 0
        if level > largest['largest']:
            largest['largest'] = level
        node['leaves'] = [node['name']]
        return [node['name']], largest
    else:
        largest_children_depth = 0
        for n in node['children']:
            plus_leaves, largest = recu_compute(n, level+1, largest)
            node['leaves'] += plus_leaves
            if n['clade_depth'] > largest_children_depth:
                largest_children_depth = n['clade_depth']
        node['clade_depth'] = largest_children_depth + 1
        return node['leaves'], largest

# init the dict of the nodes in each level
def make_level_dict(largest_level):
    result_dict = {}
    for i in range(largest_level):
        result_dict[i] = {}
    return result_dict

def make_layer_dict(nlayer):

    result_dict = {}
    for i in range(1, nlayer):
        
        result_dict[i] = {}
    return result_dict

# find layer node:leaves
def recu_layer(node, result_dict):
    depth = node['clade_depth']
    if depth in result_dict.keys():
        result_dict[depth][node['name']] = node['leaves']
    for child in node['children']:
        recu_layer(child, result_dict)

# add pseudo internal node in each layer
def to_layer_leaves(result_dict, nlayer):
    for i in range(1, nlayer-1):
        top_id = i + 1
        top_large_list = to_large_list(result_dict[top_id])
        for k, v in result_dict[i].items():
            if not has_common(v, top_large_list):
                # print("add node : {} from level {} to level {}".format(k, i, top_id))
                result_dict[top_id][k] = v

# find the node:leaves in each level
def recu_level(node, largest_level, result_dict):
    if node['level'] < largest_level:
        result_dict[node['level']][node['name']] = node['leaves']
        for n in node['children']:
            recu_level(n, largest_level, result_dict)

# add pseudo internal node in each level
def to_level_leaves(result_dict, largest_level):
    for i in range(largest_level-2):
        next_id = i + 1
        large_list = to_large_list(result_dict[next_id])
        for k, v in result_dict[i].items():
            if len(v) == 0:
                result_dict[next_id][k] = k
            elif not has_common(v, large_list):
                # print("add node : {} from level {} to level {}".format(k, i, next_id))
                result_dict[next_id][k] = v    

# check if belongs to upper level
def has_common(lower, upper):
    if len(list(set(upper).intersection(set(lower))))>0:
        return True
    else:
        return False

# check the leaves in the level
def to_large_list(level2_json):
    large_list = []
    for k, v in level2_json.items():
        large_list += v
    return large_list

# convert json tree to newick tree
def call_newick(node):
    node['taken'] = True
    if len(node['children']) == 0:
        return node['name']
    else:
        s = "("
        temp_list = []
        for n in node['children']:
            newick_tree = call_newick(n)
            temp_list.append(newick_tree)
        s += ','.join(temp_list)
        s += '){}:1'.format(node['name'])
        return s

# generate newick for the depth-limited tree
def depth_limit_newick(node, depth_limit, top_node_list = []):
    if node['clade_depth'] < depth_limit:
        if not node['taken']:
            str = call_newick(node)
            top_node_list.append(str)
    else:
        for n in node['children']:
            depth_limit_newick(n, depth_limit, top_node_list)
    return top_node_list

def limit_newick_last(top_node_list):
    s = '('
    s += ','.join(top_node_list)
    s += ')root;'
    return s

# convert children dict to json tree
def call_tree(children_dict, rid):
    node_json = {'name': "n{}".format(rid), 'children': []}
    if len(children_dict[rid]) == 0:
        return node_json
    else:
        for cid in children_dict[rid]:
            cnode = call_tree(children_dict, cid)
            node_json['children'].append(copy.deepcopy(cnode))
        return node_json


# make tree without leaves
def call_internal_tree(node, leaves_set):
    node['taken'] = True
    if len(node['children']) == 0:
        return node['name']
    else:
        children = []
        for c in node['children']:
            children.append(c['name'])
    if set(children).issubset(leaves_set):
        return node['name']
    else:
        s = "("
        temp_list = []
        for n in node['children']:
            newick_tree = call_internal_tree(n, leaves_set)
            temp_list.append(newick_tree)
        s += ','.join(temp_list)
        s += '){}:1'.format(node['name'])
        return s

def rename_node(node, rename_dict):
    if node['name'] in rename_dict.keys():
        node['name'] = rename_dict[node['name']]
    for c in node['children']:
        rename_node(c, rename_dict)

def sort_children(node, sort_rule):
    if len(node['children']) > 0:
        node['children'] = sorted(node['children'], key=sort_rule)
        for c in node['children']:
            sort_children(c, sort_rule)
