from copy import deepcopy as dc
import math

def add_parentheses(definition):
    if ',' in definition:
        if not (definition.startswith('(') and definition.endswith(')')):
            definition = '({})'.format(definition)
    return definition

def delete_optional(s):
    #print('before delete -: {}'.format(s))
    pos = s.find('-')
    while pos > -1:
        
        if pos < len(s) - 1:
            # If '-' is followed by 'K', delete the next 6 characters (which should be a KO number)

            # If '-' is followed by '(', delete until the corresponding ')'
            if s[pos+1] == '(':
                close_pos = s.find(')', pos)
                if close_pos > -1:
                    s = s[:pos] + s[close_pos+1:]
                else:
                    # If no closing ')' is found, break the loop to avoid an infinite loop.
                    break
            else:
                s = s[:pos] + s[pos+7:]
        pos = s.find('-', pos)
        #print(s)
    # Remove any double spaces that may have been created
    pos = s.find("  ")
    while pos > -1:
        s = s.replace("  ", ' ')
        pos = s.find("  ")
    s = s.strip()
    #print('after delete -: {}'.format(s))
    return s

def delete_optional_draft(s):
    # operation with - and -- : delete
    s = s.replace('--', '')
    # Remove any double spaces that may have been created
    pos = s.find("  ")
    while pos > -1:
        s = s.replace("  ", ' ')
        pos = s.find("  ")
    s = s.strip()

    pos = s.find('-')
    while pos > -1:
        if pos < len(s) - 1:
            s = s[:pos] + s[pos+7:]
        pos = s.find('-', pos)
    # Remove any double spaces that may have been created
    pos = s.find("  ")
    while pos > -1:
        s = s.replace("  ", ' ')
        pos = s.find("  ")
    s = s.strip()
    return s

def max_path(path_dicts):
    max_path = path_dicts[0]
    max_ratio = max_path['exist']/ max_path['length']
    for i in range(len(path_dicts)):
        current_path = path_dicts[i]
        current_ratio = current_path['exist']/ current_path['length']
        if current_ratio > max_ratio:
            max_path = dc(current_path)
        elif current_ratio == max_ratio:
            if current_path['length'] < max_path['length']:
                max_path = dc(current_path)
        max_ratio = max_path['exist']/ max_path['length']
    return max_path


# {module_num:, exist_num: }
def completeness(str, qlist, node_dict):
    # print('*compute completeness for {}'.format(str))
    # if only one KO
    if len(str) == 6:
        if str[0] == 'K':
            # if only one KO 
            if str in qlist:
                return {"length":1, "exist":1}
            else:
                return {"length":1, "exist":0}
        else:
            # if only one NO
            return completeness(node_dict[str], qlist, node_dict)
        
    # replace () with component id
    for k, v in node_dict.items():
        if str.find(v)>-1 and str != v:
            str = str.replace(v, k)    
    # print('**compute completeness for {}'.format(str))
    str = add_parentheses(str)
    # delete -
    # str = delete_optional(str)
    str = delete_optional_draft(str)

    #print(str)
    if str.startswith('(') and str.endswith(')'):
        # comma first and use max
        sub = str[1: -1].split(',')
        completeness_list = []
        for s in sub:
            if s.find('+')<0:
                completeness_list.append(dc(completeness(s, qlist, node_dict)))
            else:
                plus_completeness = 1
                plus_elements = s.split('+')
                for pe in plus_elements:
                    path_dict = completeness(pe, qlist, node_dict)
                    part_ratio = path_dict['exist']/path_dict['length']
                    plus_completeness *= math.floor(part_ratio)
                if plus_completeness == 1:
                    completeness_list.append({"length":1, "exist":1})
                else:
                    completeness_list.append({"length":1, "exist":0})
        #print('clist', str, completeness_list)
        max_path_record = max_path(completeness_list)
        #print('inclist', max_path_record)
    else:   
        # space first and use weight 
        sub = str.split(' ') 
        max_path_record = {"length":0, "exist":0}
        for s in sub:
            if s.find('+')<0:
                path_dict = completeness(s, qlist, node_dict)
                max_path_record['length'] += path_dict['length']
                max_path_record['exist'] += path_dict['exist']
            else:
                plus_completeness = 1
                plus_elements = s.split('+')
                for pe in plus_elements:
                    path_dict = completeness(pe, qlist, node_dict)
                    part_ratio = path_dict['exist']/path_dict['length']
                    plus_completeness *= math.floor(part_ratio)
                if plus_completeness == 1:
                    max_path_record['exist'] += 1
                    max_path_record['length'] += 1
                else:
                    max_path_record['length'] += 1
    #print(str, max_path_record)
    return dc(max_path_record)

def replace_branket(test, nid=None):
    if not nid:
        nid = 0
    node_dict = {}
    stack = []
    for c in test:
        if c == '(':
            stack.append(nid)
            # generate a new node
            k = "N%05d"% nid
            node_dict[k] = ''
            nid += 1
        for id in stack:
            k = "N%05d"% id
            node_dict[k] += c
        if c == ')':
            stack.pop()
    return node_dict, nid


def main(str, qlist, node_dict):
    for k, v in node_dict.items():
        if str.find(v) >= 0:
            str = str.replace(v, k, 1)
    return completeness(str, qlist, node_dict)


# replace module in module definition 
def replace_module(mdef_df, definition):
    mdict = {}
    idx = definition.find('M')
    while idx>-1:
        module = definition[idx: idx + 6]
        replace_str = add_parentheses(mdef_df.loc[module, 'def'])
        definition = definition.replace(module, replace_str)
        mdict[module] = replace_str
        idx = definition.find('M')
    return definition, mdict
