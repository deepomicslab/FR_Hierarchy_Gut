# funcitons used to generate files for tree style

def branch_style(node_list, species, leaf_list):
    color = '#FF0000'
    line_size = 5
    selected_line_size = line_size + 3
    template = "TREE_COLORS\nSEPARATOR TAB\nDATA\n"
    first_line = "root\tclade\t#000000\tnormal\t{}\n".format(line_size)
    line_template = "{}\t{}\t{}\tnormal\t{}\n"
    s = template + first_line
    for n in node_list:
        line = line_template.format(n, 'branch', color, selected_line_size)
        s += line
    
    node_style_template = "{}\tlabel\t#000000\tnormal\t2\n"
    for leaf in leaf_list:
        node_style_line = node_style_template.format(leaf)
        s += node_style_line

    eigen_sp_line = "{}\t{}\t{}\tbold\t3\n".format(species, 'label', color)
    s += eigen_sp_line
    return s

def branch_style_color(node_list, color):
    line_size = 5
    selected_line_size = line_size + 3
    line_template = "{}\t{}\t{}\tnormal\t{}\n"
    s = ''
    for n in node_list:
        line = line_template.format(n, 'branch', color, selected_line_size)
        s += line
    return s

def branch_style_header(leaf_list):

    line_size = 5
    template = "TREE_COLORS\nSEPARATOR TAB\nDATA\n"
    first_line = "root\tclade\t#000000\tnormal\t{}\n".format(line_size)
    s = template + first_line

    node_style_template = "{}\tlabel\t#000000\tnormal\t2\n"
    for leaf in leaf_list:
        node_style_line = node_style_template.format(leaf)
        s += node_style_line
    return s

def pie_header(color_list, label_list):
    pie_setting = "DATASET_PIECHART\nSEPARATOR TAB\nDATASET_LABEL\tPE difference\n"
    pie_setting += "COLOR\t#000000\n"
    pie_setting += "FIELD_COLORS\t{}\n".format('\t'.join(color_list))
    pie_setting += "FIELD_LABELS\t{}\n".format('\t'.join(label_list))
    pie_setting += "BORDER_WIDTH\t1\n"
    pie_setting += "BORDER_COLOR\t#000000\n"
    pie_setting += "POLAR_AREA_DIAGRAM\t1\nDATA\n"
    return pie_setting
    

def pie_ann_line(node, value_list, width=30):
    # with border
    template = "{}\t1\t{}\t{}\n".format(node, width, '\t'.join(value_list))
    return template

def std_name(leave):
    contents = leave.split('__')
