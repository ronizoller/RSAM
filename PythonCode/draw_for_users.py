import sys
sys.path.append('/anaconda3/lib/python3.7/site-packages')
from networkx.drawing.nx_agraph import graphviz_layout

import matplotlib
matplotlib.use("agg")
from matplotlib import pyplot as plt
import networkx as nx
import math
import tree_operations_v1 as tree_operations
import random
import os


def number_of_scpecies_doup(G, old_sigma, number_of_douplications):
    leafs_names = {}
    should_be_found = []
    for u in G.postorder_node_iter():
        if tree_operations.is_a_leaf(u):
            specie = old_sigma[u.taxon.label]
            if specie in leafs_names:
                leafs_names[specie] += 1
            else:
                leafs_names.update({specie:1})
    for name,score in leafs_names.items():
        if score > number_of_douplications:
            should_be_found.append(name)
            should_be_found.append(name)
    return should_be_found


def draw_new_doup(marked_nodes, sigma, new_G, G, old_sigma, lables, pattern,size, S_labels_table, x_axis, y_axis,
                  draw_marked, double_mode, lables_flag, ext, color, number_of_douplications, path):
    plt.clf()
    special_colors = []
    should_be_found = number_of_scpecies_doup(G,old_sigma,number_of_douplications)
    for i in range(0, len(should_be_found)):
        special_colors.append(hex_code_colors())
    tree_to_draw = nx.DiGraph()
    index = 1
    for u in G.postorder_node_iter():
        tree_to_draw.add_node(index, label = u.label)
        if not tree_operations.is_a_leaf(u):
            child = u.child_nodes()
            i = 0
            while i < len(child):
                if len(child) >= i + 1:
                    tree_to_draw.add_edge(index, list(G.postorder_node_iter()).index(child[i]) + 1)
                i += 1
        index += 1
    labels1 = nx.get_node_attributes(new_G, 'label')
    pos1 = graphviz_layout(tree_to_draw, prog='dot')

    plt.figure(figsize=(int(x_axis), int(y_axis)))

    nodes_color = []
    nodes_size = []
    if draw_marked:
        max_size = max([x[0] for u,x in marked_nodes.items()])
    if double_mode and draw_marked:
       max_size = max([x[1]+x[0] for u,x in marked_nodes.items()])
    exp_facor = math.log(size)

    for nd in new_G.nodes(data=True):
        flag = False
        if new_G.out_degree(nd[0]) == 0 and not new_G.in_degree(nd[0]) == 0 and color:
            temp_name = list(S_labels_table.keys())[list(S_labels_table.values()).index(sigma[nd[1]['label']])]
            for should in should_be_found:
                if not flag and(temp_name.find(should) != -1 or temp_name.replace(' ','').find(should) != -1):
                    nodes_color.append(special_colors[should_be_found.index(should)])
                    nodes_size.append(200)
                    flag = True
        if nd[1]['label'] in marked_nodes and (not flag) and draw_marked:
            nodes_color.append('blue')
            temp_size = max(marked_nodes[nd[1]['label']][0],marked_nodes[nd[1]['label']][0])
            if double_mode:
                temp_size = marked_nodes[nd[1]['label']][0]+marked_nodes[nd[1]['label']][1]
            nodes_size.append(math.exp((temp_size/max_size)*exp_facor)+500)
        elif not flag:
            nodes_color.append('#FFFFFF')
            nodes_size.append(50)
    if lables:
        for r, l in labels1.items():
            for x in G.postorder_node_iter():
                if x.taxon != None:
                    if x.label == l:
                        l = l + "\n (" + str(x.taxon) + ")"
                        l = l+'\n'+str(old_sigma[x.taxon.label])
                        labels1.update({r: l})
    nx.draw(tree_to_draw, pos1, arrows=True, node_size=nodes_size, node_color=nodes_color,
            width=1)
    if lables_flag:
        text = nx.draw_networkx_labels(tree_to_draw, pos1, labels1, font_size=7)
        for _, t in text.items():
            t.set_rotation('vertical')
    plt.savefig(path + '/figures/new_G_pattern=' + pattern + '_' + ext + '.png')
    return


def draw_new_HT(marked_nodes, colors, sigma, new_G, G,old_sigma, lables,pattern,size,S_labels_table,x_axis,y_axis,draw_marked,double_mode,lables_flag,ext, color, path):
    plt.clf()

    tree_to_draw = nx.DiGraph()
    index = 1
    for u in G.postorder_node_iter():
        tree_to_draw.add_node(index, label = u.label)
        if not tree_operations.is_a_leaf(u):
            child = u.child_nodes()
            i = 0
            while i < len(child):
                if len(child) >= i + 1:
                    tree_to_draw.add_edge(index, list(G.postorder_node_iter()).index(child[i]) + 1)
                i += 1
        index += 1
    labels1 = nx.get_node_attributes(new_G, 'label')
    pos1 = graphviz_layout(tree_to_draw, prog='dot')

    plt.figure(figsize=(int(x_axis), int(y_axis)))

    nodes_color = []
    nodes_size = []
    if draw_marked:
        if marked_nodes != {}:
            max_size = max([x[0] for u, x in marked_nodes.items()])
            if double_mode:
                max_size = max([x[0]+x[1] for u,x in marked_nodes.items()])
        else:
            max_size = 200
    exp_facor = math.log(size)

    for nd in new_G.nodes(data=True):
        if new_G.out_degree(nd[0]) == 0 and not new_G.in_degree(nd[0]) == 0 and color:
            name = list(S_labels_table.keys())[list(S_labels_table.values()).index(sigma[nd[1]['label']])]
            if colors[sigma[nd[1]['label']]] == 'red':
                nodes_color.append('red')
            elif colors[sigma[nd[1]['label']]] == 'black':
                nodes_color.append('grey')
            elif colors[sigma[nd[1]['label']]] == 'pink':
                nodes_color.append('white')
            nodes_size.append(200)
        elif nd[1]['label'] in marked_nodes:
            nodes_color.append('blue')
            temp_size = max(marked_nodes[nd[1]['label']][0],marked_nodes[nd[1]['label']][0])
            if double_mode:
                temp_size = marked_nodes[nd[1]['label']][0] + marked_nodes[nd[1]['label']][1]
            nodes_size.append(math.exp((temp_size/max_size)*exp_facor)+500)
        else:
            nodes_color.append('white')
            nodes_size.append(200)
    if lables:
        for r, l in labels1.items():
            for x in G.postorder_node_iter():
                if x.taxon:
                    if x.label == l:
                        l = l + "\n (" + str(x.taxon) + ")"
                        l = l+'\n'+str(old_sigma[x.taxon.label])
                        labels1.update({r: l})

    nx.draw(tree_to_draw, pos1, arrows=True, node_size=nodes_size, node_color=nodes_color,
            width=1)
    if lables_flag:
        text = nx.draw_networkx_labels(tree_to_draw, pos1, labels1, font_size=7)
        for _, t in text.items():
            t.set_rotation('vertical')
    plt.savefig(path + '/figures/new_G_pattern='+pattern+'_'+ext+'.png')
    return

def hex_code_colors():
    a = hex(random.randrange(0,256))
    b = hex(random.randrange(0,256))
    c = hex(random.randrange(0,256))
    a = a[2:]
    b = b[2:]
    c = c[2:]
    if len(a)<2:
        a = "0" + a
    if len(b)<2:
        b = "0" + b
    if len(c)<2:
        c = "0" + c
    z = a + b + c
    return "#" + z.upper()


def main(G,sigma,old_sigma,colors,S_labels_table,p1,p2,ext,color,lables_flag,draw_marked,x_axis,y_axis, path, res, number_of_duplications):
    p1 = (p1[0],p1[1],p1[2])
    new_G = nx.DiGraph()
    new_G = tree_operations.copy_G(G,new_G)

    double_mode = p2[0] is not None

    if p2[0] is not None:
        pattern_name = '(' + str(p1) + '_' + str(p2) + ')_Double-Mode'
    else:
        pattern_name = str(p1) + '_Single-Mode'
    try:
        input = open(path + '/saved_data/marked_RSAM_finder_' + ext + '_pattern=' + pattern_name + '.txt', 'r')
    except:
        res['error'] += '/saved_data/marked_RSAM_finder_' + ext + '_pattern=' + pattern_name + ".txt was not found."
        return
    if draw_marked:
        marked_nodes = []
        for line in input:
            marked_nodes.append(eval(line))
        marked_nodes = marked_nodes[0]
    else:
        marked_nodes = {}

    big_size = 20

    to_create = path + '/figures/'
    os.makedirs(os.path.dirname(to_create), exist_ok=True)

    if 'D' in p1[0]:
        draw_new_doup(marked_nodes, sigma, new_G, G, old_sigma,
                      lables_flag, pattern_name, big_size,S_labels_table,x_axis, y_axis,
                      draw_marked, double_mode, lables_flag, ext, color, number_of_duplications, path)

    else:
        draw_new_HT(marked_nodes, colors, sigma, new_G, G, old_sigma,
                    lables_flag, pattern_name,
                                 big_size,S_labels_table, x_axis,
                    y_axis, draw_marked, double_mode, lables_flag, ext, color, path)

    return

if __name__ == "__main__":
    main()
