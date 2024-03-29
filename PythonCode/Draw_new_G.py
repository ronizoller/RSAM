import sys
sys.path.append('/anaconda3/lib/python3.7/site-packages')
from networkx.drawing.nx_agraph import graphviz_layout
import matplotlib.pyplot as plt
import networkx as nx
import dendropy as tr
import math
import utiles
import tree_operations_v1 as tree_operations
import inits_v1 as inits
import random
import draw

p1 = (['S'], None, False)
p2 = (None,False,False)

color = True
double_mode = p2[0] is not None
lables_flag = False
draw_marked = False
draw_S = False

ext = 'proteobacteria'
ext_G = 'proteobacteria'
path = '/Users/ronizoller/Google Drive (ronizo@post.bgu.ac.il)/COGS/COG3549/'

S_colors = {}
big_size = 3000
small_size = 7
number_of_douplications = 1
x_axis = 100
y_axis = 30


def number_of_scpecies_doup(G,old_sigma):
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
    return should_be_found

def draw_new_doup(marked_nodes, colors, sigma, new_G, G,old_sigma,k,TH_compare_subtrees, path, lables, glob,spec,pattern,size,evol,compare,number_of_fields,S_labels_table):
    print('Drawing new G...')
    plt.clf()
    special_colors = []
    should_be_found = number_of_scpecies_doup(G,old_sigma)
    print()
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
                i = i + 1
        index += 1
    labels1 = nx.get_node_attributes(new_G, 'label')
    pos1 = graphviz_layout(tree_to_draw, prog='dot')

    plt.figure(figsize=(x_axis, y_axis))

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
                if temp_name.find(should) != -1 and (not flag):
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
    print('Finished drawing new G.\n')

def draw_new_HT(marked_nodes, colors, sigma, new_G, G,old_sigma,k,TH_compare_subtrees, path, lables,glob,spec,pattern,size,evol,compare,number_of_fields,S_labels_table):
    print('Drawing new G...')
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
                i = i + 1
        index += 1
    labels1 = nx.get_node_attributes(new_G, 'label')
    pos1 = graphviz_layout(tree_to_draw, prog='dot')

    plt.figure(figsize=(x_axis, y_axis))

    nodes_color = []
    nodes_size = []
    if draw_marked:
        max_size = max([x[0] for u, x in marked_nodes.items()])
        if double_mode:
            max_size = max([x[0]+x[1] for u,x in marked_nodes.items()])
    exp_facor = math.log(size)

    for nd in new_G.nodes(data=True):
        if new_G.out_degree(nd[0]) == 0 and not new_G.in_degree(nd[0]) == 0 and color:
            name = list(S_labels_table.keys())[list(S_labels_table.values()).index(sigma[nd[1]['label']])]
            if colors[sigma[nd[1]['label']]] == 'red':
                nodes_color.append('red')
            else:
                nodes_color.append('grey')
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
    plt.savefig(path+'/figures/new_G_pattern='+pattern+'_'+ext+'.png')
    print('Finished drawing new G.\n')


if p2[0] is not None:
    pattern_name = '(' + str(p1) + '_' + str(p2) + ')_Double-Mode'
else:
    pattern_name = str(p1) + '_Single-Mode'
input = open(path + '/saved_data/marked_RSAM_finder_'+ext+'_pattern='+pattern_name+'.txt', 'r')
if draw_marked:
    marked_nodes = []
    for line in input:
        marked_nodes.append(eval(line))
    marked_nodes = marked_nodes[0]
else:
    marked_nodes = {}

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

G = tr.Tree.get_from_path(path + "/GeneTree(binary,"+ext_G+")_local.txt", schema="newick")
S = tr.Tree.get_from_path(path + "/phyliptree(binary,"+ext+").phy", schema="newick")
input = open(path + '/0/sigma0.0.txt', 'r')
sigma = []
for line in input:
    sigma.append(eval(line))
sigma = sigma[0]

input = open(path + '/0/colors0.0.txt', 'r')
colors = []
for line in input:
    colors.append(eval(line))
colors = colors[0]

G.prune_taxa_with_labels(tree_operations.remove_unsigma_genes(G, sigma, True))

S = utiles.init_internal_labels(S, 'x', sigma, path)
G = utiles.init_internal_labels(G, 'u', sigma, path)

G = tree_operations.collapse_edges(G)
S = tree_operations.collapse_edges(S)


S_labels_table, G_labels_table, sigma = inits.init_taxon_to_label_table(S, G, sigma)
sigma, old_sigma = inits.update_sigma(sigma, True, S_labels_table, G_labels_table)
G.prune_taxa_with_labels(tree_operations.remove_unsigma_genes(G, sigma, False))
colors, old_colors = inits.update_colors(S, colors, True)

new_G = nx.DiGraph()
new_G = tree_operations.copy_G(G,new_G)

if draw_S:
    draw.draw_tree(S, 'S', old_sigma, colors, sigma,path,color,x_axis,y_axis,lables_flag)
    quit()

if 'D' in p1[0]:
    draw_new_doup(marked_nodes, colors, sigma, new_G, G, old_sigma, 0, 0,
                            path, True, False, ext, pattern_name,
                             big_size, ['D'], True, 1,S_labels_table)
else:
    draw_new_HT(marked_nodes, colors, sigma, new_G, G, old_sigma, 0, 0,
                            path, True, False, ext, pattern_name,
                             big_size, ['HT'], True, 1,S_labels_table)

print('Number of leafs of S: '+str(tree_operations.number_of_leafs(S,'S')))
print('Number of leafs of G: '+str(tree_operations.number_of_leafs(G,'G')))
