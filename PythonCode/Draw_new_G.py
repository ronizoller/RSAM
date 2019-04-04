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

ext = 'epsilon_delta'
S_colors = {}
big_size = 2000
small_size = 7
should_be_found = ['Geopsychrobacterelectrodiphilus','Desulforegulaconservatrix','Arcobactercibarius','Desulfobacteriumautotrophicum','Geobacterlovleyi','Geobacterdaltonii','Desulforegulaconservatrix','Sulfurimonasdenitrificans','Geobactersp.M18','Desulfobaculasp.TS','Helicobacterbizzozeronii']
path = '/Users/ronizoller/PycharmProjects/TreeReconciliation/trees/COG3550'
def draw_new_G2(marked_nodes, colors, sigma, new_G, G,old_sigma,k,TH_compare_subtrees, path, lables, glob,spec,pattern,size,evol,compare,number_of_fields,S_labels_table,special_colors):
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

    plt.figure(figsize=(70, 10))

    nodes_color = []
    nodes_size = []
    max_size = utiles.map_max(marked_nodes,number_of_fields)
    exp_facor = math.log(size)

    for nd in new_G.nodes(data=True):
        flag = False
        if new_G.out_degree(nd[0]) == 0 and not new_G.in_degree(nd[0]) == 0:
            temp_name = list(S_labels_table.keys())[list(S_labels_table.values()).index(sigma[nd[1]['label']])]
            for should in should_be_found:
                if temp_name.find(should) != -1 and (not flag):
                    nodes_color.append(special_colors[should_be_found.index(should)])
                    nodes_size.append(400)
                    flag = True
        if nd[1]['label'] in marked_nodes and (not flag):
            nodes_color.append('blue')
            temp_size = max(marked_nodes[nd[1]['label']][0][0],marked_nodes[nd[1]['label']][0][1])
            nodes_size.append(math.exp((temp_size/max_size)*exp_facor)+500)
        elif not flag:
            nodes_color.append('white')
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
    nx.draw_networkx_labels(tree_to_draw, pos1, labels1, font_size=7)
    plt.savefig(path+'/figures/new_G_k='+str(k)+'_TH_compare_subtrees='+str(TH_compare_subtrees)+'_species='+spec+'_pattern='+pattern+"_"+str(evol)+"_comapre="+str(compare)+'.png')
    print('Finished drawing new G.\n')

input = open(path + '/saved_data/marked_RSAM_finder.txt', 'r')
marked_nodes = []
for line in input:
    marked_nodes.append(eval(line))
marked_nodes = marked_nodes[0]

special_colors = ['violet','darkorange','crimson','purple','lightskyblue','yellow','gold','navy','black','green','silver']


G = tr.Tree.get_from_path(path + "/GeneTree(binary)_local.txt", schema="newick")
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
sigma, old_sigma = inits.update_sigma(S, G, 0, sigma, False, path, True, S_labels_table, G_labels_table)
G.prune_taxa_with_labels(tree_operations.remove_unsigma_genes(G, sigma, False))
colors, old_colors = inits.update_colors(S, colors, True)

new_G = nx.DiGraph()
new_G = tree_operations.copy_G(G,new_G)
draw_new_G2(marked_nodes, colors, sigma, new_G, G, old_sigma, 0, 0,
                            path, True, False, ext, 'same_color',
                             big_size, ['D'], True, 1,S_labels_table,special_colors)