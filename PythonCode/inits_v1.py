import math
import tree_operations_v1 as tree_operations
import networkx as nx
import random
import EfficiantVersion as effi
from copy import deepcopy
import utiles

def find_max_S_d(max_S_d, S_dis_matrix):
    for h, v in S_dis_matrix.items():
        if v > max_S_d:
            max_S_d = v
    return max_S_d


def init_leafs_efficient(G, H, k, H_number_of_nodes, sigma, nodes_table):
    for leaf in G.leaf_nodes():
        if not tree_operations.isolated(leaf):
            target = ''
            if leaf.label in sigma:
                target = sigma[leaf.label]
            else:
                for species in sigma:
                    if leaf.label.find(species) != -1 or species.find(leaf.label) != -1:
                        target = sigma[species]
            H.add_node(H_number_of_nodes, s=leaf.label, t=target,l=list())
            big_node = H_number_of_nodes
            H_number_of_nodes += 1
            for i in range(0,k):
                if i == 0:
                    cost_no_losses = 0
                    cost_with_losses = 0
                else:
                    cost_no_losses = math.inf
                    cost_with_losses = math.inf
                new_item = {'s':  leaf.label,  't':  target,'cost_without_losses': cost_no_losses, 'cost_with_losses': cost_with_losses, 'event':"leaf", 'list_place':  i}
                H.nodes[big_node]['l'].insert(i,new_item)
            nodes_table[target][leaf.label] = big_node
    return H,H_number_of_nodes, nodes_table


def init_distance_S(S_dis_matrix,path,spe):
    input = open(path+'/saved_data/S_edgelist.txt', 'r')
    edgelist = []
    for line in input:
        edgelist.append(eval(line))
    edgelist = edgelist[0]
    networkx_S = nx.parse_edgelist(edgelist, nodetype = str, data=(('weight',any),))
    for nd1 in networkx_S.nodes():
        for nd2 in networkx_S.nodes():
            S_dis_matrix.update(
                {(nd1, nd2): nx.shortest_path_length(networkx_S, source=nd1, target=nd2)})
    file = open(path+'/saved_data/S_dist_matrix'+'.txt', 'w')
    file.write(str(S_dis_matrix))
    file.close()
    return S_dis_matrix


def init_nodes_table(S,G,nodes_table):
    for s_nd in S.postorder_node_iter():
        nodes_table[s_nd.label] = {}
        for g_nd in G.postorder_node_iter():
            nodes_table[s_nd.label].update({g_nd.label:-1})
    return nodes_table


def init_leafs(G, H, k, H_number_of_nodes, sigma,nodes_table):
    for leaf in G.leaf_nodes():
        if not tree_operations.isolated(leaf):
            H.add_node(H_number_of_nodes,s=leaf.label,t=sigma[leaf.label],l=list())
            big_node = H_number_of_nodes
            H_number_of_nodes += 1
            i = 0
            for i in range(0,k):
                if i==0:
                    cost=0
                else:
                    cost=math.inf
                new_item = {'s':leaf.label,'t':sigma[leaf.label],'cost':cost,'event':"leaf",'list_place':i}
                H.nodes[big_node]['l'].insert(i,new_item)
            nodes_table[sigma[leaf.label]][leaf.label] = big_node
    return H,H_number_of_nodes, nodes_table


def init_taxon_to_label_table(S, G, sigma):
    S_labels_table = {}
    G_labels_table = {}

    for leaf_S in S.leaf_nodes():
        if leaf_S.taxon:
            S_labels_table.update({leaf_S.taxon.label: leaf_S.label})
            S_labels_table.update({leaf_S.taxon.label.replace(' ',''): leaf_S.label})
            S_labels_table.update({leaf_S.taxon.label.replace(' ','_'): leaf_S.label})
    for leaf_G in G.leaf_nodes():
        if leaf_G.taxon:
            G_labels_table.update({leaf_G.taxon.label.replace(' ','_'): leaf_G.label})
            G_labels_table.update({leaf_G.taxon.label.replace(' ',''): leaf_G.label})
            G_labels_table.update({leaf_G.taxon.label: leaf_G.label})
        else:
            sigma.update({leaf_G.label:'x'+leaf_G.label[1:]})
    return S_labels_table,G_labels_table,sigma


def update_sigma(sigma, S_labels_table, G_labels_table):
    old_sigma = sigma.copy()
    sigma = {}
    for u, x in old_sigma.items():
        gene = utiles.is_prefix_of(u, G_labels_table)
        species = utiles.is_prefix_of(x, S_labels_table)
        if gene and species:
            sigma = dict((u1, x1) for u1, x1 in sigma.items() if not (S_labels_table[species].replace("_", " ") == x1 and (G_labels_table[gene]).replace("_"," ") == u1))
            sigma.update({G_labels_table[gene]: S_labels_table[species]})
        else:
            sigma = dict((u1, x1) for u1, x1 in sigma.items() if not (u1 == u and x1 == x))     #remove unsigma mappings
            old_sigma = dict((u1, x1) for u1, x1 in old_sigma.items() if not (u1 == u and x1 == x))

    return sigma,old_sigma


def update_colors(S,colors):
    old_colors = colors.copy()
    for leaf_S in S.leaf_nodes():
        if not leaf_S.taxon is None:
            degel = False
            for x,color in old_colors.items():
                    if utiles.find_in_substreing_different_forms(x,leaf_S.taxon.label):
                        colors = dict((x, color) for x,color in colors.items() if not utiles.find_in_substreing_different_forms(x,leaf_S.taxon.label))
                        colors.update({leaf_S.label: color})
                        degel = True
            if not degel:
                colors = dict((x, color) for x, color in colors.items() if not x == leaf_S.taxon.label)
                #colors.update({leaf_S.label:  random.choice(['black', 'red'])})
                colors.update({leaf_S.label: 'pink'})
                old_colors.update({leaf_S.taxon.label: 'pink'})
    return colors,old_colors


def init_H_field(H, field, init,in_list, for_edge):
    if not for_edge:
        for nd in H.nodes(data=True):
            if in_list:
                for i in range(0,len(nd[1]['l'])):
                    if field in nd[1]['l'][i]:
                        nd[1][field] = init
            else:
                if field in nd[1]:
                    nd[1][field] = init
    else:
        for e in H.edges(data=True):
            if field in e[2]:
                e[2][field] = init
    return H


def init_dict_inf(H,S,G,k,nodes_table,dict,sigma,S_dis_matrix,loss_cost):
    res = {}
    for u in G.postorder_node_iter():
        res.update({u.label:{}})
        for x in S.postorder_node_iter():
            if dict == 'subtree' or dict == 'subtreeLoss':
                if effi.find_nodes_in_hypergraph(H, u.label, x.label, -1, nodes_table):
                    res[u.label].update({x.label:[effi.find_nodes_in_hypergraph(H, u.label, x.label, i, nodes_table)[0] for i in range (0,k)]})
                    S.find_node(lambda n: (n.label == sigma[u.label]))
                elif u.is_leaf() and (x in S.find_node(lambda n: (n.label == sigma[u.label])).ancestor_iter()):
                    if dict == 'subtree':
                        res[u.label].update({x.label: [effi.find_nodes_in_hypergraph(H, u.label, sigma[u.label], i, nodes_table)[0]
                                                       for i in range(0, k)]})
                    else:
                        loss_cost_new = S_dis_matrix[(x.label, sigma[u.label])] * loss_cost
                        list_of_values = [deepcopy(effi.find_nodes_in_hypergraph(H, u.label, sigma[u.label], i, nodes_table)[0])
                                          for i in range(0, k)]
                        for i in range(0, k):
                            list_of_values[i][1]['cost_with_losses'] += loss_cost_new
                        res[u.label].update({x.label: list_of_values})
                else:
                    res[u.label].update({x.label: []})
            else:
                res[u.label].update({x.label: []})
    return res

def init_dic (list_to_be_init, init):
    res = {}
    for item in list_to_be_init:
        res.update({item.label:init})
    return res
