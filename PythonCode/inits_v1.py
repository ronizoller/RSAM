import math
import tree_operations
import networkx as nx
import random
import EfficiantVersion as effi

def find_max_S_d(max_S_d, S_dis_matrix):
    for h, v in S_dis_matrix.items():
        if v > max_S_d:
            max_S_d = v
    return max_S_d

def init_leafs_efficient(S,G, H, k, H_number_of_nodes,sigma,nodes_table,subtree):
    print('Initialasing efficient hypergraph leafs...')
    for leaf in G.leaf_nodes():
        if not tree_operations.isolated(leaf):
            H.add_node(H_number_of_nodes,s=leaf.label,t=sigma[leaf.label],l=list())
            big_node = H_number_of_nodes
            H_number_of_nodes += 1
            for i in range(0,k):
                if i==0:
                    cost=0
                else:
                    cost=math.inf
                new_item = {'s':leaf.label,'t':sigma[leaf.label],'cost':cost,'event':"leaf",'list_place':i}
                H.nodes[big_node]['l'].insert(i,new_item)
            nodes_table[sigma[leaf.label]][leaf.label] = big_node
    print('Finished initialasing hypergraph leafs.\n')
    return H,H_number_of_nodes, nodes_table,subtree


def init_distance_S(S_dis_matrix, k, test, path,spe):
    print('Initing Distances...')

    if test:
        print("     Reading file 'S_dist_matrix"+".txt'...")
        input = open(path+'/saved_data/S_dist_matrix'+'.txt', 'r')
        S_disr_matrix = []
        for line in input:
            S_disr_matrix.append(eval(line))
        S_dis_matrix = S_disr_matrix[0]
        print("     Finished reading file 'S_dist_matrix.txt'.")
    else:
        input = open(path+'/saved_data/S_edgelist_'+spe+'.txt', 'r')
        edgelist = []
        for line in input:
            edgelist.append(eval(line))
        edgelist = edgelist[0]
        networkx_S = nx.parse_edgelist(edgelist, nodetype = str, data=(('weight',any),))
        for nd1 in networkx_S.nodes():
            for nd2 in networkx_S.nodes():
                #S_dis_matrix.update({(nd1,nd2):nx.shortest_path_length(networkx_S, source=nd1, target=nd2, weight='weight')})      #for weighted tree
                S_dis_matrix.update(
                    {(nd1, nd2): nx.shortest_path_length(networkx_S, source=nd1, target=nd2)})
        print('     Writing S_edge_list...')
        file = open(path+'/saved_data/S_dist_matrix'+'.txt', 'w')
        file.write(str(S_dis_matrix))
        file.close()
        print('     Finished writing S_dist_matrix.\n')
    return S_dis_matrix
    print('Finished initing Distances...')

def init_nodes_table(S,G,nodes_table):
    print('Initilasinig node table...')
    for s_nd in S.postorder_node_iter():
        nodes_table[s_nd.label] = {}
        for g_nd in G.postorder_node_iter():
            nodes_table[s_nd.label].update({g_nd.label:-1})
    return nodes_table
    print('Finished initilasinig node table.\n')

def init_leafs(G, H, k, H_number_of_nodes, sigma,nodes_table):
    print('Initialasing hypergraph leafs...')
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
    print('Finished initialasing hypergraph leafs.\n')
    return H,H_number_of_nodes, nodes_table

def init_taxon_to_label_table(S,G,sigma):
    S_labels_table = {}
    G_labels_table = {}
    for leaf_S in S.leaf_nodes():
        S_labels_table.update({leaf_S.taxon.label:leaf_S.label})
    for leaf_G in G.leaf_nodes():
        G_labels_table.update({leaf_G.taxon.label: leaf_G.label})
    return S_labels_table,G_labels_table

def update_sigma(S, G, k, sigma, test, path,exect_names,S_labels_table,G_labels_table):
    print('Updating sigma...')
    old_sigma = sigma
    if test:
        print("     Reading file 'sigma.txt'...")
        input = open(path+'/saved_data/sigma'+'.txt', 'r')
        new_sigma = []
        for line in input:
            new_sigma.append(eval(line))
        sigma = new_sigma[0]
        print("     Finished reading file 'sigma.txt'.")
    else:
        for u, x in sigma.items():
            if (u in G_labels_table) and (x in S_labels_table):
                if exect_names:
                    sigma = dict((u1, x1) for u1, x1 in sigma.items() if not (
                    (S_labels_table[x]).replace("_", " ") == x1 and (G_labels_table[u]).replace("_"," ") == u1))  # if lables are strings
                    sigma.update({G_labels_table[u]: S_labels_table[x]})
                else:
                    sigma = dict((u,x) for u,x in sigma.items() if not (((S_labels_table[x]).replace("_"," ").find(x) != -1 or x.find((S_labels_table[x]).replace("_"," ")) != -1) and (((G_labels_table[u]).replace("_"," ").find(u) != -1) or (u.find(G_labels_table[u]).replace("_"," ")) != 1)))      #if lables are strings
                    sigma.update({G_labels_table[x]: S_labels_table[u]})
            else:
                print ('        Couldnt find match for: '+u+' and '+x)
                sigma = dict((u1, x1) for u1, x1 in sigma.items() if not (u1 == u and x1 == x))     #remove unsigma mappings
                old_sigma = dict((u1, x1) for u1, x1 in old_sigma.items() if not (u1 == u and x1 == x))
        #print('     Writing sigma...')
        #file = open(path+'/saved_data/sigma.txt', 'w')
        #file.write(str(sigma))
        #file.close()
        #print('     Finished writing sigma.\n')
    print('Finished updating sigma.\n')
    return sigma,old_sigma

def update_colors(S,colors,exact_names):
    print('Updating colors...')
    old_colors = colors.copy()
    for leaf_S in S.leaf_nodes():
        degel = False
        for x,color in colors.items():
            if not leaf_S.taxon is None:
                #print(str(leaf_S.taxon.label) + '  **  ' + str(x))
                if not exact_names:
                    if x.replace("_","").find(leaf_S.taxon.label) != -1:
                        #del colors[int(leaf_S.taxon.label)]         #if lables are numbers
                        colors = dict((x, color) for x,color in colors.items() if (x.replace("_","").find(leaf_S.taxon.label) == -1 or leaf_S.taxon.label .find(x.replace("_","")) == -1))
                        colors.update({leaf_S.label: color})
                        degel = True
                else:
                    if x.replace("_","") == leaf_S.taxon.label:
                        #del colors[int(leaf_S.taxon.label)]         #if lables are numbers
                        colors = dict((x, color) for x,color in colors.items() if (x != leaf_S.taxon.label))
                        colors.update({leaf_S.label: color})
                        degel = True
        if degel == False:
            print('     No lable was assigned to: '+str(leaf_S))
            colors = dict((x, color) for x, color in colors.items() if not x == leaf_S.taxon.label)
            colors.update({leaf_S.label:  random.choice(['black', 'red'])})
            old_colors.update({leaf_S.taxon.label:  random.choice(['black', 'red'])})
    print('Finished updating colors.\n')
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

def init_dict_inf(H,S,G,k,nodes_table,dict,sigma):
    res = {}
    for u in G.postorder_node_iter():
        res.update({u.label:{}})
        for x in S.postorder_node_iter():
            if dict == 'subtree':
                if effi.find_nodes_in_hypergraph(H, u.label, x.label, -1, nodes_table) != []:
                    res[u.label].update({x.label:[effi.find_nodes_in_hypergraph(H, u.label, x.label, i, nodes_table)[0] for i in range (0,k)]})
                    S.find_node(lambda n: (n.label == sigma[u.label]))
                elif u.is_leaf() and (x in S.find_node(lambda n: (n.label == sigma[u.label])).ancestor_iter()):
                    res[u.label].update({x.label: [effi.find_nodes_in_hypergraph(H, u.label, sigma[u.label], i, nodes_table)[0]
                                                   for i in range(0, k)]})
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
