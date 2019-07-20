import sys
sys.path.append('/anaconda3/lib/python3.6/site-packages')

import tree_operations_v1 as tree_operations
import heapq
import utiles
import math
import networkx as nx
import inits_v1 as inits
import random
from copy import deepcopy



def build_hyper_garph(S, G, k,nodes_table, D_cost, S_cost,loss_cost, HT_cost, path, sigma, save_data,S_dis_matrix):
    print('Building hypergraph...')
    H = nx.MultiDiGraph()
    H.clear()

    H, H_number_of_nodes, nodes_table = inits.init_leafs_efficient(S,G, H, k, 0,sigma,nodes_table)
    incomp = inits.init_dict_inf(H, S, G, k, nodes_table, 'incomp', sigma, S_dis_matrix, loss_cost)
    subtree = inits.init_dict_inf(H, S, G, k, nodes_table, 'subtree', sigma, S_dis_matrix, loss_cost)
    subtreeLoss = inits.init_dict_inf(H, S, G, k, nodes_table, 'subtreeLoss', sigma, S_dis_matrix, loss_cost)
    for u in G.postorder_node_iter():
        if not tree_operations.is_a_leaf(u):
            for x in S.postorder_node_iter():
                key_counter = 0
                S_list = SpeciationEvent_effi(u, x, S_cost, subtreeLoss,  k,  H,  nodes_table)
                SE = list(map(lambda nd: (nd,'S'), S_list))
                D_list = DuplicationEvent_effi(H, u, x, D_cost, loss_cost, nodes_table, subtreeLoss, k)
                DE = list(map(lambda nd: (nd,'D'),D_list))
                HT_list = HTEvent_effi(u, x, HT_cost, subtreeLoss, incomp, k, H, nodes_table)
                HTE = list(map(lambda nd: (nd,'HT'),HT_list))
                Kbest = SE + DE + HTE
                random.shuffle(Kbest,random.random)
                heapq.heapify(Kbest)

                if len(Kbest) > 0:
                    H.add_node(H_number_of_nodes, s=u.label, t=x.label,l=[])
                    nodes_table[x.label][u.label] = H_number_of_nodes
                    big_node = H_number_of_nodes
                    H_number_of_nodes += 1

                for i in range(0,k):
                    if len(Kbest) > 0:
                        match = heapq.heappop(Kbest)
                        new_cost_no_losses = match[0].val[0]
                        new_cost_with_losses = match[0].val[1]
                        event = match[1]
                        match1 = match[0].val[2]
                        match2 = match[0].val[3]

                        if event == 'S':
                            H.nodes[big_node]['l'] = H.nodes[big_node]['l'] + [{'s': u.label, 't': x.label, 'cost_without_losses':new_cost_no_losses,
                                                                                'cost_with_losses': new_cost_with_losses, 'event': "S",'list_place':len(H.nodes[big_node]['l'])}]
                        elif event == 'D':
                            H.nodes[big_node]['l'] = H.nodes[big_node]['l'] + [{'s': u.label, 't': x.label, 'cost_without_losses':new_cost_no_losses,
                                                                                'cost_with_losses': new_cost_with_losses, 'event': "D",'list_place':len(H.nodes[big_node]['l'])}]
                        elif event == 'HT':
                            H.nodes[big_node]['l'] = H.nodes[big_node]['l'] + [{'s': u.label, 't': x.label, 'cost_without_losses':new_cost_no_losses,
                                                                                'cost_with_losses': new_cost_with_losses, 'event': "HT",'list_place':len(H.nodes[big_node]['l'])}]

                        H.add_edge(match1[0], big_node, key=key_counter, source=match1[1]['list_place'],target=len(H.nodes[big_node]['l'])-1,probability=0)
                        key_counter += 1
                        H.add_edge(match2[0], big_node, key=key_counter, source=match2[1]['list_place'], target=len(H.nodes[big_node]['l'])-1,probability=0)
                        key_counter += 1

                        new_node1 = find_nodes_in_hypergraph(H, match1[1]['s'], match1[1]['t'], match1[1]['list_place'] + 1, nodes_table)
                        new_node2 = find_nodes_in_hypergraph(H, match2[1]['s'], match2[1]['t'],match2[1]['list_place'] + 1, nodes_table)
                        if new_node1 and new_node2:
                            new_node1 = new_node1[0]
                            new_node2 = new_node2[0]
                            if event == 'D':
                                additional_cost = D_cost
                            elif event == 'S':
                                additional_cost = S_cost
                            else:
                                additional_cost = HT_cost
                            cost_with_losses = match1[1]['cost_with_losses'] + new_node2[1]['cost_with_losses'] + additional_cost
                            cost_without_losses = match1[1]['cost_with_losses'] + new_node2[1]['cost_without_losses'] + additional_cost
                            new_node12 = (utiles.heap_items(cost_without_losses, cost_with_losses, match1, new_node2), event)
                            new_node21 = (utiles.heap_items(cost_without_losses, cost_with_losses, match2, new_node1), event)

                            if new_node12[0] == new_node21[0]:
                                Kbest.insert(len(Kbest), new_node12)
                                Kbest.insert(len(Kbest), new_node21)
                            else:
                                heapq.heappush(Kbest, tuple(new_node12))
                                heapq.heappush(Kbest, tuple(new_node21))
                if not tree_operations.is_a_leaf(x):
                    y = x.adjacent_nodes()[0]
                    z = x.adjacent_nodes()[1]
                    subtree[u.label][x.label] += utiles.kmin_list(find_nodes_in_hypergraph(H, u.label, x.label,  -1,  nodes_table),
                                                                  subtree[u.label][y.label], subtree[u.label][z.label], H, nodes_table,'cost_without_losses')

                    list1 = [deepcopy(tocopy) for tocopy in subtreeLoss[u.label][y.label]]
                    list2 = [deepcopy(tocopy) for tocopy in subtreeLoss[u.label][z.label]]
                    for lst in (list1,list2):
                        for node in lst:
                            node[1]['cost_with_losses'] += loss_cost
                    #print('u: %s, x: %s\nsubtree[u.label][y.label]:%s\nsubtree[u.label][z.label]:%s\nc(u,x): %s\nsubtreeLoss[u.label][y.label]: %s\nsubtreeLoss[u.label][z.label]: %s\n' %
                    #      (str(u.label),str(x.label),str(subtree[u.label][y.label]),str(subtree[u.label][z.label]),
                    #       str(find_nodes_in_hypergraph(H,u.label,x.label,-1,nodes_table)),str(subtreeLoss[u.label][y.label]),subtreeLoss[u.label][z.label]))
                    subtreeLoss[u.label][x.label] += utiles.kmin_list(find_nodes_in_hypergraph(H, u.label, x.label,  -1, nodes_table),
                                                                      list1,
                                                                      list2,
                                                                      H, nodes_table,'cost_with_losses')
                else:
                    subtree[u.label][x.label] += utiles.kmin_list(find_nodes_in_hypergraph(H, u.label, x.label, -1, nodes_table), [],  [],  H, nodes_table, 'cost_without_losses')
                    subtreeLoss[u.label][x.label] += utiles.kmin_list(find_nodes_in_hypergraph(H, u.label, x.label, -1, nodes_table),[],[],  H,  nodes_table, 'cost_with_losses')
                    #print('u: %s, x: %s\nsubtree[u.label][y.label]:%s\nsubtree[u.label][z.label]:%s\nc(u,x): %s\nsubtreeLoss[u.label][y.label]: %s\nsubtreeLoss[u.label][z.label]: %s\n' %
                    #      (str(u.label),str(x.label),str(subtree[u.label][y.label]),str(subtree[u.label][z.label]),
                    #       str(find_nodes_in_hypergraph(H,u.label,x.label,-1,nodes_table)),str(subtreeLoss[u.label][y.label]),subtreeLoss[u.label][z.label]))

        for x in S.preorder_node_iter():
            if not tree_operations.is_a_leaf(x):
                y = x.adjacent_nodes()[0]
                z = x.adjacent_nodes()[1]
                incomp[u.label][y.label] += utiles.kmin_positive(incomp[u.label][x.label]+subtree[u.label][z.label],k,H,nodes_table,'cost_without_losses')
                incomp[u.label][z.label] += utiles.kmin_positive(incomp[u.label][x.label]+subtree[u.label][y.label],k,H,nodes_table,'cost_without_losses')
    if save_data:
        print('     Writing nodes...')
        file = open(path+'/saved_data/H_nodes_effi.txt', 'w')
        file.write(str(H.nodes(data=True)))
        file.close()
        print('     Finished writing nodes.\n')

        print('     Writing edges...')
        file = open(path+'/saved_data/H_edges_k='+str(k)+'.txt', 'w')
        file.write(str(H.edges(data=True)))
        file.close()
        print('     Finished writing edges.\n')

        print('     Writing nodes table...')
        file = open(path+'/saved_data/nodes_table_k='+str(k)+'.txt', 'w')
        file.write(str(nodes_table))
        file.close()
        print('     Finished writing nodes table.\n')

        file = open(path+'/saved_data/incomp.txt', 'w')
        file.write(str(incomp))
        file.close()
        print('     Finished writing nodes table.\n')

        file = open(path+'/saved_data/subtree.txt', 'w')
        file.write(str(subtree))
        file.close()
        print('     Finished writing nodes table.\n')

    #print('     No. of nodes: '+str(H_number_of_nodes * k)+'        No. on edges: '+str(len(H.edges())))
    H_root = [nd for nd in list(H.node(data=True)) if
              nd[1]['s'] == G.seed_node.label and nd[1]['t'] == S.seed_node.label]
    print(track_a_solution(H_root, H, S, G, nx.DiGraph(), 0,-1)[0].nodes(data=True))
    print('Finished building hypergraph.\n')
    return H,H_number_of_nodes, nodes_table


def find_nodes_in_hypergraph(H, s, t, lp, nodes_table):
    #print('nodes table = %s, \n s = %s, t = %s, lp = %s \nH_nodes = %s' % (str(nodes_table), str(s), str(t),str(lp), str(H.node(data=True))))
    if (lp > -1):
        if nodes_table[t][s] >= 0 and len(H.nodes(data=True)[nodes_table[t][s]]['l']) > lp:
            ans = H.nodes(data=True)[nodes_table[t][s]]['l'][lp]
            index = nodes_table[t][s]
            return([(index,ans)])
        else:
            return []
    else:
        if nodes_table[t][s] >= 0 and len(H.nodes(data=True)[nodes_table[t][s]]['l']) > lp:
            ans = H.nodes(data=True)[nodes_table[t][s]]
            index = nodes_table[t][s]
            return([(index,ans)])
        else:
            return []


def SpeciationEvent_effi (u, x, S_cost,subtreeLoss,k,H,nodes_table):
    heap = []
    if tree_operations.has_right_child(u) and tree_operations.has_left_child(u):
        w = u.adjacent_nodes()[0]
        v = u.adjacent_nodes()[1]
        if not tree_operations.is_a_leaf(x):
            if (tree_operations.has_left_child(x)) and (tree_operations.has_right_child(x)):
                y = x.adjacent_nodes()[0]
                z = x.adjacent_nodes()[1]
                list_v_to_left = utiles.kmin_positive(subtreeLoss[v.label][y.label],k,H,nodes_table,'cost_with_losses')
                list_w_to_left = utiles.kmin_positive(subtreeLoss[w.label][y.label],k,H,nodes_table,'cost_with_losses')
                list_v_to_right = utiles.kmin_positive(subtreeLoss[v.label][z.label],k,H,nodes_table,'cost_with_losses')
                list_w_to_right = utiles.kmin_positive(subtreeLoss[w.label][z.label],k,H,nodes_table,'cost_with_losses')

                for hyper_node1 in list_v_to_right:
                    for hyper_node2 in list_w_to_left:
                        cost_without_losses = hyper_node1[1]['cost_without_losses']+hyper_node2[1]['cost_without_losses'] + S_cost
                        cost_with_losses = hyper_node1[1]['cost_with_losses']+hyper_node2[1]['cost_with_losses'] + S_cost
                        heapq.heappush(heap,utiles.heap_items(cost_without_losses, cost_with_losses, hyper_node1 ,hyper_node2))

                for hyper_node1 in list_v_to_left:
                    for hyper_node2 in list_w_to_right:
                        cost_without_losses = hyper_node1[1]['cost_without_losses'] + hyper_node2[1][
                            'cost_without_losses'] + S_cost
                        cost_with_losses = hyper_node1[1]['cost_with_losses'] + hyper_node2[1][
                            'cost_with_losses'] + S_cost
                        heapq.heappush(heap,utiles.heap_items(cost_without_losses, cost_with_losses, hyper_node1, hyper_node2))

    return heap


def DuplicationEvent_effi (H, u, x, D_cost,loss_cost, nodes_table,subtreeLoss,k):
    heap = []
    if not tree_operations.is_a_leaf(x):
        list_v_to_right = []
        list_w_to_right = []
        list_v_to_left = []
        list_w_to_left = []
        if tree_operations.has_right_child(x) and tree_operations.has_right_child(u):
            y = x.adjacent_nodes()[0]
            v = u.adjacent_nodes()[0]
            list_v_to_left = utiles.kmin_positive(subtreeLoss[v.label][y.label],k,H,nodes_table,'cost_with_losses')
        if tree_operations.has_right_child(x) and tree_operations.has_left_child(u):
            y = x.adjacent_nodes()[0]
            w = u.adjacent_nodes()[1]
            list_w_to_left = utiles.kmin_positive(subtreeLoss[w.label][y.label],k,H,nodes_table,'cost_with_losses')
        if tree_operations.has_left_child(x) and tree_operations.has_right_child(u):
            z = x.adjacent_nodes()[1]
            v = u.adjacent_nodes()[0]
            list_v_to_right = utiles.kmin_positive(subtreeLoss[v.label][z.label],k,H,nodes_table,'cost_with_losses')
        if tree_operations.has_left_child(x) and tree_operations.has_left_child(u):
            z = x.adjacent_nodes()[1]
            w = u.adjacent_nodes()[1]
            list_w_to_right = utiles.kmin_positive(subtreeLoss[w.label][z.label],k,H,nodes_table,'cost_with_losses')

        for hyper_node1 in list_v_to_right:
            for hyper_node2 in list_w_to_right:
                cost_without_losses = hyper_node1[1]['cost_without_losses'] + hyper_node2[1]['cost_without_losses'] + D_cost
                cost_with_losses = hyper_node1[1]['cost_with_losses'] + hyper_node2[1]['cost_with_losses'] + D_cost + 2 * loss_cost
                heapq.heappush(heap, utiles.heap_items(cost_without_losses, cost_with_losses, hyper_node1,
                                                hyper_node2))
            hyper_node3 = find_nodes_in_hypergraph(H ,w.label, x.label, 0, nodes_table)
            if hyper_node3:
                hyper_node3 = hyper_node3[0]
                cost_without_losses = hyper_node1[1]['cost_without_losses'] + hyper_node3[1]['cost_without_losses']+ D_cost
                cost_with_losses = hyper_node1[1]['cost_with_losses'] + hyper_node3[1]['cost_with_losses'] + D_cost + loss_cost
                heapq.heappush(heap,
                           utiles.heap_items(cost_without_losses, cost_with_losses, hyper_node1,
                                             hyper_node3))
        for hyper_node1 in list_v_to_left:
            for hyper_node2 in list_w_to_left:
                cost_without_losses = hyper_node1[1]['cost_without_losses'] + hyper_node2[1]['cost_without_losses'] + D_cost
                cost_with_losses = hyper_node1[1]['cost_with_losses'] + hyper_node2[1]['cost_with_losses'] + D_cost + 2 * loss_cost
                heapq.heappush(heap, utiles.heap_items(cost_without_losses, cost_with_losses, hyper_node1,
                                                hyper_node2))
            hyper_node3 = find_nodes_in_hypergraph(H,w.label,x.label,0,nodes_table)

            if hyper_node3:
                hyper_node3 = hyper_node3[0]
                cost_without_losses = hyper_node1[1]['cost_without_losses'] + hyper_node3[1]['cost_without_losses'] + D_cost
                cost_with_losses = hyper_node1[1]['cost_with_losses'] + hyper_node3[1]['cost_with_losses'] + D_cost + loss_cost
                heapq.heappush(heap,
                           utiles.heap_items(cost_without_losses, cost_with_losses, hyper_node1,
                                             hyper_node3))

        for hyper_node1 in list_w_to_right:
            hyper_node2 = find_nodes_in_hypergraph(H, v.label, x.label, 0,nodes_table)
            if hyper_node2:
                hyper_node2 = hyper_node2[0]
                cost_without_losses = hyper_node1[1]['cost_without_losses'] + hyper_node2[1]['cost_without_losses'] + D_cost
                cost_with_losses = hyper_node1[1]['cost_with_losses'] + hyper_node2[1]['cost_with_losses'] + D_cost + loss_cost
                heapq.heappush(heap,
                           utiles.heap_items(cost_without_losses, cost_with_losses,
                                             hyper_node1,
                                             hyper_node2))
        for hyper_node1 in list_w_to_left:
            hyper_node2 = find_nodes_in_hypergraph(H, v.label, x.label, 0,nodes_table)
            if hyper_node2 != []:
                hyper_node2 = hyper_node2[0]
                cost_without_losses = hyper_node1[1]['cost_without_losses'] + hyper_node2[1]['cost_without_losses'] + D_cost
                cost_with_losses = hyper_node1[1]['cost_with_losses'] + hyper_node2[1]['cost_with_losses'] + D_cost + loss_cost
                heapq.heappush(heap,
                           utiles.heap_items(cost_without_losses, cost_with_losses,
                                             hyper_node1,
                                             hyper_node2))
        for hyper_node1 in list_v_to_left:
            for hyper_node2 in list_w_to_right:
                cost_without_losses = hyper_node1[1]['cost_without_losses'] + hyper_node2[1]['cost_without_losses'] + D_cost
                cost_with_losses = hyper_node1[1]['cost_with_losses'] + hyper_node2[1]['cost_with_losses'] + D_cost + 2 * loss_cost
                heapq.heappush(heap, utiles.heap_items(cost_without_losses, cost_with_losses, hyper_node1,
                                                hyper_node2))
        for hyper_node1 in list_w_to_left:
            for hyper_node2 in list_v_to_right:
                cost_without_losses = hyper_node1[1]['cost_without_losses'] + hyper_node2[1]['cost_without_losses'] + D_cost
                cost_with_losses = hyper_node1[1]['cost_with_losses'] + hyper_node2[1]['cost_with_losses'] + D_cost + 2 * loss_cost
                heapq.heappush(heap, utiles.heap_items(cost_without_losses, cost_with_losses, hyper_node1,
                                                hyper_node2))
        hyper_node1 = find_nodes_in_hypergraph(H, w.label, x.label, 0,nodes_table)
        hyper_node2 = find_nodes_in_hypergraph(H, v.label, x.label, 0,nodes_table)
        if hyper_node2 and hyper_node1:
            hyper_node1 = hyper_node1[0]
            hyper_node2 = hyper_node2[0]
            cost_without_losses = hyper_node1[1]['cost_without_losses'] + hyper_node2[1]['cost_without_losses'] + D_cost
            cost_with_losses = hyper_node1[1]['cost_with_losses'] + hyper_node2[1]['cost_with_losses'] + D_cost
            heapq.heappush(heap,
                       utiles.heap_items(cost_without_losses, cost_with_losses,
                                         hyper_node1,
                                         hyper_node2))
    return heap


def HTEvent_effi (u, x, HT_cost,subtreeLoss,incomp,k,H,nodes_table):
    heap = []
    v = u.adjacent_nodes()[0]                                                                                       #left child
    w = u.adjacent_nodes()[1]                                                                                       #right child

    list_v_to_subtree = utiles.kmin_positive(subtreeLoss[v.label][x.label],k,H,nodes_table,'cost_with_losses')
    list_w_to_subtree = utiles.kmin_positive(subtreeLoss[w.label][x.label],k,H,nodes_table,'cost_with_losses')

    list_v_horizontally = utiles.kmin_positive(incomp[v.label][x.label],k,H,nodes_table,'cost_without_losses')
    list_w_horizontally = utiles.kmin_positive(incomp[w.label][x.label],k,H,nodes_table,'cost_without_losses')
    for hyper_node1 in list_v_to_subtree:
        for hyper_node2 in list_w_horizontally:
            cost_without_losses = hyper_node1[1]['cost_without_losses']+hyper_node2[1]['cost_without_losses'] + HT_cost
            cost_with_losses = hyper_node1[1]['cost_with_losses']+hyper_node2[1]['cost_with_losses'] + HT_cost
            heapq.heappush(heap,utiles.heap_items(cost_without_losses, cost_with_losses , hyper_node1, hyper_node2))
    for hyper_node1 in list_w_to_subtree:
        for hyper_node2 in list_v_horizontally:
            cost_without_losses = hyper_node1[1]['cost_without_losses']+hyper_node2[1]['cost_without_losses'] + HT_cost
            cost_with_losses = hyper_node1[1]['cost_with_losses']+hyper_node2[1]['cost_with_losses'] + HT_cost
            heapq.heappush(heap, utiles.heap_items(cost_without_losses, cost_with_losses, hyper_node1, hyper_node2))
    return heap


def calculate_color_diffrence(H, G, S, colors, k, nodes_table, test, real):
    alpha0_map = {}
    alpha1_map = {}

    S_alpha0 = color_hypergraph(H, G, S, colors, k, 0, nodes_table, test, real)[1]
    for x in S_alpha0.postorder_node_iter():
        alpha0_map.update({x.label : [x.E_red, x.E_black]})

    S_alpha1 = color_hypergraph(H, G, S, colors, k, 1, nodes_table, test, real)[1]
    for u in S_alpha1.postorder_node_iter():
        alpha1_map.update({u.label : [u.E_red, u.E_black]})

    for x in S_alpha1.postorder_node_iter():
        temp = [x.E_red, x.E_black]
        #if alpha0_map[x.label] != temp:
            #print ('        %s coloring has changed from %s to %s ([red, black])' % (x.label,str(alpha0_map[x.label]),str(temp)))
            #if (temp[0] > temp[1] and alpha0_map[x.label][0] < alpha0_map[x.label][1]) or (temp[0] < temp[1] and alpha0_map[x.label][0] > alpha0_map[x.label][1]):
                #print('             and it changed the relations between the colors.')
            #print('\n')


def color_hypergraph(H, S, colors, alpha, S_colors):
    for nd in (list(nx.topological_sort(H))):
        incoming_edges = H.in_edges([nd], data=True)
        nd = H.nodes(data=True)[nd]
        if len(incoming_edges) == 0:        #is a leaf
            if colors[nd['t']] == 'red':
                S_colors.update({nd['t'] : [1,0]})
            else:
                S_colors.update({nd['t'] : [0,1]})
        else:
            for i in range(0,len(nd['l'])):
                incoming_edges_v = [e for e in incoming_edges if e[2]['target'] == i]
                prob = nd['l'][i]['probability']
                reds = 0
                blacks = 0
                if prob > 0:
                    for e in incoming_edges_v:
                        source = H.nodes(data=True)[e[0]]['l'][e[2]['source']]['t']
                        if nd['l'][i]['event'] == 'HT' and source != nd['t'] and tree_operations.is_not_ancestor(S.find_node(lambda n: (n.label == source)),S.find_node(lambda n: (n.label == nd['t']))):
                            reds = S_colors[source][0] * alpha * prob
                            blacks = S_colors[source][1] * alpha * prob
                        elif nd['l'][i]['event'] == 'HT':
                            None
                        else:
                            reds = S_colors[source][0] * prob
                            blacks = S_colors[source][1] * prob
                        if nd['t'] in S_colors:
                            temp = S_colors[nd['t']]
                        else:
                            temp = [0,0]
                        S_colors.update({nd['t']: [reds + temp[0], blacks + temp[1]]})
    return H, S, S_colors

def remove_prob_zero(H, deleted_nodes):
    #print('Removing prob. 0...')
    for nd in list(H.node(data=True)):
        if nd[1]['l'][0]['probability'] == 0:
            H.remove_node(nd[0])
            deleted_nodes.append(nd[0])
    #print("\t\tnode deleted: "+str(deleted_nodes))
    #print('Finished removing prob. 0.\n')
    return H

def assign_probabilities(S, G, H, gamma):
    #print('Assigning probs to hypergraph....')
    flag = False
    to_try = G.seed_node.label
    to_try_sp = S.seed_node.label
    max_prob = 0
    for nd in list(reversed(list(nx.topological_sort(H)))):
        outgoing_edges = H.out_edges([nd],data=True)
        incoming_edges = H.in_edges([nd],data=True)
        nd = H.nodes(data=True)[nd]
        for i in range(0,len(nd['l'])):
            outgoing_edges_v = [e for e in outgoing_edges if e[2]['source'] == i]
            incoming_edges_v = [e for e in incoming_edges if e[2]['target'] == i]
            if (nd['s'] != to_try or nd['t'] != to_try_sp):
                if outgoing_edges_v == []:
                    s = 0
                else:
                    s = 0
                    for cost in outgoing_edges_v:
                        cost = cost[2]['probability']
                        s = s+cost
                new_prob = s
                nd['l'][i].update({'probability':new_prob})
                for e in incoming_edges_v:
                    e[2]['probability'] = new_prob
            else:
                flag = True
                nd['l'] = assign_weights_to_list(nd['l'], gamma)
                nd['l'][i].update({'probability': nd['l'][i]['weight']})
                for e in incoming_edges_v:
                    e[2]['probability'] = nd['l'][i]['weight']
    if flag == False:
        print('     ** No reconciliation of '+to_try+' exists. **')
        return None, 0
    for nd in H.nodes(data = True):
        for g in range(0,len(nd[1]['l'])):
            if nd[1]['l'][g]['probability'] > max_prob:
                max_prob = nd[1]['l'][g]['probability']
    #print('Finished assigning probs to hypergraph.\n')
    return H, max_prob

def assign_weights_to_list(l, gamma):
    x_min = l[0]['cost_with_losses']
    x_max = l[len(l)-1]['cost_with_losses']
    index = len(l)
    while (x_max == math.inf):
        index = index -1
        x_max = l[index]['cost_with_losses']
    if (index < len(l)):
        index = index+1
    if (x_max>x_min):
        z = 0
        for j in range(0, index):
            s = l[j]['cost_with_losses']
            z = z+math.pow(math.e,(gamma*(x_min-s)/(x_max-x_min)))
        for j in range(0,index):
            l[j].update({'weight':(1/z)*math.pow(math.e,(gamma*(x_min-l[j]['cost_with_losses'])/(x_max-x_min)))})
        for j in range(index,len(l)):
            l[j].update({'weight': 0})
    else:
        for j in range(0, index):
            l[j].update({'weight': 1/index})
        for j in range(index,len(l)):
            l[j].update({'weight': 0})
    return l

def track_a_solution(root, H, S, G, solution, list_place, cost):
    new_nodes_table = inits.init_nodes_table(S, G, {})
    root_numbers_in_H = [root[0][0]]
    root_numbers_in_solution = [0]

    solution_number_of_nodes = 0
    if cost >= 0:
        roots = list(filter(lambda x: x['cost_with_losses'] == cost, root[0][1]['l']))
        if not roots:
            print([x['cost_with_losses'] for x in root[0][1]['l']])
        else:
            roots = [roots[0]]
    else:
        roots = [root[0][1]['l'][list_place]]

    while not roots == []:
        curr_root = roots.pop()
        curr_index_in_H = root_numbers_in_H.pop()
        index_to_edge = root_numbers_in_solution.pop()
        solution.add_node(solution_number_of_nodes, s=curr_root['s'], t=curr_root['t'],
                          l=[curr_root])
        new_nodes_table[curr_root['t']][curr_root['s']] = solution_number_of_nodes
        if not (solution_number_of_nodes == 0 and index_to_edge == 0):
            solution.add_edge(solution_number_of_nodes, index_to_edge, source=0, target=0)
        solution_number_of_nodes += 1

        incoming_edges = H.in_edges(curr_index_in_H, data=True)
        curr_incoming_edges = [e for e in incoming_edges if e[2]['target'] == curr_root['list_place']]
        if not curr_incoming_edges == []:

            right_child = H.nodes(data=True)[curr_incoming_edges[0][0]]['l'][curr_incoming_edges[0][2]['source']]

            roots.append([nd for nd in list(H.node(data=True)) if nd[1]['s'] == right_child['s'] and nd[1]['t'] == right_child['t']][0][1]['l'][right_child['list_place']])
            root_numbers_in_H.append([nd for nd in list(H.node(data=True)) if nd[1]['s'] == right_child['s'] and nd[1]['t'] == right_child['t']][0][0])

            left_child = H.nodes(data=True)[curr_incoming_edges[1][0]]['l'][curr_incoming_edges[1][2]['source']]

            roots.append([nd for nd in list(H.node(data=True)) if nd[1]['s'] == left_child['s'] and nd[1]['t'] == left_child['t']][0][1]['l'][left_child['list_place']])
            root_numbers_in_H.append([nd for nd in list(H.node(data=True)) if nd[1]['s'] == left_child['s'] and nd[1]['t'] == left_child['t']][0][0])
            root_numbers_in_solution.append(solution_number_of_nodes-1)
            root_numbers_in_solution.append(solution_number_of_nodes-1)
    #print('Finished tracking ' + str(list_place) + 'th solution...\n')
    return solution, new_nodes_table

def find_number_of_cooptimal(H,G,S):
    counter = 1
    for nd in list(reversed(list(nx.topological_sort(H)))):
        nd = H.nodes(data=True)[nd]
        if (nd['s'] == G.seed_node.label and nd['t'] == S.seed_node.label):
            optimal_cost = nd['l'][0]['cost_with_losses']
            for j in range(1, len(nd['l'])):
                if nd['l'][j]['cost_with_losses'] == optimal_cost:
                    counter += 1
                else:
                    return counter,optimal_cost
            return counter,optimal_cost

def mostly_speciation_event_in_subtree(H, nd, i):
    incoming_edges = H.in_edges([nd], data=True)
    incoming_edges_curr = [e for e in incoming_edges if e[2]['target'] == i]
    if incoming_edges_curr == []:
        return 0,0,0,0
    else:
        i0 = incoming_edges_curr[0][2]['source']
        i1 = incoming_edges_curr[1][2]['source']
        child0 = incoming_edges_curr[0][0]
        child1 = incoming_edges_curr[1][0]
        dou_res_child0,HT_res_child0,spe_res_child0,total_res_child0 = mostly_speciation_event_in_subtree(H, child0, i0)
        dou_res_child1,HT_res_child1,spe_res_child1, total_res_child1 = mostly_speciation_event_in_subtree(H, child1, i1)
        nd = H.nodes(data=True)[nd]
        if nd['l'][i]['event'] == 'S':
            return dou_res_child1 + dou_res_child0,HT_res_child0 + HT_res_child1,1 + spe_res_child0 + spe_res_child1, 1 + total_res_child0 + total_res_child1
        elif nd['l'][i]['event'] == 'HT':
            return dou_res_child1 + dou_res_child0,1 + HT_res_child0 + HT_res_child1,spe_res_child0 + spe_res_child1, 1 + total_res_child0 + total_res_child1
        else:
            return 1+dou_res_child1 + dou_res_child0,HT_res_child0 + HT_res_child1, spe_res_child0 + spe_res_child1, 1 + total_res_child0 + total_res_child1

