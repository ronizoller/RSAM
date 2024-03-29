import sys
sys.path.append('/anaconda3/lib/python3.6/site-packages')

import tree_operations_v1 as tree_operations
import heapq
import utiles
import math
import networkx as nx
import inits_v1 as inits
import random
from numpy import inf


def build_hyper_garph(S, G, test, k,nodes_table, D_cost, S_cost, HT_cost, path, alpha, sigma,save_data):
    #print('Building hypergraph...')
    H = nx.MultiDiGraph()
    H.clear()

    if test:
        print("     Reading file 'H_edges.txt'...")
        input =  open(path+'/saved_data/H_edges_k='+str(k)+'.txt','r')
        new_edges = []
        for line in input:
            new_edges.append(eval(line))
        new_edges = new_edges[0]
        print("     Finished reading file 'H_edges.txt'")

        print("     Reading file 'H_nodes.txt'...")
        input = open(path+'/saved_data/H_nodes_k='+str(k)+'_alpha='+str(alpha)+'.txt', 'r')
        new_nodes = []
        for line in input:
            new_nodes.append(eval(line))
        new_nodes = new_nodes[0]
        print("     Finished reading file 'H_nodes.txt'.")

        print("     Reading file 'nodes_table.txt'...")
        input = open(path+'/saved_data/nodes_table_k='+str(k)+'.txt', 'r')
        nodes = {}
        for line in input:
            nodes.update(eval(line))
        nodes_table = nodes
        print("     Finished reading file 'nodes_table_k="+str(k)+".txt'")


        H.add_nodes_from(new_nodes)
        H.add_edges_from(new_edges)
        H_number_of_nodes = (len(H.nodes()))
    else:
        H, H_number_of_nodes, nodes_table = inits.init_leafs(G, H, k, 0, sigma, nodes_table)
        for u in G.postorder_node_iter():
            if(not tree_operations.is_a_leaf(u)):
                for x in S.postorder_node_iter():
                    key_counter = 0
                    SE = list(map(lambda nd:(nd,'S'), SpeciationEvent(H, u, x, S_cost, nodes_table)))
                    DE = list(map(lambda nd:( nd,'D'),DuplicationEvent(H, u, x, D_cost, nodes_table)))
                    HTE = list(map(lambda nd:(nd,'HT'), HTEvent(S, H, u, x, HT_cost, nodes_table)))

                    Kbest = SE+DE+HTE
                    random.shuffle(Kbest,random.random)
                    heapq.heapify(Kbest)

                    if (len(Kbest) > 0):
                        H.add_node(H_number_of_nodes, s=u.label, t=x.label,l=[])
                        nodes_table[x.label][u.label] = H_number_of_nodes
                        big_node = H_number_of_nodes
                        H_number_of_nodes += 1


                    for i in range(0,k):
                        if (len(Kbest)>0):
                            match = heapq.heappop(Kbest)
                            new_cost = match[0].val[0]
                            event = match[1]
                            match1=match[0].val[1]
                            match2 = match[0].val[2]

                            if (event=='S'):
                                H.nodes[big_node]['l'] = H.nodes[big_node]['l'] + [{'s': u.label, 't': x.label, 'cost': new_cost, 'event': "S",'list_place':len(H.nodes[big_node]['l'])}]
                            elif(event=='D'):
                                H.nodes[big_node]['l'] = H.nodes[big_node]['l'] + [{'s': u.label, 't': x.label, 'cost': new_cost, 'event': "D",'list_place':len(H.nodes[big_node]['l'])}]
                            else:
                                H.nodes[big_node]['l'] = H.nodes[big_node]['l'] + [{'s': u.label, 't': x.label, 'cost': new_cost, 'event': "HT",'list_place':len(H.nodes[big_node]['l'])}]

                            H.add_edge(match1[0], big_node, key=key_counter, source=match1[1]['list_place'],target=len(H.nodes[big_node]['l'])-1,probability=0)
                            key_counter += 1
                            H.add_edge(match2[0], big_node, key=key_counter, source=match2[1]['list_place'], target=len(H.nodes[big_node]['l'])-1,probability=0)
                            key_counter += 1

                            new_node1 = find_nodes_in_hypergraph(H, match1[1]['s'], match1[1]['t'], match1[1]['list_place']+1, nodes_table)
                            new_node2 = find_nodes_in_hypergraph(H, match2[1]['s'], match2[1]['t'],match2[1]['list_place']+1, nodes_table)
                            if (new_node1!=[] and new_node2!=[]):
                                new_node1 = new_node1[0]
                                new_node2 = new_node2[0]
                                if event=='D':
                                    additional_cost = D_cost
                                elif event=='S':
                                    additional_cost = S_cost
                                else:
                                    additional_cost = HT_cost
                                new_node12 = (utiles.heap_items(match1[1]['cost']+new_node2[1]['cost']+additional_cost,match1, new_node2),event)
                                new_node21 = (utiles.heap_items(match2[1]['cost']+new_node1[1]['cost']+additional_cost,match2, new_node1),event)

                                if (new_node12[0]==new_node21[0]):
                                    Kbest.insert(len(Kbest), new_node12)
                                    Kbest.insert(len(Kbest), new_node21)

                                else:
                                    heapq.heappush(Kbest, tuple(new_node12))
                                    heapq.heappush(Kbest, tuple(new_node21))
        if save_data:
            print('     Writing nodes...')
            file = open(path+'/saved_data/H_nodes_naive.txt', 'w')
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

    #print('     No. of nodes: '+str(H_number_of_nodes * k)+'        No. on edges: '+str(len(H.edges())))
    #print('Finished building hypergraph.\n')
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

def SpeciationEvent (H, u, x, S_cost, nodes_table):
    heap = []
    if (not tree_operations.is_a_leaf(x)):
        v = []                                                                                     #left child
        w = []                                                                                       #right child
        list_v_to_right = []
        list_w_to_right = []
        list_v_to_left = []
        list_w_to_left = []

        if (tree_operations.has_right_child(u)):
            w = u.adjacent_nodes()[0]

        if (tree_operations.has_left_child(u)):
            v = u.adjacent_nodes()[1]

        if (tree_operations.has_left_child(x) and tree_operations.has_right_child(x)):
            for y in x.adjacent_nodes()[0].postorder_iter():                                                                #left subtree
                list_v_to_left =  list_v_to_left+find_nodes_in_hypergraph(H, v.label, y.label, 0, nodes_table)
                list_w_to_left = list_w_to_left+find_nodes_in_hypergraph(H, w.label,y.label,0, nodes_table)


            for z in x.adjacent_nodes()[1].postorder_iter():                                                                # right subtree
                list_v_to_right = list_v_to_right+find_nodes_in_hypergraph(H, v.label, z.label, 0, nodes_table)
                list_w_to_right = list_w_to_right+find_nodes_in_hypergraph(H, w.label, z.label, 0, nodes_table)

            for hyper_node1 in list_v_to_right:
                for hyper_node2 in list_w_to_left:
                    heapq.heappush(heap,utiles.heap_items(hyper_node1[1]['cost']+hyper_node2[1]['cost']+S_cost,hyper_node1,hyper_node2))

            for hyper_node1 in list_v_to_left:
                for hyper_node2 in list_w_to_right:
                    heapq.heappush(heap,utiles.heap_items(hyper_node1[1]['cost']+hyper_node2[1]['cost']+S_cost,hyper_node1,hyper_node2))
    return heap

def DuplicationEvent (H, u, x, D_cost, nodes_table):
    heap = []
    if (not tree_operations.is_a_leaf(x)):
        v = []
        w = []
        list_v_to_right = []
        list_w_to_right = []
        list_v_to_left = []
        list_w_to_left = []

        if (tree_operations.has_right_child(u)):
            v = u.adjacent_nodes()[0]

        if (tree_operations.has_left_child(u)):
            w = u.adjacent_nodes()[1]

        if (tree_operations.has_right_child(x)):
            for y in x.adjacent_nodes()[0].postorder_iter():  # left subtree
                list_v_to_left = list_v_to_left + find_nodes_in_hypergraph(H, v.label, y.label, 0, nodes_table)
                list_w_to_left = list_w_to_left + find_nodes_in_hypergraph(H, w.label, y.label, 0, nodes_table)

        if (tree_operations.has_left_child(x)):
            for z in x.adjacent_nodes()[1].postorder_iter():  # right subtree
                list_v_to_right = list_v_to_right + find_nodes_in_hypergraph(H, v.label, z.label, 0, nodes_table)
                list_w_to_right = list_w_to_right + find_nodes_in_hypergraph(H, w.label, z.label, 0, nodes_table)

        for hyper_node1 in list_v_to_right:
            for hyper_node2 in list_w_to_right:
                heapq.heappush(heap, utiles.heap_items(hyper_node1[1]['cost'] + hyper_node2[1]['cost'] + D_cost, hyper_node1,
                                                hyper_node2))

        for hyper_node1 in list_v_to_left:
            for hyper_node2 in list_w_to_left:
                heapq.heappush(heap, utiles.heap_items(hyper_node1[1]['cost'] + hyper_node2[1]['cost'] + D_cost, hyper_node1,
                                                hyper_node2))

    return heap

def HTEvent (S, H, u, x, HT_cost, nodes_table):
    heap = []
    v = u.adjacent_nodes()[0]                                                                                       #left child
    w = u.adjacent_nodes()[1]                                                                                       #right child

    list_v_to_subtree = []
    list_w_to_subtree = []
    list_v_horizontally = []
    list_w_horizontally = []

    for y in x.postorder_iter():
        list_v_to_subtree =  list_v_to_subtree + find_nodes_in_hypergraph(H, v.label,y.label,0, nodes_table)
        list_w_to_subtree = list_w_to_subtree + find_nodes_in_hypergraph(H, w.label,y.label,0, nodes_table)

    for z in S.postorder_node_iter(lambda nd: (tree_operations.is_not_ancestor(nd,x) and tree_operations.is_not_ancestor(x,nd))):
        list_v_horizontally = list_v_horizontally + find_nodes_in_hypergraph(H, v.label, z.label, 0, nodes_table)
        list_w_horizontally = list_w_horizontally + find_nodes_in_hypergraph(H, w.label, z.label, 0, nodes_table)

    for hyper_node1 in list_v_to_subtree:
        for hyper_node2 in list_w_horizontally:
            heapq.heappush(heap,utiles.heap_items(hyper_node1[1]['cost']+hyper_node2[1]['cost']+HT_cost,hyper_node1,hyper_node2))

    for hyper_node1 in list_w_to_subtree:
        for hyper_node2 in list_v_horizontally:
            heapq.heappush(heap, utiles.heap_items(hyper_node1[1]['cost']+hyper_node2[1]['cost']+HT_cost,hyper_node1, hyper_node2))
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

def color_hypergraph(H, G, S, colors, k, alpha, nodes_table, S_colors):
    #print('Coloring hypergraph...')
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
                        #print('add [%s,%s] to %s (for %s, %s event)' % (str(reds),str(blacks),str(temp),str(nd['t']),str(nd['l'][i]['event'])))
                        S_colors.update({nd['t']: [reds + temp[0], blacks + temp[1]]})
                        #print('S_colors = %s\nnd = %s\nsource nd = %s\n' % (str(S_colors), str(nd), str(source)))
    #print('Finished coloring hypergraph.\n')
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

def assign_probabilities(S, G, H, test, k, gamma, path, alpha):
    #print('Assigning probs to hypergraph....')
    flag = False
    to_try = G.seed_node.label
    to_try_sp = S.seed_node.label
    max_prob = 0
    if True:
        for nd in list(reversed(list(nx.topological_sort(H)))):
            outgoing_edges = H.out_edges([nd],data=True)
            incoming_edges = H.in_edges([nd],data=True)
            nd = H.nodes(data=True)[nd]
            for i in range(0,len(nd['l'])):
                outgoing_edges_v = [e for e in outgoing_edges if e[2]['source'] == i]
                incoming_edges_v = [e for e in incoming_edges if e[2]['target'] == i]
                if (nd['s'] != to_try or nd['t'] != to_try_sp):
                    if outgoing_edges_v==[]:
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
                    #if new_prob > 1:
                        #print ('This is bad: nd = '+str(nd)+' with prob: '+str(new_prob)+'\n')
                else:
                    flag = True
                    nd['l'] = assign_weights_to_list(nd['l'], gamma)
                    nd['l'][i].update({'probability': nd['l'][i]['weight']})
                    for e in incoming_edges_v:
                        e[2]['probability'] = nd['l'][i]['weight']
        if flag == False:
            print('     ** No reconciliation of '+to_try+' exists. **')
            return None, 0
        #print('     Writing nodes...')
        #file = open(path+'/saved_data/H_nodes_naive.txt', 'w')
        #file.write(str(H.nodes(data=True)))
        #file.close()
        #print('     Finished writing nodes.\n')
    #else :
        #print('     Probs were already assinged')
    for nd in H.nodes(data = True):
        for g in range(0,len(nd[1]['l'])):
            if nd[1]['l'][g]['probability'] > max_prob:
                max_prob = nd[1]['l'][g]['probability']
    #print('Finished assigning probs to hypergraph.\n')
    return H, max_prob

def assign_weights_to_list(l, gamma):
    x_min = l[0]['cost']
    x_max = l[len(l)-1]['cost']
    index = len(l)
    while (x_max == math.inf):
        index = index -1
        x_max = l[index]['cost']
    if (index < len(l)):
        index = index+1
    if (x_max>x_min):
        z = 0
        for j in range(0, index):
            s = l[j]['cost']
            z = z+math.pow(math.e,(gamma*(x_min-s)/(x_max-x_min)))
        for j in range(0,index):
            l[j].update({'weight':(1/z)*math.pow(math.e,(gamma*(x_min-l[j]['cost'])/(x_max-x_min)))})
        for j in range(index,len(l)):
            l[j].update({'weight': 0})
    else:
        for j in range(0, index):
            l[j].update({'weight': 1/index})
        for j in range(index,len(l)):
            l[j].update({'weight': 0})
    return l

##### def assign_weights():
    #print('Assigning weights to hypergraph....')
    #global gamma
    #for nd in H.nodes(data=True):
        #assign_weights_to_list(nd[1]['l'])

def track_a_solution(root, H, S, G, solution, list_place):
    #print('Tracking ' + str(list_place) +'th solution...')

    new_nodes_table = inits.init_nodes_table(S, G, {})
    #print('root: '+str(root))
    root_numbers_in_H = [root[0][0]]
    root_numbers_in_solution = [0]

    solution_number_of_nodes = 0
    roots = [root[0][1]['l'][list_place]]

    while not roots == []:
        #print('**\n roots = %s\n    root_numbers_in_H = %s\n    root_numbers_in_solution = %s\n' % (str(roots),str(root_numbers_in_H),str(root_numbers_in_solution)))

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
    #print('     solution number %s is: \n   %s' % (str(list_place),str(solution.nodes(data=True))))
    #print('Finished tracking ' + str(list_place) + 'th solution...\n')
    return solution, new_nodes_table

def find_number_of_cooptimal(H,G,S,k):
    #print('computing number of co-optimal solutions...')
    counter = 1
    for nd in list(reversed(list(nx.topological_sort(H)))):
        nd = H.nodes(data=True)[nd]
        if (nd['s'] == G.seed_node.label and nd['t'] == S.seed_node.label):
            optimal_cost = nd['l'][0]['cost']
            for j in range(1, len(nd['l'])):
                if nd['l'][j]['cost'] == optimal_cost:
                    counter += 1
                else:
                    quit()
            quit()
    #print('Finished computing number of co-optimal solutions.'+str(counter))
