import math
import networkx as nx
import utiles
import tree_operations_v1 as tree_operations

start = '*'
end = '*'

def identify_pattern2(G, H, k, G_nodes_to_weight, G_nodes_identified):
    print('Identifing patterns...')
    avg_diffrence = utiles.compute_avg_diff(G_nodes_to_weight)
    for u in G.postorder_node_iter():
        G_nodes_identified.update({u.label:0})

    for nd in (list(nx.topological_sort(H))):
        incoming_edges = H.in_edges([nd], data=True)
        nd = H.nodes(data=True)[nd]
        for i in range(0, k):
            curr = nd['l'][i]
            if curr['probability'] > 0:
                incoming_edges_curr = [e for e in incoming_edges if e[2]['target'] == i]
                if len(incoming_edges_curr) > 0:
                    left_egde = incoming_edges_curr[0]
                    right_egde = incoming_edges_curr[1]

                    child1_in_H = H.nodes(data=True)[left_egde[0]]['l'][left_egde[2]['source']]
                    child2_in_H = H.nodes(data=True)[right_egde[0]]['l'][right_egde[2]['source']]

                    if math.fabs(G_nodes_to_weight[(curr['s'],child1_in_H['s'])][0] - G_nodes_to_weight[(curr['s'],child2_in_H['s'])][1]) > avg_diffrence:
                        new_prob = G_nodes_identified[curr['s']] + curr['probability']
                        G_nodes_identified.update({curr['s']:new_prob})
                    if math.fabs(G_nodes_to_weight[(curr['s'],child1_in_H['s'])][1] - G_nodes_to_weight[(curr['s'],child2_in_H['s'])][0]) > avg_diffrence:
                        new_prob = G_nodes_identified[curr['s']] + curr['probability']
                        G_nodes_identified.update({curr['s']:new_prob})
    print('Finished dentifing patterns...\n')
    return G_nodes_identified

def find_signi_distance(new_G, all_vertices, TH_compare_subtrees, k, doup,compare_subtrees, TH_edges_in_subtree,max_score_TH,max_score_doup):
    marked_nodes = {}
    for u in (list(nx.topological_sort(new_G))):
        outgoing_edges = new_G.out_edges([u], data=True)
        outgoing_edges = [e for e in outgoing_edges]
        u = new_G.nodes(data=True)[u]
        if len(outgoing_edges) == 2:
            v = new_G.nodes(data=True)[outgoing_edges[0][1]]
            w = new_G.nodes(data=True)[outgoing_edges[1][1]]

            u_red_HT = 0
            u_black_HT = 0
            u_red_doup = 0
            u_black_doup = 0

            if not u['edges_in_subtree'] == 0:
                u_red_HT = u['same_HT_score'][0]
                u_black_HT = u['same_HT_score'][1]
                u_red_doup = u['same_doup_score'][0]
                u_black_doup = u['same_doup_score'][1]


            if not compare_subtrees:
                if not doup:
                    print(
                        '     %s (v = %s, w = %s) :\n        [red HT: %s ,black HT: %s] (unnormlized: [red HT: %s ,black HT: %s]), edges in subtree u: %s\n         TH_edges: %s, max_scours_HT: %s,TH_compare_subtrees: %s\n' %
                        (u['label'], str(v['label']), str(w['label']), str(u_red_HT), str(u_black_HT),
                         str(u['same_HT_score'][0]), str(u['same_HT_score'][1]),
                         str(u['edges_in_subtree']), str(TH_edges_in_subtree), str(max_score_TH),
                         str(TH_compare_subtrees)))
                    if u['edges_in_subtree'] >= TH_edges_in_subtree:
                        if u_black_HT >= TH_compare_subtrees * u_red_HT:
                            all_vertices.update({u['label']: (u_red_HT, u_black_HT)})
                            if u_black_HT in max_score_TH:
                                marked_nodes.update({u['label']: [(u_red_HT, u_black_HT), (0, 0),
                                                              'black HT']})
                        if u_red_HT >= TH_compare_subtrees * u_black_HT:
                            all_vertices.update({u['label']: (u_red_HT, u_black_HT)})
                            if u_red_HT in max_score_TH:
                                marked_nodes.update({u['label']: [(u_red_HT, u_black_HT), (0, 0),
                                                                  'red HT']})
                else:
                    #print(
                    #    '     %s (v = %s, w = %s) :\n        [red doup: %s ,black doup: %s] (unnormlized: [red doup: %s ,black doup: %s]), edges in subtree u: %s\n         TH_edges: %s, TH_pattern_in_subtree: %s,TH_compare_subtrees: %s\n' %
                    #    (u['label'], str(v['label']), str(w['label']), str(u_red_doup), str(u_black_doup),
                    #     str(u['same_doup_score'][0]), str(u['same_doup_score'][1]),
                    #     str(u['edges_in_subtree']), str(TH_edges_in_subtree), str(TH_pattern_in_subtree),
                    #     str(TH_compare_subtrees)))
                    if u['edges_in_subtree'] >= TH_edges_in_subtree:
                        if u_black_doup > TH_compare_subtrees * u_red_doup and u_black_doup in max_score_doup:
                            marked_nodes.update({u['label']: [(u_red_doup, u_black_doup), (0, 0),
                                                                'red doup']})
                        if u_red_doup > TH_compare_subtrees * u_black_doup and u_red_doup in max_score_doup:
                            marked_nodes.update({u['label']: [(u_red_doup, u_black_doup), (0, 0),
                                                                'black doup']})
        else:
            if not compare_subtrees:
                if (not doup) and (TH_edges_in_subtree == 0) and (TH_compare_subtrees == 0):
                    all_vertices.update({u['label']: (0, 0)})
                    marked_nodes.update(
                        {u['label']: [(0, 0), (0, 0), 'leaf node']})

    #print ('        marked nodes: %s' % str(marked_nodes))
    #print('     Writing marked nodes...')
    #if check_diff_sol:
    #    file = open(path+'/saved_data/marked_nodes_for_TH_'+str(index)+'.txt', 'w')
    #else:
    #    file = open(path + '/saved_data/marked.txt', 'w')
    #file.write(str(marked_nodes))
    #file.close()
    #print('     Finished writing marked nodes.\n')

    return marked_nodes,all_vertices

def find_avg_diff(G):
    sum = 0
    for u in G.nodes(data = True):
        diff = math.fabs(u[1]['same_HT_score'][0] - u[1]['same_HT_score'][1])
        sum += diff
    return sum/len(G.edges())

def compare_marked_nodes (alpha1, alpha2, path, k1, k2, both1, boTH_both, index_start, index_end, times,spe):
    if index_start > -1:        ##using only first values for all indices
        alpha_maps = {}
        for index in range(index_start,index_end):
            #print("     Reading file 'saved_data/marked_nodes_k=" + str(k1) + "_alpha=" + str(
            #    alpha1) + ".txt'...")
            input = open(path+'/saved_data/marked_nodes_k=' + str(k1) + '_alpha=' + str(alpha1) + '_both=' + str(both1) +"_(index=" + str(index)+ ')_species='+spe+'.txt', 'r')
            alpha_maps[index] = []
            for line in input:
                alpha_maps[index].append(eval(line))
            if len(alpha_maps[index]) == 0:
                alpha_maps[index] = {}
            else:
                alpha_maps[index] = alpha_maps[index][0]
            #print("     Finished reading file 'saved_data/marked_nodes_k=" + str(k1) + "_alpha=" + str(
            #    alpha1) + ".txt'")

        for index1 in range(index_start, index_end):
            for index2 in range(index_start, index_end):
                for key in alpha_maps[index1]:
                    if not key in alpha_maps[index2]:
                        print('        %s is marked in the %s solution and not in the %s solution\n' % (
                        str(key), str(index1), str(index2)))
    else:
        #print("     Reading file 'saved_data/marked_nodes_k=" + str(k1) + "_alpha=" +str(alpha1) + ".txt'...")
        input = open(path + '/saved_data/marked_nodes_k=' + str(k1) + '_alpha=' +str(alpha1) + '_both=' + str(both1) + '.txt', 'r')
        alpha1_map = []
        for line in input:
            alpha1_map.append(eval(line))
        if len(alpha1_map) == 0:
            alpha1_map = {}
        else:
            alpha1_map = alpha1_map[0]
        #print("     Finished reading file 'saved_data/marked_nodes_k=" + str(k1)  + "_alpha=" +str(alpha1) + ".txt'")

        #print("     Reading file 'saved_data/marked_nodes_k=" + str(k2) + "_alpha=" + str(alpha2) + '_both=' + str(boTH_both) + ".txt'...")
        input = open(path+ '/saved_data/marked_nodes_k=' + str(k2)  + '_alpha=' + str(alpha2) +  '_both=' + str(boTH_both) +'.txt', 'r')
        alpha2_map = []
        for line in input:
            alpha2_map.append(eval(line))
        if len(alpha2_map) == 0:
            alpha2_map = {}
        else:
            alpha2_map = alpha2_map[0]
        #print("     Finished reading file 'saved_data/marked_nodes_k=" + str(k2) + "_alpha=" + str(
        #    alpha2) + ".txt'")

        for key in alpha1_map:
            if not key in alpha2_map:
                print ('        %s is marked for k=%s, alpha=%s, both=%s and not for k=%s, alpha=%s, both=%s\n' % (str(key), str(k1),str(alpha1),str(both1),str(k2),str(alpha2),str(boTH_both)))
        for key1 in alpha2_map:
            if not key1 in alpha1_map:
                print ('        %s is marked for k=%s, alpha=%s, both=%s and not for k=%s, alpha=%s, both=%s\n' % (str(key1), str(k2),str(alpha2),str(boTH_both),str(k1),str(alpha1),str(both1)))
    quit()


