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

def find_signi_distance(new_G, all_vertices, p1,p2,max_score_p1_list,max_score_p1_and_p2_list ,check_diffrence_between_solutions):

    marked_nodes = {}
    for u in (list(nx.topological_sort(new_G))):
        outgoing_edges = new_G.out_edges([u], data=True)
        outgoing_edges = [e for e in outgoing_edges]
        u = new_G.nodes(data=True)[u]
        if len(outgoing_edges) == 2:
            v = new_G.nodes(data=True)[outgoing_edges[0][1]]
            w = new_G.nodes(data=True)[outgoing_edges[1][1]]

            u_p1_score = 0

            v_p1_score = 0
            v_p2_score = 0
            w_p1_score = 0
            w_p2_score = 0

            if not u['edges_in_subtree'] == 0:
                u_p1_score = u['p1']
                if p2[0] is not None:
                    if not v['edges_in_subtree'] == 0:
                        v_p1_score = v['p1']
                        v_p2_score = v['p2']
                    if not w['edges_in_subtree'] == 0:
                        w_p1_score = w['p1']
                        w_p2_score = w['p2']

            if p2[0] == None:
                print(
                    '(Single-mode)     %s (v = %s, w = %s) :\n        p1 score:  %s , edges in subtree u: %s\n     TH_edges: %s, max_scores p1: %s\n' %
                    (u['label'], str(v['label']), str(w['label']), str(u_p1_score),
                     str(u['edges_in_subtree']),str(p1[3]), str(max_score_p1_list)))

                if u['edges_in_subtree'] >= p1[3]:
                        if check_diffrence_between_solutions:
                            all_vertices.update({u['label']: u_p1_score})
                        elif u_p1_score in max_score_p1_list:
                                marked_nodes.update({u['label']: [u_p1_score,'Single-mode']})
            else:
                print(
                    '(Douple-mode)     %s (v = %s, w = %s) :\n        v p1 score:  %s, w p2 score: %s,w p1 score:  %s, v p2 score: %s , edges in subtree u: %s\n     TH_edges p1: %s\n      max_scores p1+p2: %s\n' %
                    (u['label'], str(v['label']), str(w['label']), str(v_p1_score), str(w_p2_score),str(w_p1_score), str(v_p2_score),
                     str(u['edges_in_subtree']),str(p1[3]), str(max_score_p1_and_p2_list)))

                if v['edges_in_subtree'] >= p1[3] and w['edges_in_subtree'] >= p1[3]:
                        if check_diffrence_between_solutions:
                            all_vertices.update({u['label']: u_p1_score})
                        elif v_p1_score + w_p2_score in max_score_p1_and_p2_list and v['edges_in_subtree'] > w['edges_in_subtree']:
                                marked_nodes.update({u['label']: [v_p1_score,w_p2_score, 'Double-mode,'+str(v['label']+'_p1')]})
                        elif w_p1_score + v_p2_score in max_score_p1_and_p2_list and w['edges_in_subtree'] > v['edges_in_subtree']:
                                marked_nodes.update({u['label']: [w_p1_score,v_p2_score, 'Double-mode,'+str(w['label']+'_p1')]})
        else:
            if check_diffrence_between_solutions and p1[3] == 0:
                all_vertices.update({u['label']: (0, 0)})
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


