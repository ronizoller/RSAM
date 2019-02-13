import math
import networkx as nx
import utiles

start = '*'
end = '*'

def identify_pattern2(G, H, k, G_nodes_to_weight, G_nodes_identified):
    print('Identifing patterns...')

    avg_diffrence = utiles.compute_avg_diff(G_nodes_to_weight)
    print('     avg_diff = %s' % avg_diffrence)

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
                        print('     %s has alot of red HT and %s has alot of black HT.' % (child1_in_H['s'],child2_in_H['s']))
                    if math.fabs(G_nodes_to_weight[(curr['s'],child1_in_H['s'])][1] - G_nodes_to_weight[(curr['s'],child2_in_H['s'])][0]) > avg_diffrence:
                        new_prob = G_nodes_identified[curr['s']] + curr['probability']
                        G_nodes_identified.update({curr['s']:new_prob})
                        print('     %s has alot of black HT and %s has alot of red HT.' % (child1_in_H['s'],child2_in_H['s']))
    print ('        G_nodes_identified = %s' % str(G_nodes_identified))
    print('Finished dentifing patterns...\n')
    return G_nodes_identified

def find_signi_distance(new_G, all_vertices, TH_compare_subtrees, TH_both, TH_pattern_in_subtree, path, k, alpha, both, G_internal_colors,index,spec,compare_subtrees, TH_edges_in_subtree,check_diff_sol):        #if both = True the pattern will be two sibs which are riched with same color HT
    marked_nodes = {}
    print('Searching for significent diffrence...')
    #print('     edges TH : %s' % str(TH_edges_in_subtree))
    #print('     max_S_d_of_HT_red = %s, max_S_d_of_HT_black = %s, max_prob = %s, threshold_red = %s, threshold_black = %s, len(red_HT_vertices_in_G) = %s, len(blacks_HT_vertices_in_G) = %s' % (str(max_S_d_of_HT[0]),str(max_S_d_of_HT[1]),str(max_prob), str(TH_compare_subtrees), str(TH_compare_subtrees), str(len(red_HT_vertices_in_G)), str(len(blacks_HT_vertices_in_G))))
    for u in (list(nx.topological_sort(new_G))):
        outgoing_edges = new_G.out_edges([u], data=True)
        outgoing_edges = [e for e in outgoing_edges]
        u = new_G.nodes(data=True)[u]
        if len(outgoing_edges) == 2:
            v = new_G.nodes(data=True)[outgoing_edges[0][1]]
            w = new_G.nodes(data=True)[outgoing_edges[1][1]]

            v_red_HT = 0
            w_red_HT = 0
            u_red_HT = 0

            v_black_HT = 0
            w_black_HT = 0
            u_black_HT = 0

            if not v['edges_in_subtree'] == 0:
                v_red_HT = v['same_HT_score'][0] / (v['edges_in_subtree'] * k)
                v_black_HT = v['same_HT_score'][1] / (v['edges_in_subtree'] * k)
            if not w['edges_in_subtree'] == 0:
                w_red_HT = w['same_HT_score'][0] / (w['edges_in_subtree'] * k)
                w_black_HT = w['same_HT_score'][1] / (w['edges_in_subtree'] * k)
            if not u['edges_in_subtree'] == 0:
                u_red_HT = u['same_HT_score'][0] / (u['edges_in_subtree'] * k)
                u_black_HT = u['same_HT_score'][1] / (u['edges_in_subtree'] * k)


            reds_under_w = G_internal_colors[w['label']][0]
            blacks_under_w = G_internal_colors[w['label']][1]
            reds_under_v = G_internal_colors[v['label']][0]
            blacks_under_v = G_internal_colors[v['label']][1]
            all_leafs_v = reds_under_v + blacks_under_v
            all_leafs_w = reds_under_w + blacks_under_w

            #print('       %s (v = %s, w = %s)\n     all_leafs_v: %s, all_leafs_w: %s' % (u['label'] ,str(v['label']),str(w['label']),str(all_leafs_v), str(all_leafs_w)))
            #print('     %s (v = %s, w = %s) :\n        [red HT v: %s ,black HT v: %s], [red HT w: %s ,black HT w: %s]\n      [red under v: %s ,black under v: %s], [red under w: %s ,black under w: %s]\n       edges in subtree u: %s, edges in subtree v: %s, edges in subtree w: %s, TH_edges: %s, TH_pattern_in_subtree: %s\n' %
            #     (u['label'] ,str(v['label']),str(w['label']),str(v_red_HT),str(v_black_HT),str(w_red_HT),str(w_black_HT),
            #       str(reds_under_v/all_leafs_v),str(blacks_under_v/all_leafs_v),str(reds_under_w/all_leafs_w),str(blacks_under_w/all_leafs_w),
            #       str(u['edges_in_subtree']),str(v['edges_in_subtree']),str(w['edges_in_subtree']),str(TH_edges_in_subtree),str(TH_pattern_in_subtree)))
            if not compare_subtrees:
                if not both:
                    if v['edges_in_subtree'] > TH_edges_in_subtree and w['edges_in_subtree'] > TH_edges_in_subtree:
                        if blacks_under_w / all_leafs_w >= TH_both:
                            if v_red_HT > TH_compare_subtrees * v_black_HT:
                                all_vertices.update({u['label']: (v_red_HT, v_black_HT)})
                                if v_red_HT >= TH_pattern_in_subtree:
                                    marked_nodes.update({u['label']: [(v_red_HT, v_black_HT), (reds_under_w, blacks_under_w),
                                                                'v red HT and blacks under w']})
                        if blacks_under_v / all_leafs_v >= TH_both:
                            if w_red_HT > TH_compare_subtrees * w_black_HT:
                                all_vertices.update({u['label']: (w_red_HT, w_black_HT)})
                                if w_red_HT >= TH_pattern_in_subtree:
                                    marked_nodes.update({u['label']: [(w_red_HT, w_black_HT), (reds_under_v, blacks_under_v),
                                                                  'w red HT and blacks under v']})
                        if reds_under_w / all_leafs_w >= TH_both:
                            if v_black_HT > TH_compare_subtrees * v_red_HT:
                                all_vertices.update({u['label']: (v_red_HT, v_black_HT)})
                                if v_black_HT >= TH_pattern_in_subtree:
                                    marked_nodes.update({u['label']: [(v_red_HT, v_black_HT), (reds_under_w, blacks_under_w),
                                                              'v blacks HT and reds under w']})
                        if reds_under_v / all_leafs_v >= TH_both:
                            if w_black_HT > TH_compare_subtrees * w_red_HT:
                                all_vertices.update({u['label']: (w_red_HT, w_black_HT)})
                                if w_black_HT >= TH_pattern_in_subtree:
                                    marked_nodes.update({u['label']: [(w_red_HT, w_black_HT), (reds_under_w, blacks_under_w),
                                                                  'w black HT and reds under v']})
            elif  compare_subtrees:
                if u['edges_in_subtree'] > TH_edges_in_subtree:
                    #print("u = %s, u_red_HT = %s, u_black_HT = %s, TH = %s" % (str(u),str(u_red_HT),str(u_black_HT),str(TH_pattern_in_subtree)))
                    if u_red_HT > TH_pattern_in_subtree:
                        marked_nodes.update({u['label']: [(u_red_HT, u_black_HT),(0,0),'u_red']})
                    if u_black_HT > TH_pattern_in_subtree:
                        marked_nodes.update({u['label']: [(u_red_HT, u_black_HT),(0,0),'u_black']})

    print ('        marked nodes: %s' % str(marked_nodes))
    print('     Writing marked nodes...')
    if check_diff_sol:
        file = open(path+'/saved_data/marked_nodes/'+str(index)+'.txt', 'w')
    else:
        file = open(path + '/saved_data/marked.txt', 'w')
    file.write(str(marked_nodes))
    file.close()
    print('     Finished writing marked nodes.\n')

    print('Finished searching for significent diffrence.\n')
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
            print("     Reading file 'saved_data/marked_nodes_k=" + str(k1) + "_alpha=" + str(
                alpha1) + ".txt'...")
            input = open(path+'/saved_data/marked_nodes_k=' + str(k1) + '_alpha=' + str(alpha1) + '_both=' + str(both1) +"_(index=" + str(index)+ ')_species='+spe+'.txt', 'r')
            alpha_maps[index] = []
            for line in input:
                alpha_maps[index].append(eval(line))
            if len(alpha_maps[index]) == 0:
                alpha_maps[index] = {}
            else:
                alpha_maps[index] = alpha_maps[index][0]
            print("     Finished reading file 'saved_data/marked_nodes_k=" + str(k1) + "_alpha=" + str(
                alpha1) + ".txt'")

        for index1 in range(index_start, index_end):
            for index2 in range(index_start, index_end):
                for key in alpha_maps[index1]:
                    if not key in alpha_maps[index2]:
                        print('        %s is marked in the %s solution and not in the %s solution\n' % (
                        str(key), str(index1), str(index2)))
    else:
        print("     Reading file 'saved_data/marked_nodes_k=" + str(k1) + "_alpha=" +str(alpha1) + ".txt'...")
        input = open(path + '/saved_data/marked_nodes_k=' + str(k1) + '_alpha=' +str(alpha1) + '_both=' + str(both1) + '.txt', 'r')
        alpha1_map = []
        for line in input:
            alpha1_map.append(eval(line))
        if len(alpha1_map) == 0:
            alpha1_map = {}
        else:
            alpha1_map = alpha1_map[0]
        print("     Finished reading file 'saved_data/marked_nodes_k=" + str(k1)  + "_alpha=" +str(alpha1) + ".txt'")

        print("     Reading file 'saved_data/marked_nodes_k=" + str(k2) + "_alpha=" + str(alpha2) + '_both=' + str(boTH_both) + ".txt'...")
        input = open(path+ '/saved_data/marked_nodes_k=' + str(k2)  + '_alpha=' + str(alpha2) +  '_both=' + str(boTH_both) +'.txt', 'r')
        alpha2_map = []
        for line in input:
            alpha2_map.append(eval(line))
        if len(alpha2_map) == 0:
            alpha2_map = {}
        else:
            alpha2_map = alpha2_map[0]
        print("     Finished reading file 'saved_data/marked_nodes_k=" + str(k2) + "_alpha=" + str(
            alpha2) + ".txt'")

        for key in alpha1_map:
            if not key in alpha2_map:
                print ('        %s is marked for k=%s, alpha=%s, both=%s and not for k=%s, alpha=%s, both=%s\n' % (str(key), str(k1),str(alpha1),str(both1),str(k2),str(alpha2),str(boTH_both)))
        for key1 in alpha2_map:
            if not key1 in alpha1_map:
                print ('        %s is marked for k=%s, alpha=%s, both=%s and not for k=%s, alpha=%s, both=%s\n' % (str(key1), str(k2),str(alpha2),str(boTH_both),str(k1),str(alpha1),str(both1)))

    quit()


