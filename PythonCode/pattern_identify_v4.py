import math
import networkx as nx
import utiles
import tree_operations_v1 as tree_operations

start = '*'
end = '*'

def identify_pattern2(G, H, k, G_nodes_to_weight, G_nodes_identified):
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

            if not p2[0]:
                if u['edges_in_subtree'] >= p1[3]:
                    if check_diffrence_between_solutions:
                        all_vertices.update({u['label']: u_p1_score})
                    elif u_p1_score in max_score_p1_list:
                        print(u)
                        marked_nodes.update({u['label']: [u_p1_score,'Single-mode']})
            else:
                if v['edges_in_subtree'] >= p1[3] and w['edges_in_subtree'] >= p1[3]:
                        if check_diffrence_between_solutions:
                            all_vertices.update({u['label']: u_p1_score})
                        elif v_p1_score + w_p2_score in max_score_p1_and_p2_list:
                                marked_nodes.update({u['label']: [v_p1_score,w_p2_score, 'Double-mode,'+str(v['label']+'_p1')]})
                        elif w_p1_score + v_p2_score in max_score_p1_and_p2_list:
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


def RSAM_finder_multithread(parameters):
    noise_in = ''

    list_of_scores_for_rand_num = {}
    max_score_p1_list = []
    max_score_p1_and_p2_list = []
    random_for_prec_curr = random_for_prec
    for rand_num in range(0, random_for_prec_curr):

        noise_level, noise_in, p1, p2 = parameters
        path_change_in = path + '/' + noise_in
        os.makedirs(os.path.dirname(path_change_in), exist_ok=True)
        path_change_in = path_change_in + '/saved_data'
        os.makedirs(os.path.dirname(path_change_in), exist_ok=True)
        path_change_in = path + '/' + noise_in

        S_dis_matrix = {}
        nodes_table = {}
        S_colors = {}
        all_vertices = {}
        new_G = nx.DiGraph()

        nCr_lookup_table = {}
        fact_lookup_table = {}

        G = tr.Tree.get_from_path(path + "/GeneTree(binary," + geneExt + ")_local.txt", schema="newick")
        S = tr.Tree.get_from_path(path + "/phyliptree(binary," + speciesTreespecification + ").phy", schema="newick")

        input = open(
            path_change_in + '/' + str(noise_level) + '/sigma' + str(noise_level) + '.' + str(rand_num) + '.txt', 'r')
        sigma = []
        for line in input:
            sigma.append(eval(line))
        sigma = sigma[0]

        input = open(
            path_change_in + '/' + str(noise_level) + '/colors' + str(noise_level) + '.' + str(rand_num) + '.txt', 'r')
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
        sigma, old_sigma = inits.update_sigma(S, G, k, sigma, test, path_change_in, exact_names, S_labels_table,
                                              G_labels_table)
        G.prune_taxa_with_labels(tree_operations.remove_unsigma_genes(G, sigma, False))
        colors, old_colors = inits.update_colors(S, colors, exact_names)

        S_dis_matrix = inits.init_distance_S(S_dis_matrix, k, test, path, speciesTreespecification)
        nodes_table = inits.init_nodes_table(S, G, nodes_table)

        H, H_number_of_nodes, nodes_table = hypergraph.build_hyper_garph(S, G, test, k,
                                                                         nodes_table, D_cost, S_cost, HT_cost,
                                                                         path_change_in, alpha,
                                                                         sigma, save_data)

        H, max_prob = hypergraph.assign_probabilities(S, G, H, gamma)
        if H is None:
            list_of_scores_for_rand_num.update({rand_num: {}})
        else:
            ##      PROBABILITIES, COLORS, PATTERN      ##

            S_colors = tree_operations.color_tree(S, 'S', S_colors, colors, sigma)

            interesting_vertices_p1, nCr_lookup_table, fact_lookup_table = find_Pattern(H, S, S_dis_matrix,
                                                                                        nCr_lookup_table,
                                                                                        fact_lookup_table, p1, S_colors)
            interesting_vertices_p2, nCr_lookup_table, fact_lookup_table = find_Pattern(H, S, S_dis_matrix,
                                                                                        nCr_lookup_table,
                                                                                        fact_lookup_table, p2, S_colors)

            max_S_d_of_HT = tree_operations.find_max_d_of_HT(S_dis_matrix, interesting_vertices_p1, p1)

            new_G = tree_operations.weight_G_based_on_same_color_HT(G, new_G, interesting_vertices_p1,
                                                                    interesting_vertices_p2, max_S_d_of_HT, p1, p2,
                                                                    False)

            new_G = tree_operations.number_of_edges_in_subtree(new_G)
            new_G = tree_operations.normlize_weights(new_G, k, p1, 'p1')
            new_G = tree_operations.normlize_weights(new_G, k, p2, 'p2')

            if p2[0] is None:
                max_score_p1_list = tree_operations.find_max_scores(new_G, number_of_planted_vertices, 'p1', p1[3])
            else:
                max_score_p1_and_p2_list = tree_operations.find_max_scores(new_G, number_of_planted_vertices, 'p2',
                                                                           p1[3])

            marked_nodes, all_vertices = pattern_identify.find_signi_distance(new_G, all_vertices, p1, p2,
                                                                              max_score_p1_list,
                                                                              max_score_p1_and_p2_list, False)

            list_of_scores_for_rand_num.update({rand_num: all_vertices})
    return (utiles.average_of_list(list_of_scores_for_rand_num, random_for_prec_curr), noise_in)


def extract_and_tarce_a_solution(parameters):
    iter, new_G, solutions, S_dis_matrix, nCr_lookup_table, fact_lookup_table, S_colors, H, S, G, TH_edges_in_subtree, k, p1, p2 = parameters

    new_G[iter] = nx.DiGraph()
    max_score_p1_list = []
    max_score_p1_and_p2_list = []

    solutions[iter] = nx.DiGraph()
    H_root = [nd for nd in list(H.node(data=True)) if
              nd[1]['s'] == G.seed_node.label and nd[1]['t'] == S.seed_node.label]
    solutions[iter], nodes_table = hypergraph.track_a_solution(H_root, H, S, G, solutions[iter],
                                                               random.choice(range(0, k)))

    solutions[iter], max_prob = hypergraph.assign_probabilities(S, G, solutions[iter], gamma)

    interesting_vertices_p1, nCr_lookup_table, fact_lookup_table = find_Pattern(solutions[iter], S, S_dis_matrix,
                                                                                nCr_lookup_table,
                                                                                fact_lookup_table, p1, S_colors)
    interesting_vertices_p2, nCr_lookup_table, fact_lookup_table = find_Pattern(solutions[iter], S, S_dis_matrix,
                                                                                nCr_lookup_table,
                                                                                fact_lookup_table, p2, S_colors)

    max_S_d_of_HT = tree_operations.find_max_d_of_HT(S_dis_matrix, interesting_vertices_p1, p1)

    new_G[iter] = tree_operations.weight_G_based_on_same_color_HT(G, new_G[iter], interesting_vertices_p1,
                                                                  interesting_vertices_p2, max_S_d_of_HT, p1, p2, False)
    new_G[iter] = tree_operations.number_of_edges_in_subtree(new_G[iter])
    new_G[iter] = tree_operations.normlize_weights(new_G[iter], 1, p1, 'p1')
    new_G[iter] = tree_operations.normlize_weights(new_G[iter], 1, p2, 'p2')
    if p2[0] is None:
        max_score_p1_list = tree_operations.find_max_scores(new_G, number_of_planted_vertices, 'p1', p1[3])
    else:
        max_score_p1_and_p2_list = tree_operations.find_max_scores(new_G, number_of_planted_vertices, 'p2', p1[3])
    all_vertices = {}
    marked_nodes, all_vertices = pattern_identify.find_signi_distance(new_G[iter], all_vertices, p1, p2,
                                                                      max_score_p1_list, max_score_p1_and_p2_list,
                                                                      False)

    return ([list(all_vertices.items()), list(marked_nodes.items()), new_G[iter], '(%s)' % (str(TH_edges_in_subtree))])
