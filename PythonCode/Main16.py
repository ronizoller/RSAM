on_lab = False
check_diffreance_between_solutions = True
real_data = False
evolutinary_event = ['HT']

if on_lab:
    if check_diffreance_between_solutions:
        path  = '/storage/DATA/users/ronizo/compare_test'
    elif real_data and 'HT' in evolutinary_event:
        path = 'PycharmProjects/RSAM/COG2602'
    elif real_data and 'D' in evolutinary_event:
        path = 'PycharmProjects/RSAM/COG3620'
    else:
        path = '/storage/DATA/users/ronizo/noise_data_500_k=100'
else:
    import sys
    sys.path.append('/PycharmProjects/RSAM_venv/lib/python3.6/site-packages/graphviz/')
    if check_diffreance_between_solutions:
        path = '/Users/ronizoller/PycharmProjects/TreeReconciliation/trees/compare_test'
    elif real_data and 'HT' in evolutinary_event:
        path = '/Users/ronizoller/PycharmProjects/TreeReconciliation/trees/COG2602'
    elif real_data and 'D' in evolutinary_event:
        path = '/Users/ronizoller/PycharmProjects/TreeReconciliation/trees/COG3620'
    else:
        path = '/Users/ronizoller/PycharmProjects/TreeReconciliation/trees/duplications_test'


import networkx as nx
import dendropy as tr
import math
import utiles
import tree_operations_v1 as tree_operations
import inits_v1 as inits
import EfficiantVersion as hypergraph
import pattern_identify_v4 as pattern_identify
from multiprocessing import Pool
from datetime import datetime
import random
import os

speciesTreespecification = 'all'
test = False                                         # if True all data will be loaded from outter files, otherwise all data will be calculated and saved
glob = False                                        # if True global alignment is used, otherwise local
compare_subtrees = False                             # if true the algorithm will look for a signi different between two children of u in G, otherwise it will look for u in G s.t. in G(u) there are alot of same color HT
dis_flag = True                                     #count the patterns and take in count the distance of the HT
k = 100
exact_names = True

pattern = "same_color"

TH_mostly_speciations = 0
min_spesiction_events = 0

HT_cost = 1
D_cost = 1
S_cost = 0
save_data = False

planted_vertices = []
number_of_planted_vertices = 5

if not real_data:
    input = open(path + '/saved_data/planted_nodes_correct_names.txt', 'r')
    for line in input:
        planted_vertices.append(eval(line))
    planted_vertices = planted_vertices[0]
random_for_prec = 1
gamma = 1                                           # factor for probability assignment
alpha = 1                                           # factor for HT counting in the coloring stage
accur = 5                                           # calculations acuuracy
noise_level_list = [5]
p = 0.05                                            #p_value

#compare several optimal solutions
if check_diffreance_between_solutions:
    iterations = 1                                  #will check for each i=0 to iteration the solution i*factor
    factor = 1
    real_data = False

big_size = 2000
small_size = 7

def find_Pattern(H, S,S_dis_matrix, nCr_lookup_table, fact_lookup_table, red_HT_vertices_in_G, black_HT_vertices_in_G, red_doup,black_doup, pattern, evolutinary_event,S_colors):
    total_red = S_colors[S.seed_node.label][0]
    total_black = S_colors[S.seed_node.label][1]
    num_of_leafs = tree_operations.number_of_leafs(S,'S')
    Pr_red = total_red / num_of_leafs
    Pr_black = total_black / num_of_leafs
    for nd in (list(nx.topological_sort(H))):
        nd_index = nd
        incoming_edges = H.in_edges([nd], data=True)
        nd = H.nodes(data=True)[nd]
        for i in range(0,len(nd['l'])):
            curr = nd['l'][i]
            if 'HT' in evolutinary_event and curr['event'] == 'HT' and curr['probability'] > 0:
                x = S.find_node(lambda n: (n.label == curr['t']))
                incoming_edges_curr = [e for e in incoming_edges if e[2]['target'] == i]
                horizontally_trans_to_option1 = S.find_node(
                    lambda n: (n.label == H.nodes(data=True)[incoming_edges_curr[0][0]]['t']))
                horizontally_trans_to_option2 = S.find_node(
                    lambda n: (n.label == H.nodes(data=True)[incoming_edges_curr[1][0]]['t']))
                if ((tree_operations.is_not_ancestor(horizontally_trans_to_option1, x)) and (
                tree_operations.is_not_ancestor(x, horizontally_trans_to_option1))):
                    HT_to = horizontally_trans_to_option1
                    HT_to_in_G = H.nodes(data=True)[incoming_edges_curr[0][0]]
                else:
                    HT_to = horizontally_trans_to_option2
                    HT_to_in_G = H.nodes(data=True)[incoming_edges_curr[1][0]]

                curr_red = S_colors[x.label][0]
                curr_blacks = S_colors[x.label][1]
                total_curr = curr_red + curr_blacks

                HT_to_reds = S_colors[HT_to.label][0]
                HT_to_blacks =  S_colors[HT_to.label][1]
                total_HT_to = HT_to_reds + HT_to_blacks

                if (pattern == "same_color" or pattern == "different_colors" or pattern == "only_red" or pattern == "only_red_to_black" or pattern == "only_black_to_red"):
                    p_value_curr_red, nCr_lookup_table, fact_lookup_table = utiles.p_value_calculation(curr_red, total_curr, nCr_lookup_table, fact_lookup_table, accur, Pr_red,
                                    Pr_black, x.label)
                    p_value_curr_black, nCr_lookup_table, fact_lookup_table = utiles.p_value_calculation(curr_blacks,
                                                                                                       total_curr,
                                                                                                       nCr_lookup_table,
                                                                                                       fact_lookup_table,
                                                                                                       accur, Pr_black,
                                                                                                        Pr_red,
                                                                                                       x.label)
                    if total_HT_to > 0:
                        p_value_HT_to_red, nCr_lookup_table, fact_lookup_table = utiles.p_value_calculation(
                            math.ceil(HT_to_reds), total_HT_to, nCr_lookup_table, fact_lookup_table, accur, Pr_red,
                            Pr_black,x.label)
                        p_value_HT_to_black, nCr_lookup_table, fact_lookup_table = utiles.p_value_calculation(
                            math.ceil(HT_to_blacks), total_HT_to, nCr_lookup_table, fact_lookup_table, accur, Pr_black,
                            Pr_red, x.label)

                    #print('        curr = %s curr_red = %s, curr_black = %s, p_value_curr_red = %s, p_value_curr_black = %s' % (
                    #str(curr), str(curr_red), str(curr_blacks),str(p_value_curr_red),str(p_value_curr_black)))
                    #print('        HT_to = %s, HT_to.E_red = %s, HT_to.E_black = %s, HT_to p_value_red = %s, HT_to p_value_black = %s\n' % (
                    #str(HT_to), str(S_colors[HT_to.label][0]), str(S_colors[HT_to.label][1]),str(p_value_HT_to_red),str(p_value_HT_to_black)))

                    if p_value_curr_red < p: #HT from red
                        if (pattern == "same_color" or pattern == "only_red") and (p_value_HT_to_red < p): #HT to red
                            red_HT_vertices_in_G.append({'curr': curr, 'HT_to_in_G': HT_to_in_G,
                                                       'probability': curr['probability'], 'distance' : S_dis_matrix[(curr['t'], HT_to.label)]})
                            print('     this is a good pattern (red to red): \n       from %s to %s (HT to %s), Total vertices in S below %s: %s, below HT to: %s\n       p value red of %s = %s, p value black = %s (Pr red = %s, Pr black = %s)\n       p value red of HT to = %s,  p value black of HT to = %s.\n       Distance = %s \n**\n'%
                               (str(curr['s']),curr['t'],str(HT_to.label),str(curr['t']),str(total_curr),str(total_HT_to),str(curr['t']),str(p_value_curr_red) ,str(p_value_curr_black) ,str(Pr_red),str(Pr_black),str(p_value_HT_to_red),str(p_value_HT_to_black),str(S_dis_matrix[(curr['t'], HT_to.label)])))
                    elif p_value_curr_black < p: #HT from black
                        if (pattern == "same_color") and (p_value_HT_to_black < p): #HT to black
                            black_HT_vertices_in_G.append({'curr': curr, 'HT_to_in_G': HT_to_in_G,
                                                       'probability': curr['probability'], 'distance' : S_dis_matrix[(curr['t'], HT_to.label)]})
                            print('     this is a good pattern (black to black): \n       from %s to %s (HT to %s), Total vertices in S below %s: %s, below HT to: %s\n       p value red of %s = %s, p value black = %s (Pr red = %s, Pr black = %s)\n       p value red of HT to = %s, p value black of HT to = %s.\n       Distance = %s \n**\n'%
                                (str(curr['s']),curr['t'],str(HT_to.label),str(curr['t']),str(total_curr),str(total_HT_to),str(curr['t']),str(p_value_curr_red),str(p_value_curr_black) ,str(Pr_red),str(Pr_black),str(p_value_HT_to_red),str(p_value_HT_to_black),str(S_dis_matrix[(curr['t'], HT_to.label)])))
                elif pattern == "any":               #mark HT (in the red list, just for convi)
                    red_HT_vertices_in_G.append({'curr': curr, 'HT_to_in_G': HT_to_in_G,
                                                 'probability': curr['probability'],
                                                 'distance': S_dis_matrix[(curr['t'], HT_to.label)]})
            if 'D' in evolutinary_event and curr['event'] == 'D' and curr['probability'] > 0:
                doupli,HT,spe,total = hypergraph.mostly_speciation_event_in_subtree(H,nd_index,i)
                print('curr = %s\n       spe: %s, HT: %s, doupli:%s,total: %s\n' %
                      (str(curr), str(spe), str(HT), str(doupli),str(total)))
                total = total -1
                if total > 0:
                    if (spe+doupli)/total >= TH_mostly_speciations and total > min_spesiction_events:
                        red_doup.append({'curr': curr,'probability': curr['probability']})

    return red_HT_vertices_in_G, black_HT_vertices_in_G,red_doup,[], nCr_lookup_table, fact_lookup_table

def init_G_nodes_to_weight(G, G_nodes_to_weight):
    for u in G.postorder_node_iter():
        for v in G.postorder_node_iter():
            G_nodes_to_weight.update({(u.label,v.label):[0,0]})
    return G_nodes_to_weight

def weight_HT_in_G(H, G, G_edges_to_weight, S_dis_matrix):       #
    #print('weighting G...')
    G_nodes_to_weight = init_G_nodes_to_weight(G, G_edges_to_weight)

    for nd in (list(nx.topological_sort(H))):
        incoming_edges = H.in_edges([nd], data=True)
        nd = H.nodes(data=True)[nd]
        for i in range(0,len(nd['l'])) :
            curr = nd['l'][i]
            if curr['event'] == 'HT' and curr['probability'] > 0:
                x = S.find_node(lambda n : (n.label == curr['t']))
                u = G.find_node(lambda n: (n.label == curr['s']))

                incoming_edges_curr = [e for e in incoming_edges if e[2]['target'] == i]
                horizontally_trans_to_option1 = S.find_node(lambda n : (n.label == H.nodes(data=True)[incoming_edges_curr[0][0]]['t']))
                horizontally_trans_to_option2 = S.find_node(lambda n : (n.label == H.nodes(data=True)[incoming_edges_curr[1][0]]['t']))
                if ((tree_operations.is_not_ancestor(horizontally_trans_to_option1,x)) and (tree_operations.is_not_ancestor(x,horizontally_trans_to_option1))):
                    HT_to = horizontally_trans_to_option1
                    HT_to_in_G = H.nodes(data=True)[incoming_edges_curr[0][0]]['s']
                else :
                    HT_to = horizontally_trans_to_option2
                    HT_to_in_G = H.nodes(data=True)[incoming_edges_curr[1][0]]['s']

                number_of_HT = G_nodes_to_weight[(u.label,HT_to_in_G)][0]
                total_weight = G_nodes_to_weight[(u.label,HT_to_in_G)][1]

                number_of_HT += 1
                total_weight += S_dis_matrix[(curr['t'],HT_to.label)]

                G_nodes_to_weight.update({(u.label,HT_to_in_G): [number_of_HT, total_weight]})
    #print('Finished weighting G...')
    return G_nodes_to_weight

def RSAM_finder_multithread(parameters):
    noise_in = ''
    real_data = parameters[0]
    noise_level = parameters[1]
    check_diffreance_between_solutions = parameters[3]

    list_of_scores_for_rand_num = {}
    max_score_TH = []
    max_score_doup = []
    random_for_prec_curr = random_for_prec
    for rand_num in range (0, random_for_prec_curr):
        if not check_diffreance_between_solutions:
            if not real_data:
                noise_in = parameters[4]
                path_change_in = path + '/' + noise_in
                os.makedirs(os.path.dirname(path_change_in), exist_ok=True)
                path_change_in = path_change_in + '/saved_data'
                os.makedirs(os.path.dirname(path_change_in), exist_ok=True)
                path_change_in = path +  '/' + noise_in
            else:
                path_change_in = path
            TH_edges_in_subtree = parameters[5]                                                    # smallest subtree that will be counted when not comparing subtrees
            TH_compare_subtrees = parameters[7]
            S_dis_matrix = {}
            nodes_table = {}
            S_colors = {}
            all_vertices = {}
            new_G = nx.DiGraph()
            red_HT_vertices_in_G = []
            black_HT_vertices_in_G = []
            red_doup = []
            black_doup = []
            nCr_lookup_table = {}
            fact_lookup_table = {}

            G = tr.Tree.get_from_path(path  + "/GeneTree(binary)_local.txt", schema="newick")
            S = tr.Tree.get_from_path(path  +"/phyliptree(binary,"+speciesTreespecification+").phy", schema="newick")

            input = open(path_change_in +'/'+ str(noise_level)+'/sigma'+str(noise_level)+'.'+str(rand_num)+'.txt', 'r')
            sigma = []
            for line in input:
                sigma.append(eval(line))
            sigma = sigma[0]

            input = open(path_change_in + '/'+ str(noise_level)+'/colors'+str(noise_level)+'.'+str(rand_num)+'.txt', 'r')
            colors = []
            for line in input:
                colors.append(eval(line))
            colors = colors[0]
            G.prune_taxa_with_labels(tree_operations.remove_unsigma_genes(G, sigma, True))

            S = utiles.init_internal_labels(S, 'x', sigma, path)
            G = utiles.init_internal_labels(G, 'u', sigma, path)

            G = tree_operations.collapse_edges(G)
            S = tree_operations.collapse_edges(S)

            S_labels_table, G_labels_table, sigma = inits.init_taxon_to_label_table(S,G,sigma)
            sigma, old_sigma = inits.update_sigma(S, G, k, sigma, test, path_change_in,exact_names,S_labels_table,G_labels_table)
            G.prune_taxa_with_labels(tree_operations.remove_unsigma_genes(G, sigma, False))
            colors,old_colors = inits.update_colors(S, colors,exact_names)


            S_dis_matrix = inits.init_distance_S(S_dis_matrix, k, test, path,speciesTreespecification)
            nodes_table = inits.init_nodes_table(S, G, nodes_table)

            H, H_number_of_nodes, nodes_table = hypergraph.build_hyper_garph(S, G, test, k,
                                                                             nodes_table, D_cost, S_cost, HT_cost, path_change_in, alpha,
                                                                             sigma, save_data)

            H, max_prob = hypergraph.assign_probabilities(S, G, H,gamma)
        else:
            S = parameters[4]
            G = parameters[5]
            S_colors = parameters[6]
            colors = parameters[7]
            deleted_nodes = parameters[8]
            S_dis_matrix = parameters[9]
            nCr_lookup_table = parameters[10]
            fact_lookup_table = parameters[11]
            red_HT_vertices_in_G = parameters[12]
            black_HT_vertices_in_G = parameters[13]
            sigma = parameters[14]
            new_G = parameters[15]
            all_vertices = parameters[17]
            TH_edges_in_subtree = parameters[19]
            H = parameters[21]
            TH_compare_subtrees = parameters[22]
            red_doup = parameters[23]
            black_doup = parameters[24]

        if H is None:
            list_of_scores_for_rand_num.update({rand_num: {}})
        else:
            ##      PROBABILITIES, COLORS, PATTERN      ##

            S_colors = tree_operations.color_tree(S, 'S', S_colors, colors, sigma)

            red_HT_vertices_in_G, black_HT_vertices_in_G,red_doup,black_doup, nCr_lookup_table, fact_lookup_table = find_Pattern(H,S,
                                                                                                                   S_dis_matrix,
                                                                                                                   nCr_lookup_table,
                                                                                                                   fact_lookup_table,
                                                                                                                   red_HT_vertices_in_G,
                                                                                                                   black_HT_vertices_in_G,red_doup,black_doup, pattern,evolutinary_event,S_colors)

            max_S_d_of_HT = tree_operations.find_max_d_of_HT(S_dis_matrix, red_HT_vertices_in_G, black_HT_vertices_in_G,evolutinary_event)

            new_G= tree_operations.weight_G_based_on_same_color_HT(G, new_G, red_HT_vertices_in_G,
                                                                    black_HT_vertices_in_G,red_doup,black_doup, max_S_d_of_HT,dis_flag,evolutinary_event,False,k)
            new_G = tree_operations.number_of_edges_in_subtree(new_G)
            new_G = tree_operations.normlize_weights(new_G,k)
            if not check_diffreance_between_solutions:
                max_score_TH,max_score_doup = tree_operations.find_max_scores(new_G, number_of_planted_vertices,TH_edges_in_subtree,TH_compare_subtrees)
            marked_nodes,all_vertices = pattern_identify.find_signi_distance(new_G, all_vertices, TH_compare_subtrees, k, 'D' in evolutinary_event,
                                                                compare_subtrees,TH_edges_in_subtree,max_score_TH,max_score_doup,check_diffreance_between_solutions,real_data)

            list_of_scores_for_rand_num.update({rand_num:all_vertices})
    return (utiles.average_of_list(list_of_scores_for_rand_num,random_for_prec_curr),noise_in)

def extract_and_tarce_a_solution(parameters):
    iter = parameters[0]
    new_G = parameters[1]
    max_dis = parameters[2]
    solutions = parameters[3]
    S_dis_matrix = parameters[4]
    nCr_lookup_table = parameters[5]
    fact_lookup_table = parameters[6]
    red_HT_vertices_in_G = parameters[7]
    black_HT_vertices_in_G = parameters[8]
    S_colors = parameters[9]
    H = parameters[11]
    S = parameters[12]
    G = parameters[13]
    TH_edges_in_subtree = parameters[14]
    TH_compare_subtrees = parameters[15]
    k = parameters[16]
    red_doup = parameters[17]
    black_doup = parameters[18]

    new_G[iter] = nx.DiGraph()
    deleted_nodes = []
    max_score_TH = []
    max_score_doup = []

    solutions[iter] = nx.DiGraph()
    H_root = [nd for nd in list(H.node(data=True)) if
              nd[1]['s'] == G.seed_node.label and nd[1]['t'] == S.seed_node.label]
    solutions[iter], nodes_table = hypergraph.track_a_solution(H_root, H, S, G, solutions[iter], k-1)

    print('     Writing nodes...')
    file = open(path + '/saved_data/H_nodes_effi.txt', 'w')
    file.write(str(solutions[iter].nodes(data=True)))
    file.close()
    print('     Finished writing nodes.\n')

    solutions[iter], max_prob = hypergraph.assign_probabilities(S, G, solutions[iter], gamma)

    red_HT_vertices_in_G, black_HT_vertices_in_G,red_doup,black_doup, nCr_lookup_table, fact_lookup_table = find_Pattern(
        solutions[iter], S,S_dis_matrix, nCr_lookup_table, fact_lookup_table, red_HT_vertices_in_G,
        black_HT_vertices_in_G,red_doup,black_doup, pattern, evolutinary_event, S_colors)
    max_S_d_of_HT = tree_operations.find_max_d_of_HT(S_dis_matrix, red_HT_vertices_in_G, black_HT_vertices_in_G,
                                                     evolutinary_event)

    new_G[iter] = tree_operations.weight_G_based_on_same_color_HT(G, new_G[iter], red_HT_vertices_in_G,
                                                                  black_HT_vertices_in_G,red_doup,black_doup, max_S_d_of_HT, dis_flag,
                                                                  evolutinary_event, check_diffreance_between_solutions,
                                                                  k)
    new_G[iter] = tree_operations.number_of_edges_in_subtree(new_G[iter])
    new_G[iter] = tree_operations.normlize_weights(new_G[iter], 1)
    if not check_diffreance_between_solutions:
        max_score_TH, max_score_doup = tree_operations.find_max_scores(new_G[iter], number_of_planted_vertices,TH_edges_in_subtree,TH_compare_subtrees)
    all_vertices = {}
    marked_nodes, all_vertices = pattern_identify.find_signi_distance(new_G[iter], all_vertices, TH_compare_subtrees,
                                                                      k, 'D' in evolutinary_event,
                                                                      compare_subtrees,
                                                                      TH_edges_in_subtree,max_score_TH,max_score_doup,check_diffreance_between_solutions,real_data)

    return([list(all_vertices.items()),list(marked_nodes.items()),new_G[iter],'(%s,%s,%s)' % (str(TH_compare_subtrees),0,str(TH_edges_in_subtree))])

##********  MAIN ***********

def main():
    global S, G, iterations, sigma, alpha, gamma, colors,planted_vertices, number_of_planted_vertices, both, path, speciesTreespecification, evolutinary_event,exact_names,noise_level_list,random_for_prec
    starting_time = datetime.now()

    all_vertices_with_index = {}
    all_RSAM_marked = {}
    all_RSAM_unmarked = {}
    list_of_scores_for_rand_num = {}
    noise_level = 0
    rand_num = 0
    random_for_prec_curr = 1
    print('                 ****     Iteration %s.%s/%s.%s     ***** \n' % (str(noise_level),str(rand_num),str(noise_level_list[len(noise_level_list)-1]),str(random_for_prec)))
    S_dis_matrix = {}
    nodes_table = {}
    S_colors = {}
    all_vertices = {}

    G = tr.Tree.get_from_path(path + "/GeneTree(binary)_local.txt", schema="newick")
    S = tr.Tree.get_from_path(path+"/phyliptree(binary,"+speciesTreespecification+").phy", schema="newick")
    if check_diffreance_between_solutions:
        input = open(path + '/colors_and_HT/' + str(noise_level) + '/sigma' + str(noise_level) + '.' + str(rand_num) + '.txt', 'r')
    else:
        input = open(path+'/'+ str(noise_level)+'/sigma'+str(noise_level)+'.'+str(rand_num)+'.txt', 'r')
    sigma = []
    for line in input:
        sigma.append(eval(line))
    sigma = sigma[0]

    if check_diffreance_between_solutions:
        input = open(path + '/colors_and_HT/' + str(noise_level) + '/colors' + str(noise_level) + '.' + str(rand_num) + '.txt', 'r')
    else:
        input = open(path + '/'+ str(noise_level)+'/colors'+str(noise_level)+'.'+str(rand_num)+'.txt', 'r')
    colors = []
    for line in input:
        colors.append(eval(line))
    colors = colors[0]
    G.prune_taxa_with_labels(tree_operations.remove_unsigma_genes(G, sigma, True))

    S = utiles.init_internal_labels(S, 'x', sigma, path)
    G = utiles.init_internal_labels(G, 'u', sigma, path)

    G = tree_operations.collapse_edges(G)
    S = tree_operations.collapse_edges(S)

    S_labels_table, G_labels_table,sigma = inits.init_taxon_to_label_table(S,G,sigma)
    sigma, old_sigma = inits.update_sigma(S, G, k, sigma, test, path,exact_names,S_labels_table,G_labels_table)
    G.prune_taxa_with_labels(tree_operations.remove_unsigma_genes(G, sigma, False))
    colors,old_colors = inits.update_colors(S, colors,exact_names)
    if 'D' in evolutinary_event:
        TH_compare_subtrees = 0
        TH_edges_in_subtree = 30
    else:
        TH_compare_subtrees = 2.5
        TH_edges_in_subtree = 50
    #if not on_lab:
        #draw.draw_S_and_G(S, G, old_sigma, colors, sigma, path, None, speciesTreespecification)
    S_dis_matrix = inits.init_distance_S(S_dis_matrix, k, test, path,speciesTreespecification)
    nodes_table = inits.init_nodes_table(S, G, nodes_table)
    H, H_number_of_nodes, nodes_table = hypergraph.build_hyper_garph(S, G, test, k,
                                                                     nodes_table, D_cost, S_cost, HT_cost, path, alpha,
                                                                     sigma,save_data)
    if check_diffreance_between_solutions:
        max_dis = tree_operations.max_dis(S_dis_matrix)
        if iterations * factor < k:
            all_marked_for_TH = {}
            all_unmarked_for_TH = {}
            H, max_prob = hypergraph.assign_probabilities(S, G, H, gamma)
            if H == None:
                quit()
            parameters = []
            p = Pool(15)
            combined = [(1,0,0)]
            for i in range(0, len(combined)):
                TH_compare_subtrees = combined[i][0]
                TH_both = combined[i][1]
                TH_edges_in_subtree = combined[i][2]
                parameters.append([real_data,noise_level_list[0], TH_both, check_diffreance_between_solutions, S, G, S_colors, colors, [], S_dis_matrix, {}, {},
                                   [], [], sigma, nx.DiGraph(), {}, all_vertices, 0, TH_edges_in_subtree,i, H,TH_compare_subtrees,[],[]])
            list_of_RSAM_results = p.map(RSAM_finder_multithread, parameters)
            p.close()
            p.join()
            ind = 0
            for res in list_of_RSAM_results:
                all_vertices_with_index.update({combined[ind]: [res[0]]})
                all_RSAM_unmarked.update({combined[ind]: utiles.find_unmarked(res[0],G,True)})
                ind += 1
            file = open(path + '/saved_data/all_marked_vertices_RSAM_finder.txt', 'w')
            file.write(str(all_vertices_with_index))
            file.close()
            file = open(path + '/saved_data/all_unmarked_vertices_RSAM_finder.txt', 'w')
            file.write(str(all_RSAM_unmarked))
            file.close()

            new_G = {}
            solutions = {}
            red_HT_vertices_in_G = []
            black_HT_vertices_in_G = []
            red_doup = []
            black_doup = []
            nCr_lookup_table = {}
            fact_lookup_table = {}
            all_vertices_with_index = {}

            S_colors = tree_operations.color_tree(S, 'S', S_colors, colors, sigma)
            p1 = Pool(15)
            parameters_list = [(x,new_G,max_dis,solutions,S_dis_matrix,nCr_lookup_table,fact_lookup_table,red_HT_vertices_in_G,
                                black_HT_vertices_in_G,S_colors,None,H,S,G,TH_edges_in_subtree,TH_compare_subtrees, k,red_doup,black_doup,None)
                               for x in range(0,iterations) for (TH_compare_subtrees,TH_both,TH_edges_in_subtree) in combined]
            list_of_results = p1.map(extract_and_tarce_a_solution, parameters_list)

            p1.close()
            p1.join()
            new_G_to_save = []
            for res in list_of_results:
                all_vertices_with_index.update({res[3]:res[0]})
                new_G_to_save.append(res[2])
                all_marked_for_TH.update({res[3]:res[1]})
                all_unmarked_for_TH.update({res[3]:utiles.find_unmarked(all_vertices_with_index[res[3]],G,False)})
            file = open(path + '/saved_data/all_vertices_TH.txt', 'w')
            file.write(str(all_vertices_with_index))
            file.close()
            file = open(path + '/saved_data/all_unmarked_nodes_for_TH.txt', 'w')
            file.write(str(all_unmarked_for_TH))
            file.close()
            quit()
        else:
            print('**   not enough solutions were calculated in order to track solution number %s   **' % str(iterations * factor))
        quit()

    else:
        new_G = nx.DiGraph()
        red_HT_vertices_in_G = []
        black_HT_vertices_in_G = []
        red_doup = []
        black_doup = []
        nCr_lookup_table = {}
        fact_lookup_table = {}
        max_score_TH = []
        max_score_doup = []
        H, max_prob = hypergraph.assign_probabilities(S, G, H, gamma)
        if H is None:
            list_of_scores_for_rand_num.update({rand_num: {}})
        else:
            S_colors = tree_operations.color_tree(S, 'S', S_colors, colors, sigma)
            red_HT_vertices_in_G, black_HT_vertices_in_G,red_doup,black_doup, nCr_lookup_table, fact_lookup_table = find_Pattern(H,S,
                                                                                                                   S_dis_matrix,
                                                                                                                   nCr_lookup_table,
                                                                                                                   fact_lookup_table,
                                                                                                                   red_HT_vertices_in_G,
                                                                                                                   black_HT_vertices_in_G,red_doup,black_doup, pattern,evolutinary_event,S_colors)
            max_S_d_of_HT = tree_operations.find_max_d_of_HT(S_dis_matrix, red_HT_vertices_in_G, black_HT_vertices_in_G,evolutinary_event)
            new_G = tree_operations.weight_G_based_on_same_color_HT(G, new_G, red_HT_vertices_in_G,
                                                                          black_HT_vertices_in_G,red_doup,black_doup, max_S_d_of_HT,dis_flag,evolutinary_event,check_diffreance_between_solutions,k)
            new_G = tree_operations.number_of_edges_in_subtree(new_G)
            new_G = tree_operations.normlize_weights(new_G, k)
            if not check_diffreance_between_solutions:
                max_score_TH,max_score_doup = tree_operations.find_max_scores(new_G, number_of_planted_vertices,TH_edges_in_subtree,TH_compare_subtrees)
            marked_nodes,all_vertices = pattern_identify.find_signi_distance(new_G, all_vertices, TH_compare_subtrees, k, 'D' in evolutinary_event,
                                                                compare_subtrees,TH_edges_in_subtree,max_score_TH,max_score_doup,check_diffreance_between_solutions,real_data)


            list_of_scores_for_rand_num.update({rand_num:marked_nodes})
            #if not on_lab:
                #draw.draw_new_G(marked_nodes, colors, sigma, new_G, G, old_sigma, k, TH_compare_subtrees,
                            #path, True, glob, speciesTreespecification, pattern,
                             #big_size, evolutinary_event, compare_subtrees, 1,S_labels_table)
    all_vertices_with_index.update({noise_level:utiles.average_of_list(list_of_scores_for_rand_num,random_for_prec_curr)})

    parameters = []
    if real_data:
        file = open(path + '/saved_data/marked_RSAM_finder_'+speciesTreespecification+'.txt', 'w')
        file.write(str(marked_nodes))
        file.close()
        for nd,x in marked_nodes.items():
            r,l = tree_operations.leaf_in_subtrees(G,'S',nd, old_sigma,False)
            print('For %s:\nlist = %s\n' % (str(nd),str(r+l)))

        print('Running time: %s\nTH_compare: %s\nk: %s\nTH_edges: %s' % (str(datetime.now()-starting_time),str(TH_compare_subtrees),str(k),str(TH_edges_in_subtree)))
    else:
        for noise_in in ['colors_and_HT','color','HT']:
            for i in range(0,len(noise_level_list)):
                parameters.append([real_data,noise_level_list[i],None,check_diffreance_between_solutions,noise_in,TH_edges_in_subtree,None,TH_compare_subtrees])
        p = Pool(15)
        list_of_results = p.map(RSAM_finder_multithread, parameters)
        p.close()
        p.join()
        temp_j = 0
        for noise_in in ['colors_and_HT', 'color', 'HT']:
            ind = 0
            for j in range(temp_j,len(noise_level_list)+temp_j):
                all_vertices_with_index.update({noise_level_list[ind]:list_of_results[j][0]})
                ind += 1
            file = open(path + '/saved_data/all_vertices_RSAM_finder_'+noise_in+'.txt', 'w')
            file.write(str(all_vertices_with_index))
            file.close()
            temp_j = len(noise_level_list)+temp_j
        print('Running time: %s\nTH_compare: %s\nk: %s\nTH_edges: %s' % (str(datetime.now()-starting_time),str(TH_compare_subtrees),str(k),str(TH_edges_in_subtree)))

if __name__ == "__main__":
    main()


