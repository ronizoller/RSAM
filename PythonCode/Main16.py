on_lab = True
check_diffreance_between_solutions = True

if on_lab:
    if check_diffreance_between_solutions:
        path  = '/users/studs/bsc/2016/ronizo/PycharmProjects/RSAM/simulator_data/comparsion'
    else:
        path = '/storage/DATA/users/ronizo/noise_data'
else:
    import sys
    sys.path.append('/anaconda3/lib/python3.6/site-packages')
    if check_diffreance_between_solutions:
        path = '/Users/ronizoller/Documents/school/Master/מחקר/DATA/example'

import networkx as nx
import dendropy as tr
import math
import utiles
import tree_operations_v1 as tree_operations
import inits_v1 as inits
import EfficiantVersion as hypergraph
import pattern_identify_v4 as pattern_identify
import EfficiantVersion as effi
from multiprocessing import Pool
from datetime import datetime


speciesTreespecification = 'all'
test = False                                         # if True all data will be loaded from outter files, otherwise all data will be calculated and saved
glob = False                                        # if True global alignment is used, otherwise local
compare_subtrees = False                             # if true the algorithm will look for a signi different between two children of u in G, otherwise it will look for u in G s.t. in G(u) there are alot of same color HT
dis_flag = True                                     #count the patterns and take in count the distance of the HT
one_enriched_on_not = False
k = 100
exact_names = True

evolutinary_event = 'HT'

HT_cost = 1
D_cost = 1
S_cost = 0
save_data = False

number_of_planted_vertices = 1
planted_vertices = []
input = open(path + '/saved_data/planted_nodes_correct_names.txt', 'r')
for line in input:
    planted_vertices.append(eval(line))
planted_vertices = planted_vertices[0]
random_for_prec = 1
gamma = 1                                           # factor for probability assignment
alpha = 1                                           # factor for HT counting in the coloring stage
both = False
accur = 5                                           # calculations acuuracy
noise_level_list = [5]
TH_both = 0.8                                      # factor for not both
p = 0.05                                            #p_value
TH_compare_subtrees = 1

#compare several optimal solutions
if check_diffreance_between_solutions:
    iterations = 5                                  #will check for each i=0 to iteration the solution i*factor
    factor = 20

####FOR HT EVOLUTNARY EVENTS###
if compare_subtrees and evolutinary_event=='HT':
    TH_compare_subtrees =  1
    pattern = "same_color"
elif not compare_subtrees and evolutinary_event == 'HT':
    TH_compare_subtrees =  1
    #pattern = "any"
    #pattern = "different_colors"
    #pattern = "only_red"
    pattern = "same_color"
    #pattern = "only_red_to_black"
    #pattern = "only_black_to_red"

###FOR DUPLICATIONS EVENTS###
elif compare_subtrees and evolutinary_event == 'D':
    TH_compare_subtrees = 1
    pattern = "any"
elif not compare_subtrees and evolutinary_event == 'D':
    pattern = "any"
    TH_pattern_in_subtree = 0.0009

big_size = 2000                                  #size of nodes
small_size = 7

def find_Pattern(H, S,S_dis_matrix, nCr_lookup_table, fact_lookup_table, red_HT_vertices_in_G, black_HT_vertices_in_G, pattern, evolutinary_event,S_colors):
    #print('Find the pattern within the tree...')
    total_red = S_colors[S.seed_node.label][0]
    total_black = S_colors[S.seed_node.label][1]
    num_of_leafs = tree_operations.number_of_leafs(S,'S')
    Pr_red = total_red / num_of_leafs
    Pr_black = total_black / num_of_leafs
    for nd in (list(nx.topological_sort(H))):
        incoming_edges = H.in_edges([nd], data=True)
        nd = H.nodes(data=True)[nd]
        for i in range(0,len(nd['l'])):
            curr = nd['l'][i]
            if evolutinary_event == 'HT' and curr['event'] == 'HT' and curr['probability'] > 0:
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

                p_value_red_HT_to = 0   #to be red
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
                        if (pattern == "different_colors" or pattern == "only_red_to_black") and (p_value_HT_to_black < p): #HT to black
                            red_HT_vertices_in_G.append({'curr': curr, 'HT_to_in_G': HT_to_in_G,
                                                         'probability': curr['probability'],
                                                         'distance': S_dis_matrix[(curr['t'], HT_to.label)]})
                            #print('     this is a good pattern (red to black): \n       from %s to %s (HT to %s), Total vertices in S below %s: %s, below HT to: %s\n       p value red of %s = %s,  p value black = %s (Pr red = %s, Pr black = %s)\n       p value red of HT to = %s, p value black of HT to = %s.\n       Distance = %s \n**\n' %
                            #    (str(curr['s']), curr['t'], str(HT_to.label), str(curr['t']), str(total_curr),
                            #     str(total_HT_to), str(curr['t']), str(p_value_curr_red), str(p_value_curr_black), str(Pr_red), str(Pr_black),
                            #     str(p_value_HT_to_red),str(p_value_HT_to_black), str(S_dis_matrix[(curr['t'], HT_to.label)])))

                        elif (pattern == "same_color" or pattern == "only_red") and (p_value_HT_to_red < p): #HT to red
                            red_HT_vertices_in_G.append({'curr': curr, 'HT_to_in_G': HT_to_in_G,
                                                       'probability': curr['probability'], 'distance' : S_dis_matrix[(curr['t'], HT_to.label)]})
                            #print('     this is a good pattern (red to red): \n       from %s to %s (HT to %s), Total vertices in S below %s: %s, below HT to: %s\n       p value red of %s = %s, p value black = %s (Pr red = %s, Pr black = %s)\n       p value red of HT to = %s,  p value black of HT to = %s.\n       Distance = %s \n**\n'%
                             #   (str(curr['s']),curr['t'],str(HT_to.label),str(curr['t']),str(total_curr),str(total_HT_to),str(curr['t']),str(p_value_curr_red) ,str(p_value_curr_black) ,str(Pr_red),str(Pr_black),str(p_value_HT_to_red),str(p_value_HT_to_black),str(S_dis_matrix[(curr['t'], HT_to.label)])))
                    elif p_value_curr_black < p: #HT from black
                        if (pattern == "different_colors" or pattern == "only_black_to_red") and (p_value_HT_to_red < p):  # HT to red
                            black_HT_vertices_in_G.append({'curr': curr, 'HT_to_in_G': HT_to_in_G,
                                                           'probability': curr['probability'],
                                                           'distance': S_dis_matrix[(curr['t'], HT_to.label)]})
                            #print(
                            #    '     this is a good pattern (black to red): \n       from %s to %s (HT to %s), Total vertices in S below %s: %s, below HT to: %s\n       p value red of %s = %s, p value black = %s (Pr red = %s, Pr black = %s)\n       p value red of HT to = %s, p value black of HT to = %s.\n       Distance = %s \n**\n' %
                            #    (str(curr['s']), curr['t'], str(HT_to.label), str(curr['t']), str(total_curr),
                            #     str(total_HT_to), str(curr['t']), str(p_value_curr_red),str(p_value_curr_black), str(Pr_red), str(Pr_black),
                            #     str(p_value_HT_to_red),str(p_value_HT_to_black), str(S_dis_matrix[(curr['t'], HT_to.label)])))

                        elif (pattern == "same_color") and (p_value_HT_to_black < p): #HT to black
                            black_HT_vertices_in_G.append({'curr': curr, 'HT_to_in_G': HT_to_in_G,
                                                       'probability': curr['probability'], 'distance' : S_dis_matrix[(curr['t'], HT_to.label)]})
                            #print('     this is a good pattern (black to black): \n       from %s to %s (HT to %s), Total vertices in S below %s: %s, below HT to: %s\n       p value red of %s = %s, p value black = %s (Pr red = %s, Pr black = %s)\n       p value red of HT to = %s, p value black of HT to = %s.\n       Distance = %s \n**\n'%
                            #    (str(curr['s']),curr['t'],str(HT_to.label),str(curr['t']),str(total_curr),str(total_HT_to),str(curr['t']),str(p_value_curr_red),str(p_value_curr_black) ,str(Pr_red),str(Pr_black),str(p_value_HT_to_red),str(p_value_HT_to_black),str(S_dis_matrix[(curr['t'], HT_to.label)])))
                elif pattern == "any":               #mark HT (in the red list, just for convi)
                    red_HT_vertices_in_G.append({'curr': curr, 'HT_to_in_G': HT_to_in_G,
                                                 'probability': curr['probability'],
                                                 'distance': S_dis_matrix[(curr['t'], HT_to.label)]})
            elif evolutinary_event == 'D' and curr['event'] == 'D' and curr['probability'] > 0:
                x = S.find_node(lambda n: (n.label == curr['t']))
                incoming_edges_curr = [e for e in incoming_edges if e[2]['target'] == i]
                child0_in_S = S.find_node(
                    lambda n: (n.label == H.nodes(data=True)[incoming_edges_curr[0][0]]['t']))
                child1_in_S = S.find_node(
                    lambda n: (n.label == H.nodes(data=True)[incoming_edges_curr[1][0]]['t']))
                child0_in_H = H.nodes(data=True)[incoming_edges_curr[0][0]]
                child1_in_H = H.nodes(data=True)[incoming_edges_curr[1][0]]

                curr_red = S_colors[x.label][0]
                curr_blacks = S_colors[x.label][1]
                total_curr = curr_red + curr_blacks

                child0_reds = S_colors[child0_in_S.label][0]
                child0_blacks = S_colors[child0_in_S.label][1]
                total_child0 = child0_reds + child0_blacks

                child1_reds = S_colors[child1_in_S.label][0]
                child1_blacks = S_colors[child1_in_S.label][1]
                total_child1 = child1_reds + child1_blacks

                p_value_curr_red, nCr_lookup_table, fact_lookup_table = utiles.p_value_calculation(curr_red, total_curr,
                                                                                                   nCr_lookup_table,
                                                                                                   fact_lookup_table,
                                                                                                   accur, Pr_red,
                                                                                                   Pr_black, x.label)
                p_value_curr_black, nCr_lookup_table, fact_lookup_table = utiles.p_value_calculation(curr_blacks,
                                                                                                     total_curr,
                                                                                                     nCr_lookup_table,
                                                                                                     fact_lookup_table,
                                                                                                     accur, Pr_black,
                                                                                                     Pr_red,
                                                                                                     x.label)
                flag = False
                if (pattern == "only_red" or pattern == "any") and p_value_curr_red < p: #duplication in red
                    red_HT_vertices_in_G.append({'curr': curr,'probability': curr['probability']})
                    #print(
                    #            '     this is a good pattern (red duplication): \n       curr = %s, total vertices under x = %s,\n       p value red = %s, p value black = %s (Pr red = %s, Pr black = %s)\n      prob = %s \n**\n' %
                    #            (str(curr), str(total_curr), str(p_value_curr_red),str(p_value_curr_black), str(Pr_red), str(Pr_black),str(curr['probability'])))
                    flag = True
                if (pattern == "only_black" or pattern == "any") and p_value_curr_black < p and not flag: #duplication in red
                    black_HT_vertices_in_G.append({'curr': curr,'probability': curr['probability']})
                    #print(
                    #            '     this is a good pattern (black duplication): \n       curr = %s, total vertices under x = %s,\n       p value red = %s, p value black = %s (Pr red = %s, Pr black = %s)\n      prob = %s \n**\n' %
                    #            (str(curr), str(total_curr), str(p_value_curr_red),str(p_value_curr_black), str(Pr_red), str(Pr_black),str(curr['probability'])))


    #print('Finished finding same-color HT.\n')
    return red_HT_vertices_in_G, black_HT_vertices_in_G, nCr_lookup_table, fact_lookup_table

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
    noise_level = parameters[0]
    TH_both = parameters[1]
    check_diffreance_between_solutions = parameters[2]

    list_of_scores_for_rand_num = {}
    random_for_prec_curr = random_for_prec
    for rand_num in range (0, random_for_prec_curr):
        print('                 ****     Iteration %s.%s/%s.%s (Multithread)    ***** \n' % (str(noise_level),str(rand_num),str(noise_level_list[len(noise_level_list)-1]),str(random_for_prec)))
        if not check_diffreance_between_solutions:
            S_dis_matrix = {}
            nodes_table = {}
            S_colors = {}
            all_vertices = {}
            iter = -1
            new_G = nx.DiGraph()
            red_HT_vertices_in_G = []
            black_HT_vertices_in_G = []
            deleted_nodes = []
            nCr_lookup_table = {}
            fact_lookup_table = {}
            G_internal_colors = {}

            G = tr.Tree.get_from_path(path + "/GeneTree(binary)_local.txt", schema="newick")
            S = tr.Tree.get_from_path(path+"/phyliptree(binary,"+speciesTreespecification+").phy", schema="newick")

            input = open(path+'/'+ str(noise_level)+'/sigma'+str(noise_level)+'.'+str(rand_num)+'.txt', 'r')
            sigma = []
            for line in input:
                sigma.append(eval(line))
            sigma = sigma[0]

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

            S_labels_table, G_labels_table, sigma = inits.init_taxon_to_label_table(S,G,sigma)
            sigma, old_sigma = inits.update_sigma(S, G, k, sigma, test, path,exact_names,S_labels_table,G_labels_table)
            G.prune_taxa_with_labels(tree_operations.remove_unsigma_genes(G, sigma, False))
            colors,old_colors = inits.update_colors(S, colors,exact_names)
            TH_edges_in_subtree = 5                                                    # smallest subtree that will be counted when not comparing subtrees
            TH_pattern_in_subtree = (TH_edges_in_subtree * 0.005)/(2*k*number_of_planted_vertices)

            S_dis_matrix = inits.init_distance_S(S_dis_matrix, k, test, path,speciesTreespecification)
            nodes_table = inits.init_nodes_table(S, G, nodes_table)

            H, H_number_of_nodes, nodes_table = hypergraph.build_hyper_garph(S, G, test, k,
                                                                             nodes_table, D_cost, S_cost, HT_cost, path, alpha,
                                                                             sigma, save_data)

            H, max_prob = hypergraph.assign_probabilities(S, G, H, test, k, gamma, path, alpha)
        else:
            S = parameters[3]
            G = parameters[4]
            S_colors = parameters[5]
            colors = parameters[6]
            deleted_nodes = parameters[7]
            S_dis_matrix = parameters[8]
            nCr_lookup_table = parameters[9]
            fact_lookup_table = parameters[10]
            red_HT_vertices_in_G = parameters[11]
            black_HT_vertices_in_G = parameters[12]
            sigma = parameters[13]
            new_G = parameters[14]
            G_internal_colors = parameters[15]
            all_vertices = parameters[16]
            TH_pattern_in_subtree = parameters[17]
            TH_edges_in_subtree = parameters[18]
            iter = parameters[19]
            H = parameters[20]
        if H is None:
            list_of_scores_for_rand_num.update({rand_num: {}})
        else:
            ##      PROBABILITIES, COLORS, PATTERN      ##

            S_colors = tree_operations.color_tree(S, 'S', S_colors, colors, sigma)

            red_HT_vertices_in_G, black_HT_vertices_in_G, nCr_lookup_table, fact_lookup_table = find_Pattern(H,S,
                                                                                                                   S_dis_matrix,
                                                                                                                   nCr_lookup_table,
                                                                                                                   fact_lookup_table,
                                                                                                                   red_HT_vertices_in_G,
                                                                                                                   black_HT_vertices_in_G, pattern,evolutinary_event,S_colors)
            H = hypergraph.remove_prob_zero(H, deleted_nodes)

            max_S_d_of_HT = tree_operations.find_max_d_of_HT(S_dis_matrix, red_HT_vertices_in_G, black_HT_vertices_in_G,evolutinary_event)

            new_G = tree_operations.weight_G_based_on_same_color_HT(G, new_G, red_HT_vertices_in_G,
                                                                          black_HT_vertices_in_G, max_S_d_of_HT,dis_flag,evolutinary_event,check_diffreance_between_solutions,k)
            new_G = tree_operations.number_of_edges_in_subtree(new_G)

            G_internal_colors = tree_operations.color_tree(G, 'G', G_internal_colors, colors, sigma)
            marked_nodes,all_vertices = pattern_identify.find_signi_distance(new_G, all_vertices, TH_compare_subtrees, TH_both, TH_pattern_in_subtree, path, k, alpha, both,
                                                                G_internal_colors, iter,speciesTreespecification,compare_subtrees,TH_edges_in_subtree,check_diffreance_between_solutions)

            list_of_scores_for_rand_num.update({rand_num:all_vertices})
    #print('         List for noise_level %s: %s' % (str(noise_level),str(list_of_scores_for_rand_num)))
    return (utiles.average_of_list(list_of_scores_for_rand_num,random_for_prec_curr))

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
    TH_both = parameters[10]
    H = parameters[11]
    S = parameters[12]
    G = parameters[13]
    TH_edges_in_subtree = parameters[14]

    print('\n                                               *** ' + str(iter + 1) + 'th iteration ***\n')

    TH_pattern_in_subtree = 2 / (max_dis * k)
    new_G[iter] = nx.DiGraph()
    deleted_nodes = []
    G_internal_colors = {}

    solutions[iter] = nx.DiGraph()
    H_root = [nd for nd in list(H.node(data=True)) if
              nd[1]['s'] == G.seed_node.label and nd[1]['t'] == S.seed_node.label]
    solutions[iter], nodes_table = hypergraph.track_a_solution(H_root, H, S, G, solutions[iter], iter * factor)
    solutions[iter], max_prob = hypergraph.assign_probabilities(S, G, solutions[iter], test, k, gamma, path, alpha)

    red_HT_vertices_in_G, black_HT_vertices_in_G, nCr_lookup_table, fact_lookup_table = find_Pattern(
        solutions[iter], S,S_dis_matrix, nCr_lookup_table, fact_lookup_table, red_HT_vertices_in_G,
        black_HT_vertices_in_G, pattern, evolutinary_event, S_colors)

    solutions[iter] = hypergraph.remove_prob_zero(solutions[iter], deleted_nodes)
    max_S_d_of_HT = tree_operations.find_max_d_of_HT(S_dis_matrix, red_HT_vertices_in_G, black_HT_vertices_in_G,
                                                     evolutinary_event)

    new_G[iter] = tree_operations.weight_G_based_on_same_color_HT(G, new_G[iter], red_HT_vertices_in_G,
                                                                  black_HT_vertices_in_G, max_S_d_of_HT, dis_flag,
                                                                  evolutinary_event, check_diffreance_between_solutions,
                                                                  k)
    new_G[iter] = tree_operations.number_of_edges_in_subtree(new_G[iter])

    G_internal_colors = tree_operations.color_tree(G, 'G', G_internal_colors, colors, sigma)
    all_vertices = {}
    marked_nodes, all_vertices = pattern_identify.find_signi_distance(new_G[iter], all_vertices, TH_compare_subtrees,
                                                                      TH_both, TH_pattern_in_subtree, path, k, alpha,
                                                                      both, G_internal_colors, iter,
                                                                      speciesTreespecification, compare_subtrees,
                                                                      TH_edges_in_subtree,
                                                                      check_diffreance_between_solutions)

    return([{iter * factor: all_vertices},list(marked_nodes.items()),new_G[iter]])

##********  MAIN ***********

def main():
    global S, G, iterations, sigma, alpha, gamma, colors, TH_compare_subtrees, TH_both,planted_vertices, TH_edges_in_subtree, number_of_planted_vertices,TH_pattern_in_subtree, both, path, speciesTreespecification, evolutinary_event,exact_names,noise_level_list,random_for_prec
    starting_time = datetime.now()

    all_vertices_with_index = {}
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

    input = open(path+'/'+ str(noise_level)+'/sigma'+str(noise_level)+'.'+str(rand_num)+'.txt', 'r')
    sigma = []
    for line in input:
        sigma.append(eval(line))
    sigma = sigma[0]

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
    TH_edges_in_subtree = 5                                                    # smallest subtree that will be counted when not comparing subtrees
    TH_pattern_in_subtree = (TH_edges_in_subtree * 0.005)/(2*k*number_of_planted_vertices)

    S_dis_matrix = inits.init_distance_S(S_dis_matrix, k, test, path,speciesTreespecification)
    nodes_table = inits.init_nodes_table(S, G, nodes_table)
    #draw.draw_S_and_G(S, G, old_sigma, colors, sigma, path, {}, 'all')

    H, H_number_of_nodes, nodes_table = hypergraph.build_hyper_garph(S, G, test, k,
                                                                     nodes_table, D_cost, S_cost, HT_cost, path, alpha,
                                                                     sigma,save_data)
    if check_diffreance_between_solutions:
        max_dis = tree_operations.max_dis(S_dis_matrix)
        if iterations * factor < k:
            all_marked_for_TH = {}
            all_unmarked_for_TH = {}
            H, max_prob = hypergraph.assign_probabilities(S, G, H, test, k, gamma, path, alpha)

            if H == None:
                quit()
            list_of_TH = utiles.frange(0,100,1)        #TH_edges_in_subtree
            parameters = []
            p = Pool(15)
            for i in range(0, len(list_of_TH)):
                parameters.append([noise_level_list[0], TH_both, check_diffreance_between_solutions, S, G, S_colors, colors, [], S_dis_matrix, {}, {},
                                   [], [], sigma, nx.DiGraph(), {}, all_vertices, TH_pattern_in_subtree, list_of_TH[i],i, H])
            list_of_RSAM_results = p.map(RSAM_finder_multithread, parameters)
            ind = 0
            for res in list_of_RSAM_results:
                all_vertices_with_index.update({list_of_TH[ind]: [res]})
                ind += 1
            file = open(path + '/saved_data/all_vertices_RSAM_finder.txt', 'w')
            file.write(str(all_vertices_with_index))
            file.close()
            for TH in list_of_TH:
                all_marked = []
                new_G = {}
                solutions = {}
                red_HT_vertices_in_G = []
                black_HT_vertices_in_G = []
                nCr_lookup_table = {}
                fact_lookup_table = {}
                S_colors = tree_operations.color_tree(S, 'S', S_colors, colors, sigma)
                p1 = Pool(15)
                parameters_list = [(x,new_G,max_dis,solutions,S_dis_matrix,nCr_lookup_table,fact_lookup_table,red_HT_vertices_in_G,
                                    black_HT_vertices_in_G,S_colors,TH_both,H,S,G,TH) for x in range(0,iterations)]
                list_of_results = p1.map(extract_and_tarce_a_solution, parameters_list)
                new_G_to_save = []
                for res in list_of_results:
                    all_vertices_with_index.update(res[0])
                    all_marked.append(res[1])
                    new_G_to_save.append(res[2])
                list_of_unmarked_TH = []
                for li in all_marked:
                    list_of_unmarked = []
                    for u in G.postorder_node_iter():
                        flag = False
                        for tup in li:
                            if u.label == tup[0]:
                                flag = True
                        if not flag:
                            list_of_unmarked.append(u.label)
                    list_of_unmarked_TH.append(list(list_of_unmarked))
                all_marked_for_TH.update({TH_both:(all_marked)})
                all_unmarked_for_TH.update({TH_both:(list_of_unmarked_TH)})
                file = open(path + '/saved_data/all_marked_nodes_for_TH.txt', 'w')
                file.write(str(all_marked_for_TH))
                file.close()
                file = open(path + '/saved_data/all_unmarked_nodes_for_TH.txt', 'w')
                file.write(str(all_unmarked_for_TH))
                if (not on_lab) and (TH_both == 0):
                    draw.draw_new_G2({}, colors, sigma, new_G_to_save[0], G, old_sigma, k, TH_compare_subtrees, TH_both,
                                     TH_pattern_in_subtree, path, both, alpha, True, glob, speciesTreespecification,
                                     pattern,
                                     big_size, evolutinary_event, compare_subtrees, 1)
            print('all_marked_for_TH: %s' % str(all_marked_for_TH))
            print('all_unmarked_for_TH: %s' % str(all_unmarked_for_TH))

            file = open(path + '/saved_data/all_marked_nodes_for_TH.txt', 'w')
            file.write(str(all_marked_for_TH))
            file.close()
            file = open(path + '/saved_data/all_unmarked_nodes_for_TH.txt', 'w')
            file.write(str(all_unmarked_for_TH))

            file.close()
            quit()
        else:
            print('**   not enough solutions were calculated in order to track solution number %s   **' % str(iterations * factor))
        quit()

    else:
        iter = -1
        new_G = nx.DiGraph()
        red_HT_vertices_in_G = []
        black_HT_vertices_in_G = []
        deleted_nodes = []
        nCr_lookup_table = {}
        fact_lookup_table = {}
        G_internal_colors = {}
        H, max_prob = hypergraph.assign_probabilities(S, G, H, test, k, gamma, path, alpha)
        if H is None:
            list_of_scores_for_rand_num.update({rand_num: {}})
        else:
            ##      PROBABILITIES, COLORS, PATTERN      ##

            S_colors = tree_operations.color_tree(S, 'S', S_colors, colors, sigma)

            red_HT_vertices_in_G, black_HT_vertices_in_G, nCr_lookup_table, fact_lookup_table = find_Pattern(H,S,
                                                                                                                   S_dis_matrix,
                                                                                                                   nCr_lookup_table,
                                                                                                                   fact_lookup_table,
                                                                                                                   red_HT_vertices_in_G,
                                                                                                                   black_HT_vertices_in_G, pattern,evolutinary_event,S_colors)
            H = hypergraph.remove_prob_zero(H, deleted_nodes)
            max_S_d_of_HT = tree_operations.find_max_d_of_HT(S_dis_matrix, red_HT_vertices_in_G, black_HT_vertices_in_G,evolutinary_event)

            new_G = tree_operations.weight_G_based_on_same_color_HT(G, new_G, red_HT_vertices_in_G,
                                                                          black_HT_vertices_in_G, max_S_d_of_HT,dis_flag,evolutinary_event,check_diffreance_between_solutions,k)
            new_G = tree_operations.number_of_edges_in_subtree(new_G)

            G_internal_colors = tree_operations.color_tree(G, 'G', G_internal_colors, colors, sigma)
            marked_nodes,all_vertices = pattern_identify.find_signi_distance(new_G, all_vertices, TH_compare_subtrees, TH_both, TH_pattern_in_subtree, path, k, alpha, both,
                                                                G_internal_colors, iter,speciesTreespecification,compare_subtrees,TH_edges_in_subtree,check_diffreance_between_solutions)

            list_of_scores_for_rand_num.update({rand_num:all_vertices})
            if not on_lab:
                draw.draw_new_G2(marked_nodes, colors, sigma, new_G, G, old_sigma, k, TH_compare_subtrees, TH_both,
                             TH_pattern_in_subtree, path, both, alpha, True, glob, speciesTreespecification, pattern,
                             big_size, evolutinary_event, compare_subtrees, 1)
    #print('         List for noise_level %s: %s' % (str(noise_level),str(list_of_scores_for_rand_num)))
    all_vertices_with_index.update({noise_level:utiles.average_of_list(list_of_scores_for_rand_num,random_for_prec_curr)})

    parameters = []
    p = Pool(15)
    for i in range(0,len(noise_level_list)):
        parameters.append([noise_level_list[i],TH_both,check_diffreance_between_solutions])
    list_of_results = p.map(RSAM_finder_multithread, parameters)
    ind = 0
    for res in list_of_results:
        all_vertices_with_index.update({noise_level_list[ind]:res})
        ind += 1
    print('all_vertices_with_index: %s' % str(all_vertices_with_index))
    file = open(path + '/saved_data/all_vertices_RSAM_finder.txt', 'w')
    file.write(str(all_vertices_with_index))
    file.close()

    #if not on_lab:
        #draw.draw_plot(all_vertices_with_index,path,marked_vertex)
    print('Running time: '+str(datetime.now()-starting_time))

if __name__ == "__main__":
    main()


