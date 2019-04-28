on_lab = False
check_diffreance_between_solutions = False
real_data = True

name = 'COG2602'

if on_lab:
        path = '/storage/DATA/users/ronizo/COGS/'
else:
    import sys
    sys.path.append('/PycharmProjects/RSAM_venv/lib/python3.6/site-packages/graphviz/')
    path = '/Users/ronizoller/Google Drive (ronizo@post.bgu.ac.il)/COGS/'

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
import draw
from utils import extract_from_FASTA_v1 as extr

speciesTreespecification = 'bacteria'
geneExt = 'bacteria'
test = False                                         # if True all data will be loaded from outter files, otherwise all data will be calculated and saved
glob = False                                        # if True global alignment is used, otherwise local
compare_subtrees = False                             # if true the algorithm will look for a signi different between two children of u in G, otherwise it will look for u in G s.t. in G(u) there are alot of same color HT
k = 50
exact_names = True

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
noise_level_list = [20]
p = 0.05                                            #p_value

#compare several optimal solutions
if check_diffreance_between_solutions:
    iterations = 1                                  #will check for each i=0 to iteration the solution i*factor
    factor = 1
    real_data = False

big_size = 2000
small_size = 7

def find_Pattern(H, S,S_dis_matrix, nCr_lookup_table, fact_lookup_table, pattern,S_colors):
    if pattern[0] is not None:
        total_red = S_colors[S.seed_node.label][0]
        total_black = S_colors[S.seed_node.label][1]
        num_of_leafs = tree_operations.number_of_leafs(S,'S')
        Pr_red = total_red / num_of_leafs
        Pr_black = total_black / num_of_leafs
        interesting_vertices = []
        for nd in (list(nx.topological_sort(H))):
            incoming_edges = H.in_edges([nd], data=True)
            nd = H.nodes(data=True)[nd]
            for i in range(0,len(nd['l'])):
                curr = nd['l'][i]

                x = S.find_node(lambda n: (n.label == curr['t']))
                curr_red = S_colors[x.label][0]
                curr_blacks = S_colors[x.label][1]
                total_curr = curr_red + curr_blacks
                if pattern[1] == 'red':
                    p_value_curr_red, nCr_lookup_table, fact_lookup_table = utiles.p_value_calculation(curr_red,
                                                                                                       total_curr,
                                                                                                       nCr_lookup_table,
                                                                                                       fact_lookup_table,
                                                                                                       accur, Pr_red,
                                                                                                       Pr_black,
                                                                                                       x.label)
                elif pattern[1] == 'black':
                    p_value_curr_black, nCr_lookup_table, fact_lookup_table = utiles.p_value_calculation(curr_blacks,
                                                                                                         total_curr,
                                                                                                         nCr_lookup_table,
                                                                                                         fact_lookup_table,
                                                                                                         accur,
                                                                                                         Pr_black,
                                                                                                         Pr_red,
                                                                                                         x.label)

                if 'HT' in pattern[0] and curr['event'] == 'HT' and curr['probability'] > 0:
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

                    HT_to_reds = S_colors[HT_to.label][0]
                    HT_to_blacks =  S_colors[HT_to.label][1]
                    total_HT_to = HT_to_reds + HT_to_blacks

                    if pattern[1] == 'red':
                        if total_HT_to > 0:
                            p_value_HT_to_red, nCr_lookup_table, fact_lookup_table = utiles.p_value_calculation(
                                math.ceil(HT_to_reds), total_HT_to, nCr_lookup_table, fact_lookup_table, accur, Pr_red,
                                Pr_black,x.label)
                            if p_value_curr_red < p and p_value_HT_to_red < p:
                                    interesting_vertices.append({'curr': curr, 'HT_to_in_G': HT_to_in_G,
                                                               'probability': curr['probability'], 'distance' : S_dis_matrix[(curr['t'], HT_to.label)]})

                    if pattern[1] == 'black':
                        if total_HT_to > 0:
                            p_value_HT_to_black, nCr_lookup_table, fact_lookup_table = utiles.p_value_calculation(
                                math.ceil(HT_to_blacks), total_HT_to, nCr_lookup_table, fact_lookup_table, accur, Pr_black,
                                Pr_red, x.label)
                            if p_value_curr_black < p and p_value_HT_to_black < p:
                                    interesting_vertices.append({'curr': curr, 'HT_to_in_G': HT_to_in_G,
                                                               'probability': curr['probability'], 'distance' : S_dis_matrix[(curr['t'], HT_to.label)]})
                    elif pattern[1] == None:
                        interesting_vertices.append({'curr': curr, 'HT_to_in_G': HT_to_in_G,
                                                     'probability': curr['probability'],
                                                     'distance': S_dis_matrix[(curr['t'], HT_to.label)]})
                if ('D' in pattern[0] and curr['event'] == 'D') or ('S' in pattern[0] and curr['event'] == 'S') and curr['probability'] > 0:
                    if pattern[1] == 'red' and p_value_curr_red < p:
                        interesting_vertices.append({'curr': curr, 'probability': curr['probability']})
                    elif pattern[1] == 'black' and p_value_curr_black < p:
                        interesting_vertices.append({'curr': curr, 'probability': curr['probability']})
                    elif pattern[1] == None:
                        interesting_vertices.append({'curr': curr, 'probability': curr['probability']})

        return interesting_vertices, nCr_lookup_table, fact_lookup_table
    return (None,None,None)

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

    check_diffreance_between_solutions = parameters[0]

    list_of_scores_for_rand_num = {}
    max_score_p1_list = []
    max_score_p1_and_p2_list = []
    random_for_prec_curr = random_for_prec
    for rand_num in range (0, random_for_prec_curr):
        if not check_diffreance_between_solutions:
            check_diffreance_between_solutions,noise_level,noise_in,p1,p2 = parameters
            path_change_in = path + '/' + noise_in
            os.makedirs(os.path.dirname(path_change_in), exist_ok=True)
            path_change_in = path_change_in + '/saved_data'
            os.makedirs(os.path.dirname(path_change_in), exist_ok=True)
            path_change_in = path +  '/' + noise_in

            S_dis_matrix = {}
            nodes_table = {}
            S_colors = {}
            all_vertices = {}
            new_G = nx.DiGraph()

            nCr_lookup_table = {}
            fact_lookup_table = {}

            G = tr.Tree.get_from_path(path  + "/GeneTree(binary,"+geneExt+")_local.txt", schema="newick")
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
            check_diffreance_between_solutions, noise_level, S, G, S_colors, colors, S_dis_matrix, nCr_lookup_table, fact_lookup_table, sigma, new_G, all_vertices, H,p1,p2 = parameters

        if H is None:
            list_of_scores_for_rand_num.update({rand_num: {}})
        else:
            ##      PROBABILITIES, COLORS, PATTERN      ##

            S_colors = tree_operations.color_tree(S, 'S', S_colors, colors, sigma)

            interesting_vertices_p1, nCr_lookup_table, fact_lookup_table = find_Pattern(H,S,S_dis_matrix,nCr_lookup_table,fact_lookup_table, p1,S_colors)
            interesting_vertices_p2, nCr_lookup_table, fact_lookup_table = find_Pattern(H,S,S_dis_matrix,nCr_lookup_table,fact_lookup_table, p2,S_colors)

            max_S_d_of_HT = tree_operations.find_max_d_of_HT(S_dis_matrix, interesting_vertices_p1,p1)

            new_G = tree_operations.weight_G_based_on_same_color_HT(G, new_G, interesting_vertices_p1,interesting_vertices_p2, max_S_d_of_HT,p1,p2,check_diffreance_between_solutions)

            new_G = tree_operations.number_of_edges_in_subtree(new_G)
            new_G = tree_operations.normlize_weights(new_G,k,p1,'p1')
            new_G = tree_operations.normlize_weights(new_G,k,p2,'p2')

            if not check_diffreance_between_solutions:
                if p2[0] is None:
                    max_score_p1_list = tree_operations.find_max_scores(new_G,number_of_planted_vertices,'p1',p1[3])
                else:
                    max_score_p1_and_p2_list = tree_operations.find_max_scores(new_G,number_of_planted_vertices,'p2',p1[3])

            marked_nodes,all_vertices = pattern_identify.find_signi_distance(new_G, all_vertices,p1, p2,max_score_p1_list, max_score_p1_and_p2_list, check_diffreance_between_solutions)

            list_of_scores_for_rand_num.update({rand_num:all_vertices})
    return (utiles.average_of_list(list_of_scores_for_rand_num,random_for_prec_curr),noise_in)

def extract_and_tarce_a_solution(parameters):
    iter, new_G, solutions, S_dis_matrix, nCr_lookup_table, fact_lookup_table, S_colors, H, S, G, TH_edges_in_subtree, k,p1,p2 = parameters

    new_G[iter] = nx.DiGraph()
    max_score_p1_list = []
    max_score_p1_and_p2_list = []

    solutions[iter] = nx.DiGraph()
    H_root = [nd for nd in list(H.node(data=True)) if
              nd[1]['s'] == G.seed_node.label and nd[1]['t'] == S.seed_node.label]
    solutions[iter], nodes_table = hypergraph.track_a_solution(H_root, H, S, G, solutions[iter], random.choice(range(0,k)))

    solutions[iter], max_prob = hypergraph.assign_probabilities(S, G, solutions[iter], gamma)

    interesting_vertices_p1, nCr_lookup_table, fact_lookup_table = find_Pattern(solutions[iter], S, S_dis_matrix, nCr_lookup_table,
                                                                                fact_lookup_table, p1, S_colors)
    interesting_vertices_p2, nCr_lookup_table, fact_lookup_table = find_Pattern(solutions[iter], S, S_dis_matrix, nCr_lookup_table,
                                                                                fact_lookup_table, p2, S_colors)

    max_S_d_of_HT = tree_operations.find_max_d_of_HT(S_dis_matrix, interesting_vertices_p1,p1)

    new_G[iter] = tree_operations.weight_G_based_on_same_color_HT(G, new_G[iter], interesting_vertices_p1,interesting_vertices_p2, max_S_d_of_HT,p1,p2,check_diffreance_between_solutions)
    new_G[iter] = tree_operations.number_of_edges_in_subtree(new_G[iter])
    new_G[iter] = tree_operations.normlize_weights(new_G[iter],1,p1,'p1')
    new_G[iter] = tree_operations.normlize_weights(new_G[iter],1,p2,'p2')
    if not check_diffreance_between_solutions:
        if p2[0] is None:
            max_score_p1_list = tree_operations.find_max_scores(new_G, number_of_planted_vertices, 'p1',p1[3])
        else:
            max_score_p1_and_p2_list = tree_operations.find_max_scores(new_G, number_of_planted_vertices, 'p2', p1[3])
    all_vertices = {}
    marked_nodes, all_vertices = pattern_identify.find_signi_distance(new_G[iter], all_vertices,p1, p2,max_score_p1_list, max_score_p1_and_p2_list, check_diffreance_between_solutions)

    return([list(all_vertices.items()),list(marked_nodes.items()),new_G[iter],'(%s)' % (str(TH_edges_in_subtree))])

def end_function(H,S,G,k,starting_time,p1,p2,marked_nodes,old_sigma,max_list_p1,max_list_p1_and_p2):
    if p2[0] is not None:
        pattern_name = '('+str(p1[0:3])+'_'+str(p2[0:3])+')_Double-Mode'
    else:
        pattern_name = str(p1[0:3]) + '_Single-Mode'
    file = open(path + '/saved_data/marked_RSAM_finder_'+speciesTreespecification+'_pattern='+pattern_name+'.txt', 'w')
    file.write(str(marked_nodes))
    file.close()

    lists = ''
    if p2[0] is not None:
        max_list = max_list_p1_and_p2
    else:
        max_list = max_list_p1
    list_of_returns = [0]*len(max_list)
    index = 0
    for nd,x in marked_nodes.items():
        if index < len(max_list):
            r,l = tree_operations.leaf_in_subtrees(G,'S',nd, old_sigma,False)
            lists += 'For %s:\nlist = %s\n' % (str(nd),str(r+l))+'\n\n'
            if p2[0] is None:
                itm = x[0]
            else:
                itm = x[0]+x[1]
            ind, list_of_returns = utiles.index_with_repeting(max_list, itm, list_of_returns)
            if p2[0] is None:
                extr.main([(r+l,'_'+speciesTreespecification)],path,speciesTreespecification,str(ind)+'th_solution',pattern_name,True)
            else:
                extr.main([(r, '_' + speciesTreespecification)], path, speciesTreespecification,
                          str(ind) + 'th_solution_right'+x[2], pattern_name, True)
                extr.main([(l, '_' + speciesTreespecification)], path, speciesTreespecification,
                          str(ind) + 'th_solution_left'+ x[2], pattern_name, True)
            index += 1

        file = open(path + '/saved_data/marked_nodes_leafs_lists_' + speciesTreespecification +'_pattern='+pattern_name+'.txt', 'w')
        file.write(str(lists))
        file.close()
    print('max_list: ' + str(max_list))

    print('COG: %s Class: %s Pattern: %s\nNumber of co-optimal out of %s solutions: %s with cost %s' % (str(name),str(speciesTreespecification),pattern_name
                                                                                                        ,str(k),str(hypergraph.find_number_of_cooptimal(H,G,S)[0]),
                                                                                                        str((hypergraph.find_number_of_cooptimal(H,G,S)[1]))))
    print('Running time: %s\nk: %s\nTH_edges: %s' % (
        str(datetime.now() - starting_time), str(k), str(p1[3])))
    quit()

##********  MAIN ***********

def main():
    global S, G, iterations, sigma, alpha, gamma, colors,planted_vertices, number_of_planted_vertices, path, speciesTreespecification,exact_names,noise_level_list,random_for_prec,k
    path = path + name +'/'
    starting_time = datetime.now()


    ### EV\subseteq {S,D,HT}, color\in {red,black,None}, distance\in {True,False}
    ### only p1 can have HT  in EV, TH_edges sould be the same
    p1 = (['HT'], 'red', True)
    p2 = (['S','D','HT'], 'black', None)

    if p2[0] is not None:
        pattern_name = '('+str(p1)+'_'+str(p2)+')_Double-Mode'
    else:
        pattern_name = str(p1) + '_Single-Mode'

    all_vertices_with_index = {}
    all_RSAM_unmarked = {}
    list_of_scores_for_rand_num = {}
    noise_level = 0
    rand_num = 0
    random_for_prec_curr = 1
    if not real_data:
        print('                 ****     Iteration %s.%s/%s.%s     ***** \n' % (str(noise_level),str(rand_num),str(noise_level_list[len(noise_level_list)-1]),str(random_for_prec)))
    else:
        print('                 ****     %s pattern %s    ***** \n' % (str(name), pattern_name))
    S_dis_matrix = {}
    nodes_table = {}
    S_colors = {}
    all_vertices = {}

    G = tr.Tree.get_from_path(path + "/GeneTree(binary,"+geneExt+")_local.txt", schema="newick")
    S = tr.Tree.get_from_path(path+"/phyliptree(binary,"+speciesTreespecification+").phy", schema="newick")
    if check_diffreance_between_solutions:
        noise_level = noise_level_list[0]
        input = open(path + '/colors_and_HT/' + str(noise_level) + '/sigma' + str(noise_level) + '.' + str(rand_num) + '.txt', 'r')
    else:
        input = open(path+'/'+ str(noise_level)+'/sigma'+str(noise_level)+'.'+str(rand_num)+'.txt', 'r')
    sigma = []
    for line in input:
        sigma.append(eval(line))
    sigma = sigma[0]

    if check_diffreance_between_solutions:
        noise_level = noise_level_list[0]
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
    if tree_operations.remove_unsigma_genes(G, sigma, False) is not []:
        G.prune_taxa_with_labels(tree_operations.remove_unsigma_genes(G, sigma, False))
    colors,old_colors = inits.update_colors(S, colors,exact_names)
    if 'D' in p1[0]:
        p1 = (p1[0],p1[1],p1[2],5)
    else:
        p1 = (p1[0], p1[1], p1[2],len(tree_operations.leaf_in_subtrees(G,'S',G.seed_node.label, old_sigma,False)[0]+tree_operations.leaf_in_subtrees(G,'S',G.seed_node.label, old_sigma,False)[1])*0.1)

    #if not on_lab:
    #    draw.draw_S_and_G(S, G, old_sigma, colors, sigma, path, None, speciesTreespecification,False)
    #tree_operations.reroot_and_save(S,'x227',path,speciesTreespecification)
    S_dis_matrix = inits.init_distance_S(S_dis_matrix, k, test, path,speciesTreespecification)
    nodes_table = inits.init_nodes_table(S, G, nodes_table)
    H, H_number_of_nodes, nodes_table = hypergraph.build_hyper_garph(S, G, test, k,
                                                                     nodes_table, D_cost, S_cost, HT_cost, path, alpha,
                                                                     sigma,save_data)
    if check_diffreance_between_solutions:
        if iterations * factor < k:
            all_marked_for_TH = {}
            all_unmarked_for_TH = {}
            H, max_prob = hypergraph.assign_probabilities(S, G, H, gamma)
            if H == None:
                quit()
            parameters = []
            p = Pool(15)
            combined = [f*10 for f in utiles.frange(7,9,0.1)+[(40,0,400)]]
            for i in range(0, len(combined)):
                TH_edges_in_subtree = combined[i]
                p1 = (p1[0], p1[1], p1[2], TH_edges_in_subtree)
                parameters.append([check_diffreance_between_solutions,noise_level_list[0] , S, G, S_colors, colors, S_dis_matrix, {}, {}, sigma, nx.DiGraph(), all_vertices, H, p1, p2])
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
            nCr_lookup_table = {}
            fact_lookup_table = {}
            all_vertices_with_index = {}

            S_colors = tree_operations.color_tree(S, 'S', S_colors, colors, sigma)
            p1 = Pool(15)
            parameters_list = [(x,new_G,solutions,S_dis_matrix,nCr_lookup_table,fact_lookup_table,S_colors,H,S,G,TH_edges_in_subtree, k,p1,p2)
                               for x in range(0,iterations) for (TH_both,TH_edges_in_subtree) in combined]

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
            end_function(H,S,G,k,starting_time,p1,p2,{},old_sigma,[],[])
        else:
            print('**   not enough solutions were calculated in order to track solution number %s   **' % str(iterations * factor))
        quit()

    else:
        new_G = nx.DiGraph()
        nCr_lookup_table = {}
        fact_lookup_table = {}
        max_score_p1_list = []
        max_score_p1_and_p2_list = []
        H, max_prob = hypergraph.assign_probabilities(S, G, H, gamma,path)
        if H is None:
            list_of_scores_for_rand_num.update({rand_num: {}})
        else:
            S_colors = tree_operations.color_tree(S, 'S', S_colors, colors, sigma)

            interesting_vertices_p1, nCr_lookup_table, fact_lookup_table = find_Pattern(H,S,S_dis_matrix,nCr_lookup_table,fact_lookup_table, p1,S_colors)
            interesting_vertices_p2, nCr_lookup_table, fact_lookup_table = find_Pattern(H,S,S_dis_matrix,nCr_lookup_table,fact_lookup_table, p2,S_colors)

            max_S_d_of_HT = tree_operations.find_max_d_of_HT(S_dis_matrix, interesting_vertices_p1,p1)

            new_G = tree_operations.weight_G_based_on_same_color_HT(G, new_G,interesting_vertices_p1,interesting_vertices_p2, max_S_d_of_HT,p1,p2,check_diffreance_between_solutions)
            new_G = tree_operations.number_of_edges_in_subtree(new_G)
            new_G = tree_operations.normlize_weights(new_G,k,p1,'p1')
            new_G = tree_operations.normlize_weights(new_G,k,p2,'p2')

            if not check_diffreance_between_solutions:
                if p2[0] is None:
                    max_score_p1_list = tree_operations.find_max_scores(new_G,number_of_planted_vertices,'p1',p1[3])
                else:
                    max_score_p1_and_p2_list = tree_operations.find_max_scores(new_G,number_of_planted_vertices,'p2',p1[3])
            marked_nodes,all_vertices = pattern_identify.find_signi_distance(new_G, all_vertices,p1, p2, max_score_p1_list, max_score_p1_and_p2_list, check_diffreance_between_solutions)

            list_of_scores_for_rand_num.update({rand_num:all_vertices})
    all_vertices_with_index.update({noise_level:utiles.average_of_list(list_of_scores_for_rand_num,random_for_prec_curr)})

    parameters = []
    if real_data:
        end_function(H,S,G,k,starting_time,p1,p2,marked_nodes,old_sigma,max_score_p1_list,max_score_p1_and_p2_list)
    else:
        for noise_in in ['colors_and_HT','color','HT']:
            for i in range(0,len(noise_level_list)):
                parameters.append([check_diffreance_between_solutions,noise_level_list[i],noise_in, p1, p2])
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
        end_function(H,S,G,k,starting_time,p1,p2,marked_nodes,old_sigma,max_score_p1_list,max_score_p1_and_p2_list)

if __name__ == "__main__":
    main()


