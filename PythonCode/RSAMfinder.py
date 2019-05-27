                ## pattern
### EV\subseteq {S,D,HT}, color\in {red,black,None}, distance\in {True,False}
### only p1 can have HT  in EV

p1 = ([], None, False)
p2 = (None, None, False)

                ### parameters
k = 50
TH_edges = 0.05
HT_cost = 1
D_cost = 1
S_cost = 0
gamma = 1                                           # factor for probability assignment
accur = 5                                           # calculations acuuracy
p = 0.05                                            #p_value
number_of_planted_vertices = 10                     # number of vertices that will be reported
save_data = False                                   #save hypergraph data
speciesTreespecification = 'bacteria'

import networkx as nx
import dendropy as tr
import math
import utiles
from utils import newick2edgelist as n2e
from utils import fix_taxa_names_in_FASTA as fix
import tree_operations_v1 as tree_operations
import inits_v1 as inits
import EfficiantVersion as hypergraph
import pattern_identify_v4 as pattern_identify
from multiprocessing import Pool
from datetime import datetime
import random
import os
import sys
from utils import extract_from_FASTA_v1 as extr
from utils import FASTA_to_taxamony_names as taxa_names
from utils import FASTA_to_gene_specie_mapping as create_sigma
from decimal import *

def find_Pattern(H, S,S_dis_matrix, nCr_lookup_table, fact_lookup_table, pattern,S_colors,p):
    print('Finding pattern '+str(pattern))
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
        print('Finished finding pattern\n')
        return interesting_vertices, nCr_lookup_table, fact_lookup_table
    return (None,None,None)

def end_function(H,S,G,k,starting_time,p1,p2,marked_nodes,old_sigma,max_list_p1,max_list_p1_and_p2):
    if p2[0] is not None:
        pattern_name = '('+str(p1[0:3])+'_'+str(p2[0:3])+')_Double-Mode'
    else:
        pattern_name = str(p1[0:3]) + '_Single-Mode'
    file = open(os.getcwd()+'/data/saved_data/marked_RSAM_finder_'+speciesTreespecification+'_pattern='+pattern_name+'.txt', 'w')
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
                extr.main([(r+l,'_'+speciesTreespecification)],os.getcwd()+'/data/',speciesTreespecification,str(ind)+'th_solution_node_'+str(nd),pattern_name,True)
            else:
                extr.main([(r, '_' + speciesTreespecification)], os.getcwd()+'/data/', speciesTreespecification,
                          str(ind) + '_node_' + str(nd) + '_'+x[2]+'_right', pattern_name, True)
                extr.main([(l, '_' + speciesTreespecification)],os.getcwd()+'data', speciesTreespecification,
                          str(ind) + '_node_'  +str(nd) + '_'+x[2]+'_left', pattern_name, True)
            index += 1

        file = open(os.getcwd()+'/data/saved_data/marked_nodes_leafs_lists_' + speciesTreespecification +'_pattern='+pattern_name+'.txt', 'w')
        file.write(str(lists))
        file.close()
    print('marked vertices: ' + str(marked_nodes))

    print('Class: %s Pattern: %s\nNumber of co-optimal out of %s solutions: %s with cost %s' % (str(speciesTreespecification),pattern_name
                                                                                                        ,str(k),str(hypergraph.find_number_of_cooptimal(H,G,S)[0]),
                                                                                                        str((hypergraph.find_number_of_cooptimal(H,G,S)[1]))))
    print('Running time: %s\nk: %s\nTH_edges: %s' % (
        str(datetime.now() - starting_time), str(k), str(p1[3])))
    quit()

def valid_pattern(ev,col,dist):
    if type(ev) == list and len(ev) <= 3:
        if col in ['red','black','None'] :
            if dist in [True,False]:
                return True
            else:
                print('distance is not valid')
        else:
            print('color is not valid')
    else:
        print('EV is not valid')

                ##********  MAIN ***********

def main (speciesTreespecification,k,TH_edges,HT_cost,D_cost,S_cost,gamma, p,number_of_planted_vertices,  p1, p2,GUI):
    starting_time = datetime.now()
    if not GUI:
        speciesTreespecification = input('Specie tree extension (will be used also for the S_edgelist extenstion)= ')
        k = int(input('k= '))
        TH_edges = int(input('Threshold Edges in Subtree= '))
        HT_cost = int(input('HT event cost= '))
        D_cost = int(input('Duplication event Cost= '))
        S_cost = int(input('Speciation event Cost= '))
        gamma = int(input('Gamma (1 recomended)= '))  # factor for probability assignment
        p = Decimal(input('p (0.05 recomended)= '))  # p_value
        number_of_planted_vertices = int(input('Number of Vertices to find= ')) # number of vertices that will be reported
        #save_data = bool(input('Save hypergraph data? (False/True) '))  # save hypergraph data
        single = input('Single mode pattern? (True/False) ') == 'True'  # save hypergraph data
        p1 = None
        p2 = None
        if single:
            while p1 == None:
                print('p1= ')
                ev = input("   EV=  (please provide at most 3 names out of 'S','D','HT', separated by space) ")  # save hypergraph data
                ev = ev.split()
                col = input("   color=  (red/black/None) ")
                dist = input("   distance=  (True/False) ") == 'True'
                if valid_pattern(ev,col,dist):
                    if col == 'None':
                        col = None
                    p1 = (ev,col,dist)
                    p2 = (None, None, False)
                else:
                    print('pattern is not valid, please try again: ')
        else:
            while p1 == None or p2 == None:
                print('p1= ')
                ev = input("   EV=  (please provide at most 3 names out of 'S','D','HT', separated by space) ") # save hypergraph data
                ev = ev.split()
                col = input("   color=  (red/black/None) ")
                dist = input("   distance=  (True/False) ") == 'True'
                if valid_pattern(ev,col,dist):
                    if col == 'None':
                        col = None
                    p1 = (ev,col,dist)
                else:
                    print('pattern is not valid, please try again: ')
                if p1 != None:
                    print('p2= ')
                    ev = input(
                        "   EV=  (please provide at most 3 names out of 'S','D','HT', separated by space) ")  # save hypergraph data
                    ev = ev.split()
                    col = input("   color=  (red/black/None) ")
                    dist = input("   distance=  (True/False) ") == 'True'
                    if valid_pattern(ev, col, dist):
                        if col == 'None':
                            col = None
                        p2 = (ev, col, dist)
                    else:
                        print('pattern is not valid, please try again: ')

    all_vertices_with_index = {}
    list_of_scores_for_rand_num = {}
    noise_level = 0
    rand_num = 0
    random_for_prec_curr = 1
    S_dis_matrix = {}
    nodes_table = {}
    S_colors = {}
    all_vertices = {}

    taxa_names.main(os.getcwd()+'/data/',[''],speciesTreespecification)
    fix.main(os.getcwd()+'/data/',[''])
    create_sigma.main(os.getcwd()+'/data/',[''],speciesTreespecification)

    G = tr.Tree.get_from_path(os.getcwd()+"/data/G.txt", schema="newick")
    S = tr.Tree.get_from_path(os.getcwd()+"/data/S_"+speciesTreespecification+".txt", schema="newick")

    input1 = open(os.getcwd()+'/data/sigma.txt', 'r')
    sigma = []
    for line in input1:
        sigma.append(eval(line))
    sigma = sigma[0]

    input1 = open(os.getcwd()+'/data/colors.txt', 'r')
    colors = []
    for line in input1:
        colors.append(eval(line))
    colors = colors[0]
    G.prune_taxa_with_labels(tree_operations.remove_unsigma_genes(G, sigma, True))

    S = utiles.init_internal_labels(S, 'x', sigma, os.getcwd()+'/data/')
    G = utiles.init_internal_labels(G, 'u', sigma, os.getcwd()+'/data/')

    G = tree_operations.collapse_edges(G)
    S = tree_operations.collapse_edges(S)

    n2e.main(os.getcwd()+'/data/',speciesTreespecification)

    S_labels_table, G_labels_table,sigma = inits.init_taxon_to_label_table(S,G,sigma)
    sigma, old_sigma = inits.update_sigma(sigma,True,S_labels_table,G_labels_table)
    if tree_operations.remove_unsigma_genes(G, sigma, False) is not []:
        G.prune_taxa_with_labels(tree_operations.remove_unsigma_genes(G, sigma, False))
    colors,old_colors = inits.update_colors(S, colors,True)

    p1 = (p1[0],p1[1],p1[2],len(tree_operations.leaf_in_subtrees(G,'S',G.seed_node.label, old_sigma,False)[0]+tree_operations.leaf_in_subtrees(G,'S',G.seed_node.label, old_sigma,False)[1])*TH_edges)

    S_dis_matrix = inits.init_distance_S(S_dis_matrix, os.getcwd()+'/data/' ,speciesTreespecification)
    nodes_table = inits.init_nodes_table(S, G, nodes_table)
    H, H_number_of_nodes, nodes_table = hypergraph.build_hyper_garph(S, G, k,
                                                                     nodes_table, D_cost, S_cost, HT_cost,os.getcwd(),
                                                                     sigma,save_data)

    new_G = nx.DiGraph()
    nCr_lookup_table = {}
    fact_lookup_table = {}
    max_score_p1_list = []
    max_score_p1_and_p2_list = []
    H, max_prob = hypergraph.assign_probabilities(S, G, H, gamma)
    if H is None:
        list_of_scores_for_rand_num.update({rand_num: {}})
    else:
        S_colors = tree_operations.color_tree(S, 'S', S_colors, colors, sigma)

        interesting_vertices_p1, nCr_lookup_table, fact_lookup_table = find_Pattern(H,S,S_dis_matrix,nCr_lookup_table,fact_lookup_table, p1,S_colors,p)
        interesting_vertices_p2, nCr_lookup_table, fact_lookup_table = find_Pattern(H,S,S_dis_matrix,nCr_lookup_table,fact_lookup_table, p2,S_colors,p)

        max_S_d_of_HT = tree_operations.find_max_d_of_HT(S_dis_matrix, interesting_vertices_p1,p1)

        new_G = tree_operations.weight_G_based_on_same_color_HT(G, new_G,interesting_vertices_p1,interesting_vertices_p2, max_S_d_of_HT,p1,p2,False)
        new_G = tree_operations.number_of_edges_in_subtree(new_G)
        new_G = tree_operations.normlize_weights(new_G,k,p1,'p1')
        new_G = tree_operations.normlize_weights(new_G,k,p2,'p2')

        if p2[0] is None:
            max_score_p1_list = tree_operations.find_max_scores(new_G,number_of_planted_vertices,'p1',p1[3])
        else:
            max_score_p1_and_p2_list = tree_operations.find_max_scores(new_G,number_of_planted_vertices,'p2',p1[3])
        marked_nodes,all_vertices = pattern_identify.find_signi_distance(new_G, all_vertices,p1, p2, max_score_p1_list, max_score_p1_and_p2_list, False)

        list_of_scores_for_rand_num.update({rand_num:all_vertices})
    all_vertices_with_index.update({noise_level:utiles.average_of_list(list_of_scores_for_rand_num,random_for_prec_curr)})

    end_function(H,S,G,k,starting_time,p1,p2,marked_nodes,old_sigma,max_score_p1_list,max_score_p1_and_p2_list)

if __name__ == "__main__":
    main()


