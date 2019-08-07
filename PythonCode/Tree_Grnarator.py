import sys
sys.path.append('/anaconda3/lib/python3.6/site-packages')
from ete3 import Tree
import random
from random import randint
import utils.newick2edgelist
import networkx as nx
import dendropy as tr
import tree_operations_v1 as tree_operations
import utiles
import inits_v1 as inits
import draw
import os
import numpy
from functools import reduce
from multiprocessing import Pool
from datetime import datetime

on_lab = False
compare = True
running_time = False
minimum_HT_under_planted = 3
number_of_leaves = 600
if on_lab:
    if compare:
        path = '/storage/DATA/users/ronizo/comparsion_600_20_planted'
    else:
        path = '/storage/DATA/users/ronizo/noise_data_500_k=100'
else:
    if compare:
        path = '/Users/ronizoller/PycharmProjects/TreeReconciliation/trees/compare_test'
    else:
        path = '/Users/ronizoller/PycharmProjects/TreeReconciliation/trees/duplications_test'

add_noise = False
number_of_planted_vertices = 20
S = Tree()
G = Tree()
k = 50
both = False
TH_both = 0
compare_subtrees = True
evolutinary_event = 'HT'
noise_level = [20]
number_of_nodes = 0
random_for_precentage = 1                             #number of different random noise for each noise %
accur = 5
p = 0.05                                                #p_value
TH_edges_in_subtree = 10                                # smallest subtree that will be counted when not comparing subtrees
TH_pattern_in_subtree = 0
if compare_subtrees and evolutinary_event=='HT':
    TH_compare_subtrees =  2
    pattern = "same_color"
sym = 'Specie'


def randome_leave_from_tree(t,number_of_leaves,dend):
    ran = (randint(1,number_of_leaves))
    i = 1
    if dend:
        for leaf in t.leaf_node_iter():
            if i == ran:
                return leaf
            else:
                i += 1
    else:
        for leaf in t.iter_leaves():
            if i == ran:
                return leaf
            else:
                i += 1


def random_vertex_in_tree(number_of_vertex,tree):
    ran = (randint(1, number_of_vertex))
    p = 0
    for u in tree.postorder_node_iter():
        if p == ran-1:
            return u
        p += 1


def randome_leave_from_node (nd):
    to_return = random.choice(nd.leaf_nodes())
    return to_return


def random_again(t,iter):
    number_of_iter = 0
    for node in t.traverse("postorder"):
        if node != t and number_of_iter < iter:
            remove = random.choice([True, False])
            if remove:
                removed_subtree = node.detach()
                print('removed_subtree: '+str(removed_subtree))
                randome_leave_from_tree(t,len(t),False).add_child(removed_subtree)
        number_of_iter += 1
    return t

def number_of_HT_needed(G,u,all_random_sources,TH,color):
    print('     Calculating number of HT under %s' % str(u))
    if color == 'red':
        to_chack_in = all_random_sources[0]
    if color == 'black':
        to_chack_in = all_random_sources[1]
    if color == 'all':
        to_chack_in = all_random_sources[0]+all_random_sources[1]
    bad_HT = 0
    u_in_G = G.find_node(lambda n: (n.label == u['label']))
    for leaf in u_in_G.leaf_nodes():
        if leaf in to_chack_in:
            bad_HT += 1
    if int(TH*bad_HT) > 0:
        print('     number of HT needed under %s is %s bad HT: %s' % (str(u), str((TH * bad_HT)+minimum_HT_under_planted), str(bad_HT)))
        return int(TH*bad_HT)+minimum_HT_under_planted
    else:
        print('     number of HT needed under %s is %s bad HT: %s' % (str(u), str(minimum_HT_under_planted), str(bad_HT)))
        return minimum_HT_under_planted

def create_good_HT(G,S,nCr_lookup_table,fact_lookup_table,number_of_HT,u,child_to_fhind_HT_in,Pr_red,Pr_black,color,max_dis,S_dis_matrix,new_G,G_internal_colors):
    HT_sources_for_subtree = [None]     #to chack if HT has already been choose
    HT_targets_for_subtree = {}
    list_of_answers = []               #to return
    i = 0
    while i < number_of_HT and i < u['edges_in_subtree']:
        print('     Start looking for HT source under %s (%sth HT)' % (str(child_to_fhind_HT_in),str(i)))
        nCr_lookup_table, fact_lookup_table,HT = find_enriched_subtree_rec(G,S, new_G, child_to_fhind_HT_in, nCr_lookup_table, fact_lookup_table,
                                       accur, Pr_red, Pr_black, p, color,
                                       G_internal_colors,HT_sources_for_subtree)
        j = 0
        while HT[1] in HT_sources_for_subtree and j < u['edges_in_subtree']:
            print('     Looking for %sth HT in %s' % (str(i),str(child_to_fhind_HT_in)))
            nCr_lookup_table, fact_lookup_table,HT = find_enriched_subtree_rec(G, S,new_G, child_to_fhind_HT_in, nCr_lookup_table, fact_lookup_table,
                                                                  accur, Pr_red, Pr_black, p, color,
                                                                  G_internal_colors,HT_sources_for_subtree)
            j += 1
        if HT[0]:
            print('     Found HT source')
            print('     Start looking for HT traget')
            HT_random_son = tree_operations.random_son(new_G, HT)
            HT_sources_for_subtree.append(HT[1])
            HT_source_in_S = S.find_node(lambda n: (n.label == 'x' + HT[1]['label'][1:]))
            HT_to = (False, None)
            for x in S.postorder_node_iter():
                dis = S_dis_matrix[((HT_source_in_S.label, x.label))]
                if not HT_to[0]:
                    if ((tree_operations.is_not_ancestor(HT_source_in_S, x)) and (
                    tree_operations.is_not_ancestor(x, HT_source_in_S))) and dis > max_dis/2 and len(x.leaf_nodes()) > minimum_HT_under_planted:
                        x = tree_operations.find_node_in_networkx_tree(new_G, 'u' + x.label[1:])
                        nCr_lookup_table, fact_lookup_table,HT_to = find_enriched_subtree_rec(G, S,new_G, x, nCr_lookup_table,
                                                          fact_lookup_table, accur, Pr_red,
                                                          Pr_black, p, color, G_internal_colors,[])
                        #print('Dis between %s and %s is %s, max: %s' % (str(HT_source_in_S.label), str(x), str(dis),str(max_dis)))
            if HT_to[0]:
                print('HT_to: '+str(HT_to[1])+' HT_targets_for_subtree: '+str(HT_targets_for_subtree))
                if HT_to[1]['label'] in HT_targets_for_subtree:                  ## it is possible to map only 2 verices to 1 species
                    if HT_targets_for_subtree[HT_to[1]['label']] == 1:
                        print('     Found HT target')
                        print('     HT: %s' % str(HT))
                        print('     HT_to: %s\n' % str(HT_to))
                        list_of_answers.append((HT_random_son['label'], HT_to[1]['label']))
                        HT_targets_for_subtree[HT_to[1]['label']] = 2
                        i += 1
                else:
                    list_of_answers.append((HT_random_son['label'], HT_to[1]['label']))
                    HT_targets_for_subtree.update({HT_to[1]['label']:1})
                    i += 1
            else:
                return nCr_lookup_table, fact_lookup_table, (False, list_of_answers)
        else:
            return nCr_lookup_table,fact_lookup_table,(False,list_of_answers)
    print('     list of HT:'+str(list_of_answers))
    return nCr_lookup_table,fact_lookup_table,(child_to_fhind_HT_in['label'],list_of_answers)

def change_color_of_vertex(u,colors,color,sigma,switch_colors):
    if switch_colors:
        if colors[sigma[u.label]] == 'red':
            color = 'black'
        elif colors[sigma[u.label]] == 'black':
            color = 'red'
        else:
            color = 'pink'
        #print ('Changing %s color from %s to %s' % (str(u),str(colors[sigma[u.label]]),str(color)))
    colors.update({sigma[u.label]:color})
    return colors

def change_colors(v, w, colors, G_internal_colors,S_colors,G,S, sigma):
    old_colors = colors.copy()
    old_S = S
    old_S_colors = S_colors.copy()
    old_internal_G_colors = G_internal_colors.copy()
    reds_under_w = G_internal_colors[w['label']][0]
    blacks_under_w = G_internal_colors[w['label']][1]
    reds_under_v = G_internal_colors[v['label']][0]
    blacks_under_v = G_internal_colors[v['label']][1]
    all_leafs_v = reds_under_v + blacks_under_v
    all_leafs_w = reds_under_w + blacks_under_w
    if (blacks_under_w / all_leafs_w) > (reds_under_v / all_leafs_v):
        color = 'black'
        to_return_color = 'red'
        color_under_father = blacks_under_w
        father_to_change = w
        to_return = v
        all_leafs_father = all_leafs_w
    else:
        color = 'red'
        to_return_color = 'black'
        color_under_father = reds_under_v
        father_to_change = v
        to_return = w
        all_leafs_father = all_leafs_v
    father_to_change_in_G = G.find_node(lambda n: (n.label == father_to_change['label']))
    for to_change in father_to_change_in_G.leaf_nodes():
        #print('     Color: %s, vertex: %s, color under/all leaf: %s, TH: %s' % (
        #str(color), str(father_to_change), str(color_under_father / all_leafs_father), str(TH_both)))
        if color_under_father/all_leafs_father < TH_both:
            colors = change_color_of_vertex(to_change,colors,color,sigma,False)
            G_internal_colors = tree_operations.color_tree(G, 'G', G_internal_colors, colors, sigma)
            if color == 'red':
                color_under_father = G_internal_colors[father_to_change['label']][0]
                other_color_under_father = G_internal_colors[father_to_change['label']][1]
            else:
                other_color_under_father = G_internal_colors[father_to_change['label']][0]
                color_under_father = G_internal_colors[father_to_change['label']][1]
            all_leafs_father = color_under_father + other_color_under_father
    S_colors = tree_operations.color_tree(S, 'S', S_colors, colors, sigma)
    G_internal_colors = tree_operations.color_tree(G, 'G', G_internal_colors, colors, sigma)
    return colors,to_return,to_return_color,G_internal_colors,old_colors,old_internal_G_colors,S,S_colors,old_S,old_S_colors

def choose_planted_vertex (S_dis_matrix,new_G,S,G,G_internal_colors,TH_edges_in_subtree,compare_subtrees,TH_compare_subtrees,sigma,k,both,TH_both,vertex_number,sol,accur,nCr_lookup_table,fact_lookup_table,all_random_sources,colors,S_colors,max_dis):
    total_red = S_colors[S.seed_node.label][0]
    total_black = S_colors[S.seed_node.label][1]
    num_of_leafs = tree_operations.number_of_leafs(S,'S')
    Pr_red = total_red / num_of_leafs
    Pr_black = total_black / num_of_leafs
    for u in list(reversed(list(nx.topological_sort(new_G)))):
        old_colors = colors.copy()
        old_internal_G_colors = G_internal_colors.copy()
        old_S = S
        old_S_colors = S_colors.copy()
        ##this is the random vertex
        outgoing_edges = new_G.out_edges([u], data=True)
        outgoing_edges = [e for e in outgoing_edges]
        u = new_G.nodes(data=True)[u]
        if len(outgoing_edges) == 2:
            if compare_subtrees:
                if u['edges_in_subtree'] > TH_edges_in_subtree and not check_if_vertex_was_chosen(vertex_number,sol,u['label']):
                    print('\n***\nVertex ' + str(u) + ' was chosen to be marked.')
                    if not both:
                        number_of_HT = number_of_HT_needed(G,u, all_random_sources, TH_compare_subtrees, 'red')
                        nCr_lookup_table, fact_lookup_table,ans = create_good_HT(G, S, nCr_lookup_table,fact_lookup_table,number_of_HT, u, u, Pr_red, Pr_black, 'red',
                                      max_dis,S_dis_matrix,new_G,G_internal_colors)

                        if ans[0] == u['label']:
                            return nCr_lookup_table,fact_lookup_table,ans,colors
                        else:
                            number_of_HT = number_of_HT_needed(G, u, all_random_sources, TH_compare_subtrees, 'black')
                            nCr_lookup_table, fact_lookup_table, ans = create_good_HT(G, S, nCr_lookup_table,
                                                                                      fact_lookup_table, number_of_HT,
                                                                                      u, u, Pr_red, Pr_black, 'black',
                                                                                      max_dis, S_dis_matrix, new_G,
                                                                                      G_internal_colors)
                        colors = old_colors
                        G_internal_colors = old_internal_G_colors
                        S = old_S
                        S_colors = old_S_colors
    print('No vertex could be planted.')
    return nCr_lookup_table,fact_lookup_table,(False,False),colors

def find_enriched_subtree(G,v,color,reds_under_v, blacks_under_v,all_leafs_v,nCr_lookup_table, fact_lookup_table,
                                                                                    accur, Pr_red,
                                                                                    Pr_black,p):
    #print('     Checking if %s is an enriched subtree...' % str(v))
    v = G.find_node(lambda n: (n.label == v['label']))
    if color == 'red':
        p_value_v_red, nCr_lookup_table, fact_lookup_table = utiles.p_value_calculation(reds_under_v, all_leafs_v,
                                                                                    nCr_lookup_table, fact_lookup_table,
                                                                                    accur, Pr_red,
                                                                                    Pr_black, v.label)

        if p_value_v_red < p:
            #print('         %s is an enriched with RED subtree.' % str(v))
            return 1,nCr_lookup_table,fact_lookup_table,p_value_v_red
        else:
            return 0, nCr_lookup_table, fact_lookup_table,0
    elif color == 'black':
        p_value_v_black, nCr_lookup_table, fact_lookup_table = utiles.p_value_calculation(blacks_under_v,
                                                                                      all_leafs_v,
                                                                                      nCr_lookup_table,
                                                                                      fact_lookup_table,
                                                                                      accur, Pr_black,
                                                                                      Pr_red,
                                                                                      v.label)
        if p_value_v_black < p:
            #print('         %s is an enriched with BLACK subtree.' % str(v))
            return 1,nCr_lookup_table,fact_lookup_table,p_value_v_black
        else:
            return 0, nCr_lookup_table, fact_lookup_table,0


def find_enriched_subtree_rec(G,S,new_G,u,nCr_lookup_table,fact_lookup_table,accur,Pr_red,Pr_black,p,color,G_internal_colors,bad_vertices):
    reds_under_u = G_internal_colors[u['label']][0]
    blacks_under_u = G_internal_colors[u['label']][1]
    all_leafs_u = reds_under_u + blacks_under_u
    u_number = u['ind']
    outgoing_edges = new_G.out_edges([u_number], data=True)
    outgoing_edges = [e for e in outgoing_edges]
    if len(outgoing_edges) > 0:
        ans_u, nCr_lookup_table, fact_lookup_table,p_value = find_enriched_subtree(G,u, color, reds_under_u, blacks_under_u,all_leafs_u, nCr_lookup_table, fact_lookup_table,accur, Pr_red,Pr_black, p)
        if ans_u == 1 and u not in bad_vertices:
            print('     FOUND!')
            print('     curr = %s\n     curr_red = %s, curr_black = %s, p_value_curr = %s\n        Pr_red = %s, Pr_black = %s' % (
                str(u), str(reds_under_u), str(blacks_under_u), str(p_value),str(Pr_red),str(Pr_black)))

            return nCr_lookup_table,fact_lookup_table,(True, u)
        elif len(outgoing_edges) == 2:
            v = new_G.nodes(data=True)[outgoing_edges[0][1]]
            w = new_G.nodes(data=True)[outgoing_edges[1][1]]
            nCr_lookup_table, fact_lookup_table,(ans_v, v) = find_enriched_subtree_rec(G,S,new_G,v,nCr_lookup_table,fact_lookup_table,accur,Pr_red,Pr_black,p,color,G_internal_colors,bad_vertices)
            nCr_lookup_table, fact_lookup_table,(ans_w,w) = find_enriched_subtree_rec(G,S,new_G,w,nCr_lookup_table,fact_lookup_table,accur,Pr_red,Pr_black,p,color,G_internal_colors,bad_vertices)
            if ans_v:
                return nCr_lookup_table,fact_lookup_table,(ans_v,v)
            elif ans_w:
                return nCr_lookup_table,fact_lookup_table,(ans_w,w)
            else:
                return nCr_lookup_table,fact_lookup_table,(False, u)
        else:
            return nCr_lookup_table,fact_lookup_table,(False, u)
    else:
        return nCr_lookup_table,fact_lookup_table,(False, u)

def random_colors(S,colors):
    for leaf in S.iter_leaves():
        col = ['black', 'red']
        colors.update({leaf.name:col[random.randint(0, 1)]})
    return colors

def count_nodes_and_update_internal_names(t):
    i = 1
    res = 0
    for node in t.traverse("postorder"):
        if node.name == "":
            node.name = "internal" + str(i)
        i += 1
        res += 1
    return res-2

def create_sigme(number_of_leaves,sigma):
    for i in range(0,number_of_leaves):
        sigma.update({"Gene"+str(i):"Specie"+str(i)})
        sigma.update({"GeneI" + str(i): "internal" + str(i)})
    return sigma

def print_tree(t,name,path):
    print('** Tree: %s **\n     Number of leaves: %s\n     Number of nodes: %s\n       %s\n%s' % (
    name, str(len(t)), number_of_nodes, t.write(format=1), str(t)))
    print('     Writing '+name)
    if name == 'S':
        file = open(path + '/phyliptree(binary,all).phy', 'w')
        file.write(str(t.write(format=1)))
        file.close()
    if name == 'G':
        print('     Writing G')
        file = open(path + '/GeneTree(binary)_local.txt', 'w')
        file.write(str(t.write(format=1)))
        file.close()

def save_edgelist(S_dis_matrix,path):
    path_curr = path + '/saved_data/S_edgelist_all.txt'
    os.makedirs(os.path.dirname(path_curr), exist_ok=True)
    input = open(path + '/saved_data/S_edgelist_all.txt', 'r')
    edgelist = []
    for line in input:
        edgelist.append(eval(line))
    edgelist = edgelist[0]
    networkx_S = nx.parse_edgelist(edgelist, nodetype=str, data=(('weight', any),))
    for nd1 in networkx_S.nodes():
        for nd2 in networkx_S.nodes():
            S_dis_matrix.update(
                {(nd1, nd2): nx.shortest_path_length(networkx_S, source=nd1, target=nd2)})
    #print('     Writing S_edge_list...')
    file = open(path + '/saved_data/S_dist_matrix.txt', 'w')
    file.write(str(S_dis_matrix))
    file.close()

def change_sigma(sigma, old_sigma, S, G, couples_list,number_of_HT_to_make,S_labels_table,G_labels_table):
    #print('Updating sigma..')
    for i in range(0,len(couples_list)):
        HT = couples_list[i][0]
        HT_to = couples_list[i][1]
        HT_in_G = G.find_node(lambda n: (n.label == HT))
        HT_to_in_S = S.find_node(lambda n: (n.label == 'x'+HT_to[1:]))

        r_leafs_HT,l_leafs_HT = tree_operations.leaf_in_subtrees(G,'G',HT_in_G,old_sigma,True)
        all_HT_leafs = r_leafs_HT+l_leafs_HT
        if len(all_HT_leafs) == 0:
            all_HT_leafs = [HT_in_G.taxon.label]
        r_leafs_HT_to, l_leafs_HT_to = tree_operations.leaf_in_subtrees(S,'S', HT_to_in_S, old_sigma,True)
        all_leafs_HT_to = r_leafs_HT_to+l_leafs_HT_to

        number_of_HT_to_make = min([number_of_HT_to_make,len(all_HT_leafs),len(all_leafs_HT_to)])
        if len(all_leafs_HT_to) == 0:
            all_leafs_HT_to = [HT_to_in_S.taxon.label]
        old_sigma_temp = old_sigma.copy()
        old_sigma = dict((u,x) for u,x in old_sigma.items() if (not u in all_HT_leafs))

        changed = 0
        for u in all_HT_leafs:
            if changed < number_of_HT_to_make:
                old_sigma.update({u:all_leafs_HT_to[changed]})
                changed += 1
            else:
                old_sigma.update({u: old_sigma_temp[u]})
        sigma, old_sigma = inits.update_sigma(S, G, 0, old_sigma, False, path,True,S_labels_table,G_labels_table)

    return sigma,old_sigma, number_of_HT_to_make

def save_data(sigma,colors,planted,noise,internal_index,compare,path):
    path_curr = path +'/'+ str(noise)+'/sigma'+str(noise)+'.'+str(internal_index)+'.txt'
    os.makedirs(os.path.dirname(path_curr), exist_ok=True)
    file = open(path_curr, 'w')
    file.write(str(sigma))
    file.close()
    path_curr = path +'/'+ str(noise)+'/colors'+str(noise)+'.'+str(internal_index)+'.txt'
    os.makedirs(os.path.dirname(path_curr), exist_ok=True)
    file = open(path_curr, 'w')
    file.write(str(colors))
    file.close()
    path_curr = path +'/'+ str(noise)+'/planted'+str(noise)+'.'+str(internal_index)+'.txt'
    os.makedirs(os.path.dirname(path_curr), exist_ok=True)
    file = open(path_curr, 'w')
    file.write(str(planted))
    file.close()

def check_if_vertex_was_chosen(j,sol,curr_sol):
    for p in range(0, j):  # check if this solution was already chosen
        if sol[p] != None:
            if sol[p]['Marked'] == curr_sol:
                return True
    return False

def return_color_to_taxon(S,colors):
    to_return = {}
    for u,color in colors.items():
        node = S.find_node(lambda n: (n.label == u))
        to_return.update({node.taxon.label:color})
    return to_return

def check_if_enriched(S,G,source, target,nCr_lookup_table,fact_lookup_table, accur,G_internal_colors,S_colors):
    reds_under_source = G_internal_colors[source.label][0]
    blacks_under_source = G_internal_colors[source.label][1]
    reds_under_target = G_internal_colors[target.label][0]
    blacks_under_target = G_internal_colors[target.label][1]

    total_source = reds_under_source + blacks_under_source
    total_target = reds_under_target + blacks_under_target
    total_red = S_colors[S.seed_node.label][0]
    total_black = S_colors[S.seed_node.label][1]
    num_of_leafs = tree_operations.number_of_leafs(S, 'S')
    Pr_red = total_red / num_of_leafs
    Pr_black = total_black / num_of_leafs

    #print('source: %s, target: %s, reds_under_source: %s, blacks_under_source: %s, reds_under_target: %s, blacks_under_target: %s' % (str(source),str(target),str(reds_under_source), str(blacks_under_source),str(reds_under_target),str(blacks_under_target)))


    p_value_source_red, nCr_lookup_table, fact_lookup_table = utiles.p_value_calculation(reds_under_source, total_source,
                                                                                       nCr_lookup_table,
                                                                                       fact_lookup_table, accur,
                                                                                       Pr_red,
                                                                                       Pr_black, source.label)
    p_value_source_black, nCr_lookup_table, fact_lookup_table = utiles.p_value_calculation(blacks_under_source,
                                                                                         total_source,
                                                                                         nCr_lookup_table,
                                                                                         fact_lookup_table,
                                                                                         accur, Pr_black,
                                                                                         Pr_red,
                                                                                         source.label)
    p_value_target_red, nCr_lookup_table, fact_lookup_table = utiles.p_value_calculation(reds_under_target,
                                                                                         total_target,
                                                                                         nCr_lookup_table,
                                                                                         fact_lookup_table, accur,
                                                                                         Pr_red,
                                                                                         Pr_black, target.label)
    p_value_target_black, nCr_lookup_table, fact_lookup_table = utiles.p_value_calculation(blacks_under_target,
                                                                                         total_target,
                                                                                         nCr_lookup_table,
                                                                                         fact_lookup_table, accur,
                                                                                         Pr_black,
                                                                                         Pr_red,
                                                                                         target.label)
    #print('     p_value_source_red: %s\n        p_value_target_red:%s' % (
    #str(p_value_source_red), str(p_value_target_red)))
    #print('     p_value_source_black: %s\n        p_value_target_black:%s' % (
    #str(p_value_source_black), str(p_value_target_black)))

    if p_value_source_red < p and p_value_target_red < p:
        return nCr_lookup_table,fact_lookup_table,'red-to-red'
    if p_value_source_black < p and p_value_target_black < p:
        return nCr_lookup_table,fact_lookup_table,'black-to-black'
    else: return nCr_lookup_table,fact_lookup_table,'nothing-to-nothing'

def create_tree_for_HT_noise(parameters):
    print('Creating noise in HT')
    noise = parameters[0]
    number_of_HT_under_planted = parameters[1]
    G_internal_colors = parameters[2]
    S_colors = parameters[3]
    nCr_lookup_table = parameters[4]
    fact_lookup_table = parameters[5]
    number_of_leaves = parameters[6]

    number_of_random_changes = number_of_leaves * (noise / 100)
    random_for_prec = random_for_precentage
    for rand_num in range(0, random_for_prec):
        input = open(path + '/0/sigma0.0' + '.txt', 'r')
        sigma = []
        for line in input:
            sigma.append(eval(line))
        sigma = sigma[0]

        input = open(path + '/0/colors0.0' + '.txt', 'r')
        colors = []
        for line in input:
            colors.append(eval(line))
        colors = colors[0]

        S = tr.Tree.get_from_path(path + "/phyliptree(binary,all).phy", schema="newick")
        G = tr.Tree.get_from_path(path + "/GeneTree(binary)_local.txt", schema="newick")

        S = utiles.init_internal_labels(S, 'x', sigma, path)
        G = utiles.init_internal_labels(G, 'u', sigma, path)

        G = tree_operations.collapse_edges(G)
        S = tree_operations.collapse_edges(S)

        S_labels_table, G_labels_table,sigma = inits.init_taxon_to_label_table(S, G, sigma)

        sigma, old_sigma = inits.update_sigma(S, G, 0, sigma, False, path, True, S_labels_table, G_labels_table)
        colors, old_colors = inits.update_colors(S, colors, True)

        i = 0
        all_random_sources_red_to_red = []
        all_random_sources_black_to_black = []
        all_random_nutral = []
        while i < number_of_random_changes:
            G_internal_colors = tree_operations.color_tree(G, 'G', G_internal_colors, colors, sigma)
            S_colors = tree_operations.color_tree(S, 'S', S_colors, colors, sigma)
            number_of_nodes = tree_operations.number_of_leafs(G, 'G')
            random_source = random_vertex_in_tree(number_of_nodes, G)
            random_target = random_vertex_in_tree(number_of_nodes, G)
            sigma, old_sigma, changed = change_sigma(sigma, old_sigma, S, G,
                                                     [(random_source.label, random_target.label)],
                                                     number_of_HT_under_planted,S_labels_table,G_labels_table)
            nCr_lookup_table, fact_lookup_table, enriched = check_if_enriched(S,G,random_source, random_target,
                                                                              nCr_lookup_table, fact_lookup_table,
                                                                              accur, G_internal_colors,S_colors)
            if enriched == 'red-to-red':
                all_random_sources_red_to_red.append(random_source)
            elif enriched == 'black_to_balck':
                all_random_sources_black_to_black.append(random_source)
            else:
                all_random_nutral.append(random_source)
            i += changed
        old_colors = return_color_to_taxon(S, colors)
        save_data(old_sigma, old_colors, {}, noise, rand_num,compare,path+'/HT')

def create_tree_for_HT_and_colors_noise(parameters):
    print('Creating noise in HT and colors')
    noise = parameters[0]
    number_of_HT_under_planted = parameters[1]
    G_internal_colors = parameters[2]
    S_colors = parameters[3]
    nCr_lookup_table = parameters[4]
    fact_lookup_table = parameters[5]
    number_of_leaves = parameters[6]

    number_of_random_changes = number_of_leaves * (noise / 100)
    random_for_prec = random_for_precentage
    for rand_num in range(0, random_for_prec):
        input = open(path + '/0/sigma0.0' + '.txt', 'r')
        sigma = []
        for line in input:
            sigma.append(eval(line))
        sigma = sigma[0]

        input = open(path + '/0/colors0.0' + '.txt', 'r')
        colors = []
        for line in input:
            colors.append(eval(line))
        colors = colors[0]

        S = tr.Tree.get_from_path(path + "/phyliptree(binary,all).phy", schema="newick")
        G = tr.Tree.get_from_path(path + "/GeneTree(binary)_local.txt", schema="newick")

        S = utiles.init_internal_labels(S, 'x', sigma, path)
        G = utiles.init_internal_labels(G, 'u', sigma, path)

        G = tree_operations.collapse_edges(G)
        S = tree_operations.collapse_edges(S)

        S_labels_table, G_labels_table,sigma = inits.init_taxon_to_label_table(S, G, sigma)

        sigma, old_sigma = inits.update_sigma(S, G, 0, sigma, False, path, True, S_labels_table, G_labels_table)
        colors, old_colors = inits.update_colors(S, colors, True)

        i = 0
        all_random_sources_red_to_red = []
        all_random_sources_black_to_black = []
        all_random_nutral = []
        while i < number_of_random_changes:
            G_internal_colors = tree_operations.color_tree(G, 'G', G_internal_colors, colors, sigma)
            S_colors = tree_operations.color_tree(S, 'S', S_colors, colors, sigma)
            number_of_nodes = tree_operations.number_of_leafs(G, 'G')
            random_source = random_vertex_in_tree(number_of_nodes, G)
            random_target = random_vertex_in_tree(number_of_nodes, G)
            random_vertex_to_change_color = randome_leave_from_tree(G, number_of_nodes, True)
            sigma, old_sigma, changed = change_sigma(sigma, old_sigma, S, G,
                                                     [(random_source.label, random_target.label)],
                                                     number_of_HT_under_planted,S_labels_table,G_labels_table)
            colors = change_color_of_vertex(random_vertex_to_change_color, colors, None, sigma, True)
            nCr_lookup_table, fact_lookup_table, enriched = check_if_enriched(S,G,random_source, random_target,
                                                                              nCr_lookup_table, fact_lookup_table,
                                                                              accur, G_internal_colors,S_colors)
            if enriched == 'red-to-red':
                all_random_sources_red_to_red.append(random_source)
            elif enriched == 'black_to_balck':
                all_random_sources_black_to_black.append(random_source)
            else:
                all_random_nutral.append(random_source)
            i += changed
        old_colors = return_color_to_taxon(S, colors)
        save_data(old_sigma, old_colors, {}, noise, rand_num,compare,path+"/colors_and_HT")

def create_tree_for_color_noise(parameters):
    print('Creating noise in colors')
    noise = parameters[0]
    G_internal_colors = parameters[2]
    S_colors = parameters[3]
    number_of_leaves = parameters[6]

    number_of_random_changes = number_of_leaves * (noise / 100)
    random_for_prec = random_for_precentage
    for rand_num in range(0, random_for_prec):
        input = open(path + '/0/sigma0.0' + '.txt', 'r')
        sigma = []
        for line in input:
            sigma.append(eval(line))
        sigma = sigma[0]
        input = open(path + '/0/colors0.0' + '.txt', 'r')
        colors = []
        for line in input:
            colors.append(eval(line))
        colors = colors[0]

        S = tr.Tree.get_from_path(path + "/phyliptree(binary,all).phy", schema="newick")
        G = tr.Tree.get_from_path(path + "/GeneTree(binary)_local.txt", schema="newick")

        S = utiles.init_internal_labels(S, 'x', sigma, path)
        G = utiles.init_internal_labels(G, 'u', sigma, path)

        G = tree_operations.collapse_edges(G)
        S = tree_operations.collapse_edges(S)

        S_labels_table, G_labels_table,sigma = inits.init_taxon_to_label_table(S, G, sigma)

        sigma, old_sigma = inits.update_sigma(S, G, 0, sigma, False, path, True, S_labels_table, G_labels_table)
        colors, old_colors = inits.update_colors(S, colors, True)

        i = 0
        while i < number_of_random_changes:
            G_internal_colors = tree_operations.color_tree(G, 'G', G_internal_colors, colors, sigma)
            S_colors = tree_operations.color_tree(S, 'S', S_colors, colors, sigma)
            number_of_nodes = tree_operations.number_of_leafs(G, 'G')
            random_vertex_to_change_color = randome_leave_from_tree(G, number_of_nodes, True)
            colors = change_color_of_vertex(random_vertex_to_change_color, colors, None, sigma, True)
            i += 1
        old_colors = return_color_to_taxon(S, colors)
        save_data(old_sigma, old_colors, {}, noise, rand_num,compare,path+'/color')

def return_planted_nodes_new_name(list_of_planted,prev_G,path):
    list_of_planted_nodes = []
    for k, l in list_of_planted.items():
        list_of_planted_nodes.append(l['Marked'])
    main_G = tr.Tree.get_from_path(path + "/GeneTree(binary)_local.txt", schema="newick")
    input = open(path + '/0/sigma0.0.txt', 'r')
    main_sigma = []
    for line in input:
        main_sigma.append(eval(line))
    main_sigma = main_sigma[0]
    main_G.prune_taxa_with_labels(tree_operations.remove_unsigma_genes(main_G, main_sigma, True))
    main_G = utiles.init_internal_labels(main_G, 'u', main_sigma, path)
    res = []
    for u in list_of_planted_nodes:
        planted_in_G = prev_G.find_node(lambda n: (n.label == u))
        to_find_leaf_nodes = list(map((lambda x:x.taxon.label),planted_in_G.leaf_nodes()))
        for v in main_G.postorder_node_iter():
            curr_leaf_nodes = list(map((lambda x:x.taxon.label),v.leaf_nodes()))
            if (len(curr_leaf_nodes) == len(to_find_leaf_nodes)) and (reduce((lambda x, y: x and y), list(numpy.in1d(to_find_leaf_nodes, curr_leaf_nodes)))):
                res.append(v.label)
    file = open(path + '/saved_data/planted_nodes_correct_names.txt', 'w')
    file.write(str(res))
    file.close()

##********  MAIN ***********

def main(S,G,number_of_leaves,path,k,running_time,number_of_planted_vertices):
    global random_for_precentage,all_edges,TH_edges_in_subtree,compare_subtrees,TH_pattern_in_subtree,TH_compare_subtrees,both,TH_both,accur
    starting_time = datetime.now()
    new_G = nx.DiGraph()
    noise = 0
    number_of_HT_under_planted = 10
    S = Tree()
    sigma = {}
    nCr_lookup_table = {}
    fact_lookup_table = {}
    colors = {}
    S_dis_matrix = {}
    names = []
    S_colors = {}
    G_internal_colors = {}
    sol = {}

    for i in range(0, number_of_leaves):
        names.append(sym + str(i))
    S.populate(number_of_leaves, names_library=names)
    count_nodes_and_update_internal_names(S)
    S = random_again(S, number_of_leaves / 4)
    colors = random_colors(S, colors)

    G = S.copy("newick")
    for leaf in G.iter_leaves():
        if leaf.name[:6] == 'Specie':
            leaf.name = "Gene" + leaf.name[6:]
        else:
            leaf.name = "GeneI" + leaf.name[8:]
    print_tree(G,'G',path)
    print_tree(S,'S',path)

    sigma = create_sigme(number_of_leaves, sigma)
    utils.newick2edgelist.main(path)
    save_edgelist(S_dis_matrix,path)

    S = tr.Tree.get_from_path(path + "/phyliptree(binary,all).phy", schema="newick")
    G = tr.Tree.get_from_path(path + "/GeneTree(binary)_local.txt", schema="newick")

    S = utiles.init_internal_labels(S, 'x', sigma, path)
    G = utiles.init_internal_labels(G, 'u', sigma, path)

    G = tree_operations.collapse_edges(G)
    S = tree_operations.collapse_edges(S)

    S_labels_table, G_labels_table,sigma = inits.init_taxon_to_label_table(S, G, sigma)

    sigma, old_sigma = inits.update_sigma(S, G, 0, sigma, False, path, True, S_labels_table, G_labels_table)
    colors, old_colors = inits.update_colors(S, colors, True)
    max_dis = tree_operations.max_dis(S_dis_matrix)

    flag = True
    j = 0
    all_random_sources_red_to_red = []
    all_random_sources_black_to_black = []
    all_random_nutral = []
    all_random_sources = (all_random_sources_red_to_red, all_random_sources_black_to_black, all_random_nutral)
    new_G = tree_operations.weight_G_based_on_same_color_HT(G, new_G, [],
                                                            [],[],[], 0, False,
                                                            'HT', False, k)
    new_G = tree_operations.number_of_edges_in_subtree(new_G)

    S_colors = tree_operations.color_tree(S, 'S', S_colors, colors, sigma)
    G_internal_colors = tree_operations.color_tree(G, 'G', G_internal_colors, colors, sigma)
    if not on_lab:
        draw.draw_S_and_G(S, G, old_sigma, colors, sigma, path, None, '_rand_before')

    if not running_time:
        while j < number_of_planted_vertices:
            print(
                '                                                                 *****         %sth vertex            ******' % str(
                    j))
            sol[j] = {}
            nCr_lookup_table, fact_lookup_table, (
                sol[j]['Marked'], sol[j]['list_of_couples']), colors = choose_planted_vertex(S_dis_matrix,new_G, S, G,
                                                                                            G_internal_colors,
                                                                                            TH_edges_in_subtree,
                                                                                            compare_subtrees,
                                                                                            TH_compare_subtrees,
                                                                                            sigma,
                                                                                            k,
                                                                                            both,
                                                                                            TH_both, j, sol, accur,
                                                                                            nCr_lookup_table,
                                                                                            fact_lookup_table,
                                                                                            all_random_sources,
                                                                                            colors,
                                                                                            S_colors, max_dis)
            if sol[j]['Marked'] == False:
                flag = flag and sol[j]['Marked']
            else:
                sigma, old_sigma, y = change_sigma(sigma, old_sigma, S, G, sol[j]['list_of_couples'],
                                                   number_of_HT_under_planted,S_labels_table,G_labels_table)
                S_colors = tree_operations.color_tree(S, 'S', S_colors, colors, sigma)
                G_internal_colors = tree_operations.color_tree(G, 'G', G_internal_colors, colors, sigma)
            j += 1
            if not flag:
                if not on_lab:
                    draw.draw_S_and_G(S, G, old_sigma, colors, sigma, path, None, '_rand')
                old_colors = return_color_to_taxon(S, colors)
                save_data(old_sigma, old_colors, sol, noise, 0, compare,path)
                if not running_time:
                    quit()
            print('Planted vertices:%s' % str(sol))
            if not on_lab:
                draw.draw_S_and_G(S, G, old_sigma, colors, sigma, path, sol, '_rand' + str(noise) + '.' + str(0))
            old_colors = return_color_to_taxon(S, colors)
            save_data(old_sigma, old_colors, sol, noise, 0,compare,path)
            return_planted_nodes_new_name(sol,G,path)
        p = Pool(15)
        parameters = [(noise_level[i],number_of_HT_under_planted,G_internal_colors,S_colors,nCr_lookup_table,fact_lookup_table,number_of_leaves) for i in range(0,len(noise_level))]
        p.map(create_tree_for_HT_and_colors_noise, parameters)
        p.map(create_tree_for_color_noise, parameters)
        p.map(create_tree_for_HT_noise, parameters)


    else:
        save_data(old_sigma, old_colors, sol, noise, 0, compare, path)
    print('Running time: %s' % str(datetime.now() - starting_time))

if __name__ == "__main__":
    main(S,G,number_of_leaves,path,k,running_time,number_of_planted_vertices)
