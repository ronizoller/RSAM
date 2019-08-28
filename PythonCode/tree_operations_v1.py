import networkx as nx
import random
import utiles
import math
import dendropy as tr

def isolated(x):
    return x.taxon is None

def is_a_leaf (u):
    return u.child_nodes() == []

def collapse_edges(tree):
    if tree:
        for nd in tree.postorder_node_iter():
            if len(nd.child_nodes()) == 1:
                if tree.seed_node == nd:
                    tree.reroot_at_node(nd.child_nodes()[0])
                else:
                    temp_parent = nd.parent_node
                    temp_child = nd.child_nodes()[0]    #its only child
                    temp_parent.remove_child(nd)
                    temp_parent.add_child(temp_child)
    return tree


def has_right_child(x):
    return len(x.adjacent_nodes()) >= 1


def has_left_child(x):
    return len(x.adjacent_nodes()) >= 2


def is_not_ancestor(nd,x):
    if nd == x:
        return False
    else:
        ans = True
        for y in nd.ancestor_iter():
            if y == x:
                ans = ans and False
            else:
                ans = ans and True
        return ans


def is_ancestor(nd,x):
    return not is_not_ancestor(nd,x)


def color_tree(tree, tree_name, tree_internal_colors, colors, sigma):
    for u in tree.postorder_node_iter():
        tree_internal_colors.update({u.label: [0, 0]})
    for u in tree.postorder_node_iter():
        if u.is_leaf():
            if tree_name == 'S':
                temp_color = colors[u.label]
            elif tree_name == 'G':
                temp_color = colors[sigma[u.label]]
            if temp_color == 'red':
                tree_internal_colors.update({u.label : [1,0]})
            elif temp_color == 'black':
                tree_internal_colors.update({u.label: [0,1]})
            elif temp_color == 'pink':
                tree_internal_colors.update({u.label: [0,0]})
        else:
            reds = tree_internal_colors[u.label][0]
            blacks = tree_internal_colors[u.label][1]
            child = u.child_nodes()
            i = 0
            while i < len(child):
                reds += tree_internal_colors[child[i].label][0]
                blacks += tree_internal_colors[child[i].label][1]
                i += 1
                tree_internal_colors.update({u.label: [reds, blacks]})
    return tree_internal_colors


def weight_G_based_on_same_color_HT (G, new_G, interesting_vertices_p1,interesting_vertices_p2,max_distance,p1,p2,copmare_solutions):
    index = 1

    for u in G.postorder_node_iter():
        new_G.add_node(index, label=u.label, p1=0, p2=0, ind=index)
        if not is_a_leaf(u):
            child = u.child_nodes()
            if has_left_child(u) and has_right_child(u):
                right_child_in_new_G = new_G.nodes()[list(G.postorder_node_iter()).index(child[1]) + 1]
                left_child_in_new_G = new_G.nodes()[list(G.postorder_node_iter()).index(child[0]) + 1]
                new_weight_p1 = edge_weight_based_on_interesting_vertices(u,right_child_in_new_G,left_child_in_new_G, interesting_vertices_p1,max_distance,p1,copmare_solutions,'p1')
                new_weight_p2 = edge_weight_based_on_interesting_vertices(u,right_child_in_new_G,left_child_in_new_G, interesting_vertices_p2,max_distance,p2,copmare_solutions,'p2')

                new_G.nodes[index]['p1'] += new_weight_p1
                new_G.nodes[index]['p2'] += new_weight_p2

                new_G.add_edge(index, list(G.postorder_node_iter()).index(child[0]) + 1, weight_p1 = left_child_in_new_G['p1'],weight_p2 = left_child_in_new_G['p2'])
                new_G.add_edge(index, list(G.postorder_node_iter()).index(child[1]) + 1, weight_p1 = right_child_in_new_G['p1'],weight_p2 = right_child_in_new_G['p2'])
            elif has_right_child(u):
                right_child_in_new_G = new_G.nodes[list(G.postorder_node_iter()).index(child[1]) + 1]
                new_weight_p1 = edge_weight_based_on_interesting_vertices(u, right_child_in_new_G, None, interesting_vertices_p1,max_distance,p1,copmare_solutions,'p1')
                new_weight_p2 = edge_weight_based_on_interesting_vertices(u, right_child_in_new_G, None, interesting_vertices_p2,max_distance,p2,copmare_solutions,'p2')

                new_G.nodes[index]['p1'] += new_weight_p1
                new_G.nodes[index]['p2'] += new_weight_p2
                new_G.add_edge(index, list(G.postorder_node_iter()).index(child[1])+1,weight_p1 = right_child_in_new_G['p1'],weight_p2 = right_child_in_new_G['p2'])

            else:
                left_child_in_new_G = new_G.nodes[list(G.postorder_node_iter()).index(child[0]) + 1]
                new_weight_p1 = edge_weight_based_on_interesting_vertices(u, left_child_in_new_G, None, interesting_vertices_p1,max_distance,p1,copmare_solutions,'p1')
                new_weight_p2 = edge_weight_based_on_interesting_vertices(u, left_child_in_new_G, None, interesting_vertices_p2,max_distance,p2,copmare_solutions,'p2')

                new_G.nodes[index]['p1'] += new_weight_p1
                new_G.nodes[index]['p2'] += new_weight_p2
                new_G.add_edge(index, list(G.postorder_node_iter()).index(child[0]) + 1,
                               weight_p1=left_child_in_new_G['p1'],weight_p2 = left_child_in_new_G['p2'])
        index += 1
    return new_G


def edge_weight_based_on_interesting_vertices(x, y, z, interesting_vertices_p,max_distance, p,compare_solutions,name):
    res_p = 0
    if p[0] is not None:
        if 'HT' in p[0]:
            for nd in interesting_vertices_p:
                if nd['curr']['event'] == 'HT':
                    if (nd['curr']['s'] == x.label and nd['HT_to_in_G']['s'] == y['label']) or (nd['curr']['s'] == x.label and nd['HT_to_in_G']['s'] == z['label']):
                        if p[2] or p[2] != 'False':
                            dis_factor = interesting_vertices_p[interesting_vertices_p.index(nd)]['distance']/max_distance
                        else :
                            dis_factor = 1
                        res_p += interesting_vertices_p[interesting_vertices_p.index(nd)]['probability'] * dis_factor
        if 'D' in p[0]:
            for nd in interesting_vertices_p:
                if nd['curr']['event'] == 'D':
                    if nd['curr']['s'] == x.label:
                        if not compare_solutions:
                            res_p += interesting_vertices_p[interesting_vertices_p.index(nd)]['probability']
                        else :
                            res_p += 1
        if 'S' in p[0]:
            for nd in interesting_vertices_p:
                if nd['curr']['event'] == 'S':
                    if nd['curr']['s'] == x.label:
                        if not compare_solutions:
                            res_p += interesting_vertices_p[interesting_vertices_p.index(nd)]['probability']
                        else :
                            res_p += 1
        if z:
            res_p += z[name]
        if y:
            res_p += y[name]
    return res_p


def random_son(t,u):
    out = t.out_edges(u[1]['ind'], data=True)
    out = [e for e in out]
    son = random.choice(out)
    return t.nodes(data=True)[son[1]]


def number_of_edges_in_subtree(G):
    for nd in (reversed(list(nx.topological_sort(G)))):
        out = G.out_edges([nd], data=True)
        out = [e for e in out]
        nd = G.nodes(data=True)[nd]
        i = 0
        nd.update({'edges_in_subtree' : 0})
        while i < len(out):
            nd['edges_in_subtree'] += 1 + G.nodes(data = True)[out[i][1]]['edges_in_subtree']
            i += 1
    return G


def number_of_possible_events_in_subtree(G):
    for nd in (reversed(list(nx.topological_sort(G)))):
        out = G.out_edges([nd], data=True)
        out = [e for e in out]
        nd = G.nodes(data=True)[nd]
        i = 0
        nd.update({'possible_events' : 0})
        while i < len(out):
            nd['possible_events'] += G.nodes(data = True)[out[i][1]]['possible_events']
            i += 1
        if len(out) > 0:
            nd['possible_events'] += 1
    return G


def find_max_d_of_HT(dis, interesting_vertices,pattern):
    if 'HT' in pattern[0]:
        maxi = 0
        for HT in interesting_vertices:
            if 'HT_to_in_G' in HT:
                HT_curr_dis = dis[(HT['curr']['t'],HT['HT_to_in_G']['t'])]
                if HT_curr_dis > maxi:
                    maxi = HT_curr_dis
        return maxi
    else:
        return


def number_of_leafs(tree, name):
    counter = 0
    for u in tree.postorder_node_iter():
        if u.is_leaf():
            counter += 1
    return counter


def remove_unsigma_genes(G,sigma,taxon):
    to_remove = []
    for u in G.postorder_node_iter():
        if u.taxon:
            if taxon:
                if u.taxon.label not in sigma:
                    to_remove.append(u.taxon.label.replace(' ','_'))
                    to_remove.append(u.taxon.label)
            else:
                if u.label not in sigma:
                    to_remove.append(u.taxon.label)
    return to_remove


def leaf_in_subtrees (tree, tree_name,u, old_sigma,gen):
    if not gen:
        u = find_node_in_tree(tree, u)
    right_subtree_leafs = []
    left_subtree_leafs = []
    u_child = u.child_nodes()
    if len(u_child) > 1:
        right_subtree_leafs_nodes = u_child[1].leaf_nodes()
        for leaf in right_subtree_leafs_nodes:
            if tree_name == 'G':
                if not gen:
                    right_subtree_leafs.append(old_sigma[leaf.taxon.label])
                else:
                    right_subtree_leafs.append(leaf.taxon.label)
            else:
                right_subtree_leafs.append(leaf.taxon.label)
    if len(u_child) > 0:
        left_subtree_leafs_nodes = u_child[0].leaf_nodes()
        for leaf in left_subtree_leafs_nodes:
            if tree_name == 'G':
                if not gen:
                    left_subtree_leafs.append(old_sigma[leaf.taxon.label])
                else:
                    left_subtree_leafs.append(leaf.taxon.label)
            else:
                left_subtree_leafs.append(leaf.taxon.label)

    return right_subtree_leafs, left_subtree_leafs


def max_dis(S_dis_matrix):
    max = 0
    for couple,dis in S_dis_matrix.items():
        if max < dis:
            max = dis
    return max


def find_node_in_tree(tree,nd):
    for u in tree.postorder_node_iter():
        if u.label == nd:
            return u


def find_node_in_networkx_tree(tree,label):
    for u in (list(nx.topological_sort(tree))):
        u = tree.nodes(data=True)[u]
        if u['label'] == label:
            return u
    return None


def normlize_weights(G, k, p, name, field):
    for nd in (reversed(list(nx.topological_sort(G)))):
        if p[0] is not None:
            if (G.nodes(data = True)[nd][field] * k) > 0:
                G.nodes(data = True)[nd][name] = G.nodes(data = True)[nd][name]/ (len(p[0]) * G.nodes(data = True)[nd][field] * k )
    return G


def find_max_scores(G, number_of_planted_vertices, name, TH):
    max_score_p = [-1]*number_of_planted_vertices
    max_score = max_score_p
    if name == 'p1':
        for nd in (reversed(list(nx.topological_sort(G)))):
            if G.nodes(data=True)[nd]['edges_in_subtree'] >= TH:
                max_score = utiles.update_top_ranking_list(G.nodes[nd][name],max_score_p)
    else:
        for nd in (reversed(list(nx.topological_sort(G)))):
            out_going = G.out_edges([nd], data=True)
            out_going = [e for e in out_going]
            if out_going != []:
                child1 = G.nodes(data=True)[out_going[0][1]]
                child2 = G.nodes(data=True)[out_going[1][1]]
                if child1['edges_in_subtree'] >= TH and child2['edges_in_subtree'] >= TH:
                    if child1['edges_in_subtree'] > child2['edges_in_subtree']:
                        max_score = utiles.update_top_ranking_list(child1['p1']+child2['p2'], max_score_p)
                    if child2['edges_in_subtree'] > child1['edges_in_subtree']:
                        max_score = utiles.update_top_ranking_list(child1['p2']+child2['p1'], max_score)
    return max_score


def copy_G(G,new_G):
    index = 1
    for u in G.postorder_node_iter():
        new_G.add_node(index, label=u.label,ind = index)      #same_HT_score[0] = reds score, same_HT_score[1] = blacks score
        if not is_a_leaf(u):
            child = u.child_nodes()
            if has_left_child(u) and has_right_child(u):
                new_G.add_edge(index, list(G.postorder_node_iter()).index(child[0]) + 1)
                new_G.add_edge(index, list(G.postorder_node_iter()).index(child[1]) + 1)
            elif has_right_child(u):
                new_G.add_edge(index, list(G.postorder_node_iter()).index(child[1])+1)
            else :
                new_G.add_edge(index, list(G.postorder_node_iter()).index(child[0]) + 1)
        index += 1
    return new_G


def reroot_and_save(tree,nd,path,ext):
    new_root = find_node_in_tree(tree,nd)
    tree.reroot_at_node(new_root)
    tree.write_to_path(path+'/rerooted_tree_'+ext, schema="newick")
    quit()