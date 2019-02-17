import networkx as nx
import random


def isolated(x):
    return x.taxon is None

def is_a_leaf (u):
    return (u.child_nodes() == [])

def collapse_edges(tree):
    print('Collapsing edges...')
    for nd in tree.postorder_node_iter():
        if len(nd.child_nodes()) == 1:
            if tree.seed_node == nd:
                tree.reroot_at_node(nd.child_nodes()[0])
            else:
                temp_parent = nd.parent_node
                temp_child = nd.child_nodes()[0]    #its only child
                temp_parent.remove_child(nd)
                temp_parent.add_child(temp_child)
    print('Finished collapsing edges.\n')
    return tree

def has_right_child(x):
    return len(x.adjacent_nodes()) >= 1

def has_left_child(x):
    return len(x.adjacent_nodes()) >= 2

def is_not_ancestor (nd,x):
    #print('     in is_not_ans, nd = %s, x = %s' % (str(nd),str(x)))
    if (nd == x):
        return False
    else:
        ans = True
        for y in nd.ancestor_iter():
            #print('     y = %s' % str(y))
            if (y == x):
                ans = ans and False
            else:
                ans = ans and True
        return ans

def color_tree(tree, tree_name, tree_internal_colors, colors, sigma):
    #print('Coloring internal vertices of S...')
    for u in tree.postorder_node_iter():
        tree_internal_colors.update({u.label: [0, 0]})

    for u in tree.postorder_node_iter():
        reds = 0
        blacks = 0
        if u.is_leaf():
            if tree_name == 'S':
                temp_color = colors[u.label]
            elif tree_name == 'G':
                temp_color = colors[sigma[u.label]]
            if temp_color == 'red':
                tree_internal_colors.update({u.label : [1,0]})
            elif temp_color == 'black':
                tree_internal_colors.update({u.label: [0, 1]})

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
    #print('Finished coloring internal vertices of S...')
    return tree_internal_colors

def weight_G_based_on_same_color_HT (G, new_G, red_HT_vertices_in_G,black_HT_vertices_in_G,max_distance,distance_flag,evol,compare_solutions,k):
    print('Weighting G...')
    #print ('red_HT_vertices_in_G = %s\nblack_HT_vertices_in_G = %s' % (str(red_HT_vertices_in_G),str(black_HT_vertices_in_G)))
    index = 1
    for u in G.postorder_node_iter():
        new_G.add_node(index, label=u.label,same_HT_score = [0,0],ind = index)      #same_HT_score[0] = reds score, same_HT_score[1] = blacks score
        if not is_a_leaf(u):
            child = u.child_nodes()
            if has_left_child(u) and has_right_child(u):
                right_child_in_new_G = new_G.nodes()[list(G.postorder_node_iter()).index(child[1]) + 1]
                left_child_in_new_G = new_G.nodes()[list(G.postorder_node_iter()).index(child[0]) + 1]
                #print('     u = %s, right_child_in_new_G = %s, left_child_in_new_G = %s' % (str(u), str(right_child_in_new_G), str(left_child_in_new_G)))

                new_weight= edge_weight_based_on_same_color_HT(u,right_child_in_new_G,left_child_in_new_G, red_HT_vertices_in_G,black_HT_vertices_in_G,max_distance,distance_flag,evol,compare_solutions,k)

                i = 0
                for i in [0,1]:
                    new_G.nodes[index]['same_HT_score'][i] += new_weight[i]

                new_G.add_edge(index, list(G.postorder_node_iter()).index(child[0]) + 1, weight = left_child_in_new_G['same_HT_score'])
                new_G.add_edge(index, list(G.postorder_node_iter()).index(child[1]) + 1, weight = right_child_in_new_G['same_HT_score'])

            elif has_right_child(u):
                right_child_in_new_G = new_G.nodes[list(G.postorder_node_iter()).index(child[1]) + 1]
                new_weight = edge_weight_based_on_same_color_HT(u, right_child_in_new_G, None, red_HT_vertices_in_G,black_HT_vertices_in_G,distance_flag,evol,compare_solutions,k)

                i = 0
                for i in [0,1]:
                    new_G.nodes[index]['same_HT_score'][i] += new_weight[i]

                new_G.add_edge(index, list(G.postorder_node_iter()).index(child[1])+1,weight = right_child_in_new_G['same_HT_score'])

            else :
                left_child_in_new_G = new_G.nodes[list(G.postorder_node_iter()).index(child[0]) + 1]
                new_weight = edge_weight_based_on_same_color_HT(u, left_child_in_new_G, None, red_HT_vertices_in_G,black_HT_vertices_in_G,distance_flag,evol,compare_solutions,k)

                for i in [0, 1]:
                    new_G.nodes[index]['same_HT_score'][i] += new_weight[i]

                new_G.add_edge(index, list(G.postorder_node_iter()).index(child[0]) + 1,
                               weight = left_child_in_new_G['same_HT_score'])
        index += 1
    print('Finished weighting G.\n')
    return new_G

def edge_weight_based_on_same_color_HT(x, y, z, red_HT_vertices_in_G, black_HT_vertices_in_G,max_distance, distance_flag,evol,compare_solutions,k):
    res = [0,0]         #res[0] = reds score, res[1] = blacks score
    if evol == 'HT':
        for nd in red_HT_vertices_in_G:
            if (nd['curr']['s'] == x.label and nd['HT_to_in_G']['s'] == y['label']) or (nd['curr']['s'] == x.label and nd['HT_to_in_G']['s'] == z['label']):
                if distance_flag:
                    dis_factor = red_HT_vertices_in_G[red_HT_vertices_in_G.index(nd)]['distance']/max_distance[0]
                else :
                    dis_factor = 1
                if not compare_solutions:
                    res[0] += red_HT_vertices_in_G[red_HT_vertices_in_G.index(nd)]['probability'] * dis_factor
                else:
                    print("dis %s for %s" % (str(dis_factor),str(x)))
                    res[0] += dis_factor
        for nd in black_HT_vertices_in_G:
            if (nd['curr']['s'] == x.label and nd['HT_to_in_G']['s'] == y['label']) or (nd['curr']['s'] == x.label and nd['HT_to_in_G']['s'] == z['label']):
                if distance_flag:
                    dis_factor = black_HT_vertices_in_G[black_HT_vertices_in_G.index(nd)]['distance']/max_distance[1]
                else :
                    dis_factor = 1
                if not compare_solutions:
                    res[1] += black_HT_vertices_in_G[black_HT_vertices_in_G.index(nd)]['probability'] * dis_factor
                else:
                    print("dis %s for %s" % (str(dis_factor), str(x)))
                    res[1] += dis_factor
    elif evol == 'D':
        for nd in red_HT_vertices_in_G:
            if nd['curr']['s'] == x.label:
                if not compare_solutions:
                    res[0] += red_HT_vertices_in_G[red_HT_vertices_in_G.index(nd)]['probability']
                else :
                    res[0] += 1
        for nd in black_HT_vertices_in_G:
            if nd['curr']['s'] == x.label:
                if not compare_solutions:
                    res[1] += black_HT_vertices_in_G[black_HT_vertices_in_G.index(nd)]['probability']
                else:
                    res[1] += 1
    if z != None:
        res[0] += z['same_HT_score'][0]
        res[1] += z['same_HT_score'][1]
    if y != None :
        res[0] += y['same_HT_score'][0]
        res[1] += y['same_HT_score'][1]
    return res

def random_son(t,u):
    out = t.out_edges(u[1]['ind'], data=True)
    out = [e for e in out]
    son = random.choice(out)
    return t.nodes(data=True)[son[1]]


def number_of_edges_in_subtree(G):
    print ('Counting edges in G...')
    for nd in (reversed(list(nx.topological_sort(G)))):
        out = G.out_edges([nd], data=True)
        out = [e for e in out]
        nd = G.nodes(data=True)[nd]
        i = 0
        nd.update({'edges_in_subtree' : 0})
        while i < len(out):
            nd['edges_in_subtree'] += 1 + G.nodes(data = True)[out[i][1]]['edges_in_subtree']
            i += 1
    print ('Finished counting edges in G')
    return G

def find_max_d_of_HT(dis, red_list, black_list,evolutanry_event):
    if evolutanry_event == 'HT':
        max_red = 0
        max_black = 0
        for HT in red_list:
            HT_curr_dis = dis[(HT['curr']['t'],HT['HT_to_in_G']['t'])]
            if HT_curr_dis > max_red:
                max_red = HT_curr_dis
        for HT in black_list:
            HT_curr_dis = dis[(HT['curr']['t'],HT['HT_to_in_G']['t'])]
            if HT_curr_dis > max_black:
                max_black = HT_curr_dis
        return [max_red,max_black]
    else :
        return

def number_of_leafs (tree, name):
    counter = 0
    for u in tree.postorder_node_iter():
        if u.is_leaf():
            counter += 1
    print('number of leafs in %s is %s' % (name,str(counter)))
    return counter

def remove_unsigma_genes(G,sigma,taxon):
    to_remove = []
    for u in G.postorder_node_iter():
        if u.is_leaf():
            if taxon:
                if not u.taxon.label in sigma:
                    to_remove.append(u.taxon.label)
            else:
                if not u.label in sigma:
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
    #print('right subtree = %s' % str(right_subtree_leafs_nodes))
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
        #print('left subtree = %s' % str(left_subtree_leafs_nodes))
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

def extract_list_of_leafs(name,tree,old_sigma):
    unmarked = leaf_in_subtrees(tree, name, old_sigma)
    print('%s leafs: %s\n' % (str(name),str(unmarked)))

def find_node_in_networkx_tree(tree,label):
    for u in (list(nx.topological_sort(tree))):
        u = tree.nodes(data=True)[u]
        if u['label'] == label:
            return u
    return None