import sys
sys.path.append('/anaconda3/lib/python3.6/site-packages')
import networkx as nx
import tree_operations
from networkx.drawing.nx_agraph import graphviz_layout
import matplotlib.pyplot as plt
import inits
import math
import utiles
import seaborn as sns

def draw_garph(H, S, G, deleted_nodes, big_size, small_size):
    print('Drawing hypergraph...')
    labels1 = nx.get_node_attributes(H, 's')
    labels2 = nx.get_node_attributes(H, 't')
    labels3 = nx.get_node_attributes(H, 'l')
    labels4 = nx.get_edge_attributes(H,'source')
    labels5 = nx.get_edge_attributes(H,'target')
    #labels6 = nx.get_node_attributes(H, 'red_colors')
    #labels7 = nx.get_node_attributes(H, 'black_colors')

    new_labels3 = {}
    for k,match in labels3.items():
        l = ""
        for to_seperate in match:
            l = l+str(to_seperate)+'\n'
        new_labels3.update({k:l})
    for k,l in labels2.items():
        for x in S.postorder_node_iter():
            if x.taxon != None:
                if x.label == l:
                    l = l+" ("+str(x.taxon)+")"
                    labels2.update({k:l})

    for k,l in labels1.items():
        for x in G.postorder_node_iter():
            if x.taxon != None:
                if x.label == l:
                    l = l+" ("+str(x.taxon)+")"
                    labels1.update({k:l})

    nodes_labels = {}
    for i in range(0,len(labels1)+len(deleted_nodes)):
        if not i in deleted_nodes:
            nodes_labels.update({i:("node_number: "+str(i)+"\ns: "+str(labels1[i])+"\nt: "+str(labels2[i])+"\n"+str(new_labels3[i])+"\nred color:")})
    edges_labels = {}
    for e in H.edges(keys=True):
        edges_labels.update({e:("source: "+str(labels4[e])+"\ntarget: "+str(labels5[e]))})
    sizes = []
    for nd in H.node(data=True):
        if((nd[1]['s']==G.seed_node.label and nd[1]['t']==S.seed_node.label)):
            sizes = sizes+[big_size]
        else:
            sizes = sizes+[small_size]

    pos = graphviz_layout(H, prog='dot', args="-Grankdir=BT")
    plt.figure(12,figsize=(20,20))                              #size of fig
    nx.draw(H, pos, arrows=True,node_size=sizes)
    nx.draw_networkx_labels(H, pos,nodes_labels, font_size=7)
    nx.draw_networkx_edges(H, pos, alpha = 0.1, width = 2, edge_color='b')
    nx.draw_networkx_edge_labels(H, pos,font_size=6,labels = edges_labels)
    plt.savefig('hyperGraph.png')
    print('Finished drawing hypergraph.\n')

def draw_new_G(G,G_nodes_identified,colors,sigma,new_G):
    print('Drawing new G...')
    plt.clf()
    tree_to_draw = nx.DiGraph()
    index = 1
    for u in G.postorder_node_iter():
        tree_to_draw.add_node(index, label=u.label)
        if not tree_operations.is_a_leaf(u):
            child = u.child_nodes()
            i = 0
            while i < len(child):
                if len(child) >= i + 1:
                    tree_to_draw.add_edge(index, list(G.postorder_node_iter()).index(child[i]) + 1)
                i = i + 1
        index += 1
    labels1 = nx.get_node_attributes(new_G, 'label')
    pos1 = graphviz_layout(tree_to_draw, prog='dot')

    plt.figure(12, figsize=(40, 40))  # size of fig

    print("G_nodes_identified = %s" % str(G_nodes_identified))
    nodes_color = []
    nodes_size = []

    for nd in new_G.nodes(data = True):
        if new_G.out_degree(nd[0]) == 0 and not new_G.in_degree(nd[0]) == 0:
            if colors[sigma[nd[1]['label']]] == 'red':
                nodes_color.append('red')
            else:
                nodes_color.append('black')
            nodes_size.append(200)
        elif G_nodes_identified[nd[1]['label']] > 0:
            nodes_color.append('blue')
            nodes_size.append(G_nodes_identified[nd[1]['label']] * 350)
        else:
            nodes_color.append('white')
            nodes_size.append(200)


    for r,l in labels1.items():
        for x in G.postorder_node_iter():
            if x.taxon != None:
                if x.label == l:
                    l = l+"\n ("+str(x.taxon)+")"
                    labels1.update({r:l})

    #edges_width = []
    edges_color = []

    for e in new_G.edges(data=True):
        if e[2]['weight'][0] == 0 and e[2]['weight'][1] == 0:
            edges_color.append('grey')
            #edges_width.append(1)
        elif e[2]['weight'][0] > e[2]['weight'][1]: #more red HT
            #edges_width.append(e[2]['weight'][0]*10)
            edges_color.append('red')
        else: #more black HT
            #edges_width.append(e[2]['weight'][1]*10)
            edges_color.append('black')
    #edges_width = utile.normlize (edges_width,20)

    print("len(nodes_size) = %s, len(nodes_color) = %s" % (len(nodes_size),len(nodes_color)))
    nx.draw(tree_to_draw, pos1, arrows=True,node_size = nodes_size, node_color=nodes_color,edge_color = edges_color, width = 1)
    nx.draw_networkx_labels(tree_to_draw, pos1, labels1, font_size=7,)

    #nx.draw_networkx_edges(tree_to_draw, pos1, )
    #nx.draw_networkx_edge_labels(tree_to_draw, pos1, font_size=6, labels=lables2)
    plt.savefig('new_G.png')
    print('Finished drawing new G.\n')

def draw_new_G2(marked_nodes, colors, sigma, new_G, G,old_sigma,k,TH_compare_subtrees, TH_both, TH_pattern_in_subtree, path, both, alpha, lables, glob,spec,pattern,size,evol,compare,number_of_fields):
    print('Drawing new G...')
    plt.clf()
    if glob:
        glob_text = '_global'
    else:
        glob_text = '_local'

    tree_to_draw = nx.DiGraph()
    index = 1
    for u in G.postorder_node_iter():
        tree_to_draw.add_node(index, label = u.label)
        if not tree_operations.is_a_leaf(u):
            child = u.child_nodes()
            i = 0
            while i < len(child):
                if len(child) >= i + 1:
                    tree_to_draw.add_edge(index, list(G.postorder_node_iter()).index(child[i]) + 1)
                i = i + 1
        index += 1
    labels1 = nx.get_node_attributes(new_G, 'label')
    pos1 = graphviz_layout(tree_to_draw, prog='dot')

    plt.figure(figsize=(80, 40))

    nodes_color = []
    nodes_size = []
    max_size = utiles.map_max(marked_nodes,number_of_fields)
    exp_facor = math.log(size)

    for nd in new_G.nodes(data=True):
        if new_G.out_degree(nd[0]) == 0 and not new_G.in_degree(nd[0]) == 0:
            if colors[sigma[nd[1]['label']]] == 'red':
                nodes_color.append('red')
            else:
                nodes_color.append('grey')
            nodes_size.append(200)
        elif nd[1]['label'] in marked_nodes:
            nodes_color.append('blue')
            temp_size = max(marked_nodes[nd[1]['label']][0][0],marked_nodes[nd[1]['label']][0][1])
            nodes_size.append(math.exp((temp_size/max_size)*exp_facor)+500)
        else:
            nodes_color.append('white')
            nodes_size.append(200)
    if lables:
        for r, l in labels1.items():
            for x in G.postorder_node_iter():
                if x.taxon != None:
                    if x.label == l:
                        l = l + "\n (" + str(x.taxon) + ")"
                        l = l+'\n'+str(old_sigma[x.taxon.label])
                        labels1.update({r: l})

    print('nodes_sizes = '+str(nodes_size))
    nx.draw(tree_to_draw, pos1, arrows=True, node_size=nodes_size, node_color=nodes_color,
            width=1)
    nx.draw_networkx_labels(tree_to_draw, pos1, labels1, font_size=7)
    plt.savefig(path+'/figures/new_G_k='+str(k)+'_TH_compare_subtrees='+str(TH_compare_subtrees)+'_TH_pattern_in_subtree='+str(TH_pattern_in_subtree)+'_species='+spec+'_pattern='+pattern+"_"+evol+"_comapre="+str(compare)+'.png')
    print('Finished drawing new G.\n')

def draw_tree(tree, name, old_sigma, colors, sigma):
    tree_to_draw = nx.DiGraph()
    index = 1
    nodes_color = []
    for u in tree.postorder_node_iter():
        tree_to_draw.add_node(index, label=u.label)
        if not tree_operations.is_a_leaf(u):
            child = u.child_nodes()
            i = 0
            while i < len(child):
                if len(child) >= i+1:
                    tree_to_draw.add_edge(index, list(tree.postorder_node_iter()).index(child[i]) + 1)
                i = i+1
        index += 1
    labels = nx.get_node_attributes(tree_to_draw, 'label')
    for k,l in labels.items():
        for x in tree.postorder_node_iter():
            if x.taxon != None:
                if x.label == l:
                    l = l+"\n ("+str(x.taxon)+")"
                    if name == 'G':
                        l = l+'\n'+str(old_sigma[x.taxon.label])
                        labels.update({k:l})
                    else:
                        labels.update({k: l})

    if name == 'S' :
        for u in tree.postorder_node_iter():
            if u.label != None and u.label in colors:
                if colors[u.label] == 'red':
                    nodes_color.append('red')
                elif colors[u.label] == 'black':
                    nodes_color.append('grey')
                else :
                    nodes_color.append('pink')
            else:
                nodes_color.append('white')
    else:
        for u in tree.postorder_node_iter():
            if tree_operations.is_a_leaf(u) and not tree_operations.isolated(u) and sigma[u.label] in colors:
                if colors[sigma[u.label]] == 'red':
                    nodes_color.append('red')
                else:
                    nodes_color.append('grey')
            else:
                nodes_color.append('white')

    postree = graphviz_layout(tree_to_draw, prog='dot')

    plt.figure(12, figsize=(80, 40))  # size of fig

    nx.draw(tree_to_draw, postree, arrows=True,node_color=nodes_color)
    nx.draw_networkx_labels(tree_to_draw, postree, labels, font_size=7,)
    nx.draw_networkx_edges(tree_to_draw, postree)
    plt.savefig(name+'.png')
    print('Drawing'+name)

def draw_S_and_G(S,G, old_sigma, colors, sigma,path,sol,ext):
    print('Sigma: %s\n old_sigma: %s\n' % (str(sigma),str(old_sigma)))
    plt.clf()
    S_to_draw = nx.DiGraph()
    G_to_draw = nx.DiGraph()
    index = 1
    nodes_color_S = []
    nodes_color_G = []

    ##FOR S
    for u in S.postorder_node_iter():
        S_to_draw.add_node(index, label=u.label)
        if not tree_operations.is_a_leaf(u):
            child = u.child_nodes()
            i = 0
            while i < len(child):
                if len(child) >= i+1:
                    S_to_draw.add_edge(index, list(S.postorder_node_iter()).index(child[i]) + 1)
                i = i+1
        index += 1
    labels_S = nx.get_node_attributes(S_to_draw, 'label')
    for k,l in labels_S.items():
        for x in S.postorder_node_iter():
            if x.taxon != None:
                if x.label == l:
                    l = l+"\n ("+str(x.taxon)+")"
                    labels_S.update({k: l})

    for u in S.postorder_node_iter():
            if u.label != None and u.label in colors:
                if colors[u.label] == 'red':
                    nodes_color_S.append('red')
                elif colors[u.label] == 'black':
                    nodes_color_S.append('grey')
                else :
                    nodes_color_S.append('pink')
            else:
                nodes_color_S.append('white')
    ## FOR G
    index = 1
    for u in G.postorder_node_iter():
        G_to_draw.add_node(index, label=u.label)
        if not tree_operations.is_a_leaf(u):
            child = u.child_nodes()
            i = 0
            while i < len(child):
                if len(child) >= i+1:
                    G_to_draw.add_edge(index, list(G.postorder_node_iter()).index(child[i]) + 1)
                i = i+1
        index += 1
    labels_G = nx.get_node_attributes(G_to_draw, 'label')
    for k,l in labels_G.items():
        for x in G.postorder_node_iter():
            if x.taxon != None:
                if x.label == l and (x.taxon.label in old_sigma):
                    l = l+"\n ("+str(x.taxon)+")"
                    l = l+'\n'+str(old_sigma[x.taxon.label])
                    labels_G.update({k:l})

    for u in G.postorder_node_iter():
        degel = False
        if sol != None:
            for int, temp_sol in sol.items():
                for p in range (0,len(temp_sol['list_of_couples'])):
                    if (u.label == temp_sol['Marked'] or u.label == temp_sol['list_of_couples'][p][0] or u.label == temp_sol['list_of_couples'][p][1]) and not degel:
                        nodes_color_G.append('blue')
                        degel = True
        if not degel and tree_operations.is_a_leaf(u) and not tree_operations.isolated(u) and sigma[u.label] in colors:
            if colors[sigma[u.label]] == 'red':
                nodes_color_G.append('red')
            elif colors[sigma[u.label]] == 'black':
                nodes_color_G.append('grey')
            else:
                nodes_color_G.append('pink')
        elif not degel:
            nodes_color_G.append('white')

    postree_S = graphviz_layout(S_to_draw, prog='dot')
    postree_G = graphviz_layout(G_to_draw, prog='dot')
    for k, v in postree_G.items():
        # Shift the x values of every node by 10 to the right
        lst = list(v)
        lst[0] = lst[0] + 100
        postree_G.update({k:tuple(lst)})

    fig, axes = plt.subplots(1,2,figsize=(200, 20))
    ax = axes.flatten()

    ax[0].set_title('Species tree',fontsize=50)
    ax[1].set_title('Gene tree',fontsize=50)

    nx.draw(S_to_draw, postree_S, arrows=True,node_color=nodes_color_S,ax=ax[0])
    nx.draw(G_to_draw, postree_G, arrows=True,node_color=nodes_color_G,ax=ax[1])

    nx.draw_networkx_labels(S_to_draw, postree_S, labels_S, font_size=7,ax=ax[0])
    nx.draw_networkx_labels(G_to_draw, postree_G, labels_G, font_size=7,ax=ax[1])

    nx.draw_networkx_edges(S_to_draw, postree_S,ax=ax[0])
    nx.draw_networkx_edges(G_to_draw, postree_G,ax=ax[1])
    plt.show()
    #plt.savefig(path + '/figures/S+G'+ext+'.png')
    print('Drawing S and G')

def draw_G_diffrent_optimal_solutions(marked_nodes, colors, sigma, old_sigma, new_G, G, k, path, both, alpha, labels, TH_compare_subtrees, TH_both, TH_pattern_in_subtree,compare_subtrees,evolutinary_event,pattern,iterations, factor,size):
    print('Drawing new G...')
    plt.clf()
    plt.figure(figsize=(80, 40))
    tree_to_draw = nx.DiGraph()
    index = 1
    for u in G.postorder_node_iter():
        tree_to_draw.add_node(index, label=u.label)
        if not tree_operations.is_a_leaf(u):
            child = u.child_nodes()
            i = 0
            while i < len(child):
                if len(child) >= i + 1:
                    tree_to_draw.add_edge(index, list(G.postorder_node_iter()).index(child[i]) + 1)
                i = i + 1
        index += 1
    labels1 = nx.get_node_attributes(new_G[0], 'label')
    pos1 = graphviz_layout(tree_to_draw, prog='dot')


    marked_counter = inits.init_dic(G.nodes(),0)
    max_counts = 0
    i = 1

    for solution in marked_nodes:
        print('solution number '+str(i)+":")
        for item in solution:
            marked_node = item[0]
            scoring = item[1]

            temp_high_score = 0
            for j in range(0,2):
                for m in range(0,2):
                    if temp_high_score < scoring[j][m]:
                        temp_high_score = scoring[j][m]
            print('         %s with score: %s' % (str(marked_node),str(temp_high_score)))
            marked_counter[marked_node] += temp_high_score
            if marked_counter[marked_node] > max_counts:
                max_counts = marked_counter[marked_node]
        i += factor

    for u,counter in marked_counter.items():
        if max_counts>0:
            marked_counter.update({u:(counter/max_counts)*size})
    #print('max_counter = %s, marked_counter = %s' % (str(max_counts),str(marked_counter)))

    nodes_color = []
    nodes_size = []

    for nd in new_G[0].nodes(data=True):
        if new_G[0].out_degree(nd[0]) == 0 and not new_G[0].in_degree(nd[0]) == 0:
            if colors[sigma[nd[1]['label']]] == 'red':
                nodes_color.append('red')
            else:
                nodes_color.append('grey')
            nodes_size.append(200)
        elif marked_counter[nd[1]['label']] > 0:
            nodes_color.append('blue')
            nodes_size.append(marked_counter[nd[1]['label']])
        else:
            nodes_color.append('white')
            nodes_size.append(200)
    if labels:
        for r, l in labels1.items():
            for x in G.postorder_node_iter():
                if x.taxon != None:
                    if x.label == l:
                        l = l + "\n (" + str(x.taxon) + ")"
                        l = l + '\n' + str(old_sigma[x.taxon.label])
                        labels1.update({r: l})

    nx.draw(tree_to_draw, pos1, arrows=True, node_size=nodes_size, node_color=nodes_color,
            width=1)
    nx.draw_networkx_labels(tree_to_draw, pos1, labels1, font_size=10)
    plt.savefig(path + '/figures/G_different_optimal. k=' + str(k) + '_TH_compare_subtrees = ' + str(TH_compare_subtrees) + '_TH_pattern_in_subtree = ' + str(TH_pattern_in_subtree) +"_pattern="+pattern+"_"+evolutinary_event+"compare_subtrees="+str(compare_subtrees)+'.png')
    print('Finished drawing new G.\n')

def connectpoints(x,y,p1,p2):
    x1, x2 = x[p1], x[p2]
    y1, y2 = y[p1], y[p2]
    plt.plot([x1,x2],[y1,y2],'k-')

def draw_plot(all_vertices_with_noise,path,marked_vertex):
    print('Drawing plot..')
    plt.clf()
    plt.figure(12, figsize=(20, 20))  # size of fig
    list_to_draw_reds = []
    list_to_draw_blacks = []
    length = []
    names = []
    to_connect = []
    length_to_connect = []
    max = 0
    max_noise_level = 0
    for noise_level, all_vertices in all_vertices_with_noise.items():
        for u, couple in all_vertices.items():
            names.append(u)
            length.append(noise_level)
            list_to_draw_reds.append(couple[0])
            list_to_draw_blacks.append(couple[1])
            if max < couple[0]:
                max = couple[0]
            if max < couple[1]:
                max = couple[1]
            if max_noise_level < noise_level:
                max_noise_level = noise_level
            if u in marked_vertex:
                to_connect.append(couple[1])
                length_to_connect.append(noise_level)
    fig, ax = plt.subplots()
    ax.plot(length, list_to_draw_reds, 'ro', length, list_to_draw_blacks, 'ro')
    for X, Y, Z in zip(length, list_to_draw_reds, names):
        # Annotate the points 5 _points_ above and to the left of the vertex
        ax.annotate('{}'.format(Z), xy=(X, Y), xytext=(5, -5), ha='left',
                    textcoords='offset points')
    for X, Y, Z in zip(length, list_to_draw_blacks, names):
        # Annotate the points 5 _points_ above and to the left of the vertex
        ax.annotate('{}'.format(Z), xy=(X, Y), xytext=(5, -5), ha='left',
                    textcoords='offset points')
    sns.set(font_scale=1.4)
    plt.xlabel('Noise', fontsize=10)
    plt.ylabel('Score', fontsize=10)
    for p in range(0, len(length_to_connect) - 1):
        connectpoints(length_to_connect, to_connect, p, p + 1)
    plt.axis([0, max_noise_level + 5, 0, max + 0.0001])
    plt.savefig(path + '/figures/plot_noise.png')

def draw_compare_k_plot(all_vertices_with_k,path):
    print('Drawing plot..')
    plt.clf()
    sns.set(font_scale=1.6)
    plt.figure(figsize=(25, 5))  # size of fig
    l = utiles.compare_dict_entries(all_vertices_with_k)
    print(l)
    to_display = []
    xtick = []
    ytick = []
    first = True
    for k, all_vertices_for_k in l.items():
        ytick.append(k)
        to_display_temp = []
        for vertex, couple in all_vertices_for_k.items():
            to_display_temp.append(max(couple))
            if first:
                xtick.append(vertex)
        first = False
        to_display = to_display + [to_display_temp]
    print(to_display)
    ax = sns.heatmap(utiles.flip_list(to_display,len(xtick)), yticklabels=xtick, xticklabels=ytick, linewidth=0.5)
    plt.yticks(rotation=0)
    plt.xticks(rotation=0)
    fig = ax.get_figure()
    fig.savefig(path + '/figures/plot_noise.png')
