


k = 50
times_effi = {}
times_naive = {}
number_of_marked_vertices = 1
speciesTreespecification = 'all'
HT_cost = 1
D_cost = 1
S_cost = 0
on_lab = True
number_of_iter_per_size = 10

if on_lab:
    path = '/users/studs/bsc/2016/ronizo/PycharmProjects/RSAM/simulator_data/running_time'
else:
    import sys
    sys.path.append('/anaconda3/lib/python3.6/site-packages')
    path = '/Users/ronizoller/Documents/school/Master/מחקר/DATA/running_times'

import utiles
import Tree_Grnarator
from datetime import datetime
import os
import EfficiantVersion as effi
import tree_operations_v1 as tree_operations
import inits_v1 as inits
import dendropy as tr
import hypergraph_v1 as hypergraph
import matplotlib.pyplot as plt


for number_of_leafs in utiles.frange(100,500,50):
    temp_time_effi = {}
    temp_time_naive = {}
    for d in range(0,number_of_iter_per_size):
        path_curr = path + '/number_of_leaves:'+str(number_of_leafs)+'/'
        os.makedirs(os.path.dirname(path_curr), exist_ok=True)
        Tree_Grnarator.main(number_of_leafs,path_curr,k,True,1)

        G = tr.Tree.get_from_path(path_curr + "/GeneTree(binary)_local.txt", schema="newick")
        S = tr.Tree.get_from_path(path_curr + "/phyliptree(binary," + speciesTreespecification + ").phy", schema="newick")

        print("     Reading file " + path_curr + "/sigma.txt'...")
        input = open(path_curr + '/0/sigma0.0.txt', 'r')
        sigma = []
        for line in input:
            sigma.append(eval(line))
        sigma = sigma[0]
        print("     Finished reading file  " + path_curr + "0/sigma0.0.txt'")
        G.prune_taxa_with_labels(tree_operations.remove_unsigma_genes(G, sigma, True))

        S = utiles.init_internal_labels(S, 'x', sigma, path_curr)
        G = utiles.init_internal_labels(G, 'u', sigma, path_curr)

        G = tree_operations.collapse_edges(G)
        S = tree_operations.collapse_edges(S)

        S_labels_table, G_labels_table = inits.init_taxon_to_label_table(S, G, sigma)
        sigma, old_sigma = inits.update_sigma(S, G, k, sigma, False, path_curr, True, S_labels_table, G_labels_table)
        G.prune_taxa_with_labels(tree_operations.remove_unsigma_genes(G, sigma, False))
        S_dis_matrix = inits.init_distance_S({}, k, False, path_curr, speciesTreespecification)
        nodes_table = inits.init_nodes_table(S, G, {})

        start = datetime.now()
        effi.build_hyper_garph(S, G, False, k,nodes_table, D_cost, S_cost, HT_cost, path_curr, 1,sigma)
        temp_time_effi.update({d:datetime.now()-start})

        start = datetime.now()
        hypergraph.build_hyper_garph(S, G, False, k,nodes_table, D_cost, S_cost, HT_cost, path_curr, 1, sigma)
        temp_time_naive.update({d:datetime.now()-start})

    times_effi.update({number_of_leafs:utiles.avg_dict_datetime(temp_time_effi)})
    times_naive.update({number_of_leafs:utiles.avg_dict_datetime(temp_time_naive)})

print('Effi times: %s\nNaive Times: %s' % (str(times_effi),str(times_naive)))
file = open(path + 'results_times_effi.txt', 'w')
file.write(str(times_effi))
file.close()
file = open(path + 'results_times_naive.txt', 'w')
file.write(str(times_naive))
file.close()

xs = []
effi = []
naive = []
names = []
for number_of_leafs,time in times_effi.items():
    xs.append(number_of_leafs)
    effi.append(time)
for number_of_leafs,time in times_naive.items():
    naive.append(time)
plt.xlabel('Number of Leafs', fontsize=10)
plt.ylabel('Running Time', fontsize=10)
fig, ax = plt.subplots()
plt.plot(xs, effi, 'ro')
plt.plot(xs, naive, 'g^')
plt.axis([0, max(xs)+50, 0, max(effi+naive)+1])
plt.show()