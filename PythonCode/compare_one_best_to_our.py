import sys
sys.path.append('/anaconda3/lib/python3.6/site-packages')
import numpy as np
import matplotlib.pyplot as plt
import utiles
on_lab = True
if on_lab:
    path  = '/users/studs/bsc/2016/ronizo/PycharmProjects/RSAM/simulator_data/comparsion'
else:
    path = '/Users/ronizoller/Documents/school/Master/XXX/DATA'
    import sys
    sys.path.append('/anaconda3/lib/python3.6/site-packages')

number_of_one_bests = 10

all_marked = {}
precentage = []                                             #{random_input_index:(false pos,false neg)}
our_precentage = []
our_sol = {}
planted_vertex = []
planted_vertex_input = open(path + '/saved_data/marked_nodes_correct_names.txt', 'r')
for line in planted_vertex_input:
    planted_vertex.append(eval(line))
all_marked = []
marked_input = open(path + '/saved_data/all_marked_nodes_for_TH.txt', 'r')
for line in marked_input:
    all_marked.append(eval(line))
all_marked = all_marked[0]

print(utiles.calculate_presentage(all_marked,planted_vertex))
