import sys
sys.path.append('/anaconda3/lib/python3.6/site-packages')
import numpy as np
import matplotlib.pyplot as plt
import utiles
on_lab = True
if on_lab:
    path  = '/users/studs/bsc/2016/ronizo/PycharmProjects/RSAM/simulator_data/comparsion'
else:
    path = '/Users/ronizoller/Documents/school/Master/מחקר/DATA/comparsion'
    import sys
    sys.path.append('/anaconda3/lib/python3.6/site-packages')

number_of_one_bests = 10

precentage = []                                             #{random_input_index:(false pos,false neg)}
our_precentage = []
our_sol = {}
planted_vertex = []
planted_vertex_input = open(path + '/saved_data/marked_nodes_correct_names.txt', 'r')
for line in planted_vertex_input:
    planted_vertex.append(eval(line))
all_marked = []
marked_input = open(path + '/saved_data/all_vertices_RSAM_finder.txt', 'r')
for line in marked_input:
    all_marked.append(eval(line))
all_marked = all_marked[0]

all_unmarked = []
unmarked_input = open(path + '/saved_data/all_unmarked_nodes_for_TH.txt', 'r')
for line in unmarked_input:
    all_unmarked.append(eval(line))
all_unmarked = all_unmarked[0]


result = utiles.calculate_presentage(all_marked,all_unmarked,planted_vertex)
xs = []
ys = []
names = []
for TH,tup in result.items():
    xs.append(tup[0])
    ys.append(tup[1])
    names.append(TH)
plt.xlabel('False Positive Rate (1-Specifity)', fontsize=10)
plt.ylabel('True Positive Rate (Sensitivity)', fontsize=10)
fig, ax = plt.subplots()
for X, Y, Z in zip(xs, ys, names):
    # Annotate the points 5 _points_ above and to the left of the vertex
    ax.annotate('{}'.format(Z), xy=(X, Y), xytext=(5, -5), ha='left',
                textcoords='offset points')
plt.plot(xs, ys, 'ro')
plt.axis([0, 1, 0, 1])
plt.show()