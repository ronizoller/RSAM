import sys
sys.path.append('/anaconda3/lib/python3.6/site-packages')
import matplotlib.pyplot as plt
from adjustText import adjust_text
from scipy import interpolate
import utiles
import draw
import operator
import ast
import math
import numpy as np

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
planted_vertex_input = open(path + '/saved_data/planted_nodes_correct_names.txt', 'r')
for line in planted_vertex_input:
    planted_vertex.append(eval(line))
planted_vertex = planted_vertex[0]

all_marked = []
marked_input = open(path + '/saved_data/all_marked_nodes_for_TH.txt', 'r')
for line in marked_input:
    all_marked.append(eval(line))
all_marked = all_marked[0]

all_unmarked = []
unmarked_input = open(path + '/saved_data/all_unmarked_nodes_for_TH.txt', 'r')
for line in unmarked_input:
    all_unmarked.append(eval(line))
all_unmarked = all_unmarked[0]

all_RSAM_marked = []
marked_RSAM_input = open(path + '/saved_data/all_marked_vertices_RSAM_finder.txt', 'r')
for line in marked_RSAM_input:
    all_RSAM_marked.append(eval(line))
all_RSAM_marked = all_RSAM_marked[0]

all_RSAM_unmarked = []
unmarked_RSAM_input = open(path + '/saved_data/all_unmarked_vertices_RSAM_finder.txt', 'r')
for line in unmarked_RSAM_input:
    all_RSAM_unmarked.append(eval(line))
all_RSAM_unmarked = all_RSAM_unmarked[0]

all_RSAM_marked = dict((u,list(x)) for (u,x) in all_RSAM_marked.items())

result_one_best = utiles.calculate_presentage(all_marked,all_unmarked,planted_vertex)
result_RSAM = utiles.calculate_presentage(all_RSAM_marked,all_RSAM_unmarked,planted_vertex)

def plot_eucs_covers(eucs,covers,text,color):
    to_draw = sorted(zip(eucs, covers), key=operator.itemgetter(0))
    plt.plot(list(zip(*to_draw))[0],list(zip(*to_draw))[1],color=color, alpha=0.5,linewidth=4)
    texts = []
    for xt, yt, s in zip(eucs, covers, text):
        texts.append(plt.text(xt, yt, s,bbox=dict(facecolor=color, alpha=0.5)))
    return texts

xs = []
ys = []
xs_RSAM = []
ys_RSAM = []

names = []
one_best_names = []
for TH,tup in result_RSAM.items():
    if TH[0] == TH[1]:
        xs_RSAM.append(tup[0])
        ys_RSAM.append(tup[1])
        names.append(TH)
for TH,tup in result_one_best.items():
    index0 = TH.find('(')
    index1 = TH.find(',')
    index2 = TH[index1+1:].find(',')+index1+1
    temp_name = (ast.literal_eval(TH[index0+1:index1]),ast.literal_eval(TH[index1+1:index2]),0)
    if temp_name in names:
        xs.append(tup[0])
        ys.append(tup[1])
        one_best_names.append((round(temp_name[0],2),round(temp_name[1],2)))

names = [(round(name[0],2),round(name[1],2)) for name in names]
fig, ax = plt.subplots()
ax.set_xlabel('False Positive Rate (1-Specifity)')
ax.set_ylabel('True Positive Rate (Sensitivity)')

plt.plot(xs, ys, 'ko')
plt.plot(xs_RSAM, ys_RSAM, 'go')

plt.axis([0, 1.2, 0, 1.2])
RSAM_text = plot_eucs_covers(xs_RSAM,ys_RSAM,names,'g')
one_best_text = plot_eucs_covers(xs,ys,one_best_names,'k')

adjust_text(RSAM_text+one_best_text)
#plt.show()
fig = ax.get_figure()
fig.savefig(path + '/plot_noise.png')