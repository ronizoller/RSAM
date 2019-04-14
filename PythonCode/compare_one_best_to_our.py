import sys
sys.path.append('/anaconda3/lib/python3.6/site-packages')
import matplotlib.pyplot as plt
from adjustText import adjust_text
import utiles
import operator
import ast

with_labels = True
on_lab = True
if on_lab:
    path  = '/users/studs/bsc/2016/ronizo/PycharmProjects/RSAM/compare_final_data'
else:
    path = '/Users/ronizoller/PycharmProjects/TreeReconciliation/trees/compare_test'
    import sys
    sys.path.append('/anaconda3/lib/python3.6/site-packages')

number_of_one_bests = 1

precentage = []                                             #{random_input_index:(false pos,false neg)}
our_precentage = []
our_sol = {}
planted_vertex = []
planted_vertex_input = open(path + '/saved_data/planted_nodes_correct_names.txt', 'r')
for line in planted_vertex_input:
    planted_vertex.append(eval(line))
planted_vertex = planted_vertex[0]

all_marked = []
marked_input = open(path + '/saved_data/all_vertices_TH.txt', 'r')
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
all_RSAM_marked_fixed = {}
for TH,x in all_RSAM_marked.items():
    all_RSAM_marked_fixed.update({TH:[]})
    for u,scores in x[0].items():
        all_RSAM_marked_fixed[TH].append((u,scores))


result_one_best = utiles.calculate_presentage(all_marked,all_unmarked,planted_vertex)
result_RSAM = utiles.calculate_presentage(all_RSAM_marked_fixed,all_RSAM_unmarked,planted_vertex)

def plot_eucs_covers(eucs,covers,text,color):
    to_draw = sorted(zip(eucs, covers), key=operator.itemgetter(0))
    plt.plot(list(zip(*to_draw))[0],list(zip(*to_draw))[1],color=color, alpha=0.5,linewidth=4)
    texts = []
    if with_labels:
        for xt, yt, s in zip(eucs, covers, text):
            texts.append(plt.text(xt, yt, s,bbox=dict(facecolor=color, alpha=0.5)))
    return texts

xs = []
ys = []
xs_RSAM = []
ys_RSAM = []

names = []
one_best_names = []
THs = [(f,0,f*10) for f in [0,0.5,1,2,2.5,3,3.1,3.2,3.3,3.5,3.6,3.7,3.8,3.9,4,4.1,4.2,4.3,4.4,4.5,5.2,7,8,9,10,20]]
for TH,tup in result_RSAM.items():
    if TH in THs:
        names.append(TH)
        xs_RSAM.append(tup[0])
        ys_RSAM.append(tup[1])
for TH,tup in result_one_best.items():
    index0 = TH.find('(')
    index1 = TH.find(',')
    index2 = TH[index1+1:].find(',')+index1+1
    to_eval1 = TH[index0+1:index1]
    to_eval2 = TH[index2+1:-1]
    temp_name = (ast.literal_eval(to_eval1),0,ast.literal_eval(to_eval2))
    if temp_name in names:
        xs.append(tup[0])
        ys.append(tup[1])
        one_best_names.append((round(temp_name[0],2),round(temp_name[2],2)))

names = [(round(name[0],2),round(name[2],2)) for name in names]
fig, ax = plt.subplots(figsize=(16,10))
ax.set_xlabel('False Positive Rate (1-Specifity)')
ax.set_ylabel('True Positive Rate (Sensitivity)')

#plt.plot(xs, ys, 'ko')
#plt.plot(xs_RSAM, ys_RSAM, 'go')

plt.axis([0, 0.06, 0, 1.1])
RSAM_text = plot_eucs_covers(xs_RSAM,ys_RSAM,names,'g')
one_best_text = plot_eucs_covers(xs,ys,one_best_names,'k')
if with_labels:
    adjust_text(one_best_text+RSAM_text)

print('One_best:\n    xs:%s\n     ys:%s' % (str(xs),str(ys)))
print('RSAM:\n    xs:%s\n     ys:%s' % (str(xs_RSAM),str(ys_RSAM)))
if with_labels:
    fig.savefig(path+'/plots/plot.png')
else:
    fig.savefig(path+'/plots/plot_no_labels.png')
