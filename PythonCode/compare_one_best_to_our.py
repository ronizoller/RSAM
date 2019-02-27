import sys
sys.path.append('/anaconda3/lib/python3.6/site-packages')
import matplotlib.pyplot as plt
import utiles
import draw
on_lab = True
if on_lab:
    path  = '/users/studs/bsc/2016/ronizo/PycharmProjects/RSAM/simulator_data/comparsion_test'
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
print(result_RSAM)

xs = []
ys = []
to_connect_one_best = []
xs_RSAM = []
ys_RSAM = []
to_connect_RSAM = []

names = []
for TH,tup in result_one_best.items():
    xs.append(tup[0])
    ys.append(tup[1])
    names.append(TH)
    to_connect_one_best.append(tup)
for TH,tup in result_RSAM.items():
    xs_RSAM.append(tup[0])
    ys_RSAM.append(tup[1])
    to_connect_RSAM.append(tup)
xs.sort()
ys.sort()
xs_RSAM.sort()
ys_RSAM.sort()
fig, ax = plt.subplots()
ax.set_xlabel('False Positive Rate (1-Specifity)')
ax.set_ylabel('True Positive Rate (Sensitivity)')
for p in range(0, len(xs) - 1):
    draw.connectpoints(xs, ys, p, p + 1)
for p in range(0, len(xs_RSAM) - 1):
    draw.connectpoints(xs_RSAM, ys_RSAM, p, p + 1)
#for X, Y, Z in zip(xs, ys, names):
#    # Annotate the points 5 _points_ above and to the left of the vertex
#    ax.annotate('{}'.format(Z), xy=(X, Y), xytext=(5, -5), ha='left',
#                textcoords='offset points')
#for X, Y, Z in zip(xs_RSAM, ys_RSAM, names):
#    # Annotate the points 5 _points_ above and to the left of the vertex
#    ax.annotate('{}'.format(Z), xy=(X, Y), xytext=(5, -5), ha='left',
#                textcoords='offset points')
plt.plot(xs, ys, 'ko')
plt.plot(xs_RSAM, ys_RSAM, 'go')
plt.axis([0, 1, 0, 1])
plt.show()
#fig = ax.get_figure()
#fig.savefig(path + '/plot_noise.png')