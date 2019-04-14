import sys
sys.path.append('/anaconda3/lib/python3.6/site-packages')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


path = '/Users/ronizoller/PycharmProjects/TreeReconciliation/trees/COG3550'
classes_marked = {}
classes_names = ['alpha','beta','gamma','delta','epsilon']

for clas in classes_names:
    planted_vertices = []

    input = open(path + '/saved_data/marked_RSAM_finder_'+clas+'.txt', 'r')
    marked_for_clas = []
    for line in input:
        marked_for_clas.append(eval(line))
    marked_for_clas = marked_for_clas[0]
    classes_marked.update({clas:marked_for_clas})


plt.clf()
classes = []
names_marked = []
names_of_marked = []
list_to_draw = []
xs = []
maxi = 0
max_score_level = 0
i = 1
width = 0.1

for clas_name, marked_vertices in classes_marked.items():
    temp_list_of_scores = []
    for u, score in marked_vertices.items():
        names_marked.append(clas_name)
        names_of_marked.append(u)
        temp_list_of_scores.append(max(score[0][0],score[0][1]))
        xs.append(i)
        if maxi < score[0][0]:
            maxi = score[0][0]
        if maxi < score[0][1]:
            maxi = score[0][1]
    list_to_draw += sorted(temp_list_of_scores)
    i += 1
fig, ax = plt.subplots(figsize=(5,11))

i = 0
bar0 = []
bar1 = []
bar2 = []
bar3 = []
bar4 = []
for y in list_to_draw:
    if i % 5 == 0:
        bar0.append(y)
    elif i % 5 == 1:
        bar1.append(y)
    elif i % 5 == 2:
        bar2.append(y)
    elif i % 5 == 3:
        bar3.append(y)
    else:
        bar4.append(y)
    i += 1
ind = np.arange(5)+0.2
#ax.plot(xs, list_to_draw,marker='o', color = 'fuchsia',linestyle = 'None')
rects0 = ax.bar(ind, bar0, width,color='fuchsia',edgecolor='k')
rects1 = ax.bar(ind+width, bar1, width,color='fuchsia',edgecolor='k')
rects2 = ax.bar(ind+2*width, bar2, width,color='fuchsia',edgecolor='k')
rects3 = ax.bar(ind+3*width, bar3, width,color='fuchsia',edgecolor='k')
rects4 = ax.bar(ind+4*width, bar4, width,color='fuchsia',edgecolor='k')




sns.set(font_scale=1.4)
plt.xlabel('Proteobacteria Classes', fontsize=10)

plt.xticks(ind+2*width,['Proteobacteria'+clas for clas in classes_names],rotation=45)

plt.ylabel('Score', fontsize=10)
plt.axis([0, len(classes_names), 0, maxi + 0.001])
ax.set_title('COG3550 Top Ranking')


plt.savefig(path+'/figures/proteobacteria_comparsion.png')