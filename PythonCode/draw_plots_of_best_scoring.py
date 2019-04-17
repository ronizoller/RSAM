import sys
sys.path.append('/anaconda3/lib/python3.6/site-packages')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


path = '/storage/DATA/users/ronizo/COGS/'
cog_marked_D = {}
cog_marked_HT = {}
marked_for_cog_D = []
marked_for_cog_HT = []

cogs_names = ['COG4679','COG3549','COG2856']
ext = 'proteobacteria'

for cog in cogs_names:
    planted_vertices = []

    input = open(path + cog +'/saved_data/marked_RSAM_finder_'+ext+'_pattern=(D).txt', 'r')
    marked_for_cog_D = []
    for line in input:
        marked_for_cog_D.append(eval(line))
    marked_for_cog_D = marked_for_cog_D[0]
    cog_marked_D.update({cog:marked_for_cog_D})

    input = open(path + cog +'/saved_data/marked_RSAM_finder_'+ext+'_pattern=(HT).txt', 'r')
    marked_for_cog_HT = []
    for line in input:
        marked_for_cog_HT.append(eval(line))
    marked_for_cog_HT = marked_for_cog_HT[0]
    cog_marked_HT.update({cog:marked_for_cog_HT})


plt.clf()
names_marked = []
names_of_marked_D = []
names_of_marked_HT = []
list_to_draw_D = []
list_to_draw_HT = []
xs = []
maxi = 0
max_score_level = 0
i = 1
width = 0.05
fig, ax = plt.subplots(figsize=(17,10))

for cog, marked_vertices in cog_marked_D.items():
    temp_list_of_scores = []
    for u, score in marked_vertices.items():
        names_marked.append(cog)
        names_of_marked_D.append(u)
        temp_list_of_scores.append(max(score[0][0],score[0][1]))
        xs.append(i)
        if maxi < score[0][0]:
            maxi = score[0][0]
        if maxi < score[0][1]:
            maxi = score[0][1]
    list_to_draw_D += sorted(temp_list_of_scores)[0:5]
    i += 1

i = 0
bar0 = []
bar1 = []
bar2 = []
bar3 = []
bar4 = []
for y in list_to_draw_D:
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

ind = np.arange(len(cogs_names))+0.1
#ax.plot(xs, list_to_draw,marker='o', color = 'fuchsia',linestyle = 'None')
rects0 = ax.bar(ind, bar0, width,color='blue',edgecolor='k')
rects1 = ax.bar(ind+width, bar1, width,color='blue',edgecolor='k')
rects2 = ax.bar(ind+2*width, bar2, width,color='blue',edgecolor='k')
rects3 = ax.bar(ind+3*width, bar3, width,color='blue',edgecolor='k')
rects4 = ax.bar(ind+4*width, bar4, width,color='blue',edgecolor='k',label='(D,none,False)')

for cog, marked_vertices in cog_marked_HT.items():
    temp_list_of_scores = []
    for u, score in marked_vertices.items():
        names_of_marked_HT.append(u)
        temp_list_of_scores.append(max(score[0][0],score[0][1]))
        xs.append(i)
        if maxi < score[0][0]:
            maxi = score[0][0]
        if maxi < score[0][1]:
            maxi = score[0][1]
    list_to_draw_HT += sorted(temp_list_of_scores)[0:5]
    i += 1
i = 0
bar0_HT = []
bar1_HT = []
bar2_HT = []
bar3_HT = []
bar4_HT = []
for y in list_to_draw_HT:
    if i % 5 == 0:
        bar0_HT.append(y)
    elif i % 5 == 1:
        bar1_HT.append(y)
    elif i % 5 == 2:
        bar2_HT.append(y)
    elif i % 5 == 3:
        bar3_HT.append(y)
    else:
        bar4_HT.append(y)
    i += 1
ind = ind+5*width
#ax.plot(xs, list_to_draw,marker='o', color = 'fuchsia',linestyle = 'None')
rects5 = ax.bar(ind, bar0_HT, width,color='pink',edgecolor='k')
rects6 = ax.bar(ind+width, bar1_HT, width,color='pink',edgecolor='k')
rects7 = ax.bar(ind+2*width, bar2_HT, width,color='pink',edgecolor='k')
rects8 = ax.bar(ind+3*width, bar3_HT, width,color='pink',edgecolor='k')
rects9 = ax.bar(ind+4*width, bar4_HT, width,color='pink',edgecolor='k',label='(HT,none,True)')


sns.set(font_scale=1.4)
plt.xlabel('COGS', fontsize=10)

plt.xticks(ind,cogs_names,rotation=45)

plt.ylabel('Score', fontsize=10)
plt.axis([0, len(cogs_names), 0, maxi + 0.001])

ax.legend()

plt.savefig('/users/studs/bsc/2016/ronizo/PycharmProjects/RSAM/COGS/figures/results.png')