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

cogs_names = ['COG4679','COG3549','COG2856','COG3550']
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
width = 0.04
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
    temp_list_of_scores = sorted(temp_list_of_scores)
    temp_list_of_scores = temp_list_of_scores[len(temp_list_of_scores)-10:]
    print(len(temp_list_of_scores))
    list_to_draw_D += sorted(temp_list_of_scores)
    i += 1

i = 0
bar0 = []
bar1 = []
bar2 = []
bar3 = []
bar4 = []
bar5 = []
bar6 = []
bar7 = []
bar8 = []
bar9 = []
for y in list_to_draw_D:
    if i % 10 == 0:
        bar0.append(y)
    elif i % 10 == 1:
        bar1.append(y)
    elif i % 10 == 2:
        bar2.append(y)
    elif i % 10 == 3:
        bar3.append(y)
    elif i % 10 == 4:
        bar4.append(y)
    elif i % 10 == 5:
        bar5.append(y)
    elif i % 10 == 6:
        bar6.append(y)
    elif i % 10 == 7:
        bar7.append(y)
    elif i % 10 == 8:
        bar8.append(y)
    else:
        bar9.append(y)
    i += 1
ind = np.arange(len(cogs_names))+0.1
#ax.plot(xs, list_to_draw,marker='o', color = 'fuchsia',linestyle = 'None')
rects0 = ax.bar(ind, bar0, width,color='blue',edgecolor='k')
rects1 = ax.bar(ind+width, bar1, width,color='blue',edgecolor='k')
rects2 = ax.bar(ind+2*width, bar2, width,color='blue',edgecolor='k')
rects3 = ax.bar(ind+3*width, bar3, width,color='blue',edgecolor='k')
rects4 = ax.bar(ind+4*width, bar4, width,color='blue',edgecolor='k')
rects5 = ax.bar(ind+5*width, bar5, width,color='blue',edgecolor='k')
rects6 = ax.bar(ind+6*width, bar6, width,color='blue',edgecolor='k')
rects7 = ax.bar(ind+7*width, bar7, width,color='blue',edgecolor='k')
rects8 = ax.bar(ind+8*width, bar8, width,color='blue',edgecolor='k')
rects9 = ax.bar(ind+9*width, bar9, width,color='blue',edgecolor='k',label='(D,none,False)')

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
    list_to_draw_HT += sorted(temp_list_of_scores)
    i += 1
i = 0
bar0_HT = []
bar1_HT = []
bar2_HT = []
bar3_HT = []
bar4_HT = []
bar5_HT = []
bar6_HT = []
bar7_HT = []
bar8_HT = []
bar9_HT = []
for y in list_to_draw_HT:
    if i % 10 == 0:
        bar0_HT.append(y)
    elif i % 10 == 1:
        bar1_HT.append(y)
    elif i % 10 == 2:
        bar2_HT.append(y)
    elif i % 10 == 3:
        bar3_HT.append(y)
    elif i % 10 == 4:
        bar4_HT.append(y)
    elif i % 10 == 5:
        bar5_HT.append(y)
    elif i % 10 == 6:
        bar6_HT.append(y)
    elif i % 10 == 7:
        bar7_HT.append(y)
    elif i % 10 == 8:
        bar8_HT.append(y)
    else:
        bar9_HT.append(y)
    i += 1
ind = ind+10*width
#ax.plot(xs, list_to_draw,marker='o', color = 'fuchsia',linestyle = 'None')
rects10 = ax.bar(ind, bar0_HT, width,color='pink',edgecolor='k')
rects11 = ax.bar(ind+width, bar1_HT, width,color='pink',edgecolor='k')
rects12 = ax.bar(ind+2*width, bar2_HT, width,color='pink',edgecolor='k')
rects13 = ax.bar(ind+3*width, bar3_HT, width,color='pink',edgecolor='k')
rects14 = ax.bar(ind+4*width, bar4_HT, width,color='pink',edgecolor='k')
rects15 = ax.bar(ind+5*width, bar5_HT, width,color='pink',edgecolor='k')
rects16 = ax.bar(ind+6*width, bar6_HT, width,color='pink',edgecolor='k')
rects17 = ax.bar(ind+7*width, bar7_HT, width,color='pink',edgecolor='k')
rects18 = ax.bar(ind+8*width, bar8_HT, width,color='pink',edgecolor='k')
rects19 = ax.bar(ind+9*width, bar9_HT, width,color='pink',edgecolor='k',label='(HT,none,True)')


sns.set(font_scale=1.4)
plt.xlabel('COGS', fontsize=10)

plt.xticks(ind,cogs_names,rotation=45)

plt.ylabel('Score', fontsize=10)
plt.axis([0, len(cogs_names), 0, maxi + 0.001])

ax.legend()

plt.savefig('/users/studs/bsc/2016/ronizo/PycharmProjects/RSAM/results.png')