import sys
sys.path.append('/anaconda3/lib/python3.6/site-packages')
import numpy as np
import matplotlib.pyplot as plt
import utiles

path = '/Users/ronizoller/PycharmProjects/TreeReconciliation/trees/Simulator/one_best_compare'
number_of_one_bests = 40

number_of_random_inputs = 4
planted_vertex = {}
all_marked = {}
precentage = []                                             #{random_input_index:(false pos,false neg)}
our_precentage = []
our_sol = {}
for random_input in range(0,number_of_random_inputs):
    all_marked.update({random_input:{}})
    planted_vertex.update({random_input:[]})
    planted_vertex_input = open(path + '/saved_data/marked_nodes_to_save/' + str(random_input) + '/planted.txt', 'r')
    for line in planted_vertex_input:
        planted_vertex[random_input].append(eval(line))
    our_sol.update({random_input:[]})
    our_sol_input = open(path + '/saved_data/marked_nodes_to_save/' + str(random_input) + '/marked.txt', 'r')
    for line in our_sol_input:
        our_sol[random_input].append(eval(line))
    for i in range(0,number_of_one_bests):
        marked_input = open(path + '/saved_data/marked_nodes_to_save/'+str(random_input)+'/'+str(i)+'.txt', 'r')
        all_marked[random_input].update({i:[]})
        for line in marked_input:
            all_marked[random_input][i].append(eval(line))
        all_marked[random_input][i] = all_marked[random_input][i][0]
    precentage.append(list(utiles.calculate_presentage(all_marked[random_input],planted_vertex[random_input])))
    our_precentage.append(list(utiles.calculate_presentage({1:our_sol[random_input][0]},planted_vertex[random_input])))

data = [j for i in zip(precentage,our_precentage) for j in i]

data = utiles.flip_list(data,2)
plt.figure(12, figsize=(7, 7))  # size of fig


#columns = tuple(j for i in zip(['%s one_best' % x for x in range(0,10)],['%s RSAM-finder' % x for x in range(0,number_of_random_inputs)]) for j in i)
columns = tuple(j for i in zip(['One Best']* 10,['RSAM'] *number_of_random_inputs) for j in i)
rows = ('False Positive','False Negative')

values = np.arange(0, 120, 20)

value_increment = 20

# Get some pastel shades for the colors
colors = ['red','grey']
n_rows = len(data)

index = np.arange(len(columns)) + 0.3
bar_width = 0.3

# Initialize the vertical-offset for the stacked bar chart.

# Plot bars and create text labels for the table
cell_text = []
plt.bar(index, data[0], bar_width, color=colors[0])
cell_text.append(list(map("{}%".format, data[0])))
plt.bar(index+bar_width, data[1], bar_width, color=colors[1])
cell_text.append(list(map("{}%".format, data[1])))

# Reverse colors and text labels to display the last value at the top.
colors = colors[::-1]
cell_text.reverse()
# Add a table at the bottom of the axes
the_table = plt.table(cellText=cell_text,
                      rowLabels=rows,
                      rowColours=colors,
                      colLabels=columns,
                      loc='bottom')
the_table.auto_set_font_size(False)
the_table.set_fontsize(10)
# Adjust layout to make room for the table:
plt.subplots_adjust(left=0.2, bottom=0.2)

plt.ylabel("Pracentage")
plt.yticks([0,20,40,60,80,100])
plt.xticks([])
plt.title('Comparison between One Best to RSAM')
plt.savefig(path + '/plot_comp.png')