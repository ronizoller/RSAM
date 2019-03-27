import sys
sys.path.append('/anaconda3/lib/python3.6/site-packages')
import matplotlib.pyplot as plt
import draw
import seaborn as sns

path = '/users/studs/bsc/2016/ronizo/PycharmProjects/RSAM/noise_final_data'

for noise_in in ['colors_and_HT']:
    planted_vertices = []
    input = open(path + '/saved_data/planted_nodes_correct_names.txt', 'r')
    for line in input:
        planted_vertices.append(eval(line))
    planted_vertices = planted_vertices[0]

    input = open(path + '/saved_data/all_vertices_RSAM_finder_'+noise_in+'.txt', 'r')
    all_vertices_with_noise = []
    for line in input:
        all_vertices_with_noise.append(eval(line))
    all_vertices_with_noise = all_vertices_with_noise[0]

    plt.clf()
    list_to_draw_reds = []
    list_to_draw_blacks = []
    length = []
    names = []
    names_marked = []
    to_connect = []
    length_to_connect = []
    maxi = 0
    max_noise_level = 0
    for noise_level, all_vertices in all_vertices_with_noise.items():
        if noise_level <= 20:
            for u, couple in all_vertices.items():
                if u != 'u540':
                    if u in planted_vertices:
                        to_connect.append(max(couple[1],couple[0]))
                        length_to_connect.append(noise_level)
                        names_marked.append(u)
                    else:
                        length.append(noise_level)
                        names.append(u)
                        list_to_draw_reds.append(couple[0])
                        list_to_draw_blacks.append(couple[1])
                    if maxi < couple[0]:
                        maxi = couple[0]
                    if maxi < couple[1]:
                        maxi = couple[1]
                    if max_noise_level < noise_level:
                        max_noise_level = noise_level
    fig, ax = plt.subplots(figsize=(8,12))
    ax.plot(length, list_to_draw_reds, 'ko', length, list_to_draw_blacks, 'ko')
    ax.plot(length_to_connect, to_connect, marker='o', color='fuchsia',linestyle = 'None')

    #for X, Y, Z in zip(length, list_to_draw_reds, names):
    #    # Annotate the points 5 _points_ above and to the left of the vertex
    #    ax.annotate('{}'.format(Z), xy=(X, Y), xytext=(5, -5), ha='left',
    #                textcoords='offset points')
    #for X, Y, Z in zip(length, list_to_draw_blacks, names):
    #    # Annotate the points 5 _points_ above and to the left of the vertex
    #    ax.annotate('{}'.format(Z), xy=(X, Y), xytext=(5, -5), ha='left',
    #                textcoords='offset points')
    #for X, Y, Z in zip(length_to_connect, to_connect, names_marked):
    #    # Annotate the points 5 _points_ above and to the left of the vertex
    #    ax.annotate('{}'.format(Z), xy=(X, Y), xytext=(5, -5), ha='left',
    #                textcoords='offset points')
    sns.set(font_scale=1.4)
    plt.xlabel('Noise', fontsize=10)
    plt.ylabel('Score', fontsize=10)
    plt.axis([0, max_noise_level + 0.5, 0, maxi + 0.0001])
    ax.set_title('Noise was added in '+noise_in)

    plt.savefig('/users/studs/bsc/2016/ronizo/PycharmProjects/RSAM/noise_final_data/plots/'+noise_in+'.png')