import os


path = '/Users/ronizoller/PycharmProjects/RSAM/PythonCode/data/COG2367(class_A)/'
res = ''

input1 = open(path + "/old_new_names.txt", 'r')
old_new_names = []
for line in input1:
    old_new_names.append(eval(line))
old_new_names = old_new_names[0]

input = open(path + '/taxa_names.txt', 'r')
for line in input:
    if line != '\n':
        # print(line,end='')
        res += line[0:len(line) - 1].replace(' ','_') + '\n'

file = open(path+'/missing_tags.txt', 'w')
file.write(res)
file.close()

if False:
    input = open(path + '/0/colors0.0.txt', 'r')
    colors = []
    for line in input:
        colors.append(eval(line))
    colors = colors[0]

    input = open(path + '/taxa_names.txt', 'r')
    input1 = open(path + "/old_new_names.txt", 'r')
    old_new_names = []
    for line in input1:
        old_new_names.append(eval(line))
    old_new_names = old_new_names[0]

    names = []
    names_with_space = {}
    for line in input:
        if line != '\n':
            # print(line,end='')
            names.append((line[0:len(line) - 1]).replace(' ', ''))
            names_with_space.update(
                {(line[0:len(line) - 1]).replace(' ', ''): (line[0:len(line) - 1]).replace(' ', '_')})
    res = ""
    for name in names:
        found = False
        # print(name)
        for spe, color in colors.items():
            if (old_new_names[name].find(spe) != -1 or spe.find(old_new_names[name]) != -1):
                found = True
        if not found:
            res = res + names_with_space[name] + '\n'
