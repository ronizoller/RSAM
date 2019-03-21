path = '/Users/ronizoller/PycharmProjects/TreeReconciliation/trees/COG2165'

input = open(path+'/0/colors0.0.txt', 'r')
colors = []
for line in input:
    colors.append(eval(line))
colors = colors[0]

input = open(path+'/taxa_names.txt', 'r')
names = []
names_with_space = {}
for line in input:
    if line != '\n':
        #print(line,end='')
        names.append((line[0:len(line)-1]).replace(' ',''))
        names_with_space.update({(line[0:len(line)-1]).replace(' ','') : (line[0:len(line)-1]).replace(' ','_')})
res = ""
for name in names:
    found = False
    #print(name)
    for spe,color in colors.items():
        if (name.find(spe) != -1 or spe.find(name) != -1):
            found = True
    if not found:
        res = res + names_with_space[name]+'\n'

file = open(path+'/missing_tags.txt', 'w')
file.write(res)
file.close()
