path = '/Users/ronizoller/PycharmProjects/TreeReconciliation/trees/real'
#first prtric then dina

patric = True

if patric:
    input = open(path+'/missing_tags.txt', 'r')
else:
    input = open(path+'/missing_tags_after_patric.txt', 'r')
missing_tags = []
for line in input:
    missing_tags.append(line[0:len(line)-1])

if patric:
    input = open(path+'/tags_from_patric.txt', 'r')
else:
    input = open('bactTaxa_Habitat.txt', 'r')

tags_names = []
tags_habitat = []
for line in input:
    start = line.find(';')+1
    end = line[start:len(line)].find(';')
    name = line[start:start+end]

    line = line[::-1]
    start = line.find(';')
    habitat = line[0:start]
    habitat = habitat[::-1]
    tags_names.append(name)
    tags_habitat.append(habitat[0:len(habitat)-1])
res = ""
not_anottated = []

alive = ["Pagrus major","Tick","Human","Termite","Acanthamoeba sp","Duck","Goose", "Animal","clams","Caenorhabditis","Ark clam","Fish", "Mouse", "Bos taurus","bear", "Mus musculus", "Apis mellifera",
         "rabbit", "rabbit"," Sponge", "Seal", "Drosophila","Cows","insect","Litopenaeus vannamei","Danio rerio","Sheep","Salvelinus namaycush","Amblyomma","Amblyomma","booklouse","Ctenocephalides",
         "Cat","animals"]
env = ["Soil","soil", "vents","groundwater","lactate","Sponge","Peltigera didactyla","Lycium barbarum","Globobulimina","Lotus japonicus","Rice","rice","thermal vent","hydrothermal vent","sludge", "Prunus persica",
       "Soybean", "clover", "water","Water", "Grape", "Sediment","sediment", "Marine", "Extreme Environments", 'offshore oil field', "Madracis decactis", "Oryza", "Acacia", "Algae", "algae", "Corn","nuts", "Sugarcane", "cotton", "Gossypium", "coccinea",
         "ambigua", "Robinia","Festuca rubra", "spring","Spring","Asterionella formosa BG1","Stereocaulon","Biserrula","Cicer arietinum","Biserrula pelecinus","chickpea","Bituminaria","pea","Piriformospora indica",
       "Mimosa pudica","rhizosphere","Aeschynomene","Acanthamoebae","peat bog","aquatic","cellulolytic culture","sanitary sewer","coral","Lake","Aquatic","sea","marine",]
#print("tags_names = %s\nmissing_tags = %s\ntags_habitat = %s" % (str(tags_names),str(missing_tags),str(tags_habitat)))

for name in missing_tags:
    i = 0
    flag = False
    for tag_name in tags_names:
        if tag_name.find(name) != -1:
            for tag in alive:
                if tags_habitat[i] != "" and (tag.find(tags_habitat[i]) != -1 or tags_habitat[i].find(tag) != -1) and not flag:
                    res = res + ", '" + name.replace('_','')+"':'red' "
                    flag = True
            for tag_env in env:
                if tags_habitat[i] != "" and (tag_env.find(tags_habitat[i]) != -1 or tags_habitat[i].find(tag_env) != -1) and not flag:
                    res = res + ", '" + name.replace('_', '') + "':'black' "
                    flag = True
        i = i+1
    if flag == False:
        j = 0
        for tag_name in tags_names:
            if tag_name.find(name) != -1:
                for tag in alive:
                    if tags_habitat[j] != "":
                        print("%s was not tagged, tag = %s" % (str(name),str(tags_habitat[j])))
                for tag_env in env:
                    if tags_habitat[j] != "":
                        print("%s was not tagged, tag = %s" % (str(name), str(tags_habitat[j])))
                        flag = True
            j = j + 1
        not_anottated.append(name)
with open(path+"/colors.txt", "a") as myfile:
    myfile.write(res+"}")

if patric:
    file = open(path+'/missing_tags_after_patric.txt', 'w')
    miss = ""
    for name in not_anottated:
        miss = miss +"\n" + name
    file.write(miss)
    file.close()
else:
    file = open(path + '/missing_tags_after_dina.txt', 'w')
    miss = ""
    for name in not_anottated:
        miss = miss + "\n" + name
    file.write(miss)
    file.close()
