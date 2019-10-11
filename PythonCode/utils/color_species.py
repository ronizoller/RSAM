path = '/Users/ronizoller/PycharmProjects/RSAM/PythonCode/data/COG2602(class_D,gram_strain,all)/'

tag_by_gram_pos_neg = True

input = open(path + '/taxa_names.txt', 'r')

missing_tags = []
for line in input:
    missing_tags.append(line[0:len(line)-1].replace(' ','_'))

input = open(path+'/tags_from_patric.txt', 'r')

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
res = "{ "
not_anottated = []

if tag_by_gram_pos_neg:
    to_be_black = ['Negative','negative','-','_']
    to_be_red = ['Positive','positive','+','gram +','Posotive']

else:
    to_be_red = ["HUMAN","Yellowfin tuna","bovine","Anomalops katoptron","Macaca","Canis","Milkfish","salmon","Rabbit","Tursiops","flounder","Aplysia","Paralichthys","cattle","Timber rattlesnake","Goat","ANIMAL","Pagrus major","Abalone","cow","Multiple",'mice',"rumen microbiome","Bird","snail","tern","Bat","Pig","Danaus",
                 "Entylia carinata","butterfly","Chicken","Nilaparvata","Aleurodicus","Crab","Coleopteran",'silkworm',"Human","Silk worm","Termite","Acanthamoeba sp","Duck","Goose", "Animal","clams","Caenorhabditis","Ark clam","Fish", "Mouse", "Bos taurus","bear", "Mus musculus", "Apis mellifera",
         "rabbit", "rabbit"," Sponge", "Seal","Rattus","Lizard","Photoblepharon palpebratus","human","Drosophila","Cows","insect","Litopenaeus vannamei","Danio rerio","Sheep","Salvelinus namaycush","booklouse","Ctenocephalides","Gastrointestinal"
         "Cat","animals"]
    to_be_black = ["PLANT","Bivalve","Mosquito","Eucalyptus","Taxus cuspidata","Elysia ornata","poplar","Japanese horseradish","Glossina","Polyopes lancifolius","Pinus","banana","spider","Cockroach","Amblyomma","Bemisia","Sea","Pine","Neogoniolithon fosliei","Apple","Pea","Acyrthosiphon pisum","pepper","Ponds","Saintpaulia","Terpios","Mimosa scabrella","Bivalvia","Wheat","Tassel Hyacinth",'Tick',"Maize","fly","Beet","Skeletonema","Insect","Haliclona","Steinernema","Heterorhabditis","Racocetra castanea","mud","Gigaspora","Scutellospora pellucida","Oyster","Trachelipus","Argopecten","Ciona intestinali","Cauliflower","Apium","cabbage","Cichorium","Pear","Brassica","Gymnodinium","Machaerium","carrota","Hedera","Lactuca","Lissodendoryx","Cassava","Manihot","Bursaphelenchus","Quercus robur","Coffee","Figaro","Fragaria","Environmental","Bean","Saccharum","Tubeworm","Solemya","clam","ENV","Soil","Epilithic","lotus","coal mine","cave","milk","chromium","air","ice","Vitis vinifera L. grapevine","mosses","salt","saltern","ocean","planet","sauce","galaga","habitat","tannery wastes","Compost","seafood","saltpan","roots","lake","Lotus corniculatus","rhizospheres","Terrestrial","Mud","soil","sewage sludge","Plant","environmental","tomato","Tomato","Stylosanthes","Least snoutbean","Solanum","Mandevilla","grass","chilli","eggplant","Anthyllis","Maytenus","Astragalus","Lessertia","Prunus","Maytenus","Lessertia","Thale cress","Cherry","Indigofera","Willow","Montipora","Centrosema pubescens","Black mimosa","Chinese silvergrass", "vents","Skeletonema marinoi","Tachigali versicolor","Coral","Poa pratensis","Lotononis carinata","Mortierella elongata","Weeping fig","Cottonwood","groundwater","Date palm","Limoniastrum monopetalum","Cuminum cyminum","Otholobium candicans","lactate","Sponge","Peltigera didactyla","Lycium barbarum","Globobulimina","Lotus corniculatus","Lotus japonicus","Rice","rice","thermal vent","hydrothermal vent","sludge", "Prunus persica",
       "Soybean","Mimosa flocculosa","pineapple", "Elysia rufescens","clover", "water","Water", "Grape", "Sediment","sediment", "Marine", "Extreme Environments", 'offshore oil field', "Madracis decactis", "Oryza", "Acacia", "Algae", "algae", "Corn","nuts", "Sugarcane", "cotton", "Gossypium", "coccinea","Acanthamoeba",
         "ambigua", "Robinia","Festuca rubra", "spring","Spring","Asterionella formosa BG1","Stereocaulon","Biserrula","Cicer arietinum","Biserrula pelecinus","chickpea","Bituminaria","pea","Piriformospora indica",
       "Mimosa pudica","rhizosphere","Pyrus","Aeschynomene","Acanthamoebae","peat bog","aquatic","cellulolytic culture","sanitary sewer","coral","Lake","Aquatic","sea","marine",
                   "rock","Plautia stali","Dinoflagellate","Macrotermes","Tunicate","Cucurbita","Rafflesia","bean","Sesbania"]


list_not_to_dis = ["Missing","Specialized","Host-Associated","Host-associated","HostAssociated","Host","Unknown","Host Associated"]

input1 = open(path + "/old_new_names.txt", 'r')
old_new_names = []
for line in input1:
    old_new_names.append(eval(line))
old_new_names = old_new_names[0]
tagged = []

for name in missing_tags:
    i = 0
    flag = False
    for tag_name in tags_names:
        if tag_name.replace(' ','').replace('_','').find(name.replace(' ','').replace('_','')) != -1 or\
                name.replace(' ', '').replace('_', '').find(tag_name.replace(' ', '').replace('_', '')) != -1:
            for tag in to_be_red:
                if tags_habitat[i].replace('\t',"") != "" and (tag.find(tags_habitat[i]) != -1 or tags_habitat[i].find(tag) != -1) and not flag :
                    res = res + ", '" + old_new_names[name.replace('_',' ')].replace(' ','')+"':'red' "
                    flag = True
            for tag_env in to_be_black:
                if tags_habitat[i].replace('\t',"") != "" and (tag_env.find(tags_habitat[i]) != -1 or tags_habitat[i].find(tag_env) != -1) and not flag:
                    res = res + ", '" + old_new_names[name.replace('_', ' ')].replace(' ','') + "':'black' "
                    flag = True
        i+=1
    if not flag:
        j = 0
        for tag_name in tags_names:
            flag2 = True
            if tag_name.find(name) != -1:
                for tag in to_be_red:
                    if tags_habitat[j].replace('\t','').replace(' ','') != "" and flag2 and (tags_habitat[j] not in list_not_to_dis):
                        print("%s was not tagged, tag = %s" % (str(name),str(tags_habitat[j])))
                        flag2 = False
                for tag_env in to_be_black:
                    if tags_habitat[j].replace('\t','').replace(' ','') != "" and flag2 and (tags_habitat[j] not in list_not_to_dis):
                        print("%s was not tagged, tag = %s" % (str(name), str(tags_habitat[j])))
                        flag2 = False
                        flag = True
            j += 1
        if name not in not_anottated:
            not_anottated.append(name)
    elif name not in tagged:
        tagged.append(name)

res += '}'
file = open(path + '' + "/colors.txt", 'w')
file.write(str(res))
file.close()

file = open(path+'/missing_tags_after_patric.txt', 'w')
miss = 'number of tagged species: '+str(len(tagged))+'/'+str(len(missing_tags))+'\n\n'
for name in not_anottated:
    miss = miss +"\n" + name
file.write(miss)
file.close()

print('number of tagged species: '+str(len(tagged))+'/'+str(len(missing_tags)))
