path = '/Users/ronizoller/PycharmProjects/TreeReconciliation/trees/COG2165'
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

alive = ["HUMAN","Yellowfin tuna","bovine","Anomalops katoptron","Macaca","Canis","Milkfish","salmon","Rabbit","Tursiops","flounder","Aplysia","Paralichthys","cattle","Timber rattlesnake","Goat","ANIMAL","Pagrus major","Abalone","cow","Multiple",'mice',"rumen microbiome","Bird","snail","tern","Bat","Pig","Danaus","Entylia carinata","butterfly","Chicken","Nilaparvata","Aleurodicus","Crab","Coleopteran","Human","Silk worm","Termite","Acanthamoeba sp","Duck","Goose", "Animal","clams","Caenorhabditis","Ark clam","Fish", "Mouse", "Bos taurus","bear", "Mus musculus", "Apis mellifera",
         "rabbit", "rabbit"," Sponge", "Seal","Rattus","Lizard","Photoblepharon palpebratus","human","Drosophila","Cows","insect","Litopenaeus vannamei","Danio rerio","Sheep","Salvelinus namaycush","booklouse","Ctenocephalides",
         "Cat","animals"]
env = ["PLANT","Bivalve","Mosquito","Eucalyptus","Taxus cuspidata","Elysia ornata","poplar","Japanese horseradish","Glossina","Polyopes lancifolius","Pinus","banana","spider","Cockroach","Amblyomma","Bemisia","Sea","Pine","Neogoniolithon fosliei","Apple","Pea","Acyrthosiphon pisum","pepper","Ponds","Saintpaulia","Mimosa scabrella","Bivalvia","Wheat","Tassel Hyacinth",'Tick',"Maize","fly","Beet","Skeletonema","Insect","Haliclona","Steinernema","Heterorhabditis","Racocetra castanea","mud","Gigaspora","Scutellospora pellucida","Oyster","Trachelipus","Argopecten","Ciona intestinali","Cauliflower","Apium","cabbage","Cichorium","Pear","Brassica","Gymnodinium","Machaerium","carrota","Hedera","Lactuca","Lissodendoryx","Cassava","Manihot","Bursaphelenchus","Quercus robur","Coffee","Figaro","Fragaria","Environmental","Bean","Saccharum","Tubeworm","Solemya","clam","ENV","Soil","Epilithic","lotus","coal mine","cave","milk","chromium","air","ice","Vitis vinifera L. grapevine","mosses","salt","saltern","ocean","planet","sauce","galaga","habitat","tannery wastes","Compost","seafood","saltpan","roots","lake","Lotus corniculatus","rhizospheres","Terrestrial","Mud","soil","sewage sludge","Plant","environmental","tomato","Tomato","Stylosanthes","Least snoutbean","Solanum","Mandevilla","grass","chilli","eggplant","Anthyllis","Maytenus","Astragalus","Lessertia","Prunus","Maytenus","Lessertia","Thale cress","Cherry","Indigofera","Willow","Montipora","Centrosema pubescens","Black mimosa","Chinese silvergrass", "vents","Skeletonema marinoi","Tachigali versicolor","Coral","Poa pratensis","Lotononis carinata","Mortierella elongata","Weeping fig","Cottonwood","groundwater","Date palm","Limoniastrum monopetalum","Cuminum cyminum","Otholobium candicans","lactate","Sponge","Peltigera didactyla","Lycium barbarum","Globobulimina","Lotus corniculatus","Lotus japonicus","Rice","rice","thermal vent","hydrothermal vent","sludge", "Prunus persica",
       "Soybean","Mimosa flocculosa","pineapple", "Elysia rufescens","clover", "water","Water", "Grape", "Sediment","sediment", "Marine", "Extreme Environments", 'offshore oil field', "Madracis decactis", "Oryza", "Acacia", "Algae", "algae", "Corn","nuts", "Sugarcane", "cotton", "Gossypium", "coccinea",
         "ambigua", "Robinia","Festuca rubra", "spring","Spring","Asterionella formosa BG1","Stereocaulon","Biserrula","Cicer arietinum","Biserrula pelecinus","chickpea","Bituminaria","pea","Piriformospora indica",
       "Mimosa pudica","rhizosphere","Pyrus","Aeschynomene","Acanthamoebae","peat bog","aquatic","cellulolytic culture","sanitary sewer","coral","Lake","Aquatic","sea","marine","rock",]


list_not_to_dis = ["Missing","Specialized","Host-Associated","Host-associated","HostAssociated","Host","Unknown","Host Associated"]

for name in missing_tags:
    i = 0
    flag = False
    for tag_name in tags_names:
        if tag_name.find(name) != -1:
            for tag in alive:
                if tags_habitat[i] != "" and (tag.find(tags_habitat[i]) != -1 or tags_habitat[i].find(tag) != -1) and not flag :
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
            flag2 = True
            if tag_name.find(name) != -1:
                for tag in alive:
                    if tags_habitat[j] != "" and flag2 and (tags_habitat[j] not in list_not_to_dis):
                        print("%s was not tagged, tag = %s" % (str(name),str(tags_habitat[j])))
                        flag2 = False
                for tag_env in env:
                    if tags_habitat[j] != "" and flag2 and (tags_habitat[j] not in list_not_to_dis):
                        print("%s was not tagged, tag = %s" % (str(name), str(tags_habitat[j])))
                        flag2 = False
                        flag = True
            j = j + 1
        not_anottated.append(name)
with open(path+"/0/colors0.0.txt", "a") as myfile:
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
