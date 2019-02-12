path = '/Users/ronizoller/PycharmProjects/TreeReconciliation/trees/real/'
list1 = ['Sulfurimonasdenitrificans', 'Synechococcussp.RS9917', 'Synechococcussp.CB0101', 'Synechococcussp.CB0101', 'Methylobacillusflagellatus', 'Methylococcuscapsulatus', 'Pseudomonasaeruginosa', 'Ruegeriapomeroyi', 'Ralstoniapickettii', 'VariovoraxparadoxusEPS', 'Parvibaculumlavamentivorans', 'Nocardiopsisdassonvillei', 'Thermomonosporacurvata', 'Streptosporangiumroseum', 'Succinatimonashippei', 'Brachyspirapilosicoli', 'Campylobacterlari', 'Campylobacterlari']
list2 = ['Acinetobacterbaumannii', 'Acinetobacterbaumannii', 'Acinetobacteroleivorans', 'Acinetobactercalcoaceticus', 'Acinetobacterradioresistens', 'Acinetobacterlwoffii']

lists = [(list1,'u62'),(list2,'u97')]
res = []
i = 0

for list,list_name in lists:
    with open(path+'/FASTA_aligned.txt','r') as fp:
        flag = False
        res.append('')
        res_temp = res[i]
        for line in fp:
            if line[0] == '>':
                flag = False
                name = line[line.find('[')+1:line.find(']')].replace(' ','')
                for leaf in list:
                    if leaf == name:
                        flag = True
            if flag:
                res_temp = res_temp + line
    with open(path + "/FASTA_aligned" + list_name + ".txt", "a") as myfile:
        myfile.write(res_temp)
    i += 1

#with open(path+'/FASTA_aligned.txt','r') as fp:
#    flag = False
#    left = ''#

#    for line in fp:
#        if line[0] == '>':
#            flag = False
#            name = line[line.find('[')+1:line.find(']')].replace(' ','')
#            for leaf in list1:
#                if leaf == name:
#                    flag = True
#        if flag:
#            left = left + line


#with open(path+"/FASTA_aligned"+list2_name+".txt", "a") as myfile:
#    myfile.write(left)