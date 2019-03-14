path = '/Users/ronizoller/PycharmProjects/TreeReconciliation/trees/real/'
list1 = ['Flexithrixdorotheae', 'VariovoraxparadoxusEPS', 'Nostocpunctiforme', 'Fischerellamuscicola', 'Spirulinasubsalsa', 'Oscillatoriaacuminata', 'Nodosilineanodulosa', 'Sediminibacteriumsp.C3', 'Actinomaduraflavalba', 'Microcoleusvaginatus', 'Nostocsp.PCC7120', 'Anabaenacylindrica', 'VerrucomicrobiaebacteriumDG1235', 'Nocardiopsiskunsanensis', 'Nocardiopsisxinjiangensis', 'Anabaenasp.PCC7108', 'Emticiciaoligotrophica', 'Arenimonasoryziterrae', 'Mesorhizobiummetallidurans', 'Myxococcusstipitatus', 'Actinomaduraatramentaria']
list2 = ['Bacillusokhensis', 'Parvibaculumlavamentivorans', 'Bacillussp.J33', 'Acinetobacterradioresistens', 'Acinetobacterbaumannii', 'Acinetobacterbaumannii', 'Sulfurimonasdenitrificans', 'Pseudomonasaeruginosa', 'Pseudomonastaeanensis', 'Nocardiopsisvalliformis', 'Nocardiopsisprasina', 'Nocardiopsisalba', 'Nocardiopsisganjiahuensis', 'Methylococcuscapsulatus', 'Porticoccushydrocarbonoclasticus', 'Methylobacillusflagellatus', 'Hahellaganghwensis', 'Nocardiopsisdassonvillei', 'Cetobacteriumsomerae', 'Thermomonosporacurvata', 'Synechococcussp.CB0101', 'Synechococcussp.RS9917', 'Synechococcussp.CB0101', 'Massiliatimonae', 'Haloferulasp.BvORR071', 'Cetobacteriumsomerae', 'Cupriavidustaiwanensis', 'RalstoniapickettiiDTP0602', 'Shewanellapiezotolerans', 'Enterococcuspallens', 'Brachyspirapilosicoli', 'Chondromycesapiculatus', 'Cupriavidussp.amp6', 'Shewanellamarina', 'Succinatimonashippei']

lists = [(list1,'u409'),(list2,'u450')]
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