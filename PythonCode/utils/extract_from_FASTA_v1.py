path = '/Users/ronizoller/PycharmProjects/TreeReconciliation/trees/COG2602'
list1 = ['Nocardiopsiskunsanensis', 'Nocardiopsisxinjiangensis', 'Nocardiopsisdassonvillei', 'Nocardiopsissp.CNT312', 'Nocardiopsisalba', 'Nocardiopsisprasina', 'Nocardiopsisganjiahuensis', 'Nocardiopsisvalliformis', 'VariovoraxparadoxusEPS', 'Verrucomicrobiumsp.BvORR034', 'Verrucomicrobiumsp.BvORR106', 'Ralstoniapickettii12J', 'Parvibaculumlavamentivorans', 'Synechococcussp.RS9917', 'Synechococcussp.CB0101', 'Synechococcussp.CB0101', 'Massiliatimonae', 'Pseudomonasaeruginosa', 'Pseudomonastaeanensis', 'Hahellaganghwensis', 'Methylobacillusflagellatus', 'Methylococcuscapsulatus', 'Porticoccushydrocarbonoclasticus', 'Ruegeriapomeroyi', 'Actinomaduraflavalba', 'Actinomaduraatramentaria', 'Spirillosporaalbida', 'Thermomonosporacurvata', 'Streptosporangiumroseum', 'Haloferulasp.BvORR071', 'Sulfurimonasdenitrificans']
list2 = ['Acinetobactersp.Ver3', 'Acinetobacterjunii', 'Acinetobactertandoii', 'Acinetobactergerneri', 'Acinetobacterbohemicus', 'Acinetobacterguillouiae', 'Acinetobacterbereziniae', 'Acinetobacterbaumannii', 'Acinetobactercalcoaceticus', 'Acinetobacteroleivorans', 'Acinetobacterbaumannii', 'Acinetobacterbaumannii', 'Acinetobacterradioresistens', 'Acinetobacterrudis', 'AcinetobacterlwoffiiSH145', 'Acinetobactervenetianus', 'Acinetobactersp.NCTC7422', 'Acinetobactergyllenbergii', 'Acinetobacterhaemolyticus', 'Commensalibactersp.MX01', 'Acinetobacterbouvetii', 'Helicobacterpametensis', 'Succinatimonashippei', 'Campylobacterlari', 'Shewanellamarina', 'Cetobacteriumsomerae', 'Brachyspirapilosicoli']

lists = [(list1,'VS_ACIN'),(list2,'ACIN')]
for list,vertex_name in lists:
    with open(path+'/FASTA_a.txt','r') as fp:
        flag = False
        res = ''
        for line in fp:
            if line[0] == '>':
                flag = False
                name = line[line.find('[')+1:line.find(']')].replace(' ','')
                for leaf in list:
                    if leaf == name:
                        flag = True
            if flag:
                res = res + line
    with open(path+"/FASTA_aligned"+vertex_name+".txt", "a") as myfile:
        myfile.write(res)
