path = '/Users/ronizoller/PycharmProjects/TreeReconciliation/trees/real/'
list1 = ['Legionelladrancourtii', 'Desulfovibriomagneticus', 'Legionellalongbeachae', 'Legionellapneumophilastr.Paris', 'Legionellapneumophilastr.Paris', 'Legionellapneumophilastr.Paris', 'Rickettsiasibirica', 'Rickettsiabellii', 'Rickettsiasibirica', 'Rickettsiafelis']
list2 = ['Aeromonashydrophila', 'Aeromonasveronii', 'Aeromonassalmonicida', 'Aeromonascaviae', 'Agrobacteriumvitis', 'Burkholderiapseudomallei', 'Ralstoniapickettii', 'Burkholderiacenocepacia', 'Burkholderiagladioli', 'Burkholderiaglumae', 'Delftiaacidovorans', 'Achromobacterpiechaudii', 'Herbaspirillumseropedicae', 'Rhodospirillumcentenum', 'AgrobacteriumtumefaciensF2', 'AgrobacteriumtumefaciensF2', 'AgrobacteriumtumefaciensF2', 'gammaproteobacteriumHdN1', 'Mesorhizobiumloti', 'Mesorhizobiumopportunistum', 'Mesorhizobiumciceri', 'Xanthobacterautotrophicus']
lists = [(list1,'u401'),(list2,'u302')]

for list,name in lists:
    with open(path+'/FASTA_aligned.txt','r') as fp:
        flag = False
        right = ''

        for line in fp:
            if line[0] == '>':
                flag = False
                name = line[line.find('[')+1:line.find(']')].replace(' ','')
                for leaf in list2:
                    if leaf == name:
                        flag = True
            if flag:
                right = right + line

with open(path+'/FASTA_aligned.txt','r') as fp:
    flag = False
    left = ''

    for line in fp:
        if line[0] == '>':
            flag = False
            name = line[line.find('[')+1:line.find(']')].replace(' ','')
            for leaf in list1:
                if leaf == name:
                    flag = True
        if flag:
            left = left + line

with open(path+"/FASTA_aligned"+list1_name+".txt", "a") as myfile:
    myfile.write(right)
with open(path+"/FASTA_aligned"+list2_name+".txt", "a") as myfile:
    myfile.write(left)