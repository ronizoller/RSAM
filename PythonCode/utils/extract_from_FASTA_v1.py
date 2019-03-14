path = '/Users/ronizoller/PycharmProjects/TreeReconciliation/trees/COG2602'
list1 = ['Thiorhodococcusdrewsii', 'Azoarcussp.KH32C', 'Pseudomonassp.BAY1663', 'PseudomonasstutzeriA1501', 'PseudomonasstutzeriRCH2', 'PseudomonasstutzeriNF13', 'PseudomonasstutzeriCCUG29243', 'Sulfurimonasgotlandica', 'Shewanellaputrefaciens', 'Thiocystisviolascens', 'Arcobactercibarius', 'Thiocapsamarina', 'Chlorobaculumparvum', 'Thiorhodovibriosp.970', 'Chlorobiumlimicola', 'alphaproteobacteriumBAL199', 'Algicolasagamiensis', 'Methylomonasmethanica', 'Thiocystisviolascens']
list2 = ['Sulfurospirillumsp.MES', 'Sulfurospirillummultivorans', 'Sulfurospirillumsp.SCADC', 'Arcobactersp.L', 'ArcobacterbutzleriRM4018', 'Helicobacterrodentium', 'Azovibriorestrictus', 'Synechocystissp.PCC6803', 'Synechococcussp.NKBG15041c', 'Thiorhodospirasibirica', 'Cellvibriomixtus', 'Pseudoalteromonasrubra', 'Colwelliapiezophila', 'Cellvibriosp.BR', 'Colwelliapsychrerythraea', 'Pseudomonaspelagia', 'Shewanellahalifaxensis', 'Colwelliapiezophila', 'Arcobactercibarius', 'Zunongwangiaprofunda', 'Idiomarinaloihiensis', 'Endozoicomonasmontiporae', 'Methylophagalonarensis', 'Sphingobacteriumspiritivorum', 'Woodsholeamaritima', 'Spirosomaspitsbergense', 'Spirosomaluteum', 'Sphingobacteriumsp.21', 'Haliscomenobacterhydrossis', 'Anditaleaandensis', 'Chryseobacteriumvrystaatense', 'Flavobacteriumhydatis', 'Pedobactersp.V48', 'Anditaleaandensis', 'Chryseobacteriumluteum', 'Mucilaginibacterpaludis', 'Shewanellaoneidensis', 'Arenibacteralgicola', 'Dyadobactertibetensis', 'Emticiciaoligotrophica', 'Arenibacteralgicola', 'Shewanellapealeana', 'Rheinheimerasp.A13L']

lists = [(list1,'u278'),(list2,'u317')]
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
