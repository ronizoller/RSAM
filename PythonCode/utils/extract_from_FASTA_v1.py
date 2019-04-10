path = '/Users/ronizoller/PycharmProjects/TreeReconciliation/trees/COG3550'
list1 = ['1121451.DESAM 20310', '177437.HRM2 44490', '398767.Glov 3690', '1121918.ARWE01000001 gene3251', '1121918.ARWE01000001 gene3234', '1121918.ARWE01000001 gene3321', '1121918.ARWE01000001 gene3339', '443143.GM18 3582', '941449.dsx2 2013', '1121918.ARWE01000001 gene319', '1232437.KL661963 gene3414', '1265505.ATUG01000002 gene1870', '243231.GSU2399', '398767.Glov 0325', '338966.Ppro 0263', '525897.Dbac 0226', '1121406.JAEX01000002 gene953', '1304885.AUEY01000093 gene1589', '1232410.KI421424 gene1848', '483219.LILAB 08680', '690850.Desaf 1416', '1121441.AUCX01000021 gene700', '641491.DND132 2215']

lists = [(list1,'u137')]
for list,vertex_name in lists:
    with open(path+'/FASTA_a.txt','r') as fp:
        flag = False
        res = ''
        for line in fp:
            if line[0] == '>':
                flag = False
                name = line[line.find('[') + 1:line.find(']')].replace(' ', '')
                gene_name = line[1:line.find('[')].replace('_',' ')
                for leaf in list:
                    if leaf+' ' == gene_name or leaf == gene_name:
                        flag = True
            if flag:
                res = res + line
    with open(path+"/FASTA_aligned"+vertex_name+".txt", "a") as myfile:
        myfile.write(res)
