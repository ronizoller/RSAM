path = '/Users/ronizoller/Google Drive (ronizo@post.bgu.ac.il)/COGS/COG2602/'
list1  = ['1144888.CM001467 gene305', '1232452.BAIB02000016 gene2604', '272951.rsib orf860', '1297569.MESS2 510063', '63737.Npun F3979', '452659.RrIowa 1453', '1319815.HMPREF0202 02107', '1123252.ATZF01000005 gene4030', '1235802.C823 04230', '1265503.KB905173 gene4231', '1340434.AXVA01000018 gene854', '1340434.AXVA01000018 gene853', '748247.AZKH 1955', '1442598.JABW01000001 gene1818', '452659.RrIowa 1454', '518637.EUBIFOR 00819', '411461.DORFOR 02711', '616991.JPOO01000003 gene1901', '1232452.BAIB02000016 gene2605', '1442598.JABW01000001 gene1817', '334545.CCMG01000009 gene781', '1192868.CAIU01000023 gene3319', '1382305.AZUC01000063 gene3300', '1123252.ATZF01000005 gene4031', '616991.JPOO01000003 gene1902', '686578.AFFX01000001 gene1409', '1038866.KB902813 gene2158', '1279017.AQYJ01000029 gene3878', '272951.rsib orf859', '32049.SYNPCC7002 A2334', '575586.HMPREF0016 02934', '1265503.KB905173 gene4229', '1265503.KB905168 gene1329', '1123008.KB905696 gene3083', '28229.ND2E 4012', '1265503.KB905173 gene4230']

gene_list = True
ext = 'bacteria_u569_2'
import os
pattern = ''
save_aligned = False

def main(lists,path,ext,nd,pattern,gene_list):
    if pattern == '':
        for list, vertex_name in lists:
            with open(path + '/FASTA_bacteria.txt', 'r') as fp:
                input = open(path + "/old_new_names.txt", 'r')
                old_new_names = []
                for line in input:
                    old_new_names.append(eval(line))
                old_new_names = old_new_names[0]
                flag = False
                res = ''
                for line in fp:
                    if line[0] == '>':
                        flag = False
                        name = line[line.find('[') + 1:line.find(']')]
                        gene_name = line[1:line.find('[')].replace('_', ' ')
                        for leaf in list:
                            if gene_list:
                                if leaf + ' ' == gene_name or leaf == gene_name:
                                    flag = True
                            else:
                                if leaf + ' ' == old_new_names[name].replace(' ', '') or leaf == old_new_names[name].replace(' ', ''):
                                    flag = True
                    if flag:
                        res = res + line
                path_curr_to_save = path + '/FASTA_' + ext + '.txt'
                os.makedirs(os.path.dirname(path_curr_to_save), exist_ok=True)
                file = open(path_curr_to_save, 'w')
                file.write(str(res))
                file.close()
    else:
        for list,vertex_name in lists:
            with open(path+'/FASTA_'+ext+'.txt','r') as fp:
                input = open(path + "/old_new_names.txt", 'r')
                old_new_names = []
                for line in input:
                    old_new_names.append(eval(line))
                old_new_names = old_new_names[0]
                flag = False
                res = ''
                for line in fp:
                    if line[0] == '>':
                        flag = False
                        name = line[line.find('[') + 1:line.find(']')]
                        gene_name = line[1:line.find('[')].replace('_',' ')
                        for leaf in list:
                            if gene_list:
                                if leaf+' ' == gene_name or leaf == gene_name:
                                    flag = True
                            else:
                                if leaf+' ' == old_new_names[name].replace(' ', '') or leaf == old_new_names[name].replace(' ', ''):
                                    flag = True
                    if flag:
                        res = res + line
                path_curr_to_save = path + '/saved_data/results/' + pattern + '/FASTA_result_' + nd + '.txt'
                os.makedirs(os.path.dirname(path_curr_to_save), exist_ok=True)
                file = open(path_curr_to_save, 'w')
                file.write(str(res))
                file.close()
            if save_aligned:
                with open(path+'/FASTA_aligned_'+ext+'.txt','r') as fp:
                    input = open(path + "/old_new_names.txt", 'r')
                    old_new_names = []
                    for line in input:
                        old_new_names.append(eval(line))
                    old_new_names = old_new_names[0]
                    flag = False
                    res = ''
                    for line in fp:
                        if line[0] == '>':
                            flag = False
                            name = line[line.find('[') + 1:line.find(']')]
                            gene_name = line[1:line.find('[')].replace('_',' ')
                            for leaf in list:
                                if gene_list:
                                    if leaf+' ' == gene_name or leaf == gene_name:
                                        flag = True
                                else:
                                    if leaf+' ' == old_new_names[name].replace(' ', '') or leaf == old_new_names[name].replace(' ', ''):
                                        flag = True
                        if flag:
                            res = res + line
                path_curr_to_save = path + '/saved_data/results/' + pattern + '/FASTA_aligned_result_for_pattern_' +  pattern + '_' + ext + '_' + nd + '.txt'
                os.makedirs(os.path.dirname(path_curr_to_save), exist_ok=True)
                file = open(path_curr_to_save, 'w')
                file.write(str(res))
                file.close()

if __name__ == "__main__":
    main([(list1,'_'+ext)],path,ext,'',pattern,gene_list)
