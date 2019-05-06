path = '/Users/ronizoller/Google Drive (ronizo@post.bgu.ac.il)/COGS/COG2856/'
list1  = ['1265502.KB905931 gene1593', '1265502.KB906002 gene3329', '1223521.BBJX01000010 gene79', '535289.Dtpsy 3235', '1223521.BBJX01000010 gene148', '596153.Alide 0472', '94624.Bpet4638', '216591.BCAL1127', '375286.mma 2741', '795666.MW7 1022', '228410.NE2538', '426114.THI 1922', '269482.Bcep1808 3274', '640512.BC1003 6099', '864073.HFRIS 020045', '522306.CAP2UW1 3865', '420662.Mpe A2374']
gene_list = True
ext = 'beta'
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
                path_curr_to_save = path + '/saved_data/results/' + pattern + '/FASTA_result_for_pattern_' + pattern + '_' + ext + '_' + nd + '.txt'
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
