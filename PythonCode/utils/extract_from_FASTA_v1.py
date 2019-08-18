path = '/Users/ronizoller/PycharmProjects/RSAM/PythonCode/data/COG3093'
list1 = ['323261.Noc 2763 ', '160488.PP 1585 ', '351348.Maqu 3198 ', '1055815.AYYA01000030 gene759 ', '1007105.PT7 0143 ', '398578.Daci 1211 ', '545264.KB898746 gene880 ', '396595.TK90 0224 ', '225937.HP15 4163 ', '351746.Pput 4192 ', '1112209.AHVZ01000006 gene1775 ', '472759.Nhal 1146 ', '863365.XHC 2235 ', '1255043.TVNIR 0522 ']

gene_list = True
ext = '_res'
import os
pattern = ''
save_aligned = False


def main(lists, path, ext, nd, pattern, gene_list):
    if pattern == '':
        for lst, vertex_name in lists:
            with open(path + '/FASTA.txt', 'r') as fp:
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
                        for leaf in lst:
                            if gene_list:

                                if leaf + ' ' == gene_name or leaf == gene_name or\
                                                        leaf + ' ' == gene_name.replace(' ','_') or\
                                                leaf == gene_name.replace(' ',''):
                                    flag = True
                                    print(leaf)
                                    print(gene_name)
                                    print('')
                            else:
                                if leaf + ' ' == old_new_names[name].replace(' ', '') or leaf == old_new_names[name].replace(' ', '') or\
                                        leaf + ' ' == old_new_names[name] or leaf == old_new_names[name]:
                                    flag = True
                    if flag:
                        res += line
                path_curr_to_save = path + '/FASTA_' + ext + '.txt'
                os.makedirs(os.path.dirname(path_curr_to_save), exist_ok=True)
                file = open(path_curr_to_save, 'w')
                file.write(str(res))
                file.close()
    else:
        for list,vertex_name in lists:
            with open(path+'/FASTA.txt','r') as fp:
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
                                if leaf+' ' == gene_name or leaf == gene_name :
                                    flag = True
                            else:
                                if leaf+' ' == old_new_names[name].replace(' ', '') or leaf == old_new_names[name].replace(' ', ''):
                                    flag = True
                    if flag:
                        res += line
                path_curr_to_save = path + '/saved_data/results/' + pattern + '/FASTA_result_' + ext + '_' + nd + '.txt'
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
