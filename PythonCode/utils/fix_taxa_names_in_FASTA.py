path = '/Users/ronizoller/Google Drive (ronizo@post.bgu.ac.il)/COGS/'
COGS_names = ['COG1396']
import os

def main(path,create_sigma_from_fasta,res):
    if create_sigma_from_fasta:
        if not os.path.isfile(path + '/old_new_names.txt'):
            ids = ''
            ids_species_list = {}
            names_map = {}
            with open(path+'/taxa_names.txt','r') as file1:
                for prev_name in file1:
                    flag = False
                    prev_name = prev_name[:len(prev_name) - 1]
                    if prev_name != '':
                        with open(os.getcwd()+'/utils/name_to_ID_from_STRING.txt','r') as fp:
                            for line in fp:
                                if not flag:
                                    line = line.replace('	', "|")
                                    index = 0
                                    i = 0
                                    while line[index].isdigit() and index < len(line):
                                        index += 1
                                    id_end = index
                                    while i < 2 and index < len(line):
                                        if line[index] == '|':
                                            i += 1
                                        index += 1
                                    old_name_index_start = index
                                    while i < 3 and index < len(line):
                                        if line[index] == '|':
                                            i += 1
                                        index += 1
                                    old_name_index_end = index
                                    old_name = line[old_name_index_start:old_name_index_end-1]
                                    new_name = line[old_name_index_end:len(line)-1]
                                    if prev_name == old_name:
                                        ids += line[:id_end]+'\n'
                                        ids_species_list.update({line[:id_end]:new_name})
                                        names_map.update({old_name:new_name})
                                        flag = True
                            if not flag:
                                res['error'] += ('There is no mapping for speice '+str(prev_name))
            file = open(path + "/NCBI_tax_ID.txt", 'w')
            file.write(str(ids))
            file.close()
            file = open(path + "/old_new_names.txt", 'w')
            file.write(str(names_map))
            file.close()
    else:
        file = open(path + "/NCBI_tax_ID.txt", 'w')
        file.write(str(''))
        file.close()
        file = open(path + "/old_new_names.txt", 'w')
        file.write(str(''))
        file.close()
    return ids_species_list

if __name__ == "__main__":
    main(path,COGS_names)