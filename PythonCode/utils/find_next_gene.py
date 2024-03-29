import extract_from_FASTA_v1 as extr


path_source = '/Users/ronizoller/PycharmProjects/RSAM/PythonCode/data/COG2602/'
path_target = '/Users/ronizoller/PycharmProjects/RSAM/PythonCode/data/COG3682_BlaI'
pattern = "((['HT', 'None', 'None'], 'red', 'True')_(['S', 'D', 'HT'], 'black', 'False'))_Double-Mode"
range_to_search = range(1,2)
old_genes  = []
file_to_examine = "FASTA_result_COG2602_1_node_u841_Double-mode,u840_p1_left.txt"
all = False

if all:
    path = path_source + 'FASTA.txt'
else:
    path = path_source + "/saved_data/results/"+pattern+'/'+file_to_examine

with open(path , 'r') as fp:
    for line in fp:
        if line[0] == '>':
            flag = False
            name = line[line.find('[') + 1:line.find(']')]
            gene_name = line[1:line.find('[')]
            if gene_name[len(gene_name)-1] == ' ':
                gene_name = gene_name [:len(gene_name)-1]
            old_genes.append(gene_name)
input = open(path_source+'/sigma.txt', 'r')
sigma = []
for line in input:
    sigma.append(eval(line))
sigma = sigma[0]
input = open(path_source + "/old_new_names.txt", 'r')
old_new_names = []
for line in input:
    old_new_names.append(eval(line))
old_new_names = old_new_names[0]
old_new_names_no_spaces = {}

for name,new_name in old_new_names.items():
    old_new_names_no_spaces.update({name.replace(' ',''):new_name})
species_genes = {}
species = []
res = []
found = {}

def find_gene_number (gene):
    gene_number = ''
    i = len(gene) - 1
    if gene[i] == ' ':
        i = i-1
    while gene[i].isdigit():
        gene_number += gene[i]
        i -= 1
    pref = ''
    while i >= 0:
        pref += gene[i]
        i -= 1
    gene_number = gene_number[::-1]
    if gene_number != '':
        return int(gene_number),pref[::-1]
    else:
        return -1,pref[::-1]

for gene in old_genes:
    if gene.replace('_',' ') in sigma:
        species_genes.update({gene:sigma[gene.replace('_',' ')]})
        species.append(sigma[gene.replace('_',' ')])

with open(path_target + '/FASTA.txt', 'r') as fp:
    flag = False
    for line in fp:
        if line[0] == '>':
            name = line[line.find('[') + 1:line.find(']')].replace(' ','')
            gene_name = line[1:line.find('[')].replace('_', ' ')
            new_gene_number, pref_new = find_gene_number(gene_name)
            if name in old_new_names_no_spaces:
                if old_new_names_no_spaces[name].replace(' ','') in species:
                    for gene,spec in species_genes.items():
                        old_gene_number,pref_old = find_gene_number(gene)
                        if pref_old == pref_new.replace(' ','_'):
                            for j in range_to_search:
                                if old_gene_number + j == new_gene_number or old_gene_number - j == new_gene_number:
                                    res.append(gene_name.replace('_',' '))
                                    found.update({gene:gene_name})
extr.main([(res, '_bacteria')], path_target, 'bacteria', '', pattern, True)
list_of_found_species = []
number_of_found = 0
for old_gene in old_genes:
    if old_gene.replace('_',' ') in sigma:
        flag = False
        for old,new in found.items():
            if old == old_gene:
                number_of_found += 1
                flag = True
        #if not flag:
        #    print('gene %s for specie %s was not found' % (str(old_gene),str(species_genes[old_gene])))
        if flag:
            list_of_found_species.append(old_gene)
            print('gene %s for specie %s was found' % (str(old_gene),str(species_genes[old_gene])))

    else:
        print(old_gene+' is not in sigma')
print(list_of_found_species)

print('\nfound %s/%s' % (str(number_of_found),str(len(old_genes))))