path_source = '/Users/ronizoller/Google Drive (ronizo@post.bgu.ac.il)/COGS/COG3549/'
path_target = '/Users/ronizoller/Google Drive (ronizo@post.bgu.ac.il)/COGS/COG3093/'
old_genes  = ['1265502.KB905931 gene1593', '1265502.KB906002 gene3329', '1223521.BBJX01000010 gene79', '535289.Dtpsy 3235', '1223521.BBJX01000010 gene148', '596153.Alide 0472', '94624.Bpet4638', '216591.BCAL1127', '375286.mma 2741', '795666.MW7 1022', '228410.NE2538', '426114.THI 1922', '269482.Bcep1808 3274', '640512.BC1003 6099', '864073.HFRIS 020045', '522306.CAP2UW1 3865', '420662.Mpe A2374']

import os

input = open(path_source+'/0/sigma0.0.txt', 'r')
sigma = []
for line in input:
    sigma.append(eval(line))
sigma = sigma[0]
input = open(path_source + "/old_new_names.txt", 'r')
old_new_names = []
for line in input:
    old_new_names.append(eval(line))
old_new_names = old_new_names[0]

species_genes = {}
species = []
res = {}
found = {}

def find_gene_number (gene):
    gene_number = ''
    i = len(gene) - 1
    while gene[i].isdigit():
        gene_number += gene[i]
    return int(gene_number)

for gene in old_genes:
    species_genes.update({gene:sigma[gene]})
    species.append(sigma[gene])

with open(path_target + '/FASTA.txt', 'r') as fp:
    flag = False
    for line in fp:
        if line[0] == '>':
            flag = False
            name = line[line.find('[') + 1:line.find(']')]
            gene_name = line[1:line.find('[')].replace('_', ' ')
            if name in species:
                for gene,spec in species_genes:
                    old_gene_number = find_gene_number(gene)
                    new_gene_number = find_gene_number(gene_name)
                    if old_gene_number + 1 == new_gene_number or old_gene_number -1 == new_gene_number:
                        res.update({spec:gene_name})
                        found.update({gene_name:gene})
print(res)
print(found)