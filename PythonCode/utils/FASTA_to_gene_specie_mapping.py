import dendropy as tr
import os
import tree_operations_v1 as tree_operations
path  = '/Users/ronizoller/Google Drive (ronizo@post.bgu.ac.il)/COGS/'

create_sigma = False
COGS_names = ['COG2602']
def main (path,create_form_fasta,res_map):
    if create_form_fasta:
        res = '{ '
        try:
            S_text = open(path + '/S.txt', 'r').read()
        except:
            res_map['error'] += "\nSpecies tree '/data/S.txt' was not found.\n" \
                            "In order to create it, please go to https://www.ncbi.nlm.nih.gov/Taxonomy/CommonTree/wwwcmt.cgi\n" \
                            "and follow the instructions:\n" \
                            " -> Choose File -> select the file 'NCBI_tax_ID.txt' from the folder data/your_data \n" \
                            "-> Add from file: -> select the subtree you want (e.i. Proteobacteria) and press Choose ->\n" \
                            "Expand All -> Save As phylip tree.\n " \
                            "Finally you need to rename the file 'S.txt' and move it to the data/your_data/ directory."
            return
        input = open(path+"/old_new_names.txt", 'r')
        old_new_names = []
        for line in input:
            old_new_names.append(eval(line))
        old_new_names = old_new_names[0]
        S_text = S_text.replace(' ','').replace('\n','')
        with open(path+'/FASTA.txt','r') as fp:
            for line in fp:
                if line[0] == '>':
                    start = 1
                    end = line.find(' ')
                    gene = ''
                    while start < end :
                        gene += line[start]
                        start += 1
                    gene = gene.replace("_"," ")
                    start = line.find('[')+1
                    end = line.find(']')
                    specie = ''
                    while start < end :
                        specie += line[start]
                        start += 1
                    if S_text.find(old_new_names[specie].replace(' ','')) != -1:
                        res += "'"+gene+"' : '"+old_new_names[specie].replace(' ','')+"',"

        res += ' }'
        file = open(path + "/sigma.txt", 'w')
        file.write(str(res))
        file.close()
    else:
        try:
            G = tr.Tree.get_from_path(path + "/G.txt", schema="newick")
        except:
            res_map['error'] += "\nGene tree '/data/G.txt' was not found.\n" \
                            "In order to create it, go to https://www.ebi.ac.uk/Tools/msa/clustalo/\n" \
                            "and follow the instructions:\n" \
                            "Choose File -> select 'FASTA.txt' from data/your_data/ directory\n" \
                            "-> select OUTPUT FORMAT: -> wait until job is done..\n" \
                            ""
            return
        res = '{ '
        for leaf in G.leaf_nodes():
            if not tree_operations.isolated(leaf):
                res += " '"+leaf.taxon.label+"' "+':'+" '"+leaf.taxon.label[0:1]+"'"+','
        res += ' }'
        file = open(path + '' + "/sigma.txt", 'w')
        file.write(str(res))
        file.close()
if __name__ == "__main__":
    main(path,create_sigma,{})
