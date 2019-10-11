cogs_names = ['COG1396']
import tree_operations_v1 as tree_operations
import dendropy as tr


create_sigma = False
path = '/Users/ronizoller/Google Drive (ronizo@post.bgu.ac.il)/COGS/'
def main(path,create_sigma,res):
    if create_sigma:
        with open(path+'/FASTA.txt','r') as fp:
            res = ''
            for line in fp:
                if line[0] == '>':
                    start = line.find('[')+1
                    end = line.find(']')
                    while start < end:
                        if line[start] == ' ':
                            res += ' '
                        else:
                            res += line[start]
                        start += 1
                    res += '\n'
        file = open(path + "/taxa_names.txt", 'w')
        file.write(str(res))
        file.close()

    else:
        try:
            S = tr.Tree.get_from_path(path + "/S.txt", schema="newick")
        except:
            res['error'] += "\nSpecies tree '/data/S.txt' was not found.\n" \
                        "In order to create it, please go to https://www.ncbi.nlm.nih.gov/Taxonomy/CommonTree/wwwcmt.cgi\n" \
                        "and follow the instructions:\n" \
                        " -> Choose File -> select the file 'NCBI_tax_ID.txt' from the folder data/your_data \n" \
                        "-> Add from file: -> select the subtree you want (e.i. Proteobacteria) and press Choose ->\n" \
                        "Expand All -> Save As phylip tree.\n " \
                        "Finally you need to rename the file 'S.txt' and move it to the data/your_data/ directory."
            return
        res = ''
        for leaf in S.leaf_nodes():
            if not tree_operations.isolated(leaf):
                res += leaf.taxon.label + ','
        file = open(path + '' + "/taxa_names.txt", 'w')
        file.write(str(res))
        file.close()
if __name__ == "__main__":
    main(path,create_sigma,{})
