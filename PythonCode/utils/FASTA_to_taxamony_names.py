cogs_names = ['COG1396']
import os
import tree_operations_v1 as tree_operations
import dendropy as tr


create_sigma = False
path = '/Users/ronizoller/Google Drive (ronizo@post.bgu.ac.il)/COGS/'
def main(path,cogs_names,exte,create_sigma):
    if create_sigma:
        for name in cogs_names:
            with open(path+name+'/FASTA.txt','r') as fp:
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
            file = open(path + name + "/taxa_names.txt", 'w')
            file.write(str(res))
            file.close()
    else:
        S = tr.Tree.get_from_path(os.getcwd() + "/data/" + exte +"/S.txt", schema="newick")
        res = ''
        for leaf in S.leaf_nodes():
            if not tree_operations.isolated(leaf):
                res += leaf.taxon.label + ','
        file = open(path + '' + "/taxa_names.txt", 'w')
        file.write(str(res))
        file.close()
if __name__ == "__main__":
    main(path,cogs_names,'',create_sigma)
