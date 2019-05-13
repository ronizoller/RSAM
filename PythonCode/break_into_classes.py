import dendropy as tr
import tree_operations_v1 as tree_operations
from utils import extract_from_FASTA_v1 as extr

path = '/Users/ronizoller/Google Drive (ronizo@post.bgu.ac.il)/COGS/'

COGS = ['COG3550']
classes = ['proteobacteria']
for COG in COGS:
    for clas in classes:
        S = tr.Tree.get_from_path(path + '/' + COG + "/phyliptree(binary," + clas + ").phy", schema="newick")
        input = open(path + '/' + COG + '/0/sigma0.0.txt', 'r')
        sigma = []
        for line in input:
            sigma.append(eval(line))
        sigma = sigma[0]
        list1 = tree_operations.leaf_in_subtrees(S, 'S', S.seed_node.label , sigma, False)[0]+tree_operations.leaf_in_subtrees(S, 'S', S.seed_node.label , sigma, False)[1]

        print('COG: %s, class: %s, \n%s\n' % (str(COG),str(clas),str(list1)))

        extr.main([(list1,'_'+clas)],path+'/'+COG,clas,'','',False)
