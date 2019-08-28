import dendropy as tr
import tree_operations_v1 as tree_operations
from utils import extract_from_FASTA_v1 as extr
import os

path = os.getcwd()+'/data/'

classes = ['_proteobacteria']
COG = 'COG2367(class_A)'

for clas in classes:
    S = tr.Tree.get_from_path(path + '/' + COG + "/S.txt", schema="newick")
    leaves = S.leaf_nodes()

    leaves = list(x.taxon.label for x in leaves)

    print('COG: %s, class: %s, \n%s\n' % (str(COG),str(clas),str(leaves)))

    extr.main([(leaves,'_'+clas)],path+'/'+COG,clas,'','',False)
