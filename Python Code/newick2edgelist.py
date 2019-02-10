import sys
sys.path.append('/anaconda3/lib/python3.6/site-packages')
import dendropy as tr
path = '/Users/ronizoller/PycharmProjects/TreeReconciliation/trees/Simulator'
exte = 'all'

def get_edgelist(tree):
    edgelist = []
    edges = tree.edges()
    for e in edges:
        if e.head_node != None and e.tail_node != None:
            edgelist.append(e.head_node.label+' '+e.tail_node.label+' '+str(e.length))
    return edgelist


def init_internal_labels (tree,char):
    print('Inisilasing internal leafs...')
    counter = 1
    for nd in tree.postorder_node_iter():
        nd.label = char+str(counter)
        counter += 1
    print('Finished inisilasing internal leafs.\n')
    return tree

def main():
    t = tr.Tree.get_from_path(path+"/phyliptree(binary,"+exte+").phy", schema="newick")
    t = init_internal_labels(t,'x')
    file = open(path+'/saved_data/S_edgelist_'+exte+'.txt', 'w')
    file.write(str(get_edgelist(t)))
    file.close()

if __name__ == "__main__":
    main()