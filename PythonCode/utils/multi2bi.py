import dendropy as tr
import random
import tree_operations_v1 as tree_operations


def split_randomly_into_two_sets(list_of_child):
    left_tree = tr.Tree()
    right_tree = tr.Tree()
    left_subtree_size = random.choice(range(0,len(list_of_child)))
    for i in range(0,len(list_of_child)):
        if i < left_subtree_size:
            left_tree.seed_node.add_child(list_of_child[i])
        else:
            right_tree.seed_node.add_child(list_of_child[i])
    if left_subtree_size == 0:
        left_tree = None
    if left_subtree_size == len(list_of_child):
        right_tree = None
    return tree_operations.collapse_edges(left_tree),tree_operations.collapse_edges(right_tree)


def rec_multi2bi(non_binary_tree):
    if not non_binary_tree or tree_operations.is_a_leaf(non_binary_tree.seed_node):
        return non_binary_tree

    seed = non_binary_tree.seed_node
    left_binary, right_binary = split_randomly_into_two_sets(seed.child_nodes())
    left_binary_tree = rec_multi2bi(left_binary)
    right_binary_tree = rec_multi2bi(right_binary)
    new_bin = tr.Tree()
    if seed.taxon or seed.label:
        if seed.taxon:
            new_bin.seed_node.taxon = non_binary_tree.taxon
        else:
            new_bin.seed_node.label = seed.label
    if (left_binary_tree and (left_binary_tree.seed_node.taxon or not tree_operations.is_a_leaf(left_binary_tree.seed_node))) and\
            (right_binary_tree and (right_binary_tree.seed_node.taxon or not tree_operations.is_a_leaf(right_binary_tree.seed_node))):
        new_bin.seed_node.set_child_nodes([left_binary_tree.seed_node,right_binary_tree.seed_node])
    elif not left_binary_tree:
        new_bin.seed_node.set_child_nodes([right_binary_tree.seed_node])
    else:
        new_bin.seed_node.set_child_nodes([left_binary_tree.seed_node])
    return new_bin


def main(non_binary_tree):
    binary_tree = rec_multi2bi(non_binary_tree)
    return binary_tree