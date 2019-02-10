import sys
sys.path.append('/anaconda3/lib/python3.6/site-packages')
from ete3 import Tree
import random
from random import randint
import utils.newick2edgelist
import networkx as nx
import dendropy as tr
import tree_operations_v1 as tree_operations
import utiles
import inits_v1 as inits
import draw
import os
import errno

S = Tree()

S.populate(10)

print(S)