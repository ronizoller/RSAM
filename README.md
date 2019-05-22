# RSAM-finder
RSAM-finder is a python application that identifies predefined patterns within a DLT reconciliation.

## Workflow
Stage 1 constructs a hypergraph H representing the k-best reconciliations of G
and S. Each supernode in H consists of k hypernodes, where each hypernode represents a partial
solution for the DLT-reconciliation problem. Stage 2 consists of assigning a probability to each partial
solution, that is, to each hypernode of H. (This is done in a top down strategy.) In Stage 3, instances
of the sought RSAM-pattern P are identified within H, and RSAM-ranking scores are assigned
accordingly to the vertices of G. Based on these scores, vertices representing putative RSAMs are
identified in G and subjected to biological interpretation.
![](https://user-images.githubusercontent.com/19167301/57970883-8dee4e00-798f-11e9-97fb-446883633860.jpg)

## Input Format
The input to the RSAM-finder consists of two trees, with in Newick format. The current version is supporting only binary trees.
In addition, one should provide a mapping between the species and the genes of the input trees, in the following format: 
```
{gene_name1 : specie_name1, gene_name2 : specie_name2,…,gene_name_n : specie_name_n}
```
This file should be saved as a sigma.txt file in the same directory as the trees.
One can also provide a coloring function, in the following format:
```
{gene_name1 : color1, gene_name2 :color2,…,gene_name_n : color_n}
```
where ![](https://latex.codecogs.com/gif.latex?color_i\in&space;\{&space;red,black&space;\})
This file should also be saved as colors.txt file, and in the same directory.
Note that the colors are optional, and depends on the desired query.

## Sample Data Sets
We provide two data sets we examined in order to demonstrate our tool. Each subfolder contains the Gene and Species trees for COG3549 and COG2602. In addition, we provide the corresponding Sigma function, which is a mapping between the leafs of the input trees, and the (optional) coloring function of the Specie tree leafs.
We also give the results of the queries discussed in the paper, in a FASTA forma.
