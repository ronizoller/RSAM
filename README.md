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

## Running RSAM-finder
In order to run the application, create a directory called /RSAM, and download the /PythonCode directory into it.
Inside /RSAM, the user should create /data directory, in which the input will be provided. For more information about the input format, see next section.
To run the application, open commandline and run the following code:
```
python3 /RSAM/PythonCode/RSAMfinder.py
```
and follow the instructiosn to specify the desiered parameters.
The results for the qeury will be saved inside /data/saved_data.

## Input Format
The input to the RSAM-finder consists of two trees, both in Newick format. The current version is supporting only binary trees. The gene tree should be named "G.txt" and the Species tree should be named "S_class.txt", where class is the Species subclass, e.g. bacteria, proteobacteria etc.
In addition, one should provide a FASTA file of all genes. This will be used in order to create the sigma function.

One can also provide a coloring function, in the following format:
```
{gene_name1 : color1, gene_name2 :color2,…,gene_name_n : color_n}
```
where ![](https://latex.codecogs.com/gif.latex?color_i\in[red,black])
This file should be saved as colors.txt file, and in the same directory. If the user does not want to use color, the file colors.txt should contain the following line:
```
{}
```

## Sample Data Sets
We provide two sample data sets. The directory /data contains the Gene and Species trees for COG3549 and COG2602. In addition, we provide the corresponding Sigma function, which is a mapping between the leafs of the input trees, and the (optional) coloring function of the Specie tree leafs.
We also give the results of the queries discussed in the paper, in a FASTA format.
