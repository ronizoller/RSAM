# RSAM
RSAM-finder is a python application that identifies

## Workflow
Stage 1 constructs a hypergraph H representing the k-best reconciliations of G
and S. Each supernode in H consists of k hypernodes, where each hypernode represents a partial
solution for the DLT-reconciliation problem. Stage 2 consists of assigning a probability to each partial
solution, that is, to each hypernode of H. (This is done in a top down strategy.) In Stage 3, instances
of the sought RSAM-pattern P are identified within H, and RSAM-ranking scores are assigned
accordingly to the vertices of G. Based on these scores, vertices representing putative RSAMs are
identified in G and subjected to biological interpretation.
![](https://user-images.githubusercontent.com/19167301/57970883-8dee4e00-798f-11e9-97fb-446883633860.jpg)
