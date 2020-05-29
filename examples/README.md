# Examples

## #1 OpenBioLink Network
Given a network, in this case the OpenBioLink network, we defined a set of source nodes acting as drugs and a set of
target nodes acting as indications/phenotypes. Each of these two types of nodes has to be passed to Drug2Ways as two
separate files where each row contains the name of a drug/indication. For example:

```
source_node_1
source_node_2
source_node_3
source_node_4
...
```

**Optimization and combination commands**: To optimize effects on multiple targets, a slightly different format for the
target nodes file is required. In this file, a second column is required (separated from the target node name with a
comma) where the desired effect is specified (i.e., 1/activate, and -1/inhibit). See the example below:

```
target_node_1,1
target_node_2,-1
target_node_3,1
target_node_4,-1
...
```

or the equivalent

```
target_node_1,activate
target_node_2,inhibit
target_node_3,activate
target_node_4,inhibit
...
```

### Data

- Network: https://github.com/drug2ways/drug2ways/blob/master/data/networks/data/openbiolink_network.tsv
- Drugs: https://github.com/drug2ways/drug2ways/blob/master/data/validation/data/source_nodes_openbiolink.tsv
- Indications/Phenotypes:
  https://github.com/drug2ways/drug2ways/blob/master/data/validation/data/target_nodes_openbiolink.tsv

This example is outlined in the folowing script:
- https://github.com/drug2ways/drug2ways/blob/master/examples/openbiolink_network_run.sh

## #2 In-House Network
Similarly to the OpenBioLink network, the same three files are required for the In-House network:

- Network: https://github.com/drug2ways/drug2ways/blob/master/data/networks/data/custom_network.tsv
- Drugs: https://github.com/drug2ways/drug2ways/blob/master/data/validation/data/source_nodes_custom.tsv
- Indications/Phenotypes:
  https://github.com/drug2ways/drug2ways/blob/master/data/validation/data/target_nodes_custom.tsv

This example is outlined in the folowing script:
- https://github.com/drug2ways/drug2ways/blob/master/examples/custom_network_run.sh

### References
1. Breit, A., *et al* (2020). OpenBioLink: A resource and benchmarking framework for large-scale
   biomedical link prediction. Bioinformatics, btaa274, https://doi.org/10.1093/bioinformatics/btaa274
 
2. Domingo-Fernández, *et al* (2019). PathMe: Merging and exploring mechanistic pathway knowledge. BMC Bioinformatics,
   20:243. https://doi.org/10.1186/s12859-019-2863-9.
