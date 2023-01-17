# Partition size and tree specification in PolyMoSim:

Our goal is to simulate the evolution of nucleotide or amino acid sequences.
In this file we specify the sizes of the alignment partitions, the root models and the trees used for the simulation.

Each line in the tree file specifies the following information:
- Scalar with which all branch lengths are multiplied. Typically, this has a value of 1.0. 
- Length of the partition.
- Root model. This is the evolutionary model assumed at the root of the tree.
- Tree sting in Newick format including branch lengths.

Example:
```
1.0  1000 Model_JC ((A:0.1,B:0.9):0.005,(C:0.1,D:0.9):0.005)
```

## The root model

In the above example, the tree with terminal taxa A, B, C, D will be used to simulate the evolution. The root sequence will be a random nucleotide sequences (since the root model is a nucleotide model) with nucleotide frequencies specified by the root model. Branches can evolve on models other than the root model if specified in the tree string. If no other models are specified for branches, the root model is used to simulate evolution on the whole tree.

Note: The data type of the simulated alignment partition is determined by the data type of the root model. The data type can differ among different partitions.

## Two partitions
```
1.0  10000 Model_JC    ((A:0.2,B:0.3):0.3,(C:0.1,D:0.5):0.4)
1.0  1000  Model_K2P_2 ((A:0.1,C:0.1):0.1,(B:0.1,D:0.1):0.1)
```

The two partitions have different lengths, different root models and evolved on different trees.

## Different models on different lineages
```
1.0  20000 Modell_JC ((A:0.8[&Modell_F81],B:0.8):0.01,(C:0.8,D:0.8):0.01[&Modell_F84])
```

Models can be specified for branches and clades. A string [&Modelname] after the branch length can be used to specify that beginning with this branch, the specified model should be used for the evolution on this branch. The new model will be used starting at the beginning of the branch for which it was specified. This model will be used up to the tips of the tree. In this example the Model Modell_F81 will be used for the branch leading to the taxon A and the model Modell_F84 will be used for the subtree (C:0.8,D:0.8):0.01.

## Different and same site rates

Each evolutionary model has its own site rates associated to it. Site rates can be inherited from other models.
By default, each model with site rate heterogeneity has newly determined invariant sites and site rates.
```
1.0  20000 Modell_JC_shape0.1_a ((A:0.8[&Modell_JC_shape0.1_b],B:0.8):0.01,(C:0.8[&Modell_JC_shape0.1_b],D:0.8):0.01)
```

Here the models Modell_JC_shape0.1_a and Modell_JC_shape0.1_b are in principle the same. They have identical site rate distributions, but in same site have different rates under the two models. The branches leading to A and C have the same site rates associated to the sites, but they differ from the site rates on all other branches. An exception is only made of the Modell_JC_shape0.1_a and Modell_JC_shape0.1_b have been explicitly defined to have the same site rates. 

Basically, branches leading to A and C share the same evolutions hotspots along the sequences.

The easiest way to think about this is: Each model comes with its own list of site rates. 

