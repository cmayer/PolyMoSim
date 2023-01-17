# Specifying models for PolyMoSim

PolyMoSim requires model definitions when simulating molecular evolution.
The model file can contain 1 or more model definitions in the following format:

Nucleotide models:
```
begin model      # (required) Indicates the start of a new nucleotide model
name: Model_JC   # (required) Give each model a unique name. A name can be composed of characters,
                 #            numbers and underscores.
modeltype: JC    # (required) The model type. A list of implemented models is given below.    
end model        # (required) Indicates the end of the model
```

Amino acid models:
```
begin aa-model   # (required) Indicates the start of a new amino acid model
name: Model_LG   # (required) Give each model a unique name. A name can be composed of characters,
                 #            numbers and underscores.
modeltype: LG    # (required) The model type. A list of implemented models is given below.    
end model        # (required) Indicates the end of the model
```

All character on a line after the first # symbol are treated as comments and are removed before the definition is interpreted.

Each model definition has to start with "begin model" and end with "end model".
The model "name" and "modeltype" are required lines in a model definition.


## Nucleotide models:

The list of available nucleotide "modeltype" is: JC, F81, K2P, F84, HKY, GTR
All other models such as Timur models, etc. be specified using the GTR model with specific
substitution rates.

Furthermore, invariant sites and gamma distributed site rates can be specified.

### Exemplary JC model with all available parameters:

```
begin model
name: Model_JC_shape_01_pinv_03
modeltype: JC
shape: 0.1        # (not required) If a shape parameter is specified, gamma distributed site rates 
                  # are assumed. If no heterogeneity is specified, site rates are homogeneous.
		  # The shape parameter must be >0.
pinv:  0.3        # proportion of invariant sites. Must be in the range: 0 <= pinv < 1.
end model
```

If no shape parameter is specified, site rates are homogeneous apart from invariant sites.
If no pinv is specified, no invariant site will occur apart from those that occur by chance.
In particular for small values of the shape parameter invariant site can occur by chance. Their proportion is not included in pinv.  

### Exemplary F81 model with all available parameters:
```
begin model
name: Model_F81
modeltype: F81
shape:     0.5
pinv:      0.3
basefreq:  0.2  0.3  0.2  0.3    # Base frequencies of A, C, G, T. Have to add up to 1.
end model
```

### Exemplary K2P model with all available parameters:
```
begin model
name: Model_K2P
modeltype: K2P
tstv:      2          # Transition transversion ratio
shape:     0.5
pinv:      0.3
end model
```

### Exemplary F84 model with all available parameters:
```
begin model
name: Model_F84
modeltype: F84
...
end model
```

### Exemplary HKY model with all available parameters:
```
begin model
name: Model_HKY
modeltype: HKY
tstv: 3         
basefreq: 0.2  0.3  0.2  0.3  # base frequencies
shape:     0.5
pinv:      0.3
end model
```


### Exemplary GTR model with all available parameters:
```
begin model
name: Model_GTR
modeltype: GTR
rrates:   2 10 3 1 11 1     # relative substitution rates in the order: 
                            # A->C,A->G,A->T,C->G,C->T,G->T
                            # These rates will be normalised automatically.
shape:    0.5
pinv:     0.3
basefreq: 0.2  0.3  0.2  0.3
end model
```

## Amino acid models:

Implemented amino acid models are: WAG, LG, JTT, DAY, WAG_OLD, WAG_STAR.
These are fixed rates models. The model parameters can be viewed in the file: mymodel.cpp

AA model definitions have to start with "begin aa-model" and end with "end model".

### Exemplary LG model with all available parameters:

```
begin model
name: Model_LG
modeltype: LG
shape:    0.5
pinv:     0.3
end model
```

If alternative amino acid frequencies shall be specified, a fully user specified amino acid model can be used:

```
begin aa-model
name: Model_jtt-rate-user-specified
modeltype: User
rrates:  0.531678 0.557967 0.827445 0.574478 0.556725 1.066681 1.740159 0.219970 0.361684 0.310007 0.369437 0.469395 0.138293 1.959599 3.887095 4.582565  0.084329 0.139492 2.924161
         0.451095 0.154899 1.019843 3.021995 0.318483 1.359652 3.210671 0.239195 0.372261 6.529255 0.431045 0.065314 0.710489 1.001551 0.650282 1.257961 0.235601 0.171995
         5.549530 0.313311 0.768834 0.578115 0.773313 4.025778 0.491003 0.137289 2.529517 0.330720 0.073481 0.121804 5.057964 2.351311 0.027700 0.700693 0.164525             
         0.105625 0.521646 7.766557 1.272434 1.032342 0.115968 0.061486 0.282466 0.190001 0.032522 0.127164 0.589268 0.425159 0.057466 0.453952 0.315261
         0.091304 0.053907 0.546389 0.724998 0.150559 0.164593 0.049009 0.409202 0.678335 0.123653 2.155331 0.469823 1.104181 2.114852 0.621323
         3.417706 0.231294 5.684080 0.078270 0.709004 2.966732 0.456901 0.045683 1.608126 0.548807 0.523825 0.172206 0.254745 0.179771
         1.115632 0.243768 0.111773 0.097485 1.731684 0.175084 0.043829 0.191994 0.312449 0.331584 0.114381 0.063452 0.465271
         0.201696 0.053769 0.069492 0.269840 0.130379 0.050212 0.208081 1.874296 0.316862 0.544180 0.052500 0.470140
         0.181788 0.540571 0.525096 0.329660 0.453428 1.141961 0.743458 0.477355 0.128193 5.848400 0.121827
         2.335139 0.202562 4.831666 0.777090 0.098580 0.405119 2.553806 0.134510 0.303445 9.533943
         0.146481 3.856906 2.500294 1.060504 0.592511 0.272514 0.530324 0.241094 1.761439                         
         0.624581 0.024521 0.216345 0.474478 0.965641 0.089134 0.087904 0.124066                                                     
         0.436181 0.164215 0.285564 2.114728 0.201334 0.189870 3.038533
         0.148483 0.943971 0.138904 0.537922 5.484236 0.593478
         2.788406 1.176961 0.069965 0.113850 0.211561
         4.777647 0.310927 0.628608 0.408532                                               
         0.080556 0.201094 1.143980
         0.747889 0.239697
         0.165473
basefreq 0.076862 0.051057 0.042546 0.051269 0.020279 0.041061 0.061820 0.074714 0.022983 0.052569 0.091111 0.059498 0.023414 0.040530 0.050532 0.068225 0.058518 0.014336 0.032303 0.066374
end model
```

## Special site rate heterogeneities

### Free rate models

The user can specify free rate models with a specific number of discrete site rate categories.
The relative rates are specified. They are normalised before being used. The final site rates will
be proportional to the specified numbers.
```
begin model
name: Modell_JC_free_catrates
modeltype: JC
cat-rates 0 1 1 1 1 9 9 9 9 6 # Ten equally large categories for which proportional rates can be specified. 
                              # to the values given.
                              # 40 % of sites have a rate of 1 while another 40% have a rate of 9.
			      # This distribution has two peaks at the unnormalised rates 1 and 9.
			      # PolyMoSim normalises the rates, so that their mean value is 1. 
end model
```

### Discrete gamma distribution

Heterogeneous site rates that are gamma distributed are drawn from a continuous gamma distribution.
Since tree reconstruction programs use discrete site rate categories it is also interesting to be able to conduct simulations using discrete site rates approximating a continuous gamma distribution as close a possible for the given number of categories. 

Example using discrete gamma distribution with 4 discrete rate categories
```
begin model
name: Modell_JC_shape_1_pinv03_ncat4
modeltype: JC
shape: 1          # (not required) If a shape parameter is specified, gamma distributed site rates 
                  # are assumed. Of no heterogeneity is specified, site rates are homogeneous.
pinv:  0.3        # proportion of invariant sites
ncat: 4
end model
```



