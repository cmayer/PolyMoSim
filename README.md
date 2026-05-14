# PolyMoSim: Nucleotide and amino acid phylogenetic sequence alignment simulator

Table of contents:

- [About the PolyMoSim program](#about-the-PolyMoSim-package)
- [Compiling and installing PolyMoSim](#compiling-and-installing)
  * [System requirements:](#system-requirements)
- [Quickstart](#quickstart)
- [Documentation](#documentation)
- [Frequently asked questions](#Frequently-aksed-questions)

## About the PolyMoSim program <a id="about-the-PolyMoSim-package"></a>

PolyMoSim can simulate the evolution of nucleotide or amino acid sequences for a specified evolutionary models, model parameters and phylogenetic trees. PolyMoSim is a fast and flexible program that supports a large number of evolutionary models, allows different models on different branches, mixture models and even a site heterogeneity that differs on different branches of the tree.

It has been designed to test phylogenetic tree reconstruction programs and to train machine learning models for phylogenetic reconstruction.

## Compiling and installing PolyMoSim <a id="compiling-and-installing"></a>
- Download the project or clone the project locally.
- On the command line go to the project folder and type "make".
- Make sure that you copy the PolyMoSim-vx.y.z program to a folder that is listed in your $PATH variable so that your system can always find it, copy it to the folder you want to use it from or specify the full path to the program.

### System requirements:  <a id="system-requirements"></a>
PolyMoSim can be compiled on all platforms, which provide a C++ compiler.
In particular this includes Windows, MacOS and Linux operating systems.
Here I will only explain how to compile it on Mac and Linux computers.

## Documentation <a id="documentation"></a>

### Required input files:
PolyMoSim requires two input files: 

i) The model file, which contains a list specified evolutionary models with model parameters that can be used for the simulation. See the [README-model-files.md.](README-model-files.md) for more details. 

ii) The tree file, which contains one or multiple lines specifying the partition size and evolutionary tree. Each line specifies the the information for one partition.
See the [README-tree-files.md](README-tree-files.md) for more details.

**Quickstart:**
The simplest way to start PolyMoSim is as follows:
```
PolyMoSim-vx.y.z -m model-file.txt -t tree-file.txt --outfile simulation-result.fas
```
or
```
PolyMoSim-vx.y.z -m model-file.txt -t tree-file.txt 1> simulation-result.fas
```

PolyMoSim has two required paramters: the model file and the tree file. 
All other parameters are optional. Most of them have default values. Parameters you might want to have a look at are:
output-file-name (-o), setting the seed for the random number generator (-s) setting the random number seed,
which is important if you do multiple simulations, and the output format (-f).
A full list of command line parameters is given below.


### Examples:
Example simulations are provided in the Examples folder. See the [Examples/README-Examples.md](Examples/README-Examples.md) for a list of examples and a brief explanation.

### List of evolutionary models avaiable for simulates:

Nucleotide and amino acid models have to be specified in the model file. Models that shall be used are specified in the tree file. There a model is specified for the whole tree. Models can be altered for each branch. Models that are changed on a branch will be used for the whole clade, unless another model is specified for a subclade.

A good starting point is to look at the examples in the Example folder.

Any number of nucleotide and amino acid substitution models can be specified in the model file.
A model specification consists of a model name, the model type, and the model parameters.

**Nucleotide models:**
- For nucleotide models the following model types are known: JC, F81, K2P, F84, HKY, GTR.

Depending on the model type, different model parameters can be specified: tstv, rrates, shape, ncat, propinv, base_frequencies, siterate-distribution. The siterate-distribution parameter allows specifying a site rate distribution function using the following operations and functions: "+", "-", "*", "/", "^", "#", "$", "Dist_Gauss", "Dist_Gamma", "Dist_Uniform", "Heavy", "Gamma", "Dist_Beta", "Dist_Cauchy", "Dist_Uniform_Interval". More documentation will be added. The mathematical expression is parsed and evaluated in math_expression_parser.h.

**Amino Acid models:**
JTT, LG, WAG_OLD, WAG, WAG_STAR, DAY, Qplant, Qpfam, Qmammal, Qinsect, Qbird, Qlg, Qyeast.
User specified models are also possible. Here the relative rates and the amino acid frequencies have to be specified in the model file. An example is provided in the Example-3_aa-models folder.

As for nucleotide models, the following model parameters can be specified: rrates, shape, ncat, propinv, base_frequencies, siterate-distribution. 

Remark: Amino acid substitution matrices are hard coded in the mymodel.cpp file. This file also contains references to the models.

**Mixture models:**
Mixture models can be specified by using multiple partitions in the tree file.






# PolyMoSim Command Line Reference

## Long Synopsis

```
PolyMoSim-vx.y.z -m <string> -t <string>
[-o <string>] [-l <string>] [-n <unsigned>]
[--print_siterate_data <string>]
[--print_siterate_histogram <string>]
[--print_ancestral_seq <string>] [--post <string>]
[--pre <string>] [-s <unsigned int>] [-f <nexus|phylip
|phylip_no_spaces|fasta|site_pattern_freq_absolute
|site_pattern_freq_relative
|site_pattern_freq_absolute_fill
|site_pattern_freq_relative_fill>] [--verbosity <integer>]
[--] [-v] [-h]  
```

## Short Synopsis

```
PolyMoSim-vx.y.z -m <string> -t <string> [OPTIONS]
```

## Required Arguments

| Flag | Description |
|------|-------------|
| `-m <string>`, `--modelfile <string>` | Model file for simulation. |
| `-t <string>`, `--treefile <string>` | Tree file for simulation. |

## Optional Arguments

### Output & Logging

| Flag | Description |
|------|-------------|
| `-o <string>`, `--outfile <string>` | Name of output file the simulated sequence data is written to. Several different file formats can be chosen. If not specified, output is printed to console (standard output). |
| `-l <string>`, `--log <string>` | File to write log information to. |
| `--verbosity <integer>` | Level of additional information given to the user. Valid values: `0`–`200`. Default: `1`. Set to `0` for less output. Values above >=100 are used for debugging. Steps for which output changes: 0,1,2,3,4,100,200|

### Simulation Control

| Flag | Description |
|------|-------------|
| `-n <unsigned>`, `--nreps <unsigned>` | Number of independent datasets generated in the simulation. Default: `1`. |
| `-s <unsigned int>`, `--seed <unsigned int>` | Seed value for the random number generator. Default: current time. **Do not rely on the default seed if you start many analyses simultaneously.** |

### Output Format

```
-f <format>, --OutputFormat <format>
```

Controls the output format for the sequence data. Default: `fasta`.

| Format | Description |
|--------|-------------|
| `nexus` | Nexus format. |
| `phylip` | Phylip format. |
| `phylip_no_spaces` | Phylip format without spaces. |
| `fasta` | FASTA format *(default)*. |
| `site_pattern_freq_absolute` | Site pattern frequencies (absolute counts, observed patterns only). |
| `site_pattern_freq_relative` | Site pattern frequencies (relative counts, observed patterns only). |
| `site_pattern_freq_absolute_fill` | Site pattern frequencies (absolute counts, all patterns including unobserved). |
| `site_pattern_freq_relative_fill` | Site pattern frequencies (relative counts, all patterns including unobserved). |

> The `site_pattern_*` formats output the site pattern frequencies instead of the alignment. The `*_fill` variants include all possible patterns; the non-fill variants list only patterns that occurred in the simulated dataset.

The `site_pattern_freq_relative_fill` is useful for machine learning training based on site pattern frequencies.

### Ancestral & Site Rate Output

| Flag | Description |
|------|-------------|
| `--print_ancestral_seq <string>` | Print the ancestral sequence to the specified file. Default: not printed. |
| `--print_siterate_data <string>` | Print site rate list for all models to the specified file. Default: not printed. |
| `--print_siterate_histogram <string>` | Print a site rate histogram for each model to the specified file. Default: not printed. |

### Dataset Wrapping

| Flag | Description |
|------|-------------|
| `--pre <string>` | Content of this file is copied into the output **before** each generated dataset. Useful in the nexus format. |
| `--post <string>` | Content of this file is copied into the output **after** each generated dataset. Useful in the nexus format. |

### General

| Flag | Description |
|------|-------------|
| `--`, `--ignore_rest` | Ignores all remaining labeled arguments after this flag. |
| `-v`, `--version` | Display version information and exit. |
| `-h`, `--help` | Display usage information and exit. |



## Frequently asked questions <a id="Frequently-aksed-questions"></a>
No questions have been asked so far.
