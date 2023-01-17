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

So PolyMoSim requires the model file and the tree file as two required parameters. 
All other parameters are optional. Most of them have default values. Parameters you might want to have a look at are:
output-file-name (-o), setting the seed for the random number generator (-s) which is important if you do multiple simulations, output format (-f).
A full list of command line parameters is given below.


### Example analyses:
See the Example analyses folder for example analyses. The [README-Examples.md](README-Examples.md) 

### Full list of command line parameters:
PolyMoSim-vx.y.z  [--verbosity <integer>] -m <string> -t <string> [-o <string>] [-l <string>] [-n <unsigned>]
[--print_siterate_data <string>]
[--print_siterate_histogram <string>]
[--print_ancestral_seq <string>] [--post <string>]
[--pre <string>] [-s <unsigned int>] [-f <nexus|phylip
|phylip_no_spaces|fasta|site_pattern_freq_absolute
|site_pattern_freq_relative
|site_pattern_freq_absolute_fill
|site_pattern_freq_relative_fill>] [--] [-v] [-h]  


Where:   
--verbosity <integer>  
(value required)  Adjust the level of additional information given to the user. Values from 0-5 are valid. Default 1. Set to 0 for less
output.  

-m <string>,  --modelfile <string>  
(Required)  (value required)  Model file for simulation.  

-t <string>,  --treefile <string>  
(Required)  (value required)  Tree file for simulation.  

-o <string>,  --outfile <string>  
(Value required)  Name of output file. If not specified, results are printed to standard output.  

-l <string>,  --log <string>  
(Value required)  File to write log information to.  

-n <unsigned>,  --nreps <unsigned>  
(Value required)  Number of independent data generated in simulation.
Default: 1  

--print_siterate_data <string>  
(Value required)  Print full site rates of all models to given file.
Default: No siterate information is printed.  

--print_siterate_histogram <string>  
(Value required)  Print site rates histogram for each model to given file. Default: No site rate information is printed.  

--print_ancestral_seq <string>  
(Value required)  With this option and by providing a file name, the ancestral sequence is printed to this file. Default: ancestral
sequence is not printed  

--post <string>  
(Value required)  File included in output after each generated data
set.  

--pre <string>  
(Value required)  File included in output before each generated data
set.  

-s <unsigned int>,  --seed <unsigned int>  
(Value required)  The seed value for the random number generator.
Default: time.  

-f <nexus|phylip|phylip_no_spaces|fasta|site_pattern_freq_absolute
|site_pattern_freq_relative|site_pattern_freq_absolute_fill
|site_pattern_freq_relative_fill>,  --outputFormat <nexus|phylip
|phylip_no_spaces|fasta|site_pattern_freq_absolute
|site_pattern_freq_relative|site_pattern_freq_absolute_fill
|site_pattern_freq_relative_fill>  
(Value required)  Output format of sequence data. Default: fasta. The
site_pattern* format list the site pattern frequencies instead of the alignment. The site_pattern*_fill formats list all patterns the non
fill formats only the site patterns that occurred in the simulated data
set.  

--,  --ignore_rest  
Ignores the rest of the labeled arguments following this flag.  

-v,  --version  
Displays version information and exits.  

-h,  --help  
Displays usage information and exits.  


## Frequently asked questions <a id="Frequently-aksed-questions"></a>
No questions have been asked so far.
