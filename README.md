# PolyMoSim: Nucleotide and amino acid phylogenetic sequence alignment simulator

Table of contents:

- [About the PolyMoSim program](#about-the-PolyMoSim-package)
- [Compiling and installing PolyMoSim](#compiling-and-installing)
  * [System requirements:](#system-requirements)
- [Quickstart](#quickstart)
- [Documentation](#documentation)
- [Frequently asked questions](#Frequently-aksed-questions)

## About the PolyMoSim program <a id="about-the-PolyMoSim-package"></a>

PolyMoSim is a program that simulates sequence evolution of nucleotide or amino acid sequences. For specified evolutionary models, model parameters and phylogenetic trees, it can simulate data sets that evolved unter this tree and model. PolyMoSim is a fast and flexible program that supports a large number of evolutionary models, allows different models on different branches, mixture models and even a site heterogeneity that differs on different branches of the tree.

It has been designed to test phylogenetic tree reconstruction programs and to train Machine learning models for phylogenetic reconstruction.

## Compiling and installing PolyMoSim <a id="compiling-and-installing"></a>
- Download the project or clone the project locally.
- On the command line go to the project folder and type "make".
- Make sure that you copy the PolyMoSim-vx.y.z program to a folder that listed in your $PATH variable, or copy it to the folder you want to use it from.

### System requirements:  <a id="system-requirements"></a>
PolyMoSim can be compiled on all platforms, which provide a C++ compiler.
In particular this includes Windows, MacOS and Linux operating systems.
Here I will only explain how to compile it on Mac and Linux computers.

## Documentation <a id="documentation"></a>
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


This program simulates the evolution of molecular (nucleotide and amino acid) sequences.  


## Quickstart <a id="quickstart"></a>
./PolyMoSim-vx.y.z. -i input-fasta-file.fas -o output-fasta.fas



## Frequently asked questions <a id="Frequently-aksed-questions"></a>
No questions have been asked so far.
