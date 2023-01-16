# PolyMoSim: Nucleotide and amino acid phylogenetic sequence alignment simulator

Table of contents:

- [About the PolyMoSim program](#about-the-PolyMoSim-package)
- [Compiling and installing PolyMoSim](#compiling-and-installing)
  * [System requirements:](#system-requirements)
- [Quickstart](#quickstart)
- [Documentation](#documentation)
- [Frequently aksed questions](#Frequently-aksed-questions)

## About the PolyMoSim program <a id="about-the-PolyMoSim-package"></a>
PolyMoSim can simulate the evolution of nucleotide and amino acid sequence alignments for a given phylogenetic tree and given evolutionary model and model parameters.
PolyMoSim is a fast and flexible simulation program with a wide range of models that allows to have different models of sequence evolution on different branches of the tree. Even the rate heterogeneity can be changed or kept in different parts of the tree.

PolymoSim can be used to test pylogenetic tree reconstuction programs and to train Machine learning models for phylogenetic reconstruction.

## Compiling and installing OliInSeq <a id="compiling-and-installing"></a>
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
(value required)  Adjust the level of additional information given to
the user. Values from 0-5 are valid. Default 1. Set to 0 for less
output.  

-m <string>,  --modelfile <string>  
(required)  (value required)  Model file for simulation.  

-t <string>,  --treefile <string>  
(required)  (value required)  Tree file for simulation.  

-o <string>,  --outfile <string>  
(value required)  Name of output file. If not specified, results are
printed to standard output.  

-l <string>,  --log <string>  
(value required)  File to write log information to.  

-n <unsigned>,  --nreps <unsigned>  
(value required)  Number of independent data generated in simulation.
Default: 1  

--print_siterate_data <string>  
(value required)  Print full site rates of all models to given file.
Default: No siterate information is printed.  

--print_siterate_histogram <string>  
(value required)  Print site rates histogram for each model to given
file. Default: No site rate information is printed.  

--print_ancestral_seq <string>  
(value required)  With this option and by providing a file name, the
ancestral sequence is printed to this file. Default: ancestral
sequence is not printed  

--post <string>  
(value required)  File included in output after each generated data
set.  

--pre <string>  
(value required)  File included in output before each generated data
set.  

-s <unsigned int>,  --seed <unsigned int>  
(value required)  The seed value for the random number generator.
Default: time.  

-f <nexus|phylip|phylip_no_spaces|fasta|site_pattern_freq_absolute
|site_pattern_freq_relative|site_pattern_freq_absolute_fill
|site_pattern_freq_relative_fill>,  --outputFormat <nexus|phylip
|phylip_no_spaces|fasta|site_pattern_freq_absolute
|site_pattern_freq_relative|site_pattern_freq_absolute_fill
|site_pattern_freq_relative_fill>  
(value required)  Output format of sequence data. Default: fasta. The
site_pattern* format list the site pattern frequencies instead of the
alignment. The site_pattern*_fill formats list all patterns the non
fill formats only the site patterns that occurred in the simulated data
set.  

--,  --ignore_rest  
Ignores the rest of the labeled arguments following this flag.  

-v,  --version  
Displays version information and exits.  

-h,  --help  
Displays usage information and exits.  


This program simulates the evolution of molecular (nucleotide and amino
acid) sequences.  


## Quickstart <a id="quickstart"></a>
./OliInSeq-vx.y.z. -i input-fasta-file.fas -o output-fasta.fas



## Frequently asked questions <a id="Frequently-aksed-questions"></a>
No questions have been asked so far.

