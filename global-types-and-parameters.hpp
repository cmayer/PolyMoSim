/***************************************************************************************************
*  This file is part of PolyMoSim.
*  
*  Copyright (c) 2004-2026, Christoph Mayer and contributors (see README.md)
*  
*  Author: Dr. Christoph Mayer
*  Affiliation: Leibniz Institute for the Analysis of Biodiversity Change (LIB), Bonn, Germany
*  
*  PolyMoSim is free software: you can redistribute it and/or modify it under the terms of
*  the GNU General Public License version 3 as published by the Free Software Foundation.
*  
*  PolyMoSim is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
*  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
*  See the GNU General Public License for more details.
*  
*  You should have received a copy of the GNU General Public License along with PolyMoSim
*  (see the file LICENSE). If not, see <https://www.gnu.org/licenses/>.
*  
*  SPDX-License-Identifier: GPL-3.0-only
*  
*  This project makes use of code coming from other projects. What follows is a complete
*  list of files which make use of external code. Please refer to the copyright within
*  these files.
*  
*  Files in tclap folder:        Copyright (c) 2003 Michael E. Smoot, MIT license.
*                                See tclap/COPYING for details.
*  discrete_gamma.cpp:           Copyright 1993-2004 by Ziheng Yang (PAML), GNU GPL v3.
*                                See copyright in this file for details.
*  CRandom.h, CRandom.cpp:       Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
*                                BSD 3-clause license. See copyright in these files for details.
***************************************************************************************************/

#ifndef GLOBAL_TYPES_AND_PARAMETERS_H
#define GLOBAL_TYPES_AND_PARAMETERS_H

#include <cstdio>
#include <iostream>
#include <vector>
#include <string>
#include <cstring>
// #include "faststring.h"
#include <climits>
#include <string>
#include <fstream>


enum enum_rg {
  rg_MT19937
};

enum enum_outputformat {
  outputformat_logfile, outputformat_nexus, outputformat_fasta, outputformat_phylip, outputformat_phylip_no_spaces, outputformat_pattern_absolute, outputformat_pattern_relative, outputformat_pattern_absolute_fill, outputformat_pattern_relative_fill
};

const char randomNumberGeneratorNames[][70] = {"Mersenne Twister MT19937 by Takuji Nishimura and Makoto Matsumoto"};

#ifndef GLOBAL_VARS
// #define GLOBAL_VARS
extern std::string                      global_tree_file;
extern std::string                      global_model_file;

extern std::string                      global_preanalysis_file;
extern std::string                      global_postanalysis_file;
extern std::string                      global_log_file;
extern bool                             global_logging;
extern std::string                      global_output_filename;

extern std::string                      global_ancestral_sequence_file;
extern std::string                      global_siterateshist_file;
extern std::string                      global_siteratesdata_file;

extern unsigned                         global_num_repetitions;
extern bool                             global_use_GUI;
extern unsigned                         global_seed_random_generator;
extern enum_rg                          global_rg;
extern enum_outputformat                global_outputformat;

extern unsigned                         global_verbosity;
#endif


void ReportErrorMessage(const char *);


/* void PrintMessage(FILE *of, const char *); */
/* void PrintMessage(FILE *of, const char *, const char *); */
/* void PrintMessage(FILE *of, const char *, const char *, const char *); */
/* void PrintMessage(FILE *of, const char *, const char *, const char *, const char *); */

void myPrint(FILE *of, const char *);
void myPrint(FILE *of, const char *, const char*);
void myPrint(FILE *of, const char *, const char*, const char *);
void myPrint(FILE *of, const char *, const char*, const char *, const char *);

void myPrint(std::ostream &os, const char *);
void myPrint(std::ostream &os, const char *, const char*);
void myPrint(std::ostream &os, const char *, const char*, const char *);
void myPrint(std::ostream &os, const char *, const char*, const char *, const char *);



#define macromax(x,y) ((x)<(y) ? (y) : (x))
#define macromin(x,y) ((x)<(y) ? (x) : (y))


/* #ifdef  DEBUG */
/* #define DEBUGOUT1(x)        fprintf(stderr, x); */
/* #define DEBUGOUT2(x,y)      fprintf(stderr, x), fprintf(stderr, y); */
/* #define DEBUGOUT3(x,y,z)    fprintf(stderr, x), fprintf(stderr, y), fprintf(stderr, z); */
/* #define DEBUGOUT4(x,y,z,w)  fprintf(stderr, x), fprintf(stderr, y), fprintf(stderr, z), fprintf(stderr, w); */
/* #define DEBUGOUT_INT(x,y)   fprintf(stderr, "%s%d\n", x, y); */
/* #define DEBUGOUT_ENDL       { fprintf(stderr, "\n"), fflush(stderr); } */
/* #else */
/* #define DEBUGOUT1(x) */
/* #define DEBUGOUT2(x,y) */
/* #define DEBUGOUT3(x,y,z) */
/* #define DEBUGOUT4(x,y,z,w) */
/* #define DEBUGOUT_INT(x,y) */
/* #define DEBUGOUT_ENDL */
/* #endif */



void init_param();
int  read_and_init_parameters(int argc, char** argv);

void print_analysis_parameters(FILE *of, const char *s);
void print_analysis_parameters(std::ofstream &of, const char *s);

#endif
