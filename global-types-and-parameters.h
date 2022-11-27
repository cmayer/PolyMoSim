/***************************************************************************************************
*  The PolyMoSim project is distributed under the following license:
*  
*  Copyright (c) 2006-2022, Christoph Mayer, Forschungsmuseum Alexander Koenig, Bonn, Germany
*  All rights reserved.
*  
*  Redistribution and use in source and binary forms, with or without
*  modification, are permitted provided that the following conditions are met:
*  1. Redistributions of source code (complete or in parts) must retain
*     the above copyright notice, this list of conditions and the following disclaimer.
*  2. Redistributions in binary form must reproduce the above copyright
*     notice, this list of conditions and the following disclaimer in the
*     documentation and/or other materials provided with the distribution.
*  3. All advertising materials mentioning features or any use of this software
*     e.g. in publications must display the following acknowledgement:
*     This product includes software developed by Christoph Mayer, Forschungsmuseum
*     Alexander Koenig, Bonn, Germany.
*  4. Neither the name of the organization nor the
*     names of its contributors may be used to endorse or promote products
*     derived from this software without specific prior written permission.
*  
*  THIS SOFTWARE IS PROVIDED BY CHRISTOPH MAYER ''AS IS'' AND ANY
*  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
*  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
*  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHTHOLDER OR ITS ORGANISATION BE LIABLE FOR ANY
*  DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
*  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
*  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
*  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
*  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
*  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*  
*  IMPORTANT (needs to be included, if code is redistributed):
*  Please not that this license is not compatible with the GNU Public License (GPL)
*  due to paragraph 3 in the copyright. It is not allowed under any
*  circumstances to use the code of this software in projects distributed under the GPL.
*  Furthermore, it is not allowed to redistribute the code in projects which are
*  distributed under a license which is incompatible with one of the 4 paragraphs above.
*  
*  This project makes use of code coming from other projects. What follows is a complete
*  list of files which make use of external code. Please refer to the copyright within
*  these files.
*  
*  Files in tclap foler:         Copyright (c) 2003 Michael E. Smoot
*                                See copyright in tclap/COPYRIGHT file for details.	
*  discrete_gamma.c:             Copyright 1993-2004 by Ziheng Yang.
*                                See copyright in this file for details.
*  CRandom.h:                    Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura
*                                See copyright in this file for details.
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
