/***************************************************************************************************
*  This file is part of PolyMoSim.
*  
*  Copyright (c) 2006-2025, Christoph Mayer
*  
*  Author: Christoph Mayer
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



#ifndef ASTERIX_H
#define ASTERIX_H

//**************************************
#define PROGNAME    "PolyMoSim"
#define VERSION     "1.1.7-devel"
#define EPS         0.00000001
//**************************************
// #define DEBUG
//**************************************

#include <cstdint>
#include "CFile/CFile2_3.h"
#include <cstdlib>
#include <cstdio>
#include "faststring3.h"

#include "global-types-and-parameters.hpp"

#include <iostream>
#include <fstream>


#ifdef  DEBUG
#define DEBUGCODE(x) x
#else
#define DEBUGCODE(x)
#endif

// Switch on a single DEBUGCODE section by adding a D in front
#define DDEBUGCODE(x) x
#define NODEBUGCODE(x)

#define SPECIAL_COMMENT -2;

#define UNUSED(x) (void(x))


const char welcome_str[] = 
            "\n\n"
            "      Welcome to " PROGNAME ", version " VERSION ",\n"
            "      a program to simulate the evolution of nucleotide and amino acid sequences.\n"
            "      Copyright (C) 2007-2025 Christoph Mayer.\n"
            "      This program is free software, distributed under the GNU General Public\n"
            "      License version 3. It comes with ABSOLUTELY NO WARRANTY.\n"
/*             "      This program is distributed directly by the author in\n" */
/*             "      form of binary executables.\n" */
/*             "      It can be used freely for academic purposes.\n" */
/*             "      Results obtained with it can be published without\n" */
/*             "      restrictions, provided the program and its author are\n" */
/*             "      acknowledged by name.\n" */
/*             "      It is not allowed to use Phobos for commercial\n" */
/*             "      purposes without permission from the author.\n" */
/*             "      This version of Phobos is still a test version.\n" */
/*             "      Phobos comes without warranty!\n\n" */
/*             "      Important note: The search mode: Search for repeats with\n" */
/*             "                      a minimum percentage perfection is not\n" */
/*             "                      fully implemented in this version.\n\n\n"; */
  "\n\n";



inline bool FileExists(const char *s)
{
  std::ifstream is(s);

  if (is.fail())
  {
    return false;	
  }
  is.close();
  return true;
}

inline void ignore_spaces(std::istream& is) {  //ignores whitespaces
  while(!is.fail()) {
    if (!std::isspace(is.get())) {
      is.unget();
      break;
    }
  }
}

inline void ignore_spaces(CFile& is) {
  char ch;
  while(!is.eof()) {
    ch = is.getchar();
    if(!std::isspace(ch)) {
      is.ungetchar();
      break;
    }
  }
}

inline void skip_comment(std::istream& is) {
  do {
  } while( !is.fail() && is.get() != ']' );
}


/* // ignores leading spaces and comments */
/* inline int mygetchar(std::istream& is, bool look_for_special_comments = false) { */
/*   int ch; */

/*   ignore_spaces(is); */
/*   ch = is.get(); */
/*   while (ch == '[' ) { */
/*     if (look_for_special_comments) */
/*     { */
/*       ch = is.peek(); */
/*       if (ch == '&') */
/* 	return SPECIAL_COMMENT; */
/*     } */
/*     skip_comment(is); */
/*     ignore_spaces(is); */
/*     ch = is.get(); */
/*   } */
/*   return ch; */
/* } */

int myother_main(int, char**);

#endif
