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

#ifndef ASTERIX_H
#define ASTERIX_H

//**************************************
#define PROGNAME    "PolyMoSim"
#define VERSION     "1.1.4"
#define EPS         0.00000001
//**************************************
// #define DEBUG
//**************************************

#include "CFile/CFile2_1.h"
#include <cstdlib>
#include <cstdio>
#include "faststring2.h"

#include "global-types-and-parameters.h"

#include <iostream>
#include <fstream>


#ifdef  DEBUG
#define DEBUGCODE(x) x
#else
#define DEBUGCODE(x)
#endif

// Switch on a single DEBUGCODE section by adding a D in front
#define DDEBUGCODE(x) x

#define SPECIAL_COMMENT -2;

#define UNUSED(x) (void(x))


const char welcome_str[] = 
            "\n\n"
            "      Welcome to " PROGNAME ", version " VERSION ",\n"
            "      a program to simulate the evolution of nucleotide and amino acid sequences.\n"
            "      Copyright (C) 2007-2021 Christoph Mayer.\n"
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
