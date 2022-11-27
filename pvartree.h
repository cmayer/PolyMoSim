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

#ifndef PvartreeH
#define PvartreeH

//aktuell

// varies the parameters of a specified startmodel over a given tree

#include <iostream>
#include <fstream>
#include <cstring>
#include <sstream>
#include <vector>
#include "faststring2.h"
#include "BasicTree.h"
#include "model_admin.h"
#include "mymodel.h"
#include "tree_admin.h"
#include "CFile/CFile2_1.h"

class pvartree {

 public:
  class error {
  };

  class readerror {
  private:
    int line;
    faststring unknown;
  public:
    readerror(int l, faststring u):line(l), unknown(u) {}
    int getLine(){return line; }
    faststring getUnknownKeyword() { return unknown; }
  };

 private:
  model_admin inmodel;
  model_admin outmodel;
  const nuc_model* startmodel;
  //  mymodel* newmodel;
  BasicTree tree;
  int name_i;

  // values for normal_random numbers
  // standard deviations
  double stdev_rAC;
  double stdev_rAG;
  double stdev_rAT;
  double stdev_rCG;
  double stdev_rCT;
  double stdev_rGT;
  double stdev_tstv;
  double stdev_shape;
  double stdev_inv;
  double stdev_A;
  double stdev_G;
  double stdev_C;
  //double mean;      

  // normal_random nrand;
  double (*ranf)(double);  // pointer to function which creates normal random numbers

 public:
  pvartree();
  //pvartree(double ());

  //  void set_startmodel(const mymodel*);

  int get_counter() const { return name_i;}
  void reset_counter() { name_i = 0; }

  double get_stdev_rAC() const;
  double get_stdev_rAG() const;
  double get_stdev_rAT() const;
  double get_stdev_rCG() const;
  double get_stdev_rCT() const;
  double get_stdev_rGT() const;
  double get_stdev_tstv() const;
  double get_stdev_shape() const;
  double get_stdev_inv() const;
  double get_stdev_A() const;
  double get_stdev_G() const;
  double get_stdev_C() const;
  //double get_mean() const;

  //void set_standard_deviation(double);
  //void set_mean(double);
  void set_random_generator(double (double));

  void vary_tree(int, istream&, const char* filename, double, const faststring& = "");
  void vary_models(BasicNode*, int);

  void read_parameters_from_file(const char*);
  //  void set_parameter(istream& = cin);

  void print_varied_tree(ostream&);
  void print_models(ostream&);

 private:
  void print_next(BasicNode*, ostream&);
  nuc_model* vary_parent_model(const nuc_model*, faststring&);

};

inline double pvartree::get_stdev_rAC() const {
  return stdev_rAC;
}

inline double pvartree::get_stdev_rAG() const {
  return stdev_rAG;
}

inline double pvartree::get_stdev_rAT() const {
  return stdev_rAT;
}

inline double pvartree::get_stdev_rCG() const {
  return stdev_rCG;
}

inline double pvartree::get_stdev_rCT() const {
  return stdev_rCT;
}

inline double pvartree::get_stdev_rGT() const {
  return stdev_rGT;
}

inline double pvartree::get_stdev_tstv() const {
  return stdev_tstv;
}

inline double pvartree::get_stdev_shape() const {
  return stdev_shape;
}

inline double pvartree::get_stdev_inv() const {
  return stdev_inv;
}

inline double pvartree::get_stdev_A() const {
  return stdev_A;
}

inline double pvartree::get_stdev_G() const {
  return stdev_G;
}

inline double pvartree::get_stdev_C() const {
  return stdev_C;
}


#endif
