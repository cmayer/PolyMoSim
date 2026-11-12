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

