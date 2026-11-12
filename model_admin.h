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

#ifndef Model_AdminH
#define Model_AdminH

#include "PolyMoSim.h"

#include <vector>
#include <map>
#include "mymodel.h"
#include "faststring3.h"


class model_admin {

 public:
  class indexerror {
  private:
    faststring errormsg;
  public:
    indexerror(const faststring &s):errormsg(s) {}
    const faststring& getReason() { return errormsg; }
    const char * getReason_c_str() { return errormsg.c_str(); }
  };

  class formaterror{
  };

  class readerror{
  };

 private:
  std::map<faststring, basic_model*> modelmap;
  std::vector<basic_model*>      modelvector;

  void               create_unique_model_name(faststring &) const;

 public:
  model_admin(){
  }

  void create(const char *);
  void print(std::ostream&, unsigned=1);
  void print(FILE *os=stdout, unsigned=1);

  void default_modelvector(unsigned);

  const basic_model* get_model(unsigned) const;
  const basic_model* get_model(faststring);    // const faststring reference fails since the modelmap has no const strings, which could be seen as an issue xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  void               append_a_copy(const nuc_model*, faststring="");
  void               append_a_copy(const aa_model*,  faststring="");

  unsigned           append_new_default_model();
  void               duplicate_model(unsigned i);
  void               delete_model(unsigned i);

  void               swap_models(unsigned i1, unsigned i2);

  // void            append_this(const nuc_model*);
  // void            append_this(const aa_model*);

  void               init_models(double(double, double), double());
  void               init_siterates(unsigned len, bool reinit);
  void               reset_siterates();

  bool               modelNameExists(faststring) const;
  basic_model*       find_model_with_name(faststring str);

  unsigned           getNumModels() const;
  faststring         getModelName(unsigned i) const;
  void               changeModelName(unsigned i, faststring s);

  void               print_relative_site_rates(std::ostream &os) const;
  void               print_site_rates_histogramm_data(std::ostream &os, faststring &) const;

  //model_admin& operator= (const model_admin&);



  void append_new_model(char   data_type,             // 'n' for nuc or 'p' for protein
                                     faststring modelname_param,
                                     faststring modeltype_param,
                                     std::vector<double> *rrates_param, // Parameters are supplied as in an upper triangular matrix.
                                     std::vector<double> *base_param,
                                     double         shape_param,
                                     double         pinv_param,
                                     unsigned       ncat_param,
                                     double         *tstv_param

  );

};


//  Die folgende Implementierung von operator= ist gefaehrlich, da die Modellzeiger kopiert werden, nicht aber die einzellnen Modelle.
//  Deshalb ist der Code auskommentiert. Gibt es einen Default Zuweisungsoperator? Wenn ja, hat der sicherlich das gleiche Problem.
//  TODO: xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
// inline model_admin& model_admin::operator= (const model_admin& x) {
//  modelmap = x.modelmap;
//  modelvector = x.modelvector;
//  return *this;
//}

//inline model_admin& model_admin::operator= (const model_admin& x) {
//  std::cerr << "operator= not yet implemented!!!!!!!!" << std::endl;
//  exit(0);
//  return *this;
// }



#endif
