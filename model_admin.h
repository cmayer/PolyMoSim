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

#ifndef Model_AdminH

#define Model_AdminH



#include "PolyMoSim.h"



#include <vector>

#include <map>

#include "mymodel.h"

#include "faststring2.h"





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

