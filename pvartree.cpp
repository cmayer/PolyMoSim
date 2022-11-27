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

#include "pvartree.h"



//Makros

#define DD1(x)        while((x) <  0.0000001)

#define DD_rates(x)   while((x) <  0.01)

#define DD2(x) while((x) <= 0)





using namespace std;



pvartree::pvartree() {

  name_i = 0;

  startmodel = NULL;

  tree.set_read_model_status(false);

}



//pvartree::pvartree(double (ran*)()) {

//  ranf = ran;

//}



//void pvartree::set_startmodel(const mymodel* start) {

//  startmodel = start;

//}



void pvartree::set_random_generator(double (*ran)(double)) { //Funktion bekommt nur Standardabweichung; neuer Wert wird extern berechnet

  ranf = ran;                                                //Mittelwert + neuer Wert

}



//double pvartree::get_mean() const {

//  return av;

//}



//void pvartree::set_mean(double a){

//  av = a;

//}



void pvartree::read_parameters_from_file(const char *filename) {

  vector<faststring> splitted;

  faststring dummy;

  CFile is;

  int number_of_tokens;



  stdev_rAC   = 0;

  stdev_rAG   = 0;

  stdev_rAT   = 0;

  stdev_rCG   = 0;

  stdev_rCT   = 0;

  stdev_rGT   = 0;

  stdev_tstv  = 0;

  stdev_shape = 0;

  stdev_inv   = 0;

  stdev_A     = 0;

  stdev_G     = 0;

  stdev_C     = 0;



  if(!FileExists(filename)) {

    cerr << "File: " << filename << " does not exist." << endl;

    throw error();

  }



  is.ffopen(filename);



  ignore_spaces(is);

  is.getline(dummy);

  dummy.replace_char(':' , ' ' );

  dummy.replace_char('=' , ' ' );

  dummy.ToLower();



  while(!is.fail()) {

    if(dummy[0] == '#') {

      ignore_spaces(is);

      is.getline(dummy);

      dummy.replace_char(':' , ' ' );

      dummy.replace_char('=' , ' ' );

      dummy.ToLower();

      continue;

    }

  

    number_of_tokens = split(splitted, dummy, " ");



    if(splitted[0] == "ac") {

      if(number_of_tokens != 2) {

	throw readerror(is.line(), "Bad standard deviation for rAC");

      }

      stdev_rAC = splitted[1].ToDouble();

    }

    else if(splitted[0] == "ag") {

      if(number_of_tokens != 2) {

	throw readerror(is.line(), "Bad standard deviation for rAG");

      }

      stdev_rAG = splitted[1].ToDouble();

    }

    else if(splitted[0] == "at") {

      if(number_of_tokens != 2) {

	throw readerror(is.line(), "Bad standard deviation for rAT");

      }

      stdev_rAT = splitted[1].ToDouble();

    }

    else if(splitted[0] == "cg") {

      if(number_of_tokens != 2) {

	throw readerror(is.line(), "Bad standard deviation for rCG");

      }

      stdev_rCG = splitted[1].ToDouble();

    }

    else if(splitted[0] == "ct") {

      if(number_of_tokens != 2) {

	throw readerror(is.line(), "Bad standard deviation for rCT");

      }

      stdev_rCT = splitted[1].ToDouble();

    }

    else if(splitted[0] == "gt") {

      if(number_of_tokens != 2) {

	throw readerror(is.line(), "Bad standard deviation for rGT");

      }

      stdev_rGT = splitted[1].ToDouble();

    }

    else if(splitted[0] == "tstv" || splitted[0] == "titv") {

      if(number_of_tokens != 2) {

	throw readerror(is.line(), "Bad standard deviation for tstv");

      }

      stdev_tstv = splitted[1].ToDouble();

    }

    else if(splitted[0] == "shape") {

      if(number_of_tokens != 2) {

	throw readerror(is.line(), "Bad standard deviation for shape");

      }

      stdev_shape = splitted[1].ToDouble();

    }

    else if(splitted[0] == "inv") {

      if(number_of_tokens != 2) {

	throw readerror(is.line(), "Bad standard deviation for inv");

      }

      stdev_inv = splitted[1].ToDouble();

    }

    else if(splitted[0] == "pia" || splitted[0] == "a") {

      if(number_of_tokens != 2) {

	throw readerror(is.line(), "Bad standard deviation for PI_A");

      }

      stdev_A = splitted[1].ToDouble();

    }

    else if(splitted[0] == "pic" || splitted[0] == "c") {

      if(number_of_tokens != 2) {

	throw readerror(is.line(), "Bad standard deviation for PI_C");

      }

      stdev_C = splitted[1].ToDouble();

    }

    else if(splitted[0] == "pig" || splitted[0] == "g") {

      if(number_of_tokens != 2) {

	throw readerror(is.line(), "Bad standard deviation for PI_G");

      }

      stdev_G = splitted[1].ToDouble();

    }

    else {

      faststring errormsg = "Umkown Parameter: " + splitted[0];

      throw readerror(is.line(), errormsg.c_str());

    }

    ignore_spaces(is);

    is.getline(dummy);

    dummy.replace_char(':' , ' ' );

    dummy.replace_char('=' , ' ' );

    dummy.ToLower();

  } //while

  is.ffclose();

}



/*

void pvartree::set_parameter(istream& is) {    //aus Datei einlesen: read_parameter_from_file

  string dummy;

  cout << "Set up the parameters needed to vary the startmodel: " << endl << endl;

  

  cout << "standard deviation for a: " ;

  is >> stdev_a;

 

  while(is.fail()) {

    cout << "error taking floating point value, please retry: " ;

    is.clear();

    is.ignore(100, '\n');

    is >> stdev_a;

  }



  cout << "standard deviation for b: " ;

  is >> stdev_b;

 

  while(is.fail()) {

    cout << "error taking floating point value, please retry: " ;

    is.clear();

    is.ignore(100, '\n');

    is >> stdev_b;

  }



  cout << "standard deviation for c: " ;

  is >> stdev_c;

 

  while(is.fail()) {

    cout << "error taking floating point value, please retry: " ;

    is.clear();

    is.ignore(100, '\n');

    is >> stdev_c;

  }



  cout << "standard deviation for d: " ;

  is >> stdev_d;

 

  while(is.fail()) {

    cout << "error taking floating point value, please retry: " ;

    is.clear();

    is.ignore(100, '\n');

    is >> stdev_d;

  }



  cout << "standard deviation for e: " ;

  is >> stdev_e;

 

  while(is.fail()) {

    cout << "error taking floating point value, please retry: " ;

    is.clear();

    is.ignore(100, '\n');

    is >> stdev_e;

  }



  cout << "standard deviation for f: " ;

  is >> stdev_f;

 

  while(is.fail()) {

    cout << "error taking floating point value, please retry: " ;

    is.clear();

    is.ignore(100, '\n');

    is >> stdev_f;

  }



  cout << "standard deviation for shape: " ;

  is >> stdev_shape;

 

  while(is.fail()) {

    cout << "error taking floating point value, please retry: " ;

    is.clear();

    is.ignore(100, '\n');

    is >> stdev_shape;

  }



  cout << "standard deviation for inv: " ;

  is >> stdev_inv;

 

  while(is.fail()) {

    cout << "error taking floating point value, please retry: " ;

    is.clear();

    is.ignore(100, '\n');

    is >> stdev_inv;

  }



  cout << "standard deviation for A: " ;

  is >> stdev_A;

 

  while(is.fail()) {

    cout << "error taking floating point value, please retry: " ;

    is.clear();

    is.ignore(100, '\n');

    is >> stdev_A;

  }



  cout << "standard deviation for G: " ;

  is >> stdev_G;

 

  while(is.fail()) {

    cout << "error taking floating point value, please retry: " ;

    is.clear();

    is.ignore(100, '\n');

    is >> stdev_G;

  }



  cout << "standard deviation for C: " ;

  is >> stdev_C;

 

  while(is.fail()) {

    cout << "error taking floating point value, please retry: " ;

    is.clear();

    is.ignore(100, '\n');

    is >> stdev_C;

  }

  is.clear();

  is.ignore(100, '\n');



  // get the random generator?!

}

*/



void pvartree::vary_tree(int tree_counter, istream& is, const char* modelfilename, double scalefactor, const faststring& defaultmodel) {

  BasicNode* node;

  faststring temp;

  faststring dummy;



  /*

  set_parameter();       // ausserhalb nur einmal setzen; für alle Bäume gleich xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  */



  inmodel.create(modelfilename);

  if(defaultmodel == "") {

    temp = "default";

  }else if(defaultmodel == "default"){

    temp = "default";

  }else{

    temp = defaultmodel;

  }

  startmodel = (const nuc_model*)inmodel.get_model(temp);



  dummy = "Model_default_";

  dummy.append(faststring(tree_counter) );

    

  tree.read_tree(is, scalefactor);  //can read from open stringstream         //was für scalenfaktor muss in dem fall mitgeschickt werden??????????????

  node = tree.get_treeroot();

  //  outmodel.append(dummy.c_str(), startmodel);



  // Need to copy start model using new name:





  nuc_model*     newmodel;



  newmodel = new nuc_model(dummy, 

			   (nuc_model::enummodeltype)startmodel->get_modeltype(),

			   startmodel->get_ratetype(),

			   startmodel->get_tstv(),

			   startmodel->get_rAC(),

			   startmodel->get_rAG(),

			   startmodel->get_rAT(),

			   startmodel->get_rCG(),

			   startmodel->get_rCT(),

			   startmodel->get_rGT(),

			   startmodel->get_shape(),

			   startmodel->get_inv(),

			   startmodel->get_ncat(),

			   startmodel->get_PI_A(),

			   startmodel->get_PI_C(),

			   startmodel->get_PI_G(),

			   startmodel->get_PI_T() );  



  outmodel.append_a_copy(newmodel, dummy);

  vary_models(node, tree_counter);



  delete newmodel;

}



void pvartree::vary_models(BasicNode* node, int tree_counter) {

  const nuc_model*  temp;

  faststring            dummy;

  nuc_model*        newmodel = NULL;

  unsigned          i;



  if(node->get_parent() == NULL){

    dummy = "Model_default_";

    dummy.append(faststring(tree_counter));



    node->set_model(startmodel);             // richtige Modell mit falschem Namen

    node->set_model_name(dummy);             // ------  überflüssig, da der Modelname im Knoten schon durch set_model gesetzt wird

  } else {

    temp = (const nuc_model*)node->get_parent()->get_model();

    

    dummy = "Model_t";

    dummy.append( faststring(tree_counter) );

    dummy.push_back('_');

    dummy.append( faststring(name_i));

    

    cout << endl;

    cout << "name: " << dummy << endl;

  

    newmodel = vary_parent_model(temp, dummy);

    // xxxxxxxxxx memory leak - is newmodel deleted?



    //    outmodel.append(dummy.c_str(), newmodel);

    outmodel.append_a_copy(newmodel, dummy);

    node->set_model(newmodel);

    //    delete newmodel;



    //    node->set_model_name(dummy);

  }

  

  for(i = 0; i < node->get_num_children(); i++) {

    name_i++;

    vary_models(node->get_childlist()[i], tree_counter);

  }

  if (newmodel)

    delete newmodel;

}





nuc_model* pvartree::vary_parent_model(const nuc_model* orig_model, faststring& modname) {    //die richtigen Fallunterscheidungen implementieren



  double nAC, nAG, nAT, nCG, nCT, nGT;

  double ntstv;

  double nshape;

  double ninv;

  double npiA, npiG, npiC, npiT;

  molecular_model<4>::enumratetype nrate;

  nuc_model::enummodeltype ntype;



  double orig_AC, orig_AG, orig_AT, orig_CG, orig_CT, orig_GT;

  nAC = orig_AC = orig_model->get_rAC();

  nAG = orig_AG = orig_model->get_rAG();

  nAT = orig_AT = orig_model->get_rAT();

  nCG = orig_CG = orig_model->get_rCG();

  nCT = orig_CT = orig_model->get_rCT();

  nGT = orig_GT = orig_model->get_rGT();



  double orig_tstv;

  ntstv = orig_tstv= orig_model->get_tstv();



  double orig_shape;

  double orig_inv;

  int    orig_ncat;



  double orig_piA;

  double orig_piC;

  double orig_piG;



  nshape = orig_shape =  orig_model->get_shape();

  ninv   = orig_inv   =  orig_model->get_inv();

  orig_ncat  =  orig_model->get_ncat();



  npiA = orig_piA = orig_model->get_PI_A();

  npiC = orig_piC = orig_model->get_PI_C();

  npiG = orig_piG = orig_model->get_PI_G();

  npiT = 1 - npiA - npiC - npiG;



  nrate = orig_model->get_ratetype();

  ntype = (nuc_model::enummodeltype)orig_model->get_modeltype();



  DEBUGCODE( cerr << "nrate: " << nrate << endl );

  DEBUGCODE( cerr << "ntype: " << ntype << endl );





  if(nrate == molecular_model<4>::ratetype_invgamma) {

    DD2( nshape = orig_shape + ranf(stdev_shape) ){}

    DEBUGCODE( cerr << "nshape: " << nshape << endl );

    DD2( ninv = orig_inv + ranf(stdev_inv) ){}

    DEBUGCODE( cerr << " ninv: " << ninv << endl );

  } 

  else if(nrate == molecular_model<4>::ratetype_propinv) {

    DD2( ninv = orig_inv + ranf(stdev_inv) ){}

    nshape = 1000000;

    DEBUGCODE( cerr << " ninv: " << ninv << endl );

  }

  else if(nrate == molecular_model<4>::ratetype_gamma) {

    ninv   = 0;

    DD2( nshape = orig_shape + ranf(stdev_shape) ){}

    DEBUGCODE( cerr << "nshape: " << nshape << endl );

  }

  else

  {

    nshape = 1000000;

    ninv   = 0;

  }



  if(ntype != nuc_model::JC) 

  {

    if(ntype != nuc_model::K2P)

    {

      do{

	DD1( npiA = orig_piA + ranf(stdev_A) ){};

	DEBUGCODE( cerr << "npiA: " << npiA << endl );

	DD1( npiC = orig_piC + ranf(stdev_C) ){};

	DEBUGCODE( cerr << "npiC: " << npiC << endl );

	DD1( npiG = orig_piG + ranf(stdev_G) ){};

	DEBUGCODE( cerr << "npiG: " << npiG << endl );

	npiT = 1 - (npiA + npiC + npiG);

	DEBUGCODE( cerr << "npiT: " << npiT << endl );

      } while(npiA < 0.1 || npiC < 0.1 ||npiG < 0.1 ||npiT < 0.1);  //kleinere Werte ergeben keinen biologischen Sinn

    }



    if(ntype != nuc_model::F81 && ntype != nuc_model::GTR) {

      DD1( ntstv = orig_tstv + ranf(stdev_tstv) ){};

      DEBUGCODE( cerr << "ntstv: " << ntstv << endl );

    }

    if(ntype == nuc_model::GTR) {

      DD_rates( nAC = orig_AC + ranf(stdev_rAC) ){}

      DEBUGCODE( cerr << "nAC: " << nAC << endl );



      DD_rates( nAG = orig_AG + ranf(stdev_rAG) ){}

      DEBUGCODE( cerr << "nAG: " << nAG << endl );



      DD_rates( nAT = orig_AT + ranf(stdev_rAT) ){}

      DEBUGCODE( cerr << "nAT: " << nAT << endl );



      DD_rates( nCG = orig_CG + ranf(stdev_rCG) ){}

      DEBUGCODE( cerr << "nCG: " << nCG << endl );



      DD_rates( nCT = orig_CT + ranf(stdev_rCT) ){}

      DEBUGCODE( cerr << "nCT: " << nCT << endl );



      DD_rates( nGT = orig_GT + ranf(stdev_rGT) ){}

      DEBUGCODE( cerr << "nGT: " << nGT << endl );

    }

  }



  cerr.setf(ios::fixed);

  cerr.precision(4);

  cerr << "old model: " << orig_model->get_modelname() << " tstv " << orig_tstv<< " r " << orig_AC<< " " << orig_AG<< " " << orig_AT<< " " << orig_CG<< " " << orig_CT<< " " << orig_GT<< " a " << orig_shape<< " i " << orig_inv<< " f " << orig_piA<< " " << orig_piC<< " " << orig_piG<< " " << (1-orig_piG-orig_piC-orig_piA) << endl;



  cerr << "new model: " << modname                     << " tstv " << ntstv<< " r " << nAC<< " " << nAG<< " " << nAT<< " " << nCG<< " " << nCT<< " " << nGT<< " a " << nshape<< " i " << ninv<< " f " << npiA<< " " << npiC<< " " << npiG<< " " << npiT << endl;



  return new nuc_model(modname, ntype, nrate, ntstv, nAC, nAG, nAT, nCG, nCT, nGT, nshape, ninv, orig_ncat, npiA, npiC, npiG, npiT);

}



void pvartree::print_varied_tree(ostream& os) {

  print_next(tree.get_treeroot(), os);

}



void pvartree::print_next(BasicNode* node, ostream& os){

  if(node->get_num_children() != 0) {

    os << "(";

  }

  for(unsigned i = 0; i < node->get_num_children(); i++) {

    print_next(node->get_childlist()[i], os);

    if(i < node->get_num_children() - 1) {

      os << "," ;

    }

  }

  if(node->get_num_children() != 0) {

    os << ")";

  }

  os << node->get_name() << ":" << node->get_branchlength();

  os << "[&" << node->get_model_name() << "]" ;

}



void pvartree::print_models(ostream& os) {

  os << "# Standard deviations used to vary the startmodel: " << endl;

  os << "# AC :    " << stdev_rAC << endl;

  os << "# AG :    " << stdev_rAG << endl;

  os << "# AT :    " << stdev_rAT << endl;

  os << "# CG :    " << stdev_rCG << endl;

  os << "# CT :    " << stdev_rCT << endl;

  os << "# GT :    " << stdev_rGT << endl;

  os << "# tstv :  " << stdev_tstv << endl;

  os << "# shape : " << stdev_shape << endl;

  os << "# inv :   " << stdev_inv << endl;

  os << "# piA :   " << stdev_A << endl;

  os << "# piC :   " << stdev_C << endl;

  os << "# piG :   " << stdev_G << endl;

  os << endl;



  outmodel.print(os);

}



