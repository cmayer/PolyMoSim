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

#include "model_admin.h"
#include <iterator>


using namespace std;

void model_admin::create(const char *filename) {
  CFile is;

  if(!FileExists(filename)){
    throw readerror();
  }

  // (Old?) Memory leak   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  // reset funktion
  modelvector.clear();
  modelmap.clear();

  is.ffopen(filename);

  faststring      start;
  basic_model     *temp;

  ignore_spaces(is);
  is.getline(start);
  start.shorten_to_first_occurrence_of('#');
  start.removeSpacesBack();

  while(!is.fail()) {
    if (start.size() == 0 ||  start[0] == '#')
    {
      ignore_spaces(is);
      is.getline(start);
      start.shorten_to_first_occurrence_of('#');
      start.removeSpacesBack();
      continue;
    }
    if (start == "begin nuc-model" ||start == "begin model" ) {
      temp = new nuc_model("");
      temp->read_next_model(is);
      ignore_spaces(is);
      append_a_copy((nuc_model*)temp);
      delete temp;
    }
    else if (start == "begin aa-model") {
      temp = new aa_model("");
      temp->read_next_model(is);
      ignore_spaces(is);
      append_a_copy((aa_model*)temp);
      delete temp;
    }
    else {
      throw formaterror();
    }
    ignore_spaces(is);
    is.getline(start);
    start.shorten_to_first_occurrence_of('#');
    start.removeSpacesBack();
  }
}




void model_admin::print(ostream& os, unsigned flag) {
  for(unsigned i = 0; i < modelvector.size(); i++) {
    modelvector[i]->print(os, flag);
  }
}

void model_admin::print(FILE *os, unsigned flag) {
  for(unsigned i = 0; i < modelvector.size(); i++) {
    modelvector[i]->print(os, flag);
  }
}

void model_admin::default_modelvector(unsigned s) {
  basic_model     *temp;
  faststring s1 = "Model_";
  faststring s2;

  // Memory leak   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  // reset funktion  // Models have to be deleted !!!!!!!!!!!!
  modelvector.clear();
  modelmap.clear();

  for(unsigned i = 0; i < s; i++) {
    s2 = s1 + faststring(i);
    temp = new nuc_model(s2);

    append_a_copy((nuc_model*)temp);
    delete temp;
  }
}

unsigned model_admin::append_new_default_model() {
  basic_model     *temp;
  faststring s1 = "Untitled";

  temp = new nuc_model(s1);
  append_a_copy((nuc_model*)temp);
  delete temp;

  return getNumModels() - 1; // Index of new model.
}

const basic_model* model_admin::get_model(unsigned i) const {
  if(i >= modelvector.size()) {
    faststring errormsg = "Internal fatal error when accessing the vector of models with index " + faststring(i) + ".";
    throw indexerror(errormsg);
  } else {
    return modelvector[i];
  }
}

const basic_model* model_admin::get_model(faststring str) {

  if (modelNameExists(str)) {
    str.ToLower();
    return modelmap[str];
  }
  else if ( str == "default" ) { //if there is no model named "default", use the first model of the file...
    return modelvector[0];
  }
  else {
    faststring errormsg = "Error: Unknown model: " + str + ".";
    throw indexerror(errormsg);
  }
}

void model_admin::init_models(double (*r_gamma)(double, double), double (*r_lfco)())
{
  if (global_verbosity > 3)
  {
    cerr << modelvector.size() << " models have been found in the model file, which we initialze now." << endl;
  }
  for(unsigned i = 0; i < modelvector.size(); ++i)
  {
    if (global_verbosity > 3)
    {
      cerr << "Initializing model: " << modelvector[i]->get_modelname() << endl;
    }
    modelvector[i]->init_model(r_gamma, r_lfco);
    modelvector[i]->set_matrices();
  }
}


void model_admin::reset_siterates()
{
  for(unsigned i = 0; i < modelvector.size(); ++i)
  {
    if (global_verbosity > 3)
    {
      cerr << "Resetting siterates for model: " << modelvector[i]->get_modelname() << endl;
    }
    modelvector[i]->reset_siterates();
  } 
}

void model_admin::init_siterates(unsigned len, bool reinit) 
{
  // reinit is not yet implemented. The variable is mainly a placeholder for the idea to do that in the future.

  if (global_verbosity > 3)
    cerr << "Call to function: init_siterates. Found " << modelvector.size() << " models to initialize." << endl;

  if (reinit) // should never be true in current implementation
  {
    if (global_verbosity > 3)
      cerr << "init_siterates: reinit mode" << endl;
  }
  else // No reinit: We will also reset the siterates if they have been set before
  {
    if (global_verbosity > 3)
      cerr << "init_siterates: first initialisation or complete reinit with reallocation of memory." << endl;

    //    reset_siterates();  // Currently this is called in the main function.
    // Since reinit is not yet implemented, this should not be called here.
  }

  // Debug code only:
  //  if (!reinit)
  {
    // Debug code: List all model and from where they inherit their siterates from:
    if (global_verbosity > 3)
    {
      cerr << "Model names and the name of the model they inherit siterates from:" << endl;
      for(unsigned i = 0; i < modelvector.size(); ++i)
      {
	cerr << modelvector[i]->get_modelname() << " -> ";
	if (modelvector[i]->get_modelname_siterates_are_inherited_from() == "")
	  cerr << "(no inheritance)" << endl;
	else
	  cerr << modelvector[i]->get_modelname_siterates_are_inherited_from() << endl;
      }
    }
  }

  //  cerr << "We are here 1." << endl; 

  // Initialisation of all site rates that do not inherit from other models
  {
    // Init all siterates that can/have to be initialized at the beginning, i.e. all
    // siterates that do not inherit from others.
    for (unsigned i = 0; i < modelvector.size(); ++i) 
    {
      // Does this model inherit (parts) of the siterates?
      // Only of not we can initialize it now.
      if (modelvector[i]->get_modelname_siterates_are_inherited_from() == "")
      {
	if (global_verbosity > 3)
	  cerr << "Direct initialization of siterates of model: " << modelvector[i]->get_modelname() << endl;
	modelvector[i]->init_siterates(len, reinit);
      }
      else
      {
	if (global_verbosity > 3)
	  cerr << "Initialization has been posponed due to inheritance of siterates for model: " << modelvector[i]->get_modelname() << endl;
      }
    }
  } // End block init all independent siterate.

  ///////////////////////////////////////////////////////

  // reinit is not yet used in the code below - since we do not make real use of it yet.
  // If reinit should be used in the future, this block should be revised.

  {
    // Connect models which inherit the siterate from other models
    unsigned i;
    unsigned j;
    unsigned unconnected=0;

    // We need a double loop since only after N*N iterations all dependencies are resolved for sure.
    // All models might depend on the last model and so forth. This is only important for setting the siterates pointer.
    for (i = 0; i < modelvector.size(); ++i)
    {
      unconnected = 0;
      for(j = 0; j < modelvector.size(); ++j) 
      {
	// Wants to be connected and is still unconnected:
	// In principle, the linking of the model does not have to be done in the double loop.
	// If models had been linked before and have been reset, we won't get into this if any more.
	if ( modelvector[j]->get_modelname_siterates_are_inherited_from() != "" && 
	     modelvector[j]->get_and_set_model_to_inherit_siterates_from()==NULL )
	{
	  const faststring& modelname_to_inherit_from = modelvector[j]->get_modelname_siterates_are_inherited_from();

	  // Connect models:
	  basic_model *tmp_basic_mod = find_model_with_name(modelname_to_inherit_from);

	  if (tmp_basic_mod == NULL)
	  {
	    cerr << "Error in model specification file: The model " << modelvector[j]->get_modelname() << endl;
	    cerr << "inherits from model " << modelvector[j]->get_modelname_siterates_are_inherited_from()
		 <<  " which could not be found in the list of madels in this model file."
	         << endl;
	    exit(-13);
	  }

	  modelvector[j]->get_and_set_model_to_inherit_siterates_from() = tmp_basic_mod;


	  if (modelvector[j]->get_and_set_model_to_inherit_siterates_from() == NULL)
	  {
	    cerr << "Internal error: Pointer to model this model inherits from is NULL which should not be the case here." << modelvector[j]->get_modelname() << endl;
	    exit(0);
	  }

	  // XXXXXXXXX Not OK yet - check that all are connected!!!!
	  if (!modelvector[j]->set_inherited_site_rates(len) )
	    ++unconnected;
	}
	else // Needs no connection of models, but might need linking siterates.
	{
	  if ( !modelvector[j]->siterates_initialized() )
	  {
	    if (!modelvector[j]->set_inherited_site_rates(len) )
	      ++unconnected;
	  }
	}
      } // End inner for loop
      if (unconnected == 0) // We break out of outer loop if no unconnected models remain.
	break;
    }

    {
      if (unconnected)
      {
	cerr << "Error in the model specifications. Inheritence of siterates could not be resolved.\n"
                "Check the following models for unresolved dependencies, e.g. in the form of circular dependencies." << endl;

	for(unsigned i = 0; i < modelvector.size(); ++i) 
	{
	  // Does this model inherit (parts) of the siterates?
	  // Only of not we can initialize it now.
	  if (modelvector[i]->get_modelname_siterates_are_inherited_from() != "" &&
	      !modelvector[i]->siterates_initialized() )
	  {
	    cerr << "Model with unresolved dependency: " << modelvector[i]->get_modelname() << endl;
	  }
	}
	exit(0);
      }
    }
    if (global_verbosity > 3)
      cerr << "Initialisation of site rates has been completed with success." << endl;
  }



}

bool model_admin::modelNameExists(faststring str) const
{
  str.ToLower();
  if (modelmap.find(str) != modelmap.end())
    return true;
  else
    return false;
}

basic_model* model_admin::find_model_with_name(faststring str)
{
  str.ToLower();
  map<faststring, basic_model*>::iterator findres = modelmap.find(str);

  if (findres != modelmap.end())
    return findres->second;
  else
    return NULL;
}


void model_admin::create_unique_model_name(faststring &s1) const
{
  faststring s2;
  int    i;

  if ( !modelNameExists(s1) )
    return;

  s2 = s1 + "_1";
  i  = 1;

  while ( modelNameExists(s2) )
  {
    ++i;
    s2 = s1 + "_" + faststring(i);
  }
  s1 = s2;
}

void model_admin::append_a_copy(const nuc_model *newmodel, faststring str) {

  if (str == "" )
    str = newmodel->get_modelname();
  create_unique_model_name(str);

  basic_model *tmp = new nuc_model(str, *newmodel);
  modelvector.push_back(tmp);
  str.ToLower();
  modelmap[str] = tmp;
}

void model_admin::append_a_copy(const aa_model *newmodel,  faststring str) {

  if (str == "" )
    str = newmodel->get_modelname();
  create_unique_model_name(str);

  basic_model *tmp = new aa_model(str, *newmodel);
  modelvector.push_back(tmp);
  str.ToLower();
  modelmap[str] = tmp;
}

void model_admin::duplicate_model(unsigned i)
{
  if (i >= modelvector.size() )
  {
    faststring errormsg = "Error: Call to duplicate model with index out of range.";
    throw indexerror(errormsg);
  }

  basic_model *orig_model = modelvector[i];
  basic_model *copy_model;
  faststring      str         = orig_model->get_modelname();

  str = str + "_copy";
  create_unique_model_name(str);

  if (orig_model->get_datatype() == 0)
  {
    copy_model = new nuc_model(str, *((nuc_model*)orig_model) );
  }
  else
  {
    copy_model = new aa_model(str, *((aa_model*)orig_model) );
  }

  modelvector.insert(modelvector.begin()+i+1, copy_model);
  str.ToLower();
  modelmap[str] = copy_model;
}

void model_admin::delete_model(unsigned i)
{
  if (i >= modelvector.size() )
  {
    faststring errormsg = "Error: Call to delete model with index out of range.";
    throw indexerror(errormsg);
  }

  faststring name = modelvector[i]->get_modelname();
  delete modelvector[i];
  modelvector.erase(modelvector.begin()+i);
  name.ToLower();
  modelmap.erase(modelmap.find(name) );
}


void model_admin::swap_models(unsigned i1, unsigned i2)
{
  if (i1 >= modelvector.size() || i2 >= modelvector.size() )
  {
    faststring errormsg = "Error: Call to swap models with index out of range.";
    throw indexerror(errormsg);
  }

  basic_model *tmp;
 
  tmp             = modelvector[i1];
  modelvector[i1] = modelvector[i2];
  modelvector[i2] = tmp;
}

unsigned model_admin::getNumModels() const {
  return modelvector.size();
}

faststring model_admin::getModelName(unsigned i) const {
  return modelvector[i]->get_modelname();
}

void model_admin::changeModelName(unsigned i, faststring str)
{
  basic_model *tmp;

  if (i >= modelvector.size() )
  {
    faststring errormsg = "Error: Call to changeModelName with index out of range.";
    throw indexerror(errormsg);
  }

  faststring str_orig = modelvector[i]->get_modelname();

  if (str == str_orig)
    return;

  faststring str_orig_lower = str_orig;
  faststring str_lower      = str;

  str_orig_lower.ToLower();
  str_lower.ToLower();

  if (str_lower != str_orig_lower)
  {
    // Copy the model with new name:
    create_unique_model_name(str);
  }

  if (modelvector[i]->get_datatype() == 0 ) // DNA
  {
    tmp = new nuc_model(str, *((nuc_model*)(modelvector[i])) );
  }
  else
  {
    tmp = new aa_model(str,  *((aa_model*)(modelvector[i])) );
  }

  faststring oldname = modelvector[i]->get_modelname();
  oldname.ToLower();
  modelmap.erase(modelmap.find(oldname) );
  delete modelvector[i];

  modelvector[i] = tmp;
  str.ToLower();
  modelmap[str]  = tmp;
}

/*
// set_model might be depricated. Does not handle use_distfunction
void model_admin::append_new_model(char   data_type,             // 'n' for nuc or 'p' for protein
                                   faststring modelname_param,
                                   faststring modeltype_param,
                                   vector<double> *rrates_param, // Parameters are supplied as in an upper triangular matrix.
                                   vector<double> *base_param,
                                   double         shape_param,
                                   double         pinv_param,
                                   unsigned       ncat_param,
                                   double         *tstv_param

)
{   basic_model *tmp;

    if (modelname_param == "" )
      modelname_param = "Model_1";
    create_unique_model_name(modelname_param);

    if (data_type == 'n')
        tmp = new nuc_model(modelname_param);
    else if (data_type == 'p')
        tmp = new aa_model(modelname_param);
    else // Handle unknown model
    {

    }

    tmp->set_model(modelname_param,
                   modeltype_param,
                   rrates_param, // Parameters are supplied as in an upper triangular matrix.
                   base_param,
                   shape_param,
                   pinv_param,
                   ncat_param,
                   tstv_param);

    modelvector.push_back(tmp);
    modelname_param.ToLower();
    modelmap[modelname_param] = tmp;
}
*/


void model_admin::print_relative_site_rates(std::ostream &os) const
{
  for(unsigned i = 0; i < modelvector.size(); ++i)
  {
    modelvector[i]->print_relative_site_rates(os);
  }

}


void model_admin::print_site_rates_histogramm_data(std::ostream &os, faststring &desc) const
{
  for(unsigned i = 0; i < modelvector.size(); ++i)
  {
    modelvector[i]->print_site_rates_histogramm_data(os, desc);
  }
}
