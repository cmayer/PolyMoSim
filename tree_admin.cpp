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

#include "tree_admin.h"     // error checking still inconsisten xxxxxxxxxxxxxxxxxxxxxxxx
#include <iterator>

using namespace std;

void tree_admin::destroy(int i) { // error checking xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  treevector.erase(treevector.begin()+i);
}

void tree_admin::add(int i, double f, unsigned y, const faststring &m, const faststring &treestring) { // error checking xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  tree_struct* newtree;
  newtree = new tree_struct(f, y, m, treestring);
  treevector.insert(treevector.begin()+i,newtree); 	 
}

void tree_admin::swap(unsigned i,unsigned y) {
  if(i >= treevector.size() || y >= treevector.size()) {
    throw indexerror();
  }	  
  tree_struct* temp;
  temp = treevector[i];
  treevector[i] = treevector[y];
  treevector[y] = temp;
}

void tree_admin::create(const char* filename) {
  CFile is;

  if(!FileExists(filename)) {
    faststring errormsg = faststring("File: ") + filename + " does not exist";
    throw readerror(0, errormsg);
  }

  treevector.clear();

  is.ffopen(filename);

  tree_struct*    temp;
  vector<faststring>  token;
  faststring          line;
  unsigned        number_of_tokens;

  ignore_spaces(is);
  is.getline(line);
  line.shorten_to_first_occurrence_of('#');
  line.removeSpacesBack();

  while(!is.fail()) {
    if (line[0] == '#')
    {
       ignore_spaces(is);
       is.getline(line);
       line.shorten_to_first_occurrence_of('#');
       line.removeSpacesBack();
       continue;
    }
    temp = new tree_struct;  // error checking xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    number_of_tokens = split_respect(token, line);
    if (number_of_tokens != 4)
    {
      throw readerror(is.line(), "Bad tree record. Wrong number of Elements. "
		      "I expect on one line:\nscale-factor alignment-length default-model newick-tree");
    }
    
    temp->scalingfactor = atof(token[0].c_str());  //evtl. mit token[0].ToDouble
    temp->partitionSize = atoi(token[1].c_str());  //evtl. mit token[1].ToInt
    temp->default_modelName = token[2];
    temp->default_modelName.ToLower();
    temp->treeString = token[3];
    treevector.push_back(temp);

    ignore_spaces(is);
    is.getline(line);
    line.shorten_to_first_occurrence_of('#');
    line.removeSpacesBack();
  }
}

void tree_admin::print(ostream& os, unsigned flag, const char * pretext) {
  if (flag == 0) // Print in the PolyMoSim tree file format
  {
    for(unsigned i = 0; i < treevector.size(); i++) {
      os << pretext
	 << treevector[i]->scalingfactor << " "
	 << treevector[i]->partitionSize << " "
	 << treevector[i]->default_modelName << " "
	 << treevector[i]->treeString << endl;
    }
  }
  else if (flag == 1) // for header in nexus or log file
  {
    os << pretext
       << "Partitions of data sets:" << endl
       << pretext
       << "------------------------" << endl;
    for(unsigned i = 0; i < treevector.size(); i++) {
      os << pretext <<   "Partition " << i+1 << ":" << endl
	 << pretext <<   "  scaling factor:      " << treevector[i]->scalingfactor << endl
	 << pretext <<   "  partition size:      " << treevector[i]->partitionSize << endl
	 << pretext <<   "  default start model: " << treevector[i]->default_modelName << endl;
      if ( treevector[i]->treeData != NULL)
      {
	os << pretext << "  tree with models:    ";
	treevector[i]->treeData->output(os, 0);
	os << endl;
      }
      else
      {
	os << pretext << "  Tree data is NULL at this point. This can be the expected status." << endl;
      }
    }
  }
}

void tree_admin::default_treevector(int s) {
  tree_struct* temp;
  treevector.clear();
  for(int i = 0; i < s; i++) {
    temp = new tree_struct;
    treevector.push_back(temp);
  }
}

// tree_struct*  tree_admin::getTree(int i) const { // get_tree // should return const ?? xxxxxxxxxxxxxxxxxxxxxx
//   if(i >= treevector.size()) {
//     throw indexerror();
//   }else{
//     return treevector[i];
//   }
// }

// const faststring*  tree_admin::getTreeString(int i) const { // get_tree // should return const ?? xxxxxxxxxxxxxxxxxxxxxx
//   if(i >= treevector.size()) {
//     throw indexerror();
//   }else{
//     return treevector[i];
//   }
// }

void tree_admin::set(int i, double f, unsigned y, const faststring &model, const faststring &treestring) {
  treevector[i]->scalingfactor = f;
  treevector[i]->partitionSize = y;
  treevector[i]->default_modelName = model;
  treevector[i]->treeString = treestring;
}
 	 
void tree_admin::ini(int s) {
  tree_struct *temp;
  for(int i = 0; i < s; ++i) {
    temp = new tree_struct;
    treevector.push_back(temp);
  }	 
}

unsigned tree_admin::getAlignmentLength() const
{
  unsigned i, length=0;

  for (i=0; i < treevector.size(); ++i)
  {
    length += treevector[i]->partitionSize;
  }
  return length;
}

unsigned tree_admin::getPartitionStartCoordinate(unsigned p) const
{
  unsigned i, length=0;

  for (i=0; i < p; ++i)
  {
    length += treevector[i]->partitionSize;
  }
  return length;
}


faststring tree_admin::get_dataTypeString() const
{
  int              i, N = treevector.size();
  vector<faststring>   dt_vec;
  int              partition_start, partition_end;
  bool             same_data_type = true;
  faststring           res;

  for (i=0; i<N; ++i)
  {
    dt_vec.push_back(treevector[i]->treeData->get_root_model()->get_datatypename());
    if (i > 0 && (dt_vec[i-1] != dt_vec[i]) )
      same_data_type = false;
  }

  if (same_data_type)
    return dt_vec[0];
 
  partition_start = 1;
  partition_end   = 0;

  res = "mixed(";
  for (i=0; i<N; ++i)
  {
    partition_end   += treevector[i]->partitionSize;

    res.append(dt_vec[i].c_str());
    res.push_back(':');
    res.append(faststring(partition_start));
    res.push_back('-');
    res.append(faststring(partition_end));

    if (i < N-1)
    {
      res.push_back(',');
      partition_start += treevector[i]->partitionSize;
    }
  }
  res.push_back(')');

  return res;
}


double tree_admin::getScalingFactor(unsigned i) const
{
  return treevector[i]->scalingfactor;
}

void tree_admin::set_itsBasicTree(unsigned i, BasicTree *p) // should 2nd parameter better be const?
{
  treevector[i]->treeData = p; // In future versions it should not be possible to overwrite a previous setting.
}

BasicTree* tree_admin::get_itsBasicTree(unsigned i)
{
  return treevector[i]->treeData;
} 

unsigned tree_admin::getPartitionSize(unsigned i) const
{
  return treevector[i]->partitionSize;
}

const faststring& tree_admin::getDefaultModelName(unsigned i)  const
{
  return treevector[i]->default_modelName;
}

const faststring& tree_admin::getTreeString(unsigned p) const
{
  return treevector[p]->treeString;
}

bool tree_admin::equal_takon_sets() const
{
  BasicTree::set_of_OTUs first_set, curr_set;
  unsigned    i, N = treevector.size();
  unsigned    size1;

  BasicTree::set_of_OTUs::iterator it1;
  BasicTree::set_of_OTUs::iterator it2;
  BasicTree::set_of_OTUs::iterator it1_end;
  BasicTree::set_of_OTUs::iterator it2_end;

  if ( treevector.empty() )
    return true;

  treevector[0]->treeData->get_set_of_OTUs( first_set );
  it1     = first_set.begin();
  it1_end = first_set.end();
  size1   = first_set.size();
  
//   cout << "Taxa in tree 0:" << endl;
//   for(it = last_set.begin();
//       it != last_set.end(); ++it)
//     cout << **it << " ";
//   cout << endl;

  for ( i=1; i < N; ++i )
  {
    treevector[i]->treeData->get_set_of_OTUs( curr_set );

    it1     = first_set.begin();
    it2     = curr_set.begin();
    it2_end = curr_set.end();

//     cout << "Taxa in tree " << i << ":" << endl;
//     for(it = curr_set.begin();
// 	it != curr_set.end(); ++it)
//       cout << **it << " ";
//     cout << endl;

    // The comparison operators of sets do not use the specified
    // comparison oject, but simply compares the values of the pointers.
    // So we have to do it for ourselves.
    if (size1 != curr_set.size() )
      return false;

    while (it1 != it1_end)
    {
      if ( **it1 != **it2)
	return false;
      ++it1;
      ++it2;
    }
  }
  return true;
}
