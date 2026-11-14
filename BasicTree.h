/***************************************************************************************************
*  This file is part of PolyMoSim.
*  
*  Copyright (c) 2004-2026, Christoph Mayer
*  Contributors: see README.md
*  
*  Author: Dr. Christoph Mayer
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

#ifndef BasicTreeH
#define BasicTreeH

#include "PolyMoSim.h"

#include <cctype>
#include <set>
#include <vector>
#include "BasicNode.h"
#include "model_admin.h"
#include "BSplit.h"
#include "split_admin.h"
#include "CTaxa.h"

struct lessThan_String_ptr
{
  bool operator()(const faststring *pa, const faststring *pb) const
  {
    faststring a = *pa, b = *pb;
    a.ToLower(); b.ToLower();
    return a < b;
  }
};

class tree_parser_assistent;

class BasicTree {

 public:
  class readerror {
  };

  typedef std::map<const faststring *, const faststring *, lessThan_String_ptr> map_of_OTUs;
  typedef std::set<const faststring *,                 lessThan_String_ptr> set_of_OTUs;
 
 private:
  double             scalefactor;
  BasicNode          *treeroot;
  unsigned           OTU;
  unsigned           nodes;
  bool               read_model_status; // Do we read models in nodes if available
  model_admin        *model_master;
  const basic_model  *root_model;       // root model determines the data type and is used if no other model is specified at the root of the tree.
  vector<BSplit>     split_vec;
  CTaxa              taxa_list;

  void               get_map_of_OTUs_intern(map_of_OTUs &, BasicNode *);
  void               get_set_of_OTUs_intern(set_of_OTUs &, BasicNode *);


 public:
  BasicTree();
  ~BasicTree();

  void               read_tree(istream&, double);
  void               link_models_to_tree(model_admin*, const basic_model*);
  void               link_next_model(BasicNode*);
  void               evolve_tree();
  void               get_tree_splits(split_admin&);

 private:
  void               evolve_next_node(BasicNode*);
  void               get_next_split(split_admin&, BasicNode*);

 public:
  void               read_model_name(tree_parser_assistent&, BasicNode*);
  void               set_read_model_status(bool);
  bool               get_read_model_status();

  void               output(ostream& = cout, unsigned = 0);
  void               output_node(BasicNode *, ostream&);
  void               destroy_subtree(BasicNode *);

  void               setScaleFactor(double sf) { scalefactor = sf; }
/*   void            set_root_model(faststring); */
  const basic_model* get_root_model();
  BasicNode*         get_treeroot() const;
  unsigned           get_num_nodes();
  unsigned           get_num_OTUs();
  void               get_map_of_OTUs(map_of_OTUs &);
  void               get_set_of_OTUs(set_of_OTUs &);

 private:
  void               push_child(BasicNode *, BasicNode *);
  BasicNode*         new_node(tree_parser_assistent&);
  BasicNode*         new_OTU(tree_parser_assistent&);
  void               get_name_and_branchlength(tree_parser_assistent&, BasicNode *);
  void               get_name(tree_parser_assistent&, BasicNode *);


};

inline void BasicTree::set_read_model_status(bool b) {
  read_model_status = b;
}

inline bool BasicTree::get_read_model_status() {
  return read_model_status;
}

inline BasicNode* BasicTree::get_treeroot() const{
  return treeroot;
}

inline unsigned BasicTree::get_num_nodes() {
  return nodes;
}

inline unsigned BasicTree::get_num_OTUs() {
  return OTU;
}

#endif
