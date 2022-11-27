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

  unsigned long      nodes;

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

