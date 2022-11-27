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

#include "BasicTree.h"

#include "PolyMoSim.h"



using namespace std;



class tree_parser_assistent {  // treat file read fail xxxxxxxxxxxxxxxxxxxxxxxx

  private:

  istream&        is;

  //  int             last_char;

  vector<faststring>  special_comments_vec;



  public:

  tree_parser_assistent(istream& pis):is(pis) {}



  int getchar(bool look_for_special_comments = false) {

    int ch;



    ignore_spaces(is);

    ch = is.get();

    while (ch == '[' ) {

      if (look_for_special_comments)

      {

	ch = is.get();

	if (ch == '&')

	{

          faststring s;



	  ignore_spaces(is);

	  ch = is.get();

	  while (ch != EOF && ch != ']' ) {

	    s.push_back(ch);

	    ch = is.get();

	  }

	  // remove trailing spaces from s xxxxxxxxxxxxxxxxxxxxxxxxxx

	  special_comments_vec.push_back(s);

	}

	else {

	  //	  is.unget();

	  skip_comment(is);

	}

      }

      else {

	skip_comment(is);

      }

      ignore_spaces(is);

      ch = is.get();

    }

    return ch;

  }



  void ungetchar()

  {

    is.unget();

  }



  int getrawchar()

  {

    return is.get();

  }



  void clear_special_comments()

  {

    special_comments_vec.clear();

  }



  faststring get_special_comment()

  {

    if (special_comments_vec.size() != 1)

      return "";

    else

      return special_comments_vec[0];

  }



  int get_number_of_special_comments()

  {

    return special_comments_vec.size();

  }



  double read_double()   // treat file read fail xxxxxxxxxxxxxxxxxxxxxxxx

  {

    double lf;



    is >> lf;

    return lf;

  }



  int read_int()

  {

    int d;



    is >> d;

    return d;

  }



};





BasicTree::BasicTree() {

  treeroot = NULL;

  OTU      = 0;

  nodes    = 0;

  read_model_status = false;

  root_model = NULL;           // "default";

}



BasicTree::~BasicTree() {

  if(treeroot)

    destroy_subtree(treeroot);

}





void BasicTree::read_tree(istream& is, double s) {  // istream can also be a stringstream

  int                     ch;

  tree_parser_assistent   tpa(is);



  scalefactor = s;

  ch = tpa.getchar();   // ignores leading spaces and comments, ignore special comments

  if(ch == '(' )

    treeroot = new_node(tpa);

}



void BasicTree::push_child(BasicNode *parent, BasicNode *child) {

  parent->push_back_child(child);

  child->set_parent(parent); 

}



// Can it parse that one: ??

// (            ## new_node: root, 2 children

//  (           ## new_node: d,    1 child

//   (          ## new_node: c,    1 child

//    (         ## new_node: b,    1 child

//     a:0.1    ## new_OTU:  a,    0 children

//    )b:0.1

//   )c:0.1

//  )d:0.1

//  ,

//  A:0.3       ## new_OTU:  A,    0 children

// )root        ## End of tree

//

//

// If an OTU has braces arround it (a:0.1) this is treated as 2 nodes

// equivalent to (a:0.1):0.



BasicNode* BasicTree::new_node(tree_parser_assistent &tpa) {



  nodes++;

  

  int         children = 0;

  BasicNode   *newnode;

  BasicNode   *newchild = NULL;

  int         ch;



  newnode = new BasicNode;



  do {

    ch = tpa.getchar();   // ignores spaces and comments in front of a node

    if (ch == '(' ) {

      ++children;

      newchild = new_node(tpa);

      push_child(newnode,newchild);

      newchild = NULL;

    } else {

      ++children;

      tpa.ungetchar();

      newchild = new_OTU(tpa);

      push_child(newnode, newchild);

      newchild = NULL;

    }

    ch = tpa.getchar();   // ignores spaces and comments behind node

  } while (ch == ',' );



  if (ch == ')' ) {

    /* ch = mygetchar(is);

       if(ch == '(') {

       is.unget();

       newchild = new_node(is);

       push_child(newnode->get_parent(), newchild);

       newchild = NULL;

       }else{*/



    get_name_and_branchlength(tpa, newnode);

    if (read_model_status) {

      read_model_name(tpa, newnode);

    }

  }

  return newnode;

}



BasicNode* BasicTree::new_OTU(tree_parser_assistent& tpa) {

  BasicNode   *newnode;



  newnode = new BasicNode;



  get_name_and_branchlength(tpa, newnode);

  if (read_model_status) {

    read_model_name(tpa, newnode);

  }



  // save OTU as a taxon in the Taxa-Class

  taxa_list.add_taxon( newnode->get_name() );



  ++OTU;

  ++nodes;   // An OTU is also a node xxxx

  return newnode;

}



void BasicTree::get_name_and_branchlength(tree_parser_assistent &tpa, BasicNode *node) {

  int     ch;



  get_name(tpa, node);

  ch = tpa.getchar(true);   // Check whether special comment has been found

  if(ch == ':' ) {

    node->set_branchlength(tpa.read_double() );

  } else {

    // default branch length is 0

    tpa.ungetchar();

  }

}





void BasicTree::get_name(tree_parser_assistent &tpa, BasicNode *node) {

  int   ch;



  ch = tpa.getchar(true);             // Check whether special comment has been found

//   if (tpa.get_number_of_special_comments() > 0 )

//   {

//     tpa.ungetchar();

//     return;

//   }

  if (ch == '"' ) {

    do {

      node->push_back_to_name(ch);

      ch = tpa.getrawchar();     // No true comments are allowed in "..."

    } while (ch != EOF && ch != '"' );

    if (ch != EOF) node->push_back_to_name(ch);

  } else if (ch == EOF) {

    tpa.ungetchar();

  } else {

    while (ch != EOF && !isspace(ch) && ch != ':' && ch != '(' && ch != '[' && ch != ',' && ch != ']' && ch != ')' && ch != ';' ) {

      node->push_back_to_name(ch);

      ch = tpa.getrawchar();    // comments and spaces are not allowed within the node name

    }

    tpa.ungetchar();

  }

}



void BasicTree::read_model_name(tree_parser_assistent &tpa, BasicNode* node) {



    tpa.getchar(true);

    tpa.ungetchar();



  if (tpa.get_number_of_special_comments() == 0)

  {

    return;

  }



  if (tpa.get_number_of_special_comments() > 1)

  {

    // Error xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    cerr << "Found multiple models on one node. Could be an internal problem. Please report this problem." << endl;

    exit(0);

    //    return;

  }



  faststring s = tpa.get_special_comment();

  tpa.clear_special_comments();



  //  if ( tolower(s[0]) == 'm')

  {

    node->get_model_name() = s;

    node->set_new_model_status(true);

  }

  //  else

  //    return;

}



void BasicTree::link_models_to_tree(model_admin *ma, const basic_model *r_model) {

  model_master = ma;

  root_model   = r_model;

  link_next_model(treeroot);

}





// When reading the tree, some nodes had been assigned a new model, which so far is only

// saved in the node in form of the model name.

// This recursive routine does the following:

// - nodes which are not marked as "new model" get the model name of their parent. So all nodes get the correct model name.

// - all nodes get a pointer their model object.

void BasicTree::link_next_model(BasicNode* node) {

  const basic_model    *temp;

        BasicNode      *parent;

        faststring      dummy;



  if (node->get_parent() == NULL) { // root

    if ( node->get_new_model_status() ) {

      dummy = node->get_model_name();

      temp  = model_master->get_model(dummy);

      // Check correct data type for each new model:

      if ( root_model->get_datatype() != temp->get_datatype() )

      {

	cerr << "Error: The data type of the model on this node is not compatible with the data type of the sequences that evolves here." << endl;

	throw 1;

      }

    }

    else

    {

      temp = root_model;

    }

  }

  else if (node->get_new_model_status() ) {

    dummy = node->get_model_name();

    temp = model_master->get_model(dummy);

    // Check correct data type for each new model:

    if ( root_model->get_datatype() != temp->get_datatype() ) {

      cerr << "Error: The data type of the model on this node is not compatible with the data type of the sequences that evolves here." << endl;

      throw 1;

    }

  }

  else

  {

    parent = node->get_parent();

    temp = parent->get_model();

  }

  node->set_model(temp);



  for (unsigned i = 0; i < node->get_num_children(); i++) {

    link_next_model(node->get_childlist()[i]);

  }

}



void BasicTree::evolve_tree() {

  evolve_next_node(treeroot);

}



void BasicTree::evolve_next_node(BasicNode* node) {

  //  char* start_pos_new;

  //  char* end_pos_new;



  //  faststring new_seq;

  double br_length;



  if (node->get_parent() == NULL) {

    DEBUGCODE( { cerr << "Wurzel: " << endl;   cerr << node->get_sequence() << endl << endl; } );

  } else {

    faststring &parent_seq = node->get_parent()->get_sequence();

    faststring &new_seq    = node->get_sequence();

    

    //    start_pos_new = node->get_sequence().begin();

    //    end_pos_new = node-> get_sequence().end();

    //    char * blub;



    DEBUGCODE( cerr << "scalefactor:              "  << scalefactor << endl );

    DEBUGCODE( cerr << "node->get_branchlength(): "  << node->get_branchlength() << endl );



    br_length = node->get_branchlength()*scalefactor;

    node->get_model()->evolve(parent_seq, new_seq, br_length);

    DEBUGCODE( if (0)

    { 

      cerr << "Evolve branch before node:" << endl;

      cerr << "Node name:  " << node->get_name() << endl;

      cerr << "Model name: " << node->get_model()->get_modelname() << endl;

    }

	       );

    DEBUGCODE({ cerr << "Node:         " << node->get_name() << endl; cerr << "Evolved seq.: " << new_seq << endl; });

  }



  for(unsigned i = 0; i < node->get_num_children(); ++i) {

    evolve_next_node(node->get_childlist()[i]);

  }

}



//This recursive routine finds the 2n-3 splits of the tree.

//It must get a split_admin object as a parameter which is then allocated with all splits of the tree.



void BasicTree::get_tree_splits(split_admin& sp_ad) {

  // Wurzel muss vorerst nicht besonders behandelt werden



  for(unsigned i = 0; i < treeroot->get_num_children(); i++) {

    get_next_split(sp_ad, treeroot->get_childlist()[i]);

  }

}



void BasicTree::get_next_split(split_admin& sp_ad, BasicNode* node) {

  BSplit *temp;

  unsigned tax_pos;



  // Neuer Split - auf Vektor-Stack

  temp = new BSplit(taxa_list.GetTaxaNum());

  split_vec.push_back(*temp);



  // Hat keine Kinder

  //    - Taxon Nummer beschaffen

  //    - Bit in allen splits auf Vektor-Stack setzten

  

  if( node->get_num_children() == 0 ) {      // Hat keine Kinder

    tax_pos = taxa_list.GetTaxonPosition( node->get_name() );



    DEBUGCODE( cerr << "neuer Taxon Name: " << node->get_name() << endl; )

    DEBUGCODE( cerr << "neue Bits: " << endl; )



    for(unsigned j = 0; j < split_vec.size(); j++) {

      split_vec[j].get_split().set(tax_pos);

      DEBUGCODE(  cerr << j << ": " << split_vec[j].get_split() << endl; )

    }

  }



  // Hat Kinder

  for(unsigned i = 0; i < node->get_num_children(); i++) {

    get_next_split(sp_ad, node->get_childlist()[i]);

  }

  sp_ad.add_split(split_vec[split_vec.size() - 1].get_split(), node->get_branchlength() );



  DEBUGCODE(  cerr << "Neuer Split auf split_admin: " << split_vec[split_vec.size() - 1].get_split() << endl; )



  split_vec.pop_back();   

}







void BasicTree::output(ostream& os, unsigned flag) {

  if (flag == 1) // Debug mode

  {

    faststring root_model_name;



    if (root_model != NULL)

      root_model_name = root_model->get_modelname();

    else

      root_model_name = "Not specified jet.";



    os << " scalefactor: "  << scalefactor

       << " OTUs: "        << OTU

       << " nodes: "       << nodes

       << " read_model_status: " << read_model_status

       << " root_model: " << root_model_name << endl;

  }

  output_node(treeroot, os);

}



void BasicTree::output_node(BasicNode *node, ostream& os) {

  if(node->get_num_children() != 0) {

    os << "(";

  }

  for(unsigned i = 0; i < node->get_num_children(); i++) {

    output_node(node->get_childlist()[i], os);

    if(i < node->get_num_children() - 1) {

      os << ",";

    }

  }

  if(node->get_num_children() != 0) {

    os << ")";

  }

  os << node->get_name() << ":" << node->get_branchlength();

  if ( read_model_status && node->get_model_name() != "")

    os << "[&" << node->get_model_name() <<  "]";

}



void BasicTree::destroy_subtree(BasicNode *node) {

  unsigned i;

    for(i=0; i < node->get_num_children() ; ++i)

      destroy_subtree(node->get_childlist()[i]);

  delete node;

}



// void BasicTree::set_root_model(faststring r_m) {

//   root_model = r_m;

// }



const basic_model* BasicTree::get_root_model() {

  return root_model;

}



void  BasicTree::get_map_of_OTUs(map_of_OTUs &m) {

  get_map_of_OTUs_intern(m, treeroot);

}



void  BasicTree::get_map_of_OTUs_intern(map_of_OTUs &m, BasicNode *node) {

  if(node->get_num_children() == 0) {

    m.insert(make_pair(&node->get_name(), &node->get_sequence()) );

  }

  else for (unsigned i = 0; i < node->get_num_children(); ++i) {

      get_map_of_OTUs_intern(m, node->get_childlist()[i]);

    }

}



void  BasicTree::get_set_of_OTUs(set_of_OTUs &m) {

  get_set_of_OTUs_intern(m, treeroot);

}



void  BasicTree::get_set_of_OTUs_intern(set_of_OTUs &m, BasicNode *node) {

  if(node->get_num_children() == 0) {

    m.insert( &node->get_name() );

  }

  else for (unsigned i = 0; i < node->get_num_children(); ++i) {

      get_set_of_OTUs_intern(m, node->get_childlist()[i]);

    }

}





