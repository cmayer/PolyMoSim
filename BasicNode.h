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

#ifndef BasicNodeH
#define BasicNodeH

#include <iostream>
#include <cstdlib>
#include "faststring3.h"
#include "mymodel.h"
#include <map>
#include <vector>

class BasicNode{

 protected:
  int        children;
  faststring     name;
  double     branchlength;
  std::vector<BasicNode*> childlist;
  BasicNode* my_parent;

  // modelnode components
  const basic_model    *model;  //the model used by the node (can be the same at two nodes)
  bool       new_model; 
  faststring     model_name;
  faststring     sequence;


 public:
  BasicNode();
  virtual ~BasicNode(){}
 
  void                        push_back_child(BasicNode*);
  BasicNode*                  get_parent();
  void                        set_parent(BasicNode*);
  unsigned                    get_num_children() const;
  const std::vector<BasicNode*>&   get_childlist();
  const faststring&               get_name() const;
  void                        set_name(const faststring&);
  void                        push_back_to_name(char);
  double                      get_branchlength();
  void                        set_branchlength(double);

  //modelnode components
  
  const                       basic_model* get_model();
  void                        set_model(const basic_model*);
  bool                        get_new_model_status();
  void                        set_new_model_status(bool);
  faststring&                     get_model_name();
  void                        set_model_name(const faststring&);
  void                        set_sequence(const faststring&);
  // void set_sequence(const char *, const char *);
  faststring&                     get_sequence();

};

#endif
