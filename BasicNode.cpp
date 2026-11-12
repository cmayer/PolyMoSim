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

#include "BasicNode.h"



using namespace std;



BasicNode::BasicNode(){

  children = 0;

  name = "";

  branchlength = 0.0;

  model_name = "";

  new_model = false;

  sequence = "";

  model = NULL;

  my_parent = NULL;

}



void BasicNode::push_back_child(BasicNode *child){

  childlist.push_back(child);

  children++;

}



BasicNode* BasicNode::get_parent() {

  return my_parent;

}



void BasicNode::set_parent(BasicNode *parent){

  my_parent = parent;

}



unsigned BasicNode::get_num_children() const {

  return children;

}



const std::vector<BasicNode*>& BasicNode::get_childlist(){

  return childlist;

}



const faststring& BasicNode::get_name() const {

  return name;

}



void BasicNode::set_name(const faststring& str){

  name = str;

}



void BasicNode::push_back_to_name(char c){

  name.push_back(c);

}



double BasicNode::get_branchlength(){

  return branchlength;

}



void BasicNode::set_branchlength(double b){

  branchlength = b;

}



const basic_model* BasicNode::get_model() {

  return model;

}



void BasicNode::set_model(const basic_model* newmodel) {

  model = newmodel;

  set_model_name(model->get_modelname());

}



bool BasicNode::get_new_model_status() {

  return new_model;

}



void BasicNode::set_new_model_status(bool b) {

  new_model = b;

}



faststring& BasicNode::get_model_name() {

  return model_name;

}



void BasicNode::set_model_name(const faststring& str) {

  model_name = str;

}



void BasicNode::set_sequence(const faststring& str) {

  sequence = str;

}



// void BasicNode::append_sequence(const char *b, const char *e) {

//   sequence.append(b, e);

// }



faststring& BasicNode::get_sequence() {

  return sequence;

}

