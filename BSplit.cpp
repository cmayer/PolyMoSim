/***************************************************************************************************
*  This file is part of PolyMoSim.
*  
*  Copyright (c) 2004-2026, Christoph Mayer and contributors (see README.md)
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

#include "BSplit.h"



BSplit::BSplit(const fast_dynamic_bitset& b) {

  split = b;

}



BSplit::BSplit(const BSplit& s) {

  split = s.split;

}



// Neue Reihenfolge der Argumente

BSplit::BSplit(unsigned long len, const vector<unsigned long>& vec):split(len) {

  vector<unsigned long>::const_iterator it;

  it = vec.begin();

  while( it != vec.end() ) {

    split.set(*it);       //indexerror

    ++it;

  }

}



// Neue Reihenfolge der Argumente

BSplit::BSplit(unsigned long len, unsigned long *u, unsigned long array_len):split(len) {

  split.clear();

  unsigned i = 0;

  while(i < array_len) {

    split.set(u[i]);      //indexerror

    ++i;

  }

}



bool BSplit::compatible_with(const BSplit& s_b) const {

  BSplit s_a_flipped(split);

  BSplit s_b_flipped(s_b.split);

  s_a_flipped.split.flip();

  s_b_flipped.split.flip();



  if( ((split & s_b.split).none())             || ((split & s_b_flipped.split).none())            || 

      ((s_a_flipped.split & s_b.split).none()) || ((s_a_flipped.split & s_b_flipped.split).none()) ) {

    return true;

  } else {

    return false;

  }  

}



void BSplit::print(ostream& os) const {

  os << split;

  os << endl;

}

