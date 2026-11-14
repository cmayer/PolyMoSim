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

#include "split_admin.h"



void split_admin::add_split(const BSplit& bs, double br_len) {

  split_vec.push_back(bs);

  split_map[bs] = br_len;

  num_splits++;

}



const BSplit& split_admin::get_split_vec(unsigned i) const 

{

  if(i >= split_vec.size()) 

  {

    faststring errormsg = "Error when accessing the vector of splits with index " + faststring(i);

    throw IndexError(errormsg);

  }

  else

  {

    return split_vec[i];

  }

}



double split_admin::get_branchlength(BSplit& bs) {

  map<BSplit, double>::iterator it;

  it = split_map.find(bs); 

  if(it != split_map.end())

    return split_map[bs];

  else{

    faststring errormsg = "Error: No such split in bitset map.";

    throw IndexError(errormsg);

  }

}



void split_admin::print_splits() {

  for(unsigned i = 0; i < split_vec.size(); i++) {

    split_vec[i].print();

  }

}

