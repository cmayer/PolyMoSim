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

#ifndef Split_adminH
#define Split_adminH

#include <map>
#include <vector>
#include <iostream>
#include "fast-dynamic-bitset.h"
#include "BSplit.h"
#include "faststring3.h"

using namespace std;

class split_admin {

 public:
  class IndexError{
  private:
    faststring errormsg;
  public:
    IndexError(const faststring& s):errormsg(s) {}
    const faststring& getReason() { return errormsg; }
  };

 private:

  vector<BSplit>          split_vec;
  map<BSplit, double>     split_map;
  unsigned long num_splits;
 
 public:

  split_admin() {
  }

  void            add_split(const BSplit&, double = 1.0);

  const BSplit&   get_split_vec(unsigned) const;
  double          get_branchlength(BSplit&);

  unsigned long   get_num_splits();
  
  void            print_splits();
 
};

inline unsigned long split_admin::get_num_splits() {
  return num_splits;
}

#endif
