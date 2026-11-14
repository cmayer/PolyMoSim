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

#ifndef BSplitH
#define BSplitH

#include "fast-dynamic-bitset.h"
#include "faststring3.h"
#include <iostream>
#include <vector>

using namespace std;

class BSplit {

 private:
  fast_dynamic_bitset split;

 public:
  BSplit(unsigned long u = 0):split(u) {
    split.clear();
  }

  BSplit(const fast_dynamic_bitset&);
  BSplit(const BSplit&);
  BSplit(unsigned long, const vector<unsigned long>&);
  BSplit(unsigned long, unsigned long*, unsigned long);

  friend bool operator ==(const BSplit&, const BSplit&);
  friend bool operator <(const BSplit&, const BSplit&);

  bool compatible_with(const BSplit&) const;
  fast_dynamic_bitset& get_split();

  void print(ostream& = cout) const;

};

inline bool operator ==(const BSplit& split_a, const BSplit& split_b) {
  if( split_a.split == split_b.split ) {
    return true;
  } else {
    BSplit split_b_flipped(split_b);
    split_b_flipped.split.flip();
    if( split_a.split == split_b_flipped.split ) {
      return true;
    } else {
      return false;
    }
  }
}

inline bool operator <(const BSplit& split_a, const BSplit& split_b) {
  if( split_a.split < split_b.split ) {
    return true;
  }
  else {
    return false;
  }
}

inline fast_dynamic_bitset& BSplit::get_split() {
  return split;
}

#endif
