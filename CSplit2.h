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

#ifndef CSplitH
#define CSplitH

#include "fast-dynamic-bitset.h"
#include <iostream>
#include <vector>
#include <fstream>

class CSplit
{
 private:
  fast_dynamic_bitset split;

 public:
  // Bitset with a specfied number of bits, here taxa
  CSplit(unsigned long size):split(size) {
    split.clear();
  }

 CSplit(const fast_dynamic_bitset& b):split(b){}

 CSplit(const CSplit& b):split(b.split){}
    //  CSplit(unsigned long, const vector<unsigned long>&);
    //  CSplit(unsigned long, unsigned long*, unsigned long);

  friend bool operator ==(const CSplit&, const CSplit&);
  friend bool operator <(const CSplit&, const CSplit&);
  friend bool lessThan_mod_flip_friend(const CSplit&, const CSplit&);

  fast_dynamic_bitset& get_split()
  {
    return split;
  }

  void reset()
  {
    split.reset();
  }

  void set(unsigned i)
  {
    split.set(i);
  }

  bool test(unsigned i) const
  {
    return split.test(i);
  }

  void set(std::vector<int> vec)
  {
    int n=(int)vec.size();

    for (int i=0; i<n; ++i)
    {
      split.set(vec[i]);
    }
  }

  void flip()
  {
    split.flip();
  }

  unsigned count_taxa_in_ingroup()
  {
    return (unsigned) split.count();
  }

  unsigned size()
  {
    return (unsigned) split.size();
  }

  bool compatible_with(const CSplit& s_b) const
  {
    CSplit s_a_flipped(split);
    CSplit s_b_flipped(s_b.split);
    s_a_flipped.split.flip();
    s_b_flipped.split.flip();

    if( ((split & s_b.split).none())             || ((split & s_b_flipped.split).none())            || 
	((s_a_flipped.split & s_b.split).none()) || ((s_a_flipped.split & s_b_flipped.split).none()) )
      return true;
    else
      return false;
  }

  void print(std::ostream& os) const
  {
    os << split;
  }

  const std::string as_string_lowBitsFirst() const
  {
    return split.as_string_lowBitsFirst();
  }

  const std::string as_string_highBitsFirst() const
  {
    return split.as_string_highBitsFirst();
  }

};


/* // Neue Reihenfolge der Argumente */
/* CSplit::CSplit(unsigned long len, const vector<unsigned long>& vec):split(len) { */
/*   vector<unsigned long>::const_iterator it; */
/*   it = vec.begin(); */
/*   while( it != vec.end() ) { */
/*     split.set(*it);       //indexerror */
/*     ++it; */
/*   } */
/* } */

/* // Neue Reihenfolge der Argumente */
/* CSplit::CSplit(unsigned long len, unsigned long *u, unsigned long array_len):split(len) { */
/*   split.clear(); */
/*   int i = 0; */
/*   while(i < array_len) { */
/*     split.set(u[i]);      //indexerror */
/*     i++; */
/*     *u++; */
/*   } */
/* } */





inline bool operator ==(const CSplit& split_a, const CSplit& split_b) {
  if ( split_a.split == split_b.split ) {
    return true;
  } else {
    CSplit split_b_flipped(split_b);
    split_b_flipped.split.flip();
    if( split_a.split == split_b_flipped.split ) {
      return true;
    } else {
      return false;
    }
  }
}

inline bool operator <(const CSplit& split_a, const CSplit& split_b)
{
  if ( split_a.split < split_b.split ) {
    return true;
  }
  else {
    return false;
  }
}

/* inline bool lessThan_mod_flip_friend(const CSplit& split_a, const CSplit& split_b) */
/* { */
/*   CSplit s_a_flipped(split_a); */
/*   CSplit s_b_flipped(split_b); */

/*   s_a_flipped.split.flip(); */
/*   s_b_flipped.split.flip(); */

/*   //  return (split_a < split_b) || (s_a_flipped < s_b_flipped); */
/*   return (s_a_flipped < s_b_flipped); */
/* } */

/* struct lessThan_mod_flip */
/* { */
/*   bool operator()(const CSplit& split_a, const CSplit& split_b) const */
/*   { */
/*     return lessThan_mod_flip_friend(split_a, split_b); */
/*   } */
/* }; */

//
// Split with branch length:
//
class CSplit_l : public CSplit
{
  double b_length;

 public:

 CSplit_l(unsigned long u):CSplit(u), b_length(-1)
 {}

 CSplit_l(const fast_dynamic_bitset& b, double bl = -1):CSplit(b), b_length(bl) {}

 CSplit_l(const CSplit& b,  double bl = -1):CSplit(b), b_length(bl){}

 CSplit_l(const CSplit_l& b):CSplit(b), b_length(b.b_length){}

  

  void set_b_length(double b_l)
  {
    b_length = b_l;
  }

  double get_b_length() const
  {
    return b_length;
  }


};

#endif
