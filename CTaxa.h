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

#ifndef CTAXA_H
#define CTAXA_H

#include <map>
#include <vector>
#include "faststring3.h"

using namespace std;

class CTaxa
{
private:
  map<faststring, unsigned>   taxaNames_map;
  vector<faststring>          taxaNames_vec;

public:
  unsigned       GetTaxaNum();
  const faststring&  GetTaxonLabel(unsigned) const;
  unsigned       GetTaxonPosition(const faststring&);
  void           add_taxon(const faststring&);

  unsigned       GetMaxTaxonLabelLength() const;
  //  void           GetNamesOfTaxaInSplit(const mybitset<BITSETSIZE>&,
  //					     faststring, faststring&, faststring& ) const;

};

inline unsigned CTaxa::GetTaxaNum() {
  return (unsigned)taxaNames_vec.size();
}

inline const faststring&  CTaxa::GetTaxonLabel(unsigned i) const {
  return taxaNames_vec[i]; //index abfragen
}

inline unsigned CTaxa::GetTaxonPosition(const faststring& s) {
  return taxaNames_map[s];  // ist faststring (bzw. taxonname) vorhanden?
}

inline void CTaxa::add_taxon(const faststring& s) {
  taxaNames_map[s] = (unsigned)taxaNames_vec.size();
  taxaNames_vec.push_back(s);
}


inline unsigned CTaxa::GetMaxTaxonLabelLength() const
{
  unsigned maxlen = 0, ti, len;

  for (ti = 0; ti < ( taxaNames_vec.size() - 1 ) ; ++ti)
    if ( (len = (unsigned)taxaNames_vec[ti].length())>maxlen )
      maxlen = len;

  return maxlen;
}

#endif
