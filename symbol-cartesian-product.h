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

#include <iostream>
#include <string>
#include <vector>
#include "faststring3.h"

using namespace std;

typedef faststring mystring;


class cartesian_product
{
 private:
  

  void cartesian_product_for_2strings(const vector<mystring> &a, const vector<mystring> &b, vector<mystring> &res)
  {
    res.clear();
    for (unsigned i=0; i<a.size(); ++i)
    {
      for (unsigned j=0; j<b.size(); ++j)
      {
	mystring tmp = a[i] + b[j];
	res.push_back(tmp);
      }
    }
  }

  void nth_cartesian_product_of_symbolset(unsigned replicates, const mystring symbols, vector<mystring> &res)
  {
    res.clear();

    vector<mystring> start, tmp_res, tmp_res2;
    for (unsigned i=0; i < symbols.size(); ++i)
    {
      // string:
      //    mystring tmp = mystring(1u,symbols[i]);
      // faststring:
      mystring tmp = mystring((size_t)1, (char)symbols[i]);
      
      start.push_back(tmp);
    }
    
    //    cout << "DEBUG: start" << endl;
    //    print_container(cout, start, "", ",","");
    //    cout << "DEBUG: start -- END" << endl;
    
    tmp_res = start;
    for (unsigned j=1; j<replicates; ++j)
    {
      cartesian_product_for_2strings(tmp_res, start, tmp_res2);
      tmp_res = tmp_res2;
    }

    res = tmp_res2;
  }


 public:
  void operator()(unsigned replicates, const mystring symbols, vector<mystring> &res)
  {
    nth_cartesian_product_of_symbolset(replicates, symbols, res);
  }
};



