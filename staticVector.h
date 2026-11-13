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

#ifndef STATICVECTORH

#define STATICVECTORH



#include <iostream>



template<int N>

class staticVector {



  template<int> friend class staticSquareMatrix;



 public:

  class indexerror {

  };



 private:

  double bvector[N];



 public:

  staticVector() {};

  staticVector(double);

  staticVector(double *);



  void     assign(double);

  void     assign(double *);



  void     print(std::ostream& = std::cout, double = 4);

  double&  operator [](int);

  const double&  operator [](int) const;

};





template<int N>

staticVector<N>::staticVector(double x) {

  int j;

  for (j = 0; j < N; ++j) {

    bvector[j] = x;

  }

}



template<int N>

staticVector<N>::staticVector(double *x) {

  std::memcpy(bvector, x, sizeof(bvector));

}



template<int N>

inline void staticVector<N>::assign(double x) {

  int j;

  for (j = 0; j < N; ++j) {

    bvector[j] = x;

  }

}



template<int N>

inline void staticVector<N>::assign(double *x) {

  std::memcpy(bvector, x, sizeof(bvector));

}



template<int N>

inline void staticVector<N>::print(std::ostream& os, double pres) {

  int j;



  std::cout.setf(std::ios::fixed);

  std::cout.precision(pres);

  for (j = 0; j < N-1; ++j) {

    os << bvector[j] << " ";

  }

  os << bvector[j];

}



template<int N>

inline double& staticVector<N>::operator[](int i) {

  return bvector[i];

}



template<int N>

inline const double& staticVector<N>::operator[](int i) const {

  return bvector[i];

}



#endif

