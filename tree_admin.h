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

#ifndef Tree_adminH   // error checking still inconsisten xxxxxxxxxxxxxxxxxxxxxxxx
#define Tree_adminH

#include "PolyMoSim.h"

#include <vector>
#include "BasicTree.h"

struct tree_struct {
        double     scalingfactor;
        unsigned   partitionSize;
        faststring     default_modelName;
        faststring     treeString;
        BasicTree  *treeData;

  tree_struct(double pscalingfactor = 1, unsigned ppartitionSize = 0,
              faststring pdefault_modelName = "default",
              faststring ptreeString = ""):
    scalingfactor(pscalingfactor), partitionSize(ppartitionSize),
      default_modelName(pdefault_modelName), treeString(ptreeString), treeData(NULL) //, BasicTree* treeptr = NULL)
  {}
};

class tree_admin {

 public:
  class indexerror {
  };

  class readerror {
  private:
    int    line;
    faststring unknown;
  public:
  readerror(int l, faststring u):line(l),unknown(u) {}
    int    getLine() {return line; }
    faststring getUnknwonKeyword() {return unknown; }
  };

 private:
  vector<tree_struct*> treevector;

 public: 
  void               destroy(int);
  void               add(int, double, unsigned, const faststring&, const faststring &);
  void               swap(unsigned, unsigned);

  void               create(const char*);
  void               print(ostream& = cout, unsigned = 0, const char * = "");

  void               default_treevector(int);

  
  void               set(int, double, unsigned, const faststring &, const faststring &);
  void               set_BasicTree(unsigned);
  void               set_itsBasicTree(unsigned, BasicTree *);
  BasicTree*         get_itsBasicTree(unsigned);

  void               ini(int);

  int                getNumTrees() const;
  unsigned           getAlignmentLength() const;
  unsigned           getPartitionStartCoordinate(unsigned) const;
  faststring             get_dataTypeString() const;

  double             getScalingFactor(unsigned) const;
  unsigned           getPartitionSize(unsigned) const;
  const faststring&      getDefaultModelName(unsigned)  const;
  const faststring&      getTreeString(unsigned) const;
  
  bool               equal_takon_sets() const;
};

inline int tree_admin::getNumTrees() const {
  return (int) treevector.size();
}



#endif
