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

#ifndef CDISCRETIZEDDISTRIBUTION_H
#define CDISCRETIZEDDISTRIBUTION_H

#include <iostream>
#include "math_expression_parser.h"
#include <vector>
#include "faststring3.h"
#include <algorithm>
#include "CRandom.h"


class DiscretizedDistribution
{
  double   a;
  double   b;
  unsigned N;
  double   integral;
  double   step;
  
  double *cdf; // Data of the cumulative distribution.

  //  std::vector<double> cdf_vec;

  double *debug1, *debug2;
  
  DiscretizedDistribution();

 public:
 DiscretizedDistribution(double param_a, double param_b, unsigned param_N,
			 Cexpression_evaluator  *exp): a(param_a), b(param_b), N(param_N)
  {
    //    cdf_vec.reserve(N);
    cdf = new double [N];

    unsigned i;
    double   x=a;
    step= (b-a)/N;
    double val;
    double sum = 0;

    // Start in middle of for first interval.

    for (i=0, x=a+step/2; i<N; ++i, x+=step )
    {
      val = exp->compute_for_indeterminate_x(x);
      sum += val;
      cdf[i] = sum;
      // cdf_vec.push_back(sum); // adds element i
    }

    integral = sum;
  }

  ~DiscretizedDistribution()
  {
    //    std::cout << "Entering DiscretizedDistribution destructor" << std::endl;
    delete [] cdf;
    //    std::cout << "Exiting DiscretizedDistribution destructor" << std::endl;
  }


  DiscretizedDistribution(const DiscretizedDistribution &);
  DiscretizedDistribution& operator=(const DiscretizedDistribution &);

  double get_integral()
  {
    return integral;
  }

  double get_cdf(unsigned i)
  {
    if (i<N)
      //      return cdf_vec[i];
      return cdf[i];
    else
    {
      std::cerr << "Internal error in function: CDiscrtizedDistribution.h:get_cdf. Parameter is out of range.\n";
      std::exit(0);
    }
  }

  double get_pdf(unsigned i)
  {
    if (i==0)
      return cdf[0];
    //      return cdf_vec[0];

    if (i<N)
      //      return cdf_vec[i]-cdf_vec[i-1];
    return cdf[i]-cdf[i-1];
    else
    {
      std::cerr << "Internal error in function: CDiscrtizedDistribution.h:get_cdf. Parameter is out of range.\n";
      std::exit(0);
    }
  }

  double get_coordinate(unsigned long i)
  {
    return a+step*(i+0.5);
  }

  double get_coordinate_of_cdf_proportion(double x)
  {
    double prop = x*integral;

    double *lb = std::lower_bound(cdf, cdf+N, prop);
    //    std::vector<double>::iterator it_lb = std::lower_bound(cdf_vec.begin(), cdf_vec.end(), prop);
    

/*     if ( it_lb != cdf_vec.end() ) */
/*       return get_coordinate(std::distance(cdf_vec.begin(), it_lb)); */
/*     else */
/*       return b+1; */


    if ( lb != cdf+N )
      return get_coordinate(lb-cdf);
    else
      return b+1;

  }

  void print_cdf(std::ostream &os)
  {
    unsigned i;

    for (i=0; i<N; ++i)
    {
      os << get_coordinate(i) << " " << get_cdf(i) << std::endl;
    }
  }
};



class Random_DiscretizedDistribution : public DiscretizedDistribution
{
  double (*random_lf_co)();

 public:
 Random_DiscretizedDistribution(double param_a, double param_b, unsigned param_N,
				Cexpression_evaluator  *param_exp, double (*random_lf_co_param)() ):
  DiscretizedDistribution(param_a, param_b, param_N, param_exp), random_lf_co(random_lf_co_param)
  {}

  double next_random_value()
  {
    return get_coordinate_of_cdf_proportion(random_lf_co());
  }

};


#endif
