/***************************************************************************************************
*  The PolyMoSim project is distributed under the following license:
*  
*  Copyright (c) 2006-2022, Christoph Mayer, Forschungsmuseum Alexander Koenig, Bonn, Germany
*  All rights reserved.
*  
*  Redistribution and use in source and binary forms, with or without
*  modification, are permitted provided that the following conditions are met:
*  1. Redistributions of source code (complete or in parts) must retain
*     the above copyright notice, this list of conditions and the following disclaimer.
*  2. Redistributions in binary form must reproduce the above copyright
*     notice, this list of conditions and the following disclaimer in the
*     documentation and/or other materials provided with the distribution.
*  3. All advertising materials mentioning features or any use of this software
*     e.g. in publications must display the following acknowledgement:
*     This product includes software developed by Christoph Mayer, Forschungsmuseum
*     Alexander Koenig, Bonn, Germany.
*  4. Neither the name of the organization nor the
*     names of its contributors may be used to endorse or promote products
*     derived from this software without specific prior written permission.
*  
*  THIS SOFTWARE IS PROVIDED BY CHRISTOPH MAYER ''AS IS'' AND ANY
*  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
*  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
*  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHTHOLDER OR ITS ORGANISATION BE LIABLE FOR ANY
*  DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
*  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
*  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
*  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
*  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
*  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*  
*  IMPORTANT (needs to be included, if code is redistributed):
*  Please not that this license is not compatible with the GNU Public License (GPL)
*  due to paragraph 3 in the copyright. It is not allowed under any
*  circumstances to use the code of this software in projects distributed under the GPL.
*  Furthermore, it is not allowed to redistribute the code in projects which are
*  distributed under a license which is incompatible with one of the 4 paragraphs above.
*  
*  This project makes use of code coming from other projects. What follows is a complete
*  list of files which make use of external code. Please refer to the copyright within
*  these files.
*  
*  Files in tclap foler:         Copyright (c) 2003 Michael E. Smoot
*                                See copyright in tclap/COPYRIGHT file for details.	
*  discrete_gamma.c:             Copyright 1993-2004 by Ziheng Yang.
*                                See copyright in this file for details.
*  CRandom.h:                    Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura
*                                See copyright in this file for details.
***************************************************************************************************/

#ifndef CDISCRETIZEDDISTRIBUTION_H
#define CDISCRETIZEDDISTRIBUTION_H

#include <iostream>
#include "math_expression_parser.h"
#include <vector>
#include "faststring2.h"
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

  double get_coordinate(unsigned i)
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
