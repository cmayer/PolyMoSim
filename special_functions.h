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

#ifndef SPECIAL_FUNCTIONS_H
#define SPECIAL_FUNCTIONS_H

#include <cmath>

//static const double EPS      = 0.00000001;
//static const double ONE_EPS  = 1-EPS;

inline double distribution_uniform(double x)
{
  if (x >= 0 && x < 1)
    return 1;
  else
    return 0;
}

inline double distribution_uniform_interval(double x, double a, double b)
{
  if (x > a && x < b)  // x==a is not allowed. It automatically prevents division by 0 if a==b.
    return (1.0/(b-a));
  else
    return 0.0;
}


inline  double distribution_gauss(double x, double mu, double sd)
{
  // 1/sqrt(2 M_PI) replaced by 0.39894228
  double tmp = (x-mu)/sd;
  return 0.39894228 /sd *exp(-tmp*tmp/2);
}


// From numerical recipes in C, 1992
// returns ln(gamma(x))
inline double ln_gamma(double x)
{
 int j;
 double y,tmp,ser;
 static const double cof[6]={76.18009172947146,-86.50532032941677,24.01409824083091,-1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5};
 
 y=x;
 tmp=x+5.5;
 tmp -= (x+0.5)*log(tmp);
 ser=1.000000000190015;
 for (j=0;j<6;j++) ser += cof[j]/++y;
 return -tmp+log(2.5066282746310005*ser/x);
}

inline double local_gamma(double x)
{
  return exp(ln_gamma(x));
}


inline  double distribution_gamma(double x, double a)
{
  // We use that g(x) = pow(x, a-1)*exp(-x)/gamma(a) 
  //                  = exp((a-1)*ln(x)-x -ln(gamma(a)))
  return exp((a-1)*log(x) - x - ln_gamma(a));
}

// provides:
// beta^alpha/gamma(alpha)*x(alpha-1)exp(-beta*x)
// should have mean m=alpha/beta and variance v=m/beta
// Also: alpha=m*m/V and beta=m/V
// Other implementations sometimes use: theta=1/beta instead of beta.
// This is the form usually found when beta is used as parameter,
// but the literature is also not consistent in this respect.

inline  double distribution_gamma(double x, double a, double b)
{
  // We use that g(x) = pow(x, a-1)*exp(-x)/gamma(a) 
  //                  = exp((a-1)*ln(x)-x -ln(gamma(a)))
  return distribution_gamma(x*b,a)*b;
}

// Beta - distribution:
// pdf = Gamma(a+b)/Gamma(a)/Gamma(b)(1-x)^(b-1)x^(a-1)
//     = exp(lnGamma(a+b)-lnGamma(a)-lnGamma(b)+(b-1)ln(1-x)+(a-1)ln(x))   for x!=0, x!=1
// Defined for a>0, b>0, 0 <= x <= 1
// For some a,b it is not defined for x==0 or x==1
// Sometimes the convention pdf(x) = 0 is used for x<0 or x>1.
// Since for x==0 and x==1 pdf is either 0 or undefined, we set it to 0
// Here we will use pdf(x)=0 if (x<=0 and x>=1)
inline double distribution_beta(double x, double a, double b)
{
  if (x <= 0  || x >= 1 )
  {
    return 0;
  }
  if (a <= 0 || b <=0)
  {
    return 0;
  }
  return exp( ln_gamma(a+b)-ln_gamma(a)-ln_gamma(b) +(b-1)*log(1-x)+(a-1)*log(x) );
}

inline double heaviside(double x)
{
  if (x<0)
    return 0;
  else
    return 1;
}

inline double distribution_cauchy(double x, double s, double t)
{
  if (s<=0)
    return 0;

  return 1/3.141592653589793238*s/(s*s+(x-t)*(x-t)); 
}



#endif
