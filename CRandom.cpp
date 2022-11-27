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

//*****************************************************************************
//
// Copyright and references:
//
// This class hierarchy provides an easy access to different random number
// generators.
// Random number generators that are currently implemented: (more will be added if I need them)
// - Mersenne twister MT19937 by Takuji Nishimura and Makoto Matsumoto for 32 bit random numbers.
//
// Furthermore, this class hierachy provides access to random numbers from different
// distributions.
// Distributions currently implemented: (more will be added if I need them)
// - gamma distribution
// - normal distribution
//
// Copyright for C++ class hierarchy:
//    Copyright (C) 2004-2011, Christoph Mayer
//    All rights reserved. 
//
//  As far as external code or published algorithms have been used,
//    license terms and references are given here: 

// ************************************************************************
// Mersenne twister: MT19937
//
// This C++ code for the MT19937 is based on the c-code provided by the authors of the
// Mersenne Twister, Takuji Nishimura and Makoto Matsumoto, as given on their web page:
// http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/CODES/mt19937ar.c
// Functions and source code taken from this source are indicated below.
//
//-------------------------------------------------------------
// License and copyright of the origianl Mersenne Twister:
//-------------------------------------------------------------
// A C-program for MT19937, with initialization improved 2002/1/26.
// Coded by Takuji Nishimura and Makoto Matsumoto.
//
// Copyright for random number generator:
//    Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
//    All rights reserved.                          
//
//    Redistribution and use in source and binary forms, with or without
//    modification, are permitted provided that the following conditions
//    are met:
//
//      1. Redistributions of source code must retain the above copyright
//         notice, this list of conditions and the following disclaimer.
//
//      2. Redistributions in binary form must reproduce the above copyright
//         notice, this list of conditions and the following disclaimer in the
//         documentation and/or other materials provided with the distribution.
//
//      3. The names of its contributors may not be used to endorse or promote 
//         products derived from this software without specific prior written 
//         permission.
//
//    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
//    A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
//    CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
//    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
//    PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
//    PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
//    LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
//    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//   Any feedback is very welcome.
//   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
//   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove spaces)
//   For C++ adaption: Christoph.Mayer @ rub.de (remove spaces)
//
//*****************************************************************************

//*****************************************************************************
// Random numbers from distributions:
//
// Random number drawn from gamma distribution - References:
// Chapter 3.4.1 in "The art of computer programming" by D. Knuth, volume 2, 3rd edition
// Ahrens,J.H.,Dieter,U. (1974), Computer methods for sampling from gamma, beta, Poisson and binomial distributions. Computing 12, 223-246.
// Best, D.J. (1983), A note on gamma variate generators with shape parameter less than unity, Computing 30, 185-188. 
//
// Random number drawn from Gauss or normal distribution - References:
// A convenient method for generating normal variables, G. Marsaglia and T. A. Bray, SIAM Rev. 6, 260-264, 1964




#include "CRandom.h"
#include <cmath>

using namespace std;

//***************************
// CRandom functions:
//***************************

void CRandom::init_set_max_rand(unsigned long m)
{
  max_r         = m;
  double_max_r  = m;
  one_over_max_r    = 1.0/double_max_r;
  one_over_max_r_05 = 1.0/(double_max_r+0.5);
  one_over_max_r_1  = 1.0/(double_max_r+1.0);
}


//*****************************************************************************
// Mersenne Twister:  see copyright below
//*****************************************************************************

/* The source code for this function is taken from the c-code provided by */
/* Takuji Nishimura and Makoto Matsumoto in their MT19937 program.        */
/* See copyright and license at top of this file.                         */

/* initializes mt[N] with a seed */
void CRandom_MT19937::seed(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; ++mti) {
        mt[mti] = 
	    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
} // exits with mti == N


/* The source code for this function is taken from the c-code provided by */
/* Takuji Nishimura and Makoto Matsumoto in their MT19937 program.        */
/* See copyright and license at top of this file.                         */

/***********************************************/
/* Comment in original code:                   */
/* initialize by an array with array-length    */
/* init_key is the array for initializing keys */
/* key_length is its length                    */
/* slight change for C++, 2004/2/26            */
/***********************************************/

/*** Original function header: void init_by_array(unsigned long init_key[], int key_length) ***/
/* Has been changed since I intend to add more random number generators all of which should */
/* shoud be subclasses of CRandom and all should have the same interface.                   */
void CRandom_MT19937::seed(unsigned long init_key[], unsigned key_length)
{
    unsigned i, j, k;
    seed(19650218UL);
    i=1; j=0;
    k = (N>key_length ? N : key_length);
    for (; k; --k) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        ++i; ++j;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=N-1; k; --k) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        ++i;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
    }

    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 
} // exits with mti == N


/* The source code for this function is taken from the c-code provided by */
/* Takuji Nishimura and Makoto Matsumoto in their MT19937 program.        */
/* See copyright and license at top of this file.                         */

/* The code here was taken from the function with the header:                     */ 
/* unsigned long genrand_int32(void)                                              */
/* The rest of this function is implemented in the method with the name:          */
/* virtual unsigned long random_ul()                                              */
/* Only from there we call this function.                                         */

void CRandom_MT19937::generate_N_words()
{
  unsigned      kk;
  unsigned long y;

  static unsigned long mag01[2]={0x0UL, MATRIX_A};
  /* mag01[x] = x * MATRIX_A  for x=0,1 */

  for (kk=0;kk<N-M;++kk) {
    y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
    mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
  }
  for (;kk<N-1;++kk) {
    y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
    mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
  }
  y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
  mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];
  
  mti = 0;
}


//******************************
// Random number distributions:
//******************************

//************************************************
// Ranom numbers from gamma distribution:
//************************************************
// See chapter 3.4.1 in "The art of computer programming" by D. Knuth,
// volume 2, 3rd edition, page 134.
// A very good presentation of all parts of this algorithm is given in
// Ahrens,J.H.,Dieter,U.:Computer methods for sampling from gamma, beta, Poisson and binomial distributions. Computing 12, 223-246 (1974).

//**************************************************************
// Random number from gamma distribution for 0 < a < infinity
//**************************************************************
// This function is described in the two references given above.
double CRandom::random_gamma(double a, double b) 
{
  unsigned long ui_a = (unsigned long) floor(a);

  if (a == ui_a) // For integer values:
    return b * random_gamma_int(ui_a);
  else if (ui_a == 0) // For the range a>0 and a<1
    //    return b * random_gamma_interval_01_Ahrens1(a);
    return b * random_gamma_interval_01_Ahrens2(a);
    // return b * random_gamma_interval_01_Best(a);
    //    return b * random_gamma_old_01(a);
  else // For a >= 1 we make use of: If X and Y are gamma distrubuted of order a and b, X+Y is gamma distributed of order a+b, see Knuth, section 3.4.1
    return b * (random_gamma_int(ui_a) + random_gamma_interval_01_Ahrens2(a - ui_a));
}

//**************************************************************
// Random number from gamma distribution for a being an integer:
//**************************************************************
// If a is an integer >= 1 we again use:
// If X and Y are gamma distributed of order a and b, X+Y is gamma distributed of order a+b, see Knuth, section 3.4.1
// Furthermore, for a==1, X can be computed using X = -\mu ln (U) where U is uniformly distributed, see Knuth, section 3.4.1
// So we simply add a times Xi, where Xi are gamma distributed with order 1.
// We do this by computing X <- -ln(U1*...*Ua), where U1 to Ua are uniformly distributed with Ui>0.
// This algorithm is described in both references given above.
// The potential pitfal of obtaining an underflow is mentioned in the book by Knuth.
double CRandom::random_gamma_int (unsigned a)
{
  // For large a there is the risk that we get an underflow for the product U1*..*Ua.
  // But ln(0) is not defined, so we need to make sure this can not happen.

  // We use the following thoughts:
  // The smallest random number we can get is assumed to be 2^-32. ( = 1/max_rand)
  // The smallest double is assumed to be 10^-307.
  // Then 2^(-32*15) is approximately equal to 3.2*10^-145, which is well in the double precision range.

  if (a < 15)
  {
    unsigned i;
    double   res = 1;
    
    for (i = 0; i < a; ++i)
    {
      res *= random_lf_oo();
    }

    return -log (res);
  }
  else
  {
    return random_gamma_large((double) a);
  }
}


double CRandom::random_gamma_large (double a)
{
  /* Requies a > 1, and is better if a is large.

     Has been taken from Knuth, chapter 3.4.1. (Algorithm A).
     There, it is referenced to:
     J.H. Ahrens, Ann. Inst. Stat. Math. 13 (1962), pp 231-237.

     Other interesting references for gamma distributed random numbers for large a
     are given in Knuth's book.

     It is also mentioned that a faster method has been published by:
     J.H.Ahrens and U Dieter, CACM 25 (1982), pp 47-54.

     The algorithm used in this function has also been mentioned in
     Ahrens,J.H.,Dieter,U.: Computer methods for sampling
     from gamma, beta, Poisson and binomial distributions. Computing 12, 223-246 (1974).

     Several modificatins exist in the literature.
  */

  double X, Y, V;
  static double sqrt_2a_1 = sqrt (2.0 * a - 1.0);
  for (;;)
  {
    for (;;)
    {
      Y = tan (M_PI * random_lf_co());
      X = sqrt_2a_1 * Y + a - 1.0;
      if (X > 0) break; // Accept??
    }

    V = random_lf_co();

    if (V <= ((1.0+Y*Y)*exp((a-1.0)*log (X/(a-1.0))-sqrt_2a_1*Y)))
      break;
  }
  return X;
}

// The algorithm given by in J.H.Ahrens and U Dieter, CACM 25 (1982), pp 47-54.
// should be faster. It could be implemented in a future version of this random class.
// Here we state the algorithm as given in the refernce:
// The text has been copied from the OCR-Pdf of the original paper:

// 0. Preset a' = 0 and a" = 0 (at compilation time).
// 1. If a != a' set a'	= a, s2 = a-1/2, s=sqrt(s2), and d=4*sqrt(2)- 12s = 5.656 854 249 492 38 - 12s. 
// 2. Generate T (standard normal deviate). Set X=s+T/2. If T>=0 return X*X. 
// 3. Generate U [(0, 1)-uniform deviate]. If d*U <= T*T*T return X*X. 
// 4. If a != a" set a" = a and calculate q0, b, sigma, and c as follows. 
//    q0 = SUM q_k a^(-k) (instead of ln(2*PI) - InGamma(a)- s2 + s2 ln s2);
//    For 1<=a<=3.686:              b     = 0.463 + s + 0.178*s2
//                                  sigma = 1.235
//                                  c     = 0.195/s - 0.079 + 0.16s
//    For 3.686 < a <= 13.022:      b     = 1.654   + 0.0076 s2
//                                  sigma = 1.68/s  + 0.275
//                                  c     = 0.062/s + 0.024
//    For 13.022 < a < infinity:    b     = 1.77
//                                  sigma = 0.75
//                                  c     = 0.1515/s
// 5. If X <= 0 goto Step 8 
// 6. Set V = T/(s + s) 
//    If abs(V) >  1/4: Q = q0 - sT + T*T/4 + (s2 + s2)ln(1+V).
//    If abs(V) <= 1/4: Q = q0 + (T*T/2) SUM a_k V^k. 
// 7. If ln(1-U)<=Q return X*X. 
// 8. Generate E (standard exponential deviate) and U [(0, 1)-uniform deviate].
//    Set U=U+U-1 and T=b+E*sigma*sign(U). 
// 9. If T <= tau_1 = -0.718 744 837 717 19 go to Step 8. 
// 10. Set V = T/(s + s) and calculate Q as in Step 6.
// 11. If Q<=0 or if c*abs(U)>(exp(Q)-1)*exp(E-T*T/2) goto Step 8.
//     (If Q <= 1/2 the factor exp Q -	1 is calculated as SUM e_k Q^k.)
// 12. Set X = s + T/2 and return X*X.



// The following code has been adapted from the book: Knuth, The art of computer programming vol 2, 3rd edition,
// Excersice 16, section 3.4.1, page 140 and solution on page 586.
// This book also mentiones some modifications which I did not implement since I found it
// difficult to understand why they are equivalent to the original.
// This algorithm is attributed to Ahrens.
double CRandom::random_gamma_interval_01_Ahrens1 (double a)
{
  double U, V, X;
  double p, q;
  p = M_E/(a+M_E); // Probability to use G1

  for (;;)
  {
    // Generate uniform deviates U and V
    U = random_lf_co();
    V = random_lf_oo(); // two sided open interval => V > 0
    
    if (U < p)
    {
      X = exp ((1.0 / a) * log (V)); // This is equal to V^(1/a), which is mentioned in Knuth
      q = exp (-X);
    }
    else
    {
      X = 1 - log (V);
      q = exp ((a - 1.0) * log (X));
    }
    if (random_lf_co() < q) break;
  }
  return X;
}


// A modified version of the algorithm above.
// It has some similarity to the above mentioned modification described in the book of Knuth
// which I found hard to understand.
// The algorithm has been adapted from
// Ahrens,J.H.,Dieter,U. (1974) Computer methods for sampling from gamma, beta, Poisson and binomial distributions. Computing 12, 223-246.
// There it is referred to as "Algorithm GS".

// Algorithm GS (0 < a <= 1): Pseudo code: (The variable w below was called b in the original paper.
// Since b is also used for the second parameter of the gamma distribution, we chose to call it w here.)
// 1: Generate U. Set w=(e+a)/e and P=wU. If P>1 goto 3 elso goto 2
// 2: (Case X<=1) Set X = P^(1/a). Generate U*. If U* > e^(-X) go back to 1 (rejection). Otherwise deliver X.
// 3: (Case X>1)  Set X = -ln ((w-P)/a). Generate U*. If U* > X^(a-1) go back to 1 (rejection). Otherwise deliver X.
//
// We make use of the simple formula: x^y = exp(y*ln(x)) in the hope the these are faster on our computer than computing the power.

double CRandom::random_gamma_interval_01_Ahrens2 (double a) // Double checked: YES (13.7.2012)
{
  double w, P, X, one_over_a = 1.0/a;

  w = (a + M_E)/M_E;

  for (;;)
  {
    // Generate uniform deviate U and compute P immediately. U won't be used any more to we have dropped the varialble.
    P = w*random_lf_co();

    if (P > 1) // 3:
    {
      X = -log ((w-P)*one_over_a);
      if (random_lf_oo() <= exp((a-1.0)*log(X)) )
	return X;
    }
    else // 2:
    {
      X = exp (one_over_a * log (P));

      if (random_lf_oo() <= exp(-X))
	return X;
    }
  }
}


// The following algorithm has been adapted from
// Best, D.J. (1983), A note on gamma variate generators with shape parameter less than unity, Computing 30, 185-188.
// The algorithm is referred to as Algorithm RGS (0 < a < 1): 
//
// Pseudo code:
// (The variable r below was called b in the original paper. Since b is also used for the second parameter of the gamma
// distribution, we chose to use r here.)
//
// 0: Initialize z = 0.07 + 0.75(1 -a)^(1/2), r = 1 + a/z*e^(-z).
// 1: Generate U and set P=rU. If P>1 go to 4 else goto 2
// 2: Set X=zP^(1/a). Generate U*. If U*<=(2-X)/(2+X), deliver X.
// 3: If U* > e^(-x) go to 1, otherwise deliver X.
// 4: Set X = -ln(z(r-P)/a), Y=X/z. Generate U*. If U*(a + Y - aY) < 1, deliver X.
// 5: If U*>Y^(a-1) go to 1, otherwise deliver X.

double CRandom::random_gamma_interval_01_Best (double a)
{
  double z, r;
  double P, X, Y, U_star;
  double one_over_a = 1.0/a;

  z = 0.07 + 0.75*sqrt(1.0-a); // In principle only has to be calculated once for every a
  r = 1.0 + a/z* exp(-z);         // In principle only has to be calculated once for every a

  for (;;)
  {
    // We compute P immediately and drop the variable U.
    P = r*random_lf_co();
    if (P>1) // (IV)
    {
      X = -log(z*(r-P)*one_over_a);
      Y = X/z;
      U_star = random_lf_oo();
      if (U_star*(a+Y-a*Y) < 1.0)
	return X;
      if (U_star <= exp((a-1.0)*log(Y)) )
	return X;
    }
    else // (II)
    {
      X = z*exp((one_over_a)*log(P));
      U_star = random_lf_oo();
      if (U_star <= (2.0-X)/(2.0+X))
	return X;
      if (U_star <= exp(-X) )
	return X;
    }
  }
}

//****************************************
// Gauss or normal distribution:
//****************************************
// Reimplementation of the well known Box-Mueller transform
// This polar form is attributed to Marsaglia and Bray:
// Reference: A convenient method for generating normal variables, G. Marsaglia and T. A. Bray, SIAM Rev. 6, 260-264, 1964
//
double CRandom::random_gauss_polar_box_mueller(double sd)
{
  double x1, x2, w, y1;
  static double y2;
  static bool   use_last = false;

  if (use_last) 
  {
    use_last = false;
    return y2 * sd;
  }
  else 
  {
    do 
    {
      x1 = 2.0 * random_lf_oo() - 1.0;
      x2 = 2.0 * random_lf_oo() - 1.0;
      w = x1 * x1 + x2 * x2;
    } while ( w >= 1.0 || w == 0 );

    w = sqrt((-2.0 * log(w))/w);
    y1 = x1 * w;
    y2 = x2 * w;
    use_last = true;
    return y1 * sd;
  }
}


//*****************************************************************************
// Global stuff
//*****************************************************************************

//*****************************************************************************
// Global object from which random numbers are drawn.
//*****************************************************************************


CRandom_MT19937 __intern_rand_object;



//*****************************************************************************
