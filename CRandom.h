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

#ifndef CRANDOM_H
#define CRANDOM_H

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
//    Copyright (C) 2004-2008, Christoph Mayer
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



class CRandom
{
 protected:
  unsigned long  max_r;
  double         double_max_r;
  double         one_over_max_r;
  double         one_over_max_r_05;
  double         one_over_max_r_1;

  void           init_set_max_rand(unsigned long m);

  // The following constructor will be called by the subclass.
  // This is way it is protected.
                 CRandom(unsigned long m) { init_set_max_rand(m); }


 public:
  //  CRandom():max_r(4294967295) {}
  virtual  ~CRandom() {}

  virtual void           seed(unsigned long s) = 0;

  // random number in interval [0,0xffffffff] (32 bit)
  virtual unsigned long  random_ul()           = 0;
          unsigned long  max_rand() { return max_r; }
	  double         double_max_rand() { return double_max_r; }

  //**************************************
  // Special random number functions - in intervals / floating point format
  //**************************************

  // random number in interval [0,0x7fffffff] (31 bit) - non negative int
  long     random_l()     { return (long)(random_ul()>>1); }

  // random number in interval [0,1]
  double   random_lf_cc() { return ( random_ul())*one_over_max_r; }

  // random number in interval (0,1)
  double   random_lf_oo() { return ( random_ul()+0.5)*one_over_max_r_1; }

  // random number in interval [0,1)
  double   random_lf_co() { return ( random_ul())*one_over_max_r_05; }
	
  // random number in interval (0,1]
  double   random_lf_oc() { return ( random_ul()+0.5)*one_over_max_r_05; }

  //******************************
  // Random number distributions:
  //******************************
  double random_gamma(double a, double b);
  double random_gauss_polar_box_mueller(double sd);

 

  //******************************
  // internal functions:
  //******************************

 private:
  double random_gamma_large (double);
  double random_gamma_interval_01_Ahrens1  (double);
  double random_gamma_interval_01_Ahrens2  (double);
  double random_gamma_interval_01_Best  (double);
  double random_gamma_int    (unsigned);
  //  double random_gamma_old_01  (double);

};


//*****************************************************************************
// Mersenne Twister:  see copyright notice above.
// The class CRandom_MT19937 uses code from the C-program provided by
// Takuji Nishimura and Makoto Matsumoto.
//*****************************************************************************

class CRandom_MT19937 : public CRandom
{
 private:

  // Code taken from Takuji Nishimura and Makoto Matsumoto
  enum { N=624, M=397,
	 MATRIX_A=0x9908b0dfUL,    /* constant vector a         */
	 UPPER_MASK=0x80000000UL,  /* most significant w-r bits */
	 LOWER_MASK=0x7fffffffUL   /* least significant r bits  */ };
  
  unsigned long mt[N];             /* the array for the state vector          */
  unsigned      mti;               /* mti==N+1 means mt[N] is not initialized */

  void generate_N_words();


 public:
  virtual void seed(unsigned long); 
          void seed(unsigned long init_key[], unsigned key_length);

  virtual unsigned long  random_ul() 
  {
  // Code taken from Takuji Nishimura and Makoto Matsumoto
    unsigned long y;

    if (mti == N)  /* generate N words at one time */
      generate_N_words();

    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
  }

  CRandom_MT19937(unsigned long s = 5489UL):CRandom(4294967295UL)
  {
    seed(s);
  }

  CRandom_MT19937(unsigned long init_key[], int key_length):CRandom(4294967295UL)
  {
    seed(init_key, key_length);
  }
};



//*****************************************************************************
// Global functions to obtain random numbers from a single static Object.
//*****************************************************************************

//****************************
// Caution and warnng:
//****************************
// Take caution if you use these functions !!!!!!!!!!!!!!!!!!
// If included in different source files the same series of random numbers will
// be generated in each source file. This might not be intended.

static CRandom_MT19937 MT19937_rand_object;

inline void srandom_MT19937 (unsigned long s)
{
  MT19937_rand_object.seed(s);
}

inline unsigned long next_random_ul_MT19937()
{
  return MT19937_rand_object.random_ul();
}

inline long next_random_l_MT19937()
{
  return MT19937_rand_object.random_l();
}

inline double next_random_lf_cc_MT19937()
{
  return MT19937_rand_object.random_lf_cc();
}

inline double next_random_lf_oo_MT19937()
{
  return MT19937_rand_object.random_lf_oo();
}

inline double next_random_lf_co_MT19937()
{
  return MT19937_rand_object.random_lf_co();
}

inline double next_random_lf_oc_MT19937()
{
  return MT19937_rand_object.random_lf_oc();
}

inline double next_random_gamma_MT19937(const double a, const double b)
{
  return MT19937_rand_object.random_gamma(a, b);
}

inline double next_random_gauss_polar_box_muelle_MT19937(const double sd)
{
  return MT19937_rand_object.random_gauss_polar_box_mueller(sd);
}



#endif
