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

#ifndef MymodelH
#define MymodelH

// TODO: Check what needs to be private, protected, public!!!!!!!!

#include "faststring2.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <map>
// #include "mymatrix.h"
// #include "myvector.h"
#include "PolyMoSim.h"
#include "staticSquareMatrix.h"
#include "staticVector.h"

class basic_model { /**/ // These could be moved to this basic class
  public:
  enum enumDataType {DNA, Protein};
  protected:
  const static char datatypenames[][8];

  enumDataType dataType;

 basic_model(enumDataType d):dataType(d){};

 public:
  enum enumratetype {
    ratetype_equal,
    ratetype_gamma,
    ratetype_invgamma,
    ratetype_propinv,
    ratetype_cat_rates,
    ratetype_inv_cat_rates,
    ratetype_distfunction,
    ratetype_distfunction_inv
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
  
  class setmodelerror {
  private:
    faststring reason;
  public:  
    setmodelerror(faststring r=""):reason(r) {}
    faststring getReason() {return reason; }
  };

  class error          {};

  faststring get_datatypename() const { return datatypenames[dataType]; } 
  int    get_datatype()     const { return dataType; }



  virtual ~basic_model(){};
  virtual faststring get_modelname() const=0;        /**/ // Should be moved into this class in the future 
  //  virtual void   set_modelname(const faststring &s); /**/ // Should be moved into this class in the future 
  virtual void   read_next_model(CFile&)=0;
  virtual void   set_model(faststring modelname_param, 
			   faststring modeltype_param,
			   std::vector<double> *rrates_param, // Parameters are supplied as in an upper triangular matrix.
			   std::vector<double> *base_param,
			   double         shape_param,
			   double         pinv_param,
			   unsigned       ncat_param,
			   double         *tstv_param)=0;

  virtual void   print(std::ostream&, unsigned)=0;
  virtual void   print(FILE *, unsigned)=0;

  virtual void   init_model(double(double, double), double())=0; /**/ // Should be moved into this class in the future 
  virtual void   set_matrices()=0;
  virtual void   init_siterates(unsigned, bool)=0; /**/ // Should be moved into this class in the future 
  virtual void   reset_siterates()=0; /**/ // Should be moved into this class in the future 

  virtual bool   siterates_initialized() const =0; /**/ // Should be moved into this class in the future 

  virtual void   evolve(const faststring &, faststring &, double) const=0;
  virtual void   get_random_sequence(faststring&, unsigned) const =0; 
  virtual double get_mean_siterate_value() const =0;
  virtual const faststring& get_modelname_siterates_are_inherited_from() const = 0;

  virtual basic_model*&   get_and_set_model_to_inherit_siterates_from() =0;

  virtual bool                set_inherited_site_rates(unsigned)=0;
  virtual double*             get_siterates()=0;
  virtual double              get_shape() const=0;
  virtual unsigned            get_ncat() const=0;
  virtual double              get_inv() const=0;
  virtual enumratetype        get_ratetype() const=0;
  virtual std::vector<double> get_cat_rates() const=0;

  virtual void   print_relative_site_rates(std::ostream &) const =0;
  virtual void   print_site_rates_histogramm_data(std::ostream &, faststring &)const =0;
};

template <int N>   // nucleotide: N=4, aa: N=20
class molecular_model : public basic_model {
 public:

  // Several of these variables should be moved to the basic_model:

 protected:
  unsigned               seq_length;
  faststring             modelname;

  double                 shape;           // gamma_alpha value
  faststring             distribution_math_expression;

  double                 distribution_a;
  double                 distribution_b;
  unsigned               distribution_N;

  enumratetype           ratetype;

  unsigned               ncat;            // Number of rate categories. ncat == 0 is continues.
  double                 inv;             // invariant positions
  int                    modeltype;  

  double                 *siterates;
  faststring             modelname_inherit_siterates_from;
  basic_model*           model_to_inherit_siterates_from;

  std::vector<double>    cat_rates;

  char                   lower_or_upper_matrix_in_input; // allowed values: 'l', and 'u'

  staticSquareMatrix<N>  relRates;
  staticSquareMatrix<N>  T, T_inv;
  staticVector<N>        eigenvalues;
  staticVector<N>        pi;
  //  staticVector<N>        summed_pi;


  // Random number generators needed by the model:
  double (*random_gamma)(double, double);
  double (*random_lf_co)();



 protected:

  // Constructors:
  molecular_model(const faststring &, enumDataType);
  molecular_model(const faststring &, const molecular_model&);

  virtual                      ~molecular_model() { if (siterates) { delete [] siterates; siterates = NULL;}  }

  virtual void                 print(std::ostream&, unsigned)=0;
  virtual void                 print(FILE *, unsigned)=0;

  virtual void                 evolve(const faststring &, faststring &, double) const;
  virtual void                 read_next_model(CFile& );  //--------
  virtual void                 set_model(faststring modelname_param, 
					 faststring modeltype_param,
					 std::vector<double> *rrates_param, // Parameters are supplied as in an upper triangular matrix.
					 std::vector<double> *base_param,
					 double         shape_param,
					 double         pinv_param,
					 unsigned       ncat_param,
					 double         *tstv_param);  //--------

          void                 complete_relRateMatrix();
	  void                 normalize_rrates();
	  bool                 normalize_basefreq(bool warn);
  //  void                     calculate_transition_matrix(double, int, mymatrix &) const;

  virtual const unsigned char* get_index_to_symbol() const=0;
  virtual const unsigned char* get_symbol_to_index() const=0;
  virtual void                 set_model_specific_parameters(CFile *,
							     faststring& input_modeltypename,
							     double  input_tstv,
							     bool    specified_tstv,
							     bool    specified_rrates, bool specified_base_frequencies)=0;

 public:
          faststring         get_modelname() const;
	  ////          void           set_modelname(const faststring &s); // Dangeros with model_map

  virtual int            get_modeltype() const=0;
  virtual faststring         get_modeltypename() const=0;
  virtual void           get_random_sequence(faststring&, unsigned) const;
//  const    get_modeltypenames() const=0;

  void                   init_model(double(double, double), double());
  void                   init_siterates(unsigned, bool);
  void                   reset_siterates();
  bool                   siterates_initialized() const;
  
  double                 get_mean_siterate_value() const;
  const faststring&      get_modelname_siterates_are_inherited_from() const;
  basic_model*&          get_and_set_model_to_inherit_siterates_from() { return model_to_inherit_siterates_from; }

  // Returns true if siterates are indeed initialized now.
  bool                   set_inherited_site_rates(unsigned len)
  {
    if (len == 0)
    {
      std::cerr << "Error: Sequence length is specified to be 0, which does not make sense in function set_inherited_site_rates. For model " << get_modelname() << std::endl;
      exit(0);;
    }

    if (model_to_inherit_siterates_from)
    {
      seq_length = len;
      siterates  = model_to_inherit_siterates_from->get_siterates();
      shape      = model_to_inherit_siterates_from->get_shape();
      ratetype   = model_to_inherit_siterates_from->get_ratetype();

      ncat       = model_to_inherit_siterates_from->get_ncat();
      inv        = model_to_inherit_siterates_from->get_inv();
      cat_rates  = model_to_inherit_siterates_from->get_cat_rates();
      //							if (siterates)
      //							  siterates_initialized = true;
    }
    else
    {
      std::cerr << "Error: Call to function set_inherited_site_rates(unsigned len) was not succcessful." << std::endl;
    }
    return siterates != NULL;
                                                    } 
  double*                get_siterates() { return siterates; }

  double                 compute_true_propinv() const;

  //sets matrices for calculation of transition matrix
  void                   set_matrices();
  double                 get_shape() const;
  unsigned               get_ncat() const;
  double                 get_inv() const;
//  enumratetype           get_ratetype() const;
  enumratetype           get_ratetype() const {return ratetype;}

  std::vector<double>    get_cat_rates() const;
  

  // Debug routines:
  void                   debug_print_model_symbols() const;
  void                   print_relative_site_rates(std::ostream &os) const;
  void                   print_site_rates_histogramm_data(std::ostream &os, faststring &) const;
}; // End class molecular_model : public basic_model


//*************************************************
// nuc_model
//*************************************************
class nuc_model : public molecular_model<4> {
 public:
  enum         enummodeltype  {  JC,   F81,   K2P,   F84,   HKY,   GTR };
 private:
  enum         sym_enum       {  nA, nC, nG, nT};

  const static int              number_of_known_models = 6;
  const static char             modeltypenames[][6];
  const static unsigned char    index_to_symbol[4];
  const static unsigned char    symbol_to_index[256];

  double                        tstv;

  void                          set_rates_nst12();
  double                        compute_tstv() const;

  virtual const unsigned char*  get_index_to_symbol() const { return index_to_symbol; }
  virtual const unsigned char*  get_symbol_to_index() const { return symbol_to_index; }

  virtual void         set_model_specific_parameters(CFile *,
						     faststring& input_modeltypename,
						     double  input_tstv,
						     bool    specified_tstv,
						     bool    specified_rrates, bool specified_base_frequencies);

 public:
 nuc_model(const faststring &s):molecular_model<4>(s, DNA){tstv = -1;}
 nuc_model(const faststring &s, const nuc_model &m):molecular_model<4>(s, m){tstv = m.tstv;}
 nuc_model(const faststring &s, enummodeltype model_type, enumratetype rate_type, double v_tstv, double v_rAC, double v_rAG, double v_rAT, double v_rCG, double v_rCT, double v_rGT, double v_shape, double v_inv, unsigned v_ncat, double v_PI_A, double v_PI_C, double v_PI_G, double v_PI_T):molecular_model<4>(s, DNA)
  {
    
   if(rate_type == ratetype_equal        || rate_type == ratetype_gamma ||
      rate_type == ratetype_invgamma     || rate_type == ratetype_propinv ||
      rate_type == ratetype_distfunction || rate_type == ratetype_distfunction_inv)
   {
     ratetype = rate_type;
   }
   else
   {
     throw setmodelerror();
   }

   if (model_type == JC   || model_type == F81 ||
       model_type == K2P  || model_type == F84 ||
       model_type == HKY  || model_type == GTR) {
     modeltype = model_type;
   }
   else 
   {
    throw setmodelerror();
   }

   shape = v_shape;
   inv = v_inv;
   ncat = v_ncat;   //wird nicht verändert...

  pi[0] = v_PI_A; 
  pi[1] = v_PI_C; 
  pi[2] = v_PI_G; 
  pi[3] = v_PI_T; 

  tstv = v_tstv;

  relRates(0,1) = v_rAC;
  relRates(0,2) = v_rAG;
  relRates(0,3) = v_rAT;
  relRates(1,2) = v_rCG;
  relRates(1,3) = v_rCT;
  relRates(2,3) = v_rGT;

  complete_relRateMatrix();
  normalize_rrates();       // wirklich hier oder nur bei GTR?????????????????
  normalize_basefreq(false);


  if ( modeltype == JC  || modeltype == F81 ||
       modeltype == K2P || modeltype == F84 ||
       modeltype == HKY) {
    set_rates_nst12();
  }
  else if ( modeltype == GTR ) {
    tstv = compute_tstv();
  }

  //notwendig an dieser Stelle????????????????
  double PI_all = pi[0] + pi[1] + pi[2] + pi[3];

  if(PI_all < 1-EPS || PI_all > 1+EPS)
    throw setmodelerror();

  }


  virtual int           get_modeltype() const;
  virtual faststring        get_modeltypename() const;
  virtual void          print(std::ostream&, unsigned);
  virtual void          print(FILE *, unsigned);

  double    get_rAC() const;
  double    get_rAG() const;
  double    get_rAT() const;
  double    get_rCG() const;
  double    get_rCT() const;
  double    get_rGT() const;
  double    get_tstv() const;
  double    get_PI_A() const;
  double    get_PI_G() const;
  double    get_PI_T() const;
  double    get_PI_C() const;


};


//*************************************************
// aa_model
//*************************************************
class aa_model : public molecular_model<20> {
 private:
  enum         enummodeltype  { USER, JTT, LG, WAG_OLD, WAG, WAG_STAR, DAY};
  enum         sym_enum   { aaA, aaR, aaN, aaD, aaC, aaQ, aaE, aaG, aaH, aaI,
			    aaL, aaK, aaM, aaF, aaP, aaS, aaT, aaW, aaY, aaV };

  const static int              number_of_known_models = 3;
  const static char             modeltypenames[][9];
  const static unsigned char    index_to_symbol[20];
  const static unsigned char    symbol_to_index[256];

  virtual const unsigned char*  get_index_to_symbol() const { return index_to_symbol; }
  virtual const unsigned char*  get_symbol_to_index() const { return symbol_to_index; }

  virtual void                  set_model_specific_parameters(CFile *,
							      faststring& input_modeltypename,
							      double  input_tstv,
							      bool    specified_tstv,
							      bool    specified_rrates, bool specified_base_frequencies);

 public:
  aa_model(const faststring &s):molecular_model<20>(s, Protein){};
  aa_model(const faststring &s, const aa_model &m):molecular_model<20>(s, m){};

  virtual int           get_modeltype() const;
  virtual faststring        get_modeltypename() const;
  virtual void          print(std::ostream&, unsigned);
  virtual void          print(FILE *, unsigned);



};







/* class mymodel */
/* { */



/*   void      set_model(const faststring&, enummodeltype, enumratetype, double, double, double, double, double, double, double, */
/* 		      double, double, double, double, double, double); */




/* }; */


//************************************************************************

#endif

// ToDo
// - Konstruktor
// - set_model
// - read model



// order of funciton calls: normalize rates / complete_relRateMatrix
