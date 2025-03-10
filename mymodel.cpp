/***************************************************************************************************
*  The PolyMoSim project is distributed under the following license:
*  
*  Copyright (c) 2006-2025, Christoph Mayer, Leibniz Institute for the Analysis of Biodiversity Change,
*  Bonn, Germany
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

#include "mymodel.h"
#include "CFile/CFile2_3.h"
#include "discrete_gamma.hpp"
#include <cassert>
#include "CHistogram.h"
#include "CDiscretizedDistribution.h"
#include "math_expression_parser.h"

using namespace std;

#define REPORT_PROBABILITY_DIFFERENCE

//*********************************************
// static variables basic_model
//********************************************
const          char basic_model::datatypenames[][8] = { "DNA", "Protein"};

//*********************************************
// static variables nuc_model
//*********************************************
const          char nuc_model::modeltypenames[][6] = { "JC", "F81", "K2P", "F84", "HKY", "GTR"};
const unsigned char nuc_model::index_to_symbol[4]          = {'A', 'C', 'G', 'T'};
const unsigned char nuc_model::symbol_to_index[256]       = {
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255,   0, 255,   1, 255, 255,
  255,   2, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255,   3, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255  };

//*********************************************
// static variables aa_model
//*********************************************
// A	Ala	Alanine         00   65
// R	Arg	Arginine        01   82
// N	Asn	Asparagine      02   78
// D	Asp	Aspartic acid   03   68
// C	Cys	Cysteine        04   67
// Q	Gln	Glutamine       05   81
// E	Glu	Glutamic acid   06   69
// G	Gly	Glycine         07   71
// H	His	Histidine       08   72
// I	Ile	Isoleucine      09   73
// L	Leu	Leucine         10   76
// K	Lys	Lysine          11   75 
// M	Met	Methionine      12   77
// F	Phe	Phenylalanine   13   70
// P	Pro	Proline         14   80
// S	Ser	Serine          15   83
// T	Thr	Threonine       16   84
// W	Trp	Tryptophan      17   87
// Y	Tyr	Tyrosine        18   89
// V	Val	Valine          19   86

const          char aa_model::modeltypenames[][9] = { "USER", "JTT", "LG", "WAG_OLD", "WAG", "WAG_STAR", "DAY"};
const unsigned char aa_model::index_to_symbol[20] =
{'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
  'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'};
const unsigned char aa_model::symbol_to_index[256] =
{ 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255,   0, 255,   4,   3,   6,
  13,   7,   8,   9, 255,  11,  10,  12,   2, 255,
  14,   5,   1,  15,  16, 255,  19,  17, 255,  18,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255  };


//---------------------------------------------------
// JTT model for amino acid evolution
// D.T. Jones, W.R. Taylor, and J.M. Thornton
// The rapid generation of mutation data matrices from protein sequences
// CABIOS  vol. 8 no. 3 1992 pp. 275-282
//---------------------------------------------------
/*static double jttRelativeRates_1[] = {
 0.531678, 0.557967, 0.827445, 0.574478, 0.556725, 1.066681, 1.740159, 0.219970, 0.361684, 0.310007, 0.369437, 0.469395, 0.138293, 1.959599, 3.887095, 4.582565, 0.084329, 0.139492, 2.924161,
 0.451095, 0.154899, 1.019843, 3.021995, 0.318483, 1.359652, 3.210671, 0.239195, 0.372261, 6.529255, 0.431045, 0.065314, 0.710489, 1.001551, 0.650282, 1.257961, 0.235601, 0.171995,
 5.549530, 0.313311, 0.768834, 0.578115, 0.773313, 4.025778, 0.491003, 0.137289, 2.529517, 0.330720, 0.073481, 0.121804, 5.057964, 2.351311, 0.027700, 0.700693, 0.164525,
 0.105625, 0.521646, 7.766557, 1.272434, 1.032342, 0.115968, 0.061486, 0.282466, 0.190001, 0.032522, 0.127164, 0.589268, 0.425159, 0.057466, 0.453952, 0.315261,
 0.091304, 0.053907, 0.546389, 0.724998, 0.150559, 0.164593, 0.049009, 0.409202, 0.678335, 0.123653, 2.155331, 0.469823, 1.104181, 2.114852, 0.621323,
 3.417706, 0.231294, 5.684080, 0.078270, 0.709004, 2.966732, 0.456901, 0.045683, 1.608126, 0.548807, 0.523825, 0.172206, 0.254745, 0.179771,
 1.115632, 0.243768, 0.111773, 0.097485, 1.731684, 0.175084, 0.043829, 0.191994, 0.312449, 0.331584, 0.114381, 0.063452, 0.465271,
 0.201696, 0.053769, 0.069492, 0.269840, 0.130379, 0.050212, 0.208081, 1.874296, 0.316862, 0.544180, 0.052500, 0.470140,
 0.181788, 0.540571, 0.525096, 0.329660, 0.453428, 1.141961, 0.743458, 0.477355, 0.128193, 5.848400, 0.121827,
 2.335139, 0.202562, 4.831666, 0.777090, 0.098580, 0.405119, 2.553806, 0.134510, 0.303445, 9.533943,
 0.146481, 3.856906, 2.500294, 1.060504, 0.592511, 0.272514, 0.530324, 0.241094, 1.761439,
 0.624581, 0.024521, 0.216345, 0.474478, 0.965641, 0.089134, 0.087904, 0.124066,
 0.436181, 0.164215, 0.285564, 2.114728, 0.201334, 0.189870, 3.038533,
 0.148483, 0.943971, 0.138904, 0.537922, 5.484236, 0.593478,
 2.788406, 1.176961, 0.069965, 0.113850, 0.211561,
 4.777647, 0.310927, 0.628608, 0.408532,
 0.080556, 0.201094, 1.143980,
 0.747889, 0.239697,
 0.165473
 };


 static double jttFrequencies_1[] = {
 0.076862, 0.051057, 0.042546, 0.051269, 0.020279, 0.041061, 0.061820, 0.074714, 0.022983, 0.052569, 0.091111, 0.059498, 0.023414, 0.040530, 0.050532, 0.068225, 0.058518, 0.014336, 0.032303, 0.066374
 };
 */

static double jttRelativeRates[] = {
  0.531678,
  0.557967, 0.451095,
  0.827445, 0.154899, 5.549530,
  0.574478, 1.019843, 0.313311, 0.105625,
  0.556725, 3.021995, 0.768834, 0.521646, 0.091304,
  1.066681, 0.318483, 0.578115, 7.766557, 0.053907, 3.417706,
  1.740159, 1.359652, 0.773313, 1.272434, 0.546389, 0.231294, 1.115632,
  0.219970, 3.210671, 4.025778, 1.032342, 0.724998, 5.684080, 0.243768, 0.201696,
  0.361684, 0.239195, 0.491003, 0.115968, 0.150559, 0.078270, 0.111773, 0.053769, 0.181788,
  0.310007, 0.372261, 0.137289, 0.061486, 0.164593, 0.709004, 0.097485, 0.069492, 0.540571, 2.335139,
  0.369437, 6.529255, 2.529517, 0.282466, 0.049009, 2.966732, 1.731684, 0.269840, 0.525096, 0.202562, 0.146481,
  0.469395, 0.431045, 0.330720, 0.190001, 0.409202, 0.456901, 0.175084, 0.130379, 0.329660, 4.831666, 3.856906, 0.624581,
  0.138293, 0.065314, 0.073481, 0.032522, 0.678335, 0.045683, 0.043829, 0.050212, 0.453428, 0.777090, 2.500294, 0.024521, 0.436181,
  1.959599, 0.710489, 0.121804, 0.127164, 0.123653, 1.608126, 0.191994, 0.208081, 1.141961, 0.098580, 1.060504, 0.216345, 0.164215, 0.148483,
  3.887095, 1.001551, 5.057964, 0.589268, 2.155331, 0.548807, 0.312449, 1.874296, 0.743458, 0.405119, 0.592511, 0.474478, 0.285564, 0.943971, 2.788406,
  4.582565, 0.650282, 2.351311, 0.425159, 0.469823, 0.523825, 0.331584, 0.316862, 0.477355, 2.553806, 0.272514, 0.965641, 2.114728, 0.138904, 1.176961, 4.777647,
  0.084329, 1.257961, 0.027700, 0.057466, 1.104181, 0.172206, 0.114381, 0.544180, 0.128193, 0.134510, 0.530324, 0.089134, 0.201334, 0.537922, 0.069965, 0.310927, 0.080556,
  0.139492, 0.235601, 0.700693, 0.453952, 2.114852, 0.254745, 0.063452, 0.052500, 5.848400, 0.303445, 0.241094, 0.087904, 0.189870, 5.484236, 0.113850, 0.628608, 0.201094, 0.747889,
  2.924161, 0.171995, 0.164525, 0.315261, 0.621323, 0.179771, 0.465271, 0.470140, 0.121827, 9.533943, 1.761439, 0.124066, 3.038533, 0.593478, 0.211561, 0.408532, 1.143980, 0.239697, 0.165473,
};

static double jttFrequencies[] = {
  0.076862, 0.051057, 0.042546, 0.051269, 0.020279, 0.041061, 0.061820, 0.074714, 0.022983, 0.052569, 0.091111, 0.059498, 0.023414, 0.040530, 0.050532, 0.068225, 0.058518, 0.014336, 0.032303, 0.066374
};


//---------------------------------------------------
// WAG model of amino acid evolution
// Whelan, S. and Goldman, N. (2001) A general empirical model of protein
// evolution derived from multiple protein families using a maximum-likelihood
// approach. Mol. Biol. Evol. 18, 691-699.
//---------------------------------------------------
// What is the difference between wag old and new????
static double wag_oldRelativeRates[] = {
  0.610810, 0.569079, 0.821500, 1.141050, 1.011980, 1.756410, 1.572160, 0.354813, 0.219023, 0.443935, 1.005440, 0.989475, 0.233492, 1.594890, 3.733380, 2.349220, 0.125227, 0.268987, 2.221870,
  0.711690, 0.165074, 0.585809, 3.360330, 0.488649, 0.650469, 2.362040, 0.206722, 0.551450, 5.925170, 0.758446, 0.116821, 0.753467, 1.357640, 0.613776, 1.294610, 0.423612, 0.280336,
  6.013660, 0.296524, 1.716740, 1.056790, 1.253910, 4.378930, 0.615636, 0.147156, 3.334390, 0.224747, 0.110793, 0.217538, 4.394450, 2.257930, 0.078463, 1.208560, 0.221176,
  0.033379, 0.691268, 6.833400, 0.961142, 1.032910, 0.043523, 0.093930, 0.533362, 0.116813, 0.052004, 0.472601, 1.192810, 0.417372, 0.146348, 0.363243, 0.169417,
  0.109261, 0.023920, 0.341086, 0.275403, 0.189890, 0.428414, 0.083649, 0.437393, 0.441300, 0.122303, 1.560590, 0.570186, 0.795736, 0.604634, 1.114570,
  6.048790, 0.366510, 4.749460, 0.131046, 0.964886, 4.308310, 1.705070, 0.110744, 1.036370, 1.141210, 0.954144, 0.243615, 0.252457, 0.333890,
  0.630832, 0.635025, 0.141320, 0.172579, 2.867580, 0.353912, 0.092310, 0.755791, 0.782467, 0.914814, 0.172682, 0.217549, 0.655045,
  0.276379, 0.034151, 0.068651, 0.415992, 0.194220, 0.055288, 0.273149, 1.486700, 0.251477, 0.374321, 0.114187, 0.209108,
  0.152215, 0.555096, 0.992083, 0.450867, 0.756080, 0.771387, 0.822459, 0.525511, 0.289998, 4.290350, 0.131869,
  3.517820, 0.360574, 4.714220, 1.177640, 0.111502, 0.353443, 1.615050, 0.234326, 0.468951, 8.659740,
  0.287583, 5.375250, 2.348200, 0.462018, 0.382421, 0.364222, 0.740259, 0.443205, 1.997370,
  1.032220, 0.098843, 0.619503, 1.073780, 1.537920, 0.152232, 0.147411, 0.342012,
  1.320870, 0.194864, 0.556353, 1.681970, 0.570369, 0.473810, 2.282020,
  0.179896, 0.606814, 0.191467, 1.699780, 7.154480, 0.725096,
  1.786490, 0.885349, 0.156619, 0.239607, 0.351250,
  4.847130, 0.578784, 0.872519, 0.258861,
  0.126678, 0.325490, 1.547670,
  2.763540, 0.409817,
  0.347826
};

static double wag_oldFrequencies[20] = {
  0.0866, 0.0440, 0.0391, 0.0570, 0.0193, 0.0367, 0.0581, 0.0833, 0.0244, 0.0485, 0.0862, 0.0620, 0.0195, 0.0384, 0.0458, 0.0695, 0.0610, 0.0144, 0.0353, 0.0709
};


static double wagRelativeRates[] = {
  0.551571,
  0.509848, 0.635346,
  0.738998, 0.147304, 5.429420,
  1.027040, 0.528191, 0.265256, 0.0302949,
  0.908598, 3.035500, 1.543640, 0.616783, 0.0988179,
  1.582850, 0.439157, 0.947198, 6.174160, 0.021352, 5.469470,
  1.416720, 0.584665, 1.125560, 0.865584, 0.306674, 0.330052, 0.567717,
  0.316954, 2.137150, 3.956290, 0.930676, 0.248972, 4.294110, 0.570025, 0.249410,
  0.193335, 0.186979, 0.554236, 0.039437, 0.170135, 0.113917, 0.127395, 0.0304501, 0.138190,
  0.397915, 0.497671, 0.131528, 0.0848047, 0.384287, 0.869489, 0.154263, 0.0613037, 0.499462, 3.170970,
  0.906265, 5.351420, 3.012010, 0.479855, 0.0740339, 3.894900, 2.584430, 0.373558, 0.890432, 0.323832, 0.257555,
  0.893496, 0.683162, 0.198221, 0.103754, 0.390482, 1.545260, 0.315124, 0.174100, 0.404141, 4.257460, 4.854020, 0.934276,
  0.210494, 0.102711, 0.0961621, 0.0467304, 0.398020, 0.0999208, 0.0811339, 0.049931, 0.679371, 1.059470, 2.115170, 0.088836, 1.190630,
  1.438550, 0.679489, 0.195081, 0.423984, 0.109404, 0.933372, 0.682355, 0.243570, 0.696198, 0.0999288, 0.415844, 0.556896, 0.171329, 0.161444,
  3.370790, 1.224190, 3.974230, 1.071760, 1.407660, 1.028870, 0.704939, 1.341820, 0.740169, 0.319440, 0.344739, 0.967130, 0.493905, 0.545931, 1.613280,
  2.121110, 0.554413, 2.030060, 0.374866, 0.512984, 0.857928, 0.822765, 0.225833, 0.473307, 1.458160, 0.326622, 1.386980, 1.516120, 0.171903, 0.795384, 4.378020,
  0.113133, 1.163920, 0.0719167, 0.129767, 0.717070, 0.215737, 0.156557, 0.336983, 0.262569, 0.212483, 0.665309, 0.137505, 0.515706, 1.529640, 0.139405, 0.523742, 0.110864,
  0.240735, 0.381533, 1.086000, 0.325711, 0.543833, 0.227710, 0.196303, 0.103604, 3.873440, 0.420170, 0.398618, 0.133264, 0.428437, 6.454280, 0.216046, 0.786993, 0.291148, 2.485390,
  2.006010, 0.251849, 0.196246, 0.152335, 1.002140, 0.301281, 0.588731, 0.187247, 0.118358, 7.821300, 1.800340, 0.305434, 2.058450, 0.649892, 0.314887, 0.232739, 1.388230, 0.365369, 0.314730
};

static double wagFrequencies[20] = {
  0.0866279, 0.043972, 0.0390894, 0.0570451, 0.0193078, 0.0367281, 0.0580589, 0.0832518, 0.0244313, 0.048466, 0.086209, 0.0620286, 0.0195027, 0.0384319, 0.0457631, 0.0695179, 0.0610127, 0.0143859, 0.0352742, 0.0708956
};


static double wagstarRelativeRates[] = {
  0.589718,
  0.514347, 0.67416,
  0.731152, 0.159054, 5.30821,
  1.21324,  0.568449, 0.233527,  0.0379056,
  1.03344,  3.02808,  1.62299,   0.657364,  0.0999068,
  1.55788,  0.443685, 1.00122,   6.04299,   0.0284956, 5.6037,
  1.41993,  0.629768, 1.12717,   0.88357,   0.312544,  0.346823, 0.588609,
  0.317684, 2.31211,  3.9337,    0.958529,  0.341479,  4.87366,  0.599188,  0.279542,
  0.214596, 0.187262, 0.527321,  0.0390513, 0.198958,  0.125999, 0.124553,  0.0310522, 0.162975,
  0.400822, 0.51821,  0.144354,  0.0869637, 0.451124,  0.873266, 0.154936,  0.067443,  0.508952, 3.1554,
  0.881639, 5.74119,  2.88102,   0.480308,  0.0719929, 4.19125,  2.45392,   0.381514,  0.854485, 0.320597, 0.255092,
  0.887458, 0.660816, 0.198404,  0.0992829, 0.428648,  1.64018,  0.294481,  0.184545,  0.40117,  3.94646,  4.81956,  0.877057,
  0.213179, 0.122792, 0.0848492, 0.0458258, 0.485001,  0.109241, 0.0873936, 0.0552962, 0.631713, 1.06458,  2.10414,  0.0832422, 1.14516,
  1.51861,  0.711498, 0.204905,  0.444152,  0.109081,  0.913179, 0.720567,  0.254626,  0.722123, 0.111722, 0.422851, 0.588203,  0.179858, 0.165205,
  3.52499,  1.35611,  3.90127,   1.09965,   1.35221,   0.87908,  0.822025,  1.33618,   0.876688, 0.321774, 0.351913, 1.05314,   0.554077, 0.563999, 1.54694,
  2.24161,  0.594177, 2.06787,   0.395176,  0.522957,  0.829315, 0.889765,  0.236489,  0.54992,  1.48876,  0.351564, 1.45173,   1.56873,  0.188237, 0.802531, 4.02507,
  0.135395, 1.24086,  0.0746093, 0.142159,  0.728065,  0.208163, 0.176397,  0.366467,  0.261223, 0.259584, 0.706082, 0.159261,  0.565299, 1.58681,  0.135024, 0.528249, 0.118584,
  0.270321, 0.386714, 1.05269,   0.326191,  0.481954,  0.210494, 0.209621,  0.108982,  4.31772,  0.44009,  0.427718, 0.155623,  0.437069, 6.49269,  0.212945, 0.742154, 0.286443, 2.42261,
  1.92496,  0.282892, 0.193323,  0.155419,  1.10899,   0.32893,  0.588443,  0.190095,  0.119749, 7.48376,  1.82105,  0.300343,  2.03324,  0.653015, 0.325745, 0.23769,  1.4088,   0.396884, 0.353358
};

static double wagstarFrequencies[20] = {
  0.0866279, 0.043972, 0.0390894, 0.0570451, 0.0193078, 0.0367281, 0.0580589, 0.0832518, 0.0244313, 0.048466, 0.086209, 0.0620286, 0.0195027, 0.0384319, 0.0457631, 0.0695179, 0.0610127, 0.0143859, 0.0352742, 0.0708956
};

//---------------------------------------------------
// LG model of amino acid evolution:
// Le, S. Q., and O. Gascuel. 2008. An improved general amino acid replacement matrix. Mol. Biol. Evol. 25:1307-1320.
//---------------------------------------------------
/*static double lgRelativeRates_1[] = {
 0.425093, 0.276818, 0.395144, 2.489084, 0.969894, 1.038545, 2.066040, 0.358858, 0.149830, 0.395337, 0.536518, 1.124035, 0.253701, 1.177651, 4.727182, 2.139501, 0.180717, 0.218959, 2.547870,
 0.751878, 0.123954, 0.534551, 2.807908, 0.363970, 0.390192, 2.426601, 0.126991, 0.301848, 6.326067, 0.484133, 0.052722, 0.332533, 0.858151, 0.578987, 0.593607, 0.314440, 0.170887,
 5.076149, 0.528768, 1.695752, 0.541712, 1.437645, 4.509238, 0.191503, 0.068427, 2.145078, 0.371004, 0.089525, 0.161787, 4.008358, 2.000679, 0.045376, 0.612025, 0.083688,
 0.062556, 0.523386, 5.243870, 0.844926, 0.927114, 0.010690, 0.015076, 0.282959, 0.025548, 0.017416, 0.394456, 1.240275, 0.425860, 0.029890, 0.135107, 0.037967,
 0.084808, 0.003499, 0.569265, 0.640543, 0.320627, 0.594007, 0.013266, 0.893680, 1.105251, 0.075382, 2.784478, 1.143480, 0.670128, 1.165532, 1.959291,
 4.128591, 0.267959, 4.813505, 0.072854, 0.582457, 3.234294, 1.672569, 0.035855, 0.624294, 1.223828, 1.080136, 0.236199, 0.257336, 0.210332,
 0.348847, 0.423881, 0.044265, 0.069673, 1.807177, 0.173735, 0.018811, 0.419409, 0.611973, 0.604545, 0.077852, 0.120037, 0.245034,
 0.311484, 0.008705, 0.044261, 0.296636, 0.139538, 0.089586, 0.196961, 1.739990, 0.129836, 0.268491, 0.054679, 0.076701,
 0.108882, 0.366317, 0.697264, 0.442472, 0.682139, 0.508851, 0.990012, 0.584262, 0.597054, 5.306834, 0.119013,
 4.145067, 0.159069, 4.273607, 1.112727, 0.078281, 0.064105, 1.033739, 0.111660, 0.232523, 10.649107,
 0.137500, 6.312358, 2.592692, 0.249060, 0.182287, 0.302936, 0.619632, 0.299648, 1.702745,
 0.656604, 0.023918, 0.390322, 0.748683, 1.136863, 0.049906, 0.131932, 0.185202,
 1.798853, 0.099849, 0.346960, 2.020366, 0.696175, 0.481306, 1.898718,
 0.094464, 0.361819, 0.165001, 2.457121, 7.803902, 0.654683,
 1.338132, 0.571468, 0.095131, 0.089613, 0.296501,
 6.472279, 0.248862, 0.400547, 0.098369,
 0.140825, 0.245841, 2.188158,
 3.151815, 0.189510,
 0.249313
 };

 static double lgFrequencies_1[20] = {
 0.079066, 0.055941, 0.041977, 0.053052, 0.012937, 0.040767, 0.071586, 0.057337, 0.022355, 0.062157, 0.099081, 0.064600, 0.022951, 0.042302, 0.044040, 0.061197, 0.053287, 0.012066, 0.034155, 0.069147
 };
 */

static double lgRelativeRates[] = {
  0.425093,
  0.276818, 0.751878,
  0.395144, 0.123954, 5.076149,
  2.489084, 0.534551, 0.528768, 0.062556,
  0.969894, 2.807908, 1.695752, 0.523386, 0.084808,
  1.038545, 0.363970, 0.541712, 5.243870, 0.003499, 4.128591,
  2.066040, 0.390192, 1.437645, 0.844926, 0.569265, 0.267959, 0.348847,
  0.358858, 2.426601, 4.509238, 0.927114, 0.640543, 4.813505, 0.423881, 0.311484,
  0.149830, 0.126991, 0.191503, 0.010690, 0.320627, 0.072854, 0.044265, 0.008705, 0.108882,
  0.395337, 0.301848, 0.068427, 0.015076, 0.594007, 0.582457, 0.069673, 0.044261, 0.366317, 4.145067,
  0.536518, 6.326067, 2.145078, 0.282959, 0.013266, 3.234294, 1.807177, 0.296636, 0.697264, 0.159069, 0.137500,
  1.124035, 0.484133, 0.371004, 0.025548, 0.893680, 1.672569, 0.173735, 0.139538, 0.442472, 4.273607, 6.312358, 0.656604,
  0.253701, 0.052722, 0.089525, 0.017416, 1.105251, 0.035855, 0.018811, 0.089586, 0.682139, 1.112727, 2.592692, 0.023918, 1.798853,
  1.177651, 0.332533, 0.161787, 0.394456, 0.075382, 0.624294, 0.419409, 0.196961, 0.508851, 0.078281, 0.249060, 0.390322, 0.099849, 0.094464,
  4.727182, 0.858151, 4.008358, 1.240275, 2.784478, 1.223828, 0.611973, 1.739990, 0.990012, 0.064105, 0.182287, 0.748683, 0.346960, 0.361819, 1.338132,
  2.139501, 0.578987, 2.000679, 0.425860, 1.143480, 1.080136, 0.604545, 0.129836, 0.584262, 1.033739, 0.302936, 1.136863, 2.020366, 0.165001, 0.571468, 6.472279,
  0.180717, 0.593607, 0.045376, 0.029890, 0.670128, 0.236199, 0.077852, 0.268491, 0.597054, 0.111660, 0.619632, 0.049906, 0.696175, 2.457121, 0.095131, 0.248862, 0.140825,
  0.218959, 0.314440, 0.612025, 0.135107, 1.165532, 0.257336, 0.120037, 0.054679, 5.306834, 0.232523, 0.299648, 0.131932, 0.481306, 7.803902, 0.089613, 0.400547, 0.245841, 3.151815,
  2.547870, 0.170887, 0.083688, 0.037967, 1.959291, 0.210332, 0.245034, 0.076701, 0.119013, 10.649107, 1.702745, 0.185202, 1.898718, 0.654683, 0.296501, 0.098369, 2.188158, 0.189510, 0.249313
};

static double lgFrequencies[20] = {
  0.079066, 0.055941, 0.041977, 0.053052, 0.012937, 0.040767, 0.071586, 0.057337, 0.022355, 0.062157, 0.099081, 0.064600, 0.022951, 0.042302, 0.044040, 0.061197, 0.053287, 0.012066, 0.034155, 0.069147
};


//---------------------------------------------------
// Kosiol, C., and Goldman, N. 2005. Different versions of the Dayhoff rate matrix. 
// Molecular Biology and Evolution 22:193-199.
//
// See also http://www.ebi.ac.uk/goldman/dayhoff
//---------------------------------------------------
static double dayRelativeRates[] = {
  0.267828,
  0.984474, 0.327059,
  1.199805, 0.000000, 8.931515,
  0.360016, 0.232374, 0.000000, 0.000000,
  0.887753, 2.439939, 1.028509, 1.348551, 0.000000,
  1.961167, 0.000000, 1.493409, 11.388659, 0.000000, 7.086022,
  2.386111, 0.087791, 1.385352, 1.240981, 0.107278, 0.281581, 0.811907,
  0.228116, 2.383148, 5.290024, 0.868241, 0.282729, 6.011613, 0.439469, 0.106802,
  0.653416, 0.632629, 0.768024, 0.239248, 0.438074, 0.180393, 0.609526, 0.000000, 0.076981,
  0.406431, 0.154924, 0.341113, 0.000000, 0.000000, 0.730772, 0.112880, 0.071514, 0.443504, 2.556685,
  0.258635, 4.610124, 3.148371, 0.716913, 0.000000, 1.519078, 0.830078, 0.267683, 0.270475, 0.460857, 0.180629,
  0.717840, 0.896321, 0.000000, 0.000000, 0.000000, 1.127499, 0.304803, 0.170372, 0.000000, 3.332732, 5.230115, 2.411739,
  0.183641, 0.136906, 0.138503, 0.000000, 0.000000, 0.000000, 0.000000, 0.153478, 0.475927, 1.951951, 1.565160, 0.000000, 0.921860,
  2.485920, 1.028313, 0.419244, 0.133940, 0.187550, 1.526188, 0.507003, 0.347153, 0.933709, 0.119152, 0.316258, 0.335419, 0.170205, 0.110506,
  4.051870, 1.531590, 4.885892, 0.956097, 1.598356, 0.561828, 0.793999, 2.322243, 0.353643, 0.247955, 0.171432, 0.954557, 0.619951, 0.459901, 2.427202,
  3.680365, 0.265745, 2.271697, 0.660930, 0.162366, 0.525651, 0.340156, 0.306662, 0.226333, 1.900739, 0.331090, 1.350599, 1.031534, 0.136655, 0.782857, 5.436674,
  0.000000, 2.001375, 0.224968, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.270564, 0.000000, 0.461776, 0.000000, 0.000000, 0.762354, 0.000000, 0.740819, 0.000000,
  0.244139, 0.078012, 0.946940, 0.000000, 0.953164, 0.000000, 0.214717, 0.000000, 1.265400, 0.374834, 0.286572, 0.132142, 0.000000, 6.952629, 0.000000, 0.336289, 0.417839, 0.608070,
  2.059564, 0.240368, 0.158067, 0.178316, 0.484678, 0.346983, 0.367250, 0.538165, 0.438715, 8.810038, 1.745156, 0.103850, 2.565955, 0.123606, 0.485026, 0.303836, 1.561997, 0.000000, 0.279379
};

static double dayFrequencies[20] = {
  0.087127, 0.040904, 0.040432, 0.046872, 0.033474, 0.038255, 0.049530, 0.088612, 0.033619, 0.036886, 0.085357, 0.080481, 0.014753, 0.039772, 0.050680, 0.069577, 0.058542, 0.010494, 0.029916, 0.064718
};

template <typename T>
void print_array(T beg, T end, ostream &os)
{
  while (beg != end)
  {
    os << *beg << " ";
    ++beg;
  }
}

template <typename T>
void print_array(T beg, T end, FILE *os)
{
  while (beg != end)
  {
    myPrint(os, faststring(*beg).c_str(), " ");
    ++beg;
  }
}


//***************************************************
//** molecular_model
//***************************************************

// default constructor creates JC;
template <int N>
molecular_model<N>::molecular_model(const faststring &s, enumDataType d):
basic_model(d), seq_length(0), modelname(s), shape(0), ratetype(ratetype_equal), ncat(0), inv(0), modeltype(0), siterates(NULL), model_to_inherit_siterates_from(NULL) // , siterates_initialized(false)
{
  // xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  //   set_model("JC", mymodel::JC, ratetype_equal, 0.5,
  // 	    1, 1, 1, 1, 1, 1, 0, 0, 0.25, 0.25, 0.25, 0.25);
}


template <int N>
molecular_model<N>::molecular_model(const faststring &s, const molecular_model& a):basic_model(a.dataType)
{
  *this = a;  // Copy everything

  modelname = s;
  if (seq_length > 0 && siterates != NULL)
  {
    siterates = new double[seq_length];
    memcpy(siterates, a.siterates, seq_length*sizeof(double) );
    //    siterates_initialized = true;
  }
  else
  {
    siterates = NULL;
    //    siterates_initialized = false;
  }
}

//Requires a completely specified relRates Matrix
template <int N>
void molecular_model<N>::normalize_rrates()
{
  double x;
  staticSquareMatrix<N> A, pi_diag, absolute;

  pi_diag.assign_diagonal(pi);  // Copy frequencies from vector to matrix so we can use matrix algebra here
  A.setToProductOf(pi_diag, relRates);  // A = pi * relRates; // A is a dummy variable
  absolute.setToProductOf(A, pi_diag);
  x = -absolute.trace();

  if (global_verbosity >= 100)
  {
    cerr << "Value of relative rate correction for model " << modelname << ": " << x << endl;
  }

  relRates *= 1/x;
}


//Requires a completely specified relRates Matrix
template <int N>
bool molecular_model<N>::normalize_basefreq(bool warn)
{
  double   sum=0;
  unsigned i;

  for (i=0; i < N; ++i)
  {
    sum     += pi[i];
  }

  double diff_from_1 = fabs(1 - sum);
  
  // Large deviation:
  if (diff_from_1 > 0.0001)
  {
    cerr << "Sum of base frequencies in model " << modelname << " is " << sum << endl;
    return false;
  }
  else
  {
    if (diff_from_1 > 0.5*EPS) // Very small deviation
    {
      for (i=0; i<N; ++i)
        pi[i] /= sum;
      if (warn)
      {
        cerr << "For the model " << modelname << " base frequencies have been corrected slighly. Correction factor:  " << 1.0/sum << endl;
      }
    }
  }
  return true;
}


template <int N>
void molecular_model<N>::set_matrices()
{
  // the relRates matrix must be fill completely before calling this function

  staticSquareMatrix<N> D, D_sqrt, D_sqrt_inv, DsqrtBDsqrt, DsqrtB, U, U_inv;
  staticVector<N>       pi_sqrt, pi_sqrt_inv;

  int i;
  for (i=0; i < N; ++i)
  {
    pi_sqrt[i]     =   std::sqrt(pi[i]);
    pi_sqrt_inv[i] = 1/std::sqrt(pi[i]);
  }

  if (global_verbosity >= 200)
  {
    cerr << endl << "relRates" << endl;
    relRates.print();
  }

  D.assign_diagonal(pi);
  D_sqrt.assign_diagonal(pi_sqrt);
  D_sqrt_inv.assign_diagonal(pi_sqrt_inv);

  if (global_verbosity >= 200)
  {
    cerr << endl  << "D" << endl;
    D.print();
    cerr << endl  << "D_sqrt" << endl;
    D_sqrt.print();
    cerr << endl  << "D_sqrt_inv" << endl;
    D_sqrt_inv.print();
  }

  DsqrtB.setToProductOf(D_sqrt, relRates);

  if (global_verbosity >= 200)
  {
    cerr << endl << "DsqrtB" << endl;
    DsqrtB.print();
  }

  DsqrtBDsqrt.setToProductOf(DsqrtB, D_sqrt);

  if (global_verbosity >= 200)
  {
    cerr << endl  << "DsqrtBDsqrt" << endl;
    DsqrtBDsqrt.print();
  }

  DsqrtBDsqrt.EigenVectorsValues_JCM(U, eigenvalues);

  if (global_verbosity >= 200)
  {
    cerr << endl  << "U" << endl;
    U.print();
    cerr << endl  << "eigenvalues" << endl;
    eigenvalues.print();
  }

  if (! U.orthogonal())
  {
    cerr << "The matrix of eigenvectors is not orthogonal.\nPlease make sure that all relative rates and frequencies are positive numbers (>0)."
         << endl;
    exit(1);
  }

  U_inv = U;
  U_inv.transpose();

  if (global_verbosity >= 200)
  {
    cerr << endl  << "U_inv" << endl;
    U_inv.print();
    cerr << endl;
    staticSquareMatrix<N> U_e;
    U_e.setToProductOf(U, U_inv);
    cerr << endl << "U_e" << endl;
    U_e.print();
  }

  T.setToProductOf(D_sqrt, U);
  T_inv.setToProductOf(U_inv, D_sqrt_inv);   // (D_sqrt*U^-1) = U^-1*D_sqrt^-1

  if (global_verbosity >= 200)
  {
    cerr << endl  << "T" << endl; T.print();
    cerr << endl  << "T_inv" << endl; T_inv.print();
    staticSquareMatrix<N> T_e;
    T_e.setToProductOf(T, T_inv);
    cerr << endl << "T_e" << endl;
    T_e.print();
  }
}


// Old order of nucleotides:
//  W = (rAG * PI_G + rAC * PI_C + rAT * PI_T) / PI_A;
//  X = (rAG * PI_A + rCG * PI_C + rGT * PI_T) / PI_G;
//  Y = (rAC * PI_A + rCG * PI_G + rCT * PI_T) / PI_C;
//  Z = (rAT * PI_A + rGT * PI_G + rCT * PI_C) / PI_T;

template <int N>
void molecular_model<N>::complete_relRateMatrix()
{
  int    i,j;
  double tmp;

  // Set lower triangle
  for (i=0; i<N; ++i) {                // for all rows i
    for (j=i+1; j<N; ++j) {            // for all cols j   (j > i)
      //      cerr << "Copying: " << modelname << "i, j " << i << ", " << j << " " << relRates(i,j) << endl;
      relRates(j,i) = relRates(i,j);   // relRates_ij -> relRates_ji
    }
  }

  // Set diagonal
  for (j=0; j<N; ++j)
  {
    tmp = 0;
    for (i=0; i<N; ++i)
    {
      if (i != j)
      {
        tmp += relRates(i,j) * pi[i];
      }
    }
    relRates(j,j) = -tmp/pi[j];
  }
}


template <int N>
void molecular_model<N>::init_model(double(*rgamma)(double, double), double(*rlfco)())
{
  random_gamma = rgamma;
  random_lf_co = rlfco;
}


template <int N>
bool molecular_model<N>::siterates_initialized() const
{
  return (siterates != NULL); // True of != NULL, i.e. if they are initialized
}

template <int N>
void molecular_model<N>::reset_siterates() 
{
  if (siterates && !model_to_inherit_siterates_from)
    delete [] siterates;
  siterates = NULL;
}


template <int N>
void molecular_model<N>::init_siterates(unsigned len, bool reinit) 
{
  (void) reinit;
  // reinit site rates is not yet implementetd, since the number of occasions this can be used are rare!!!!!!!

  // Nothing to be done if siterates are already initialized.
  if (siterates_initialized())
  {
    cerr << "Internal notice: Call to init_siterates even though they have already been initialized. Modelname: " << get_modelname() << endl;
    return;
  }

  if (len == 0)
  {
    cerr << "Error: Sequence length is specified to be 0, which does not make sense in function init_siterates. For model " << get_modelname() << endl;
    exit(0);;
  }

  bool      use_gamma_rates = (ratetype == ratetype_gamma)   || (ratetype == ratetype_invgamma);
  bool      use_inv_sites   = (ratetype == ratetype_propinv) || (ratetype == ratetype_invgamma) || (ratetype == ratetype_inv_cat_rates);
  bool      use_explicit_rates = (ratetype == ratetype_cat_rates) || (ratetype == ratetype_inv_cat_rates);
  bool      use_distfunction   = (ratetype == ratetype_distfunction) || (ratetype == ratetype_distfunction_inv);

  if ( use_gamma_rates )
  {
    if (shape <= 0)
    {
      cerr << "Error: Shape parameter <= 0 while trying to use a gamma distribution." << endl;
      exit(-100);
    }
  }
  //   else
  //   {
  //     shape = 0;
  //   }

  if ( use_inv_sites )
  {
    if (inv < 0 || inv >= 1)
    {
      cerr << "Internal error: Inv must have a value in the range  0 <= inv < 1." << endl;
      exit(-100);
    }
  }
  //   else
  //   {
  //     inv = 0;
  //   }

  unsigned  i;
  double    tmp_rate;
  double    one_over_shape;
  double    invar_correction;
  double    *dummy = NULL, *siterates_in_categories = NULL;
  double    rate_sum = 0;

  if (use_gamma_rates && ncat > 0)
  {
    dummy                   = new double [ncat];
    siterates_in_categories = new double [ncat];

    if (global_verbosity >= 100)
      cerr << "Calling DiscreteGamma" << endl;
    DiscreteGamma(dummy, siterates_in_categories, shape, shape, ncat, 0);
  }

  if (use_explicit_rates)
  {
    siterates_in_categories = new double [ncat];

    // We simply copy the array, so we do not have to change existing code.
    unsigned i;

    if (global_verbosity >= 100)
      cerr << "Initializing explicit rate categories." << endl;

    for (i=0; i<ncat; ++i)
      siterates_in_categories[i] = cat_rates[i];
  }

  if (use_gamma_rates)
    one_over_shape = 1.0/shape;
  else
    one_over_shape = 0; // For homogeneious data: shape is infinity, one_over_shape is 0.

  if (use_inv_sites)
    invar_correction = 1.0/(1-inv);
  else
    invar_correction = 1;

  seq_length = len;

  //   if (!reinit) // If we reinit siterates we do not have to allocate the memory again
  {
    if (siterates != NULL)
      delete [] siterates;
    siterates  = new double[len];
  }


  Random_DiscretizedDistribution *dd=NULL;

  // Interprete the expression for the distribution
  // Initialize the discretized distribution
  if (use_distfunction)
  {
    math_expression_parser math(distribution_math_expression);
    math.tokenize();
    math.infix2postfix();
    math.postfix2Cexpression_evaluator();
    Cexpression_evaluator *p_exp_eval = math.get_exp_eval_pointer(); // Object we point to will be destroyed in math destructor.
    // We only retain a pointer to the DiscretizedDistribution object.
    // It will be deleted after the siterates of this model have been initialized.
    dd = new Random_DiscretizedDistribution(distribution_a,
                                            distribution_b,
                                            distribution_N,
                                            p_exp_eval,
                                            random_lf_co
                                            );
    // The math object will be destroyed here and with it the object p_exp_eval is pointing to.
    // This object is only used in the DiscretizedDistribution constuctor and
    // after that it is not used any more.
  }

  if (ratetype == ratetype_equal) // equal rates
  {
    for (i=0; i < len; ++i) {
      siterates[i] = 1;
      // Debug
      // cerr << siterates[i] << " ";
    }
  }
  else // non-equal rates
  {
    for (i=0; i < len; ++i) // For all sites in the sequence:
    {
      tmp_rate = invar_correction;
      
      if (use_inv_sites && random_lf_co() < inv) // Do we have an invariant site
      {
        tmp_rate = -1;  // invariant site
      }
      else if ( use_gamma_rates ) //shape = alpha value // gamma distributed site rates:
      {
        if (ncat > 0) // rate categories
        {
          tmp_rate = invar_correction * siterates_in_categories[(int)(random_lf_co()*ncat)];
        }
        else // continuous gamma rates:
        {
          tmp_rate = invar_correction * random_gamma(shape, one_over_shape);
        }
      }
      else if (use_explicit_rates)
      {
        tmp_rate = ncat * invar_correction * siterates_in_categories[(int)(random_lf_co()*ncat)];
      }
      else if (use_distfunction)
      {
        tmp_rate = dd->next_random_value();
        // TEST
        cerr << "Rate: " << tmp_rate << endl;
        rate_sum += tmp_rate; // Adds nothing in case of invariant site. OK
      }

      siterates[i] = tmp_rate;
      // Debug
      //cerr << siterates[i] << " ";
    }
  }
  // Debug

  // Normalize in case of use_distfunction
  // rate_sum coulb be 0 in very rare cases. All sites are invariant and nothing needs to be done.
  if (use_distfunction && rate_sum > 0)
  {
    double correction= len/rate_sum;

    // TEST
    cerr << "correction: " << correction << endl;
    for (i=0; i < len; ++i)
    {
      siterates[i] *= correction;
    }
  }

  if (global_verbosity >= 100)
  {
    cerr.setf(ios::fixed);
    cerr.precision(6);
    cerr << "  Mean siterate and true propinv value: " << get_mean_siterate_value() << " " << compute_true_propinv() << endl;
  }

  // Clean up:
  if (dummy)
    delete [] dummy;

  if (siterates_in_categories)
  {
    if (global_verbosity >= 100)
    {
      cerr << "  rates_cats of Model: " << get_modelname() << ": ";
      print_array(siterates_in_categories, siterates_in_categories+ncat, cerr);
      cerr << endl;
    }
    delete [] siterates_in_categories;
  }

  if (use_distfunction)
  {
    delete dd;
  }

  //  siterates_initialized = true;
} // init_siterates


template <int N>
double molecular_model<N>::get_mean_siterate_value() const
{
  if (!siterates_initialized())
  {
    cerr << "Internal error: Site rates not initialized when trying to compute a mean siterate value for model: " << get_modelname() << ": ";
    exit(0);
  }

  unsigned  i;
  double    m = 0;
  for (i=0; i<seq_length; ++i)
    m += macromax(0.0, siterates[i]); // Invariant sites have a value of -1 in siterates array so that the value can be distinguished from other 0 rates.
  return m/seq_length;
}

template <int N>
const faststring&  molecular_model<N>::get_modelname_siterates_are_inherited_from() const
{
  return modelname_inherit_siterates_from;
}

template <int N>
double molecular_model<N>::compute_true_propinv() const
{
  unsigned  i;
  double    m = 0;

  if (!siterates_initialized())
  {
    cerr << "Internal error: Site rates not initialized when trying to compute the proportion of invariant sizes for model: " << get_modelname() << ": ";
    exit(0);
  }

  for (i=0; i < seq_length; ++i)
    if (siterates[i] < 0)
      m += 1.0;
  return m/seq_length;
}


// New nucleotide and aa order as in my model.h
// ------ old ------- Nucleotide order in evolve: A,G, C, T
template <int N>
void molecular_model<N>::evolve(const faststring &parent_seq, faststring &new_seq, double branchlength) const
{
  double         *siterates_pos;
  double         ran;
  double         probability_sum;
  unsigned char  sym_index;

  const unsigned char*    symbol_to_index = get_symbol_to_index();
  const unsigned char*    index_to_symbol = get_index_to_symbol();

  staticSquareMatrix<N>   P;

  const unsigned char *it     = (unsigned char *) parent_seq.begin();
  const unsigned char *it_end = (unsigned char *) parent_seq.end();

  new_seq.reserve(parent_seq.length()+1);

  siterates_pos = siterates;
  while(it != it_end) {
    sym_index = symbol_to_index[*it];
    assert(sym_index < 255);

    // Compute P:
    {
      staticSquareMatrix<N>  e_lambda, tmp;
      double                 gamma_rate = macromax(0.0, *siterates_pos);
      e_lambda.assign_diagonal(0.0);
      for (int i=0; i<N; ++i)
      {
        e_lambda(i,i) = std::exp(eigenvalues[i]*branchlength*gamma_rate);
      }
      tmp.setToProductOf(T, e_lambda);
      P.setToProductOf(tmp, T_inv);
    }

    if (global_verbosity >= 100)
    {
      cerr << endl  << "P for branch length: " << branchlength << ", gamma rate: " << *siterates_pos << endl;
      P.print();
    }

    // Lets roll the dice:
    ran = random_lf_co();

    // Compute probablity sum from P matrix. They have to be normalised if !=1.
    // For nucleotides we found that it can deviate e.g. for the K2P model
    // by 1e-9. This is expected due to the many numerical operations we
    // used to compute the P matrix.
    // This is important below.

    double prob_sum=0;
    for (int i=0; i<N; ++i)
      prob_sum += P(i,sym_index);

    if (abs(1-prob_sum) > 1e-14) // Correct if difference to 1 is too large. A smaller difference than is tolerated.
    {
      double correction_factor = 1.0/prob_sum; // This is the correction factor
      if (global_verbosity >= 100)
      {
        unsigned long pos = it - (unsigned char *) parent_seq.begin();
        cerr << "MINOR WARNING: Probability sum != 1 for site: " << pos  << " 1+" << prob_sum-1 << endl;
        cerr << "Correction needed. Correction factor: 1+" << correction_factor-1 << endl;
      }
      prob_sum=0;
      for (int i=0; i<N; ++i)
        prob_sum += P(i,sym_index)*correction_factor;

      if (global_verbosity >= 100)
      {
        unsigned long pos = it - (unsigned char *) parent_seq.begin();
        cerr << "...Probability sum after correction: " << pos  << " 1+" << prob_sum-1 << endl;
        if (abs(1-prob_sum) > 1e-14)
        {
          unsigned long pos = it - (unsigned char *) parent_seq.begin();
          cerr << "MAJOR  WARNING: Probability_sum is still != 1 for site: " << pos  << " 1+" << prob_sum-1 << endl;
        }
      }

      probability_sum=0;
      for (int i=0; i<N; ++i)
      {
        // sym_index contains the index of the column corresponding to the origianl symbol.
        probability_sum += P(i,sym_index)*correction_factor;
        if (ran <= probability_sum || i == N-1)
        {
          new_seq.push_back(index_to_symbol[i]);
          break;
        }
      }
    }
    else
    {
      probability_sum=0;
      for (int i=0; i<N; ++i)
      {
        // sym_index contains the index of the column corresponding to the origianl symbol.
        probability_sum += P(i,sym_index);
        if (ran <= probability_sum || i == N-1)
        {
          new_seq.push_back(index_to_symbol[i]);
          break;
        }
      }
    }
    ++it;
    ++siterates_pos;
  }
}

template <int N>
void molecular_model<N>::debug_print_model_symbols() const
{
  int i;
  const unsigned char*    symbol_to_index = get_symbol_to_index();
  const unsigned char*    index_to_symbol = get_index_to_symbol();

  for (i=0; i<256; ++i)
    if (symbol_to_index[i] < 255) {
      cerr << "i: " << i << " (unsigned char)i: " << (unsigned char)i << " (int)symbol_to_index[i]: " << (int)symbol_to_index[i]
      << " index_to_symbol[symbol_to_index[i]]: " << index_to_symbol[symbol_to_index[i]] << endl;
    }
}


template <int N>
void molecular_model<N>::print_relative_site_rates(ostream &os) const
{
  unsigned i;

  if ( !siterates_initialized() )
  {
    os << "Printint relative site rates failed since they are not initialized yet.\n:";
  }

  os << "Relative site rates for model:" << modelname << endl;
  for (i=0; i<seq_length; ++i)
  {
    if (siterates[i] == -1)
      os << i << "\t" << "00" << endl;
    else
      os << i << "\t" << siterates[i] << endl;
  }
}

template <int N>
void molecular_model<N>::print_site_rates_histogramm_data(ostream &os, faststring &desc) const
{
  if ( !siterates_initialized() )
  {
    os << "Printint absolute siterates as histogram failed since they are not initialized yet.\n:";
    exit(-1);
  }

  CHistogram *hist;

  double    *histrates;
  //  unsigned  numrates;

  //    if (global_verbosity >= 200)
  //      os << "Computing histogram for ncat=" << ncat << " and model: " << modelname << endl;

  // Invariant sites are not represented well in this function.
  if (inv != 0)
  {
    histrates = new double [seq_length];

    unsigned i;
    memcpy (histrates, siterates, seq_length*sizeof(double));
    for (i=0; i<seq_length; ++i)
      if (histrates[i] < 0 )
        histrates[i] = 0;
  }
  else
  {
    histrates = siterates;
  }

  if (ncat > 0) // There is only a specified number of different rates. They will be used to produce the histogram.
  {
    hist = new CHistogram((double *)histrates, (double *)(histrates+seq_length), -10);
  }
  else // The distribution must be continuous. So lets use an automatic binning:
  {
    if (shape != 0)
    {
      hist = new CHistogram((double *)histrates, (double *)(histrates+seq_length), 31);
      //       hist = new CHistogram(-0.00000000001, 15, 31);
      //       if (!hist->add((double *)histrates, (double *)(histrates+seq_length)))
      //       {
      // 	cerr << "Internal error: Out of range when adding number to histogram data." << endl;
      // 	exit(-33);
      //       }
    }
    else
    {
      hist = new CHistogram((double *)histrates, (double *)(histrates+seq_length), 0);
    }
  }

  if (inv != 0)
  {
    delete [] histrates;
  }

  //    if (global_verbosity >= 200)
  //      os << "First five site rates: " << siterates[0]
  //  	 << " "  << siterates[1]
  //  	 << " "  << siterates[2]
  //  	 << " "  << siterates[3]
  //  	 << " "  << siterates[4]
  //  	 << endl;


  std::vector<double>   bin_coords;
  std::vector<unsigned> bin_data;

  hist->get_bin_coords(bin_coords);
  bin_data = hist->get_histogram_data();

  unsigned i, n=(unsigned)bin_data.size();

  os << "Site rates histogram data for model:    " << modelname << " for " << desc << endl;

  os << "Number of site rate bins in this model: " << bin_data.size() << endl;
  os << "site-rate number-of-site-rates-in-bin" << endl;
  //  os << "Number of coordinates: " << bin_coords.size() << endl;

  if (ncat > 0)
  {
    //    os << "For case ncat >= 0" << endl;

    for (i=0; i<n; ++i)
    {
      os << bin_coords[i] << " " << bin_data[i] << endl;
    }
  }
  else
  {
    //    os << "For case ncat == 0" << endl;

    for (i=0; i<n; ++i)
    {
      //      os << bin_coords[i] << " " << bin_coords[i+1] << " "  << (bin_coords[i] + bin_coords[i+1])/2.0 << " " << bin_data[i] << endl;
      os << bin_coords[i]  << " " << bin_data[i] << endl;
    }
  }
  os << endl;

  delete hist;
}




template <int N>
faststring molecular_model<N>::get_modelname() const {
  return modelname;
}

template <int N>
double molecular_model<N>::get_shape() const {
  return shape;
}

template <int N>
unsigned molecular_model<N>::get_ncat() const {
  return ncat;
}

template <int N>
double molecular_model<N>::get_inv() const {
  return inv;
}

template <int N>
std::vector<double> molecular_model<N>::get_cat_rates() const {
  return cat_rates;
}


template <int N>
void molecular_model<N>::read_next_model(CFile& is) {
  vector<faststring> splitted;
  faststring         dummy;
  faststring         dummy_lower;
  faststring         end = "end model";
  unsigned           number_of_tokens;
  //  unsigned       i;

  //temporary variables
  //  faststring        tmodelname;                                     // required parameter
  //  staticSquareMatrix<N>                    inputRates;
  //  staticVector<N>                          inputFrequencies;
  //  double    tshape=0, tinv=0;                               // default 0,0

  faststring        input_modeltypename;
  double        input_tstv=1;

  bool          specified_name              = false;
  bool          specified_modeltype         = false;
  bool          specified_rrates            = false;
  bool          specified_tstv              = false;
  bool          specified_shape             = false;
  bool          specified_ncat              = false;
  bool          specified_cat_rates         = false;
  bool          specified_propinv           = false;
  bool          specified_base_frequencies  = false;
  bool          specified_distfunction      = false;


  relRates.assign(1);
  pi.assign(0.25);

  ignore_spaces(is);
  is.getline(dummy);
  dummy.shorten_to_first_occurrence_of('#');
  dummy.removeSpacesBack();
  dummy.replace_char(';' , ' ' );
  dummy.replace_char(':' , ' ' );
  dummy.replace_char('=' , ' ' );

  dummy_lower = dummy;
  dummy_lower.ToLower();
  dummy_lower.replace_char(',' , ' ' );

  while(!is.fail() && dummy_lower != end ) {
    //cerr << "dummy: " <<  dummy << endl;
    if (dummy.size() == 0 || dummy[0] == '#' )
    {
      ignore_spaces(is);
      is.getline(dummy);
      dummy.shorten_to_first_occurrence_of('#');
      dummy.removeSpacesBack();
      dummy.replace_char(';' , ' ' );
      dummy.replace_char(':' , ' ' );
      dummy.replace_char('=' , ' ' );

      dummy_lower = dummy;
      dummy_lower.ToLower();
      dummy_lower.replace_char(',' , ' ' );
      continue;
    }

    number_of_tokens = (unsigned)split(splitted,dummy_lower," \t");
    //cerr << splitted[0] << " " << splitted[1] << endl;

    if(splitted[0] == "name" ) {
      if (specified_name)
        throw readerror(is.line(), "Name specified twice");
      if (number_of_tokens != 2)
        throw readerror(is.line(), "Bad model name");
      modelname = splitted[1];
      specified_name = true;
    }
    else if(splitted[0] == "modeltype" ) {
      if (specified_modeltype)
        throw readerror(is.line(), "modeltype specified twice");
      if (number_of_tokens != 2)
        throw readerror(is.line(), "Bad modeltype");
      specified_modeltype = true;

      // Handle model type in model
      splitted[1].ToUpper();
      input_modeltypename = splitted[1];
    }
    //    else if(splitted[0] == "rates" ) {}  // ignore this for now xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    else if(splitted[0] == "rrates" || splitted[0] == "rrates_upper" || splitted[0] == "rrates_lower")
    {
      unsigned expected_rates = ((N-1)*N)/2;// nuc: N=4, (N-1)*N/2=6 rates, aa: N=20, (N-1)*N/2=190 rates
      unsigned i, j, k;

      if (specified_rrates)
        throw readerror(is.line(), "Rates specified twice");
      if (number_of_tokens < expected_rates ||       // < 6
          number_of_tokens > expected_rates +1)      // > 7
      {
        cerr << "Error: The number of rates specified by rrates is not the as expected for this model:" << endl;
        cerr << "Number of letters in the alphabet: " << N << endl;
        cerr << "Number of rates expected: " << expected_rates+1 << endl;
        cerr << "If only " << expected_rates << " are specified, the last rate is assumed to be 1" << endl;
        cerr << "Number of rates found in model descirption: " << number_of_tokens << endl;
        throw readerror(is.line(), "Wrong number of rates");
      }

      if (splitted[0] == "rrates" || splitted[0] == "rrates_upper")  // rates, row by row, from an upper triangular matrix
      {
        lower_or_upper_matrix_in_input = 'u';
        k=1;
        for (i=0; i<N; ++i) {           // for all rows i  // We fill the upper triangualr matrix with: relRates_ij
          for (j=i+1; j<N; ++j, ++k) {  // for all cols j  //   -  r_01=a r_02=b r_03=c
            if (k == number_of_tokens )                    //         -   r_12=d r_13=e
              break;                                       //                -   r_23=f
            relRates(i,j) = splitted[k].ToDouble();        //                       -
            //	    cerr << "relRates(i,j) for " << i << " " << j << " set to " <<  relRates(i,j) << endl;
          }
        }
      }
      else  // rates, row by row, from a lower triangular matrix are stored in the upper triangle (since matrix symmetric)
      {
        lower_or_upper_matrix_in_input = 'l';
        k=1;
        for (i=0; i<N; ++i) {           // for all rows i  // We fill the upper triangualr matrix with:
          for (j=0; j<i; ++j, ++k) {    // for all cols j  //    relRates_ij (lower) -> relRates_ji (upper)
            if (k == number_of_tokens )                    //    Reading r_ij in order r_10, r_20, r_21, r_30, r_31, r_32, ...
              break;                                       //    Storing r_ij in r_ji
            relRates(j,i) = splitted[k].ToDouble();        //
          }
        }
      }

      // If the last rate is not specified, the default rate of 1 comes from the initialization.
      specified_rrates = true;
    }
    else if (splitted[0] == "tstv"   || splitted[0] == "titv" ||
             splitted[0] == "tratio") {
      if (specified_tstv)
        throw readerror(is.line(), "Tstv specified twice");
      if (number_of_tokens != 2)
        throw readerror(is.line(), "Bad tstv value");
      input_tstv = splitted[1].ToDouble();
      specified_tstv = true;
    }
    else if(splitted[0] == "shape" ) {
      if (specified_shape)
        throw readerror(is.line(), "Shape specified twice");
      if (number_of_tokens != 2)
        throw readerror(is.line(), "Bad shape value");
      shape = splitted[1].ToDouble();
      if (shape <= 0)
        throw readerror(is.line(), "Bad shape value: Shape must be >0.");
      specified_shape = true;
    }
    else if(splitted[0] == "ncat" ) {
      unsigned tmp_ncat;
      if (specified_ncat)
        throw readerror(is.line(), "Ncat specified twice");
      if (number_of_tokens != 2)
        throw readerror(is.line(), "Bad ncat value");
      tmp_ncat = splitted[1].ToInt();
      if (specified_cat_rates && tmp_ncat != ncat) // Check consistency:
      {
        throw readerror(is.line(), "The number of rates found in cat_rates does not match the number of rate categories specified by ncat.");
      }
      ncat = tmp_ncat;
      if (ncat <= 0)
        throw readerror(is.line(), "Bad ncat value: Ncat must be >0.");
      specified_ncat = true;
    }
    else if (splitted[0] == "cat-rates" || splitted[0] == "cat_rates" ) {
      if (specified_cat_rates)
        throw readerror(is.line(), "Cat_rates specified twice");

      // This option implicitly specifies ncat. We allow ncat to be specified indiviually without error message, if the ncat values are consitent.
      if (number_of_tokens < 2)
        throw readerror(is.line(), "No cat values specified in cat_rates");
      if (specified_ncat && (number_of_tokens-1) != ncat)
        throw readerror(is.line(), "The number of rates specified in cat_rates does not correspond to the number specified with ncat.");
      ncat = number_of_tokens-1;

      cat_rates.clear();

      unsigned i;
      double   sum=0;
      double   val;

      for (i=0; i < ncat; ++i)
      {
        val = splitted[i+1].ToDouble();
        cat_rates.push_back(val);
        sum += val;
      }

      // Normalize values to that the mean rate is 1:
      for (i=0; i < ncat; ++i)
      {
        cat_rates[i] /= sum;
      }
      specified_cat_rates = true;
    }
    else if (splitted[0] == "siterate-distribution")
    {
      // We need to use the string with chars in upper and lower case - and the split needs to respect quotes
      number_of_tokens = (unsigned)split_respect(splitted,dummy," \t");

      if (number_of_tokens != 5)
        throw readerror(is.line(), "Model keyword: siterate-distribution: Wrong number of options.");

      distribution_a = splitted[1].ToDouble();
      distribution_b = splitted[2].ToDouble();
      distribution_N = splitted[3].ToUnsigned();
      distribution_math_expression = splitted[4];
      distribution_math_expression.unquote();
      specified_distfunction = true;

      if (distribution_a >= distribution_b)
        throw readerror(is.line(), "Model keyword: siterate-distribution: Options: a b N \"Math expression for distribution function\".");

      cerr << "Found siterate-distribution functionn in model: " << modelname << endl;
      cerr << distribution_a << " " << distribution_b << " " << distribution_N << " "
      << distribution_math_expression << endl;
    }
    else if(splitted[0] == "inv" || splitted[0] == "pinv") {
      if (specified_propinv)
        throw readerror(is.line(), "Inv specified twice");
      if (number_of_tokens != 2)
        throw readerror(is.line(), "Bad inv value");
      inv = splitted[1].ToDouble();
      if (inv < 0 || inv >= 1)
        throw readerror(is.line(), "Bad inv value: Inv must be in the range 0 <= inv < 1.");
      specified_propinv = true;
    }
    else if(splitted[0] == "basefreq" ) {
      unsigned    i;
      double      sum       = 0;

      if (specified_base_frequencies)
        throw readerror(is.line(), "Base frequencies specified twice");
      if (number_of_tokens < N || number_of_tokens > N+1)
        throw readerror(is.line(), "Wrong number of base frequency values");

      for (i=1; i < number_of_tokens; ++i)
      {
        pi[i-1]  = splitted[i].ToDouble();
        sum     += pi[i-1];
      }
      if (number_of_tokens < N+1)
      {
        pi[N-1] = 1 - sum;
        sum += pi[N-1];

        if (global_verbosity >= 100)
        {
          cerr << "For model " << modelname << " the last base frequency was not specified and has been determined to be: " << pi[N-1] << endl;
        }
      }
      else
      {
        // Normalization will be done further down for all models
        //	if (!normalize_basefreq())
        //	  throw readerror(is.line(), "Base frequencies do not add up to 1.");
      }

      //       if (global_verbosity > 200)
      //       {
      // 	cerr << "For model " << modelname << " the sum of specified base frequency is: " << sum << endl;
      //       }

      specified_base_frequencies = true;
    }
    else if(splitted[0] == "inherit-siterates" )
    {
      modelname_inherit_siterates_from = splitted[1];
    }
    else
    {
      //cerr << "Fehler bei end model (" << dummy <<")"<< endl;
      faststring errormsg = "Unknown keyword: " + splitted[0];
      throw readerror(is.line(), errormsg.c_str());
    }
    ignore_spaces(is);
    is.getline(dummy);
    dummy.shorten_to_first_occurrence_of('#');
    dummy.removeSpacesBack();
    dummy.replace_char( ';' , ' ' );
    dummy.replace_char( ':' , ' ' );
    dummy.replace_char( '=' , ' ' );

    dummy_lower = dummy;
    dummy_lower.ToLower();
    dummy_lower.replace_char( ',' , ' ' );
  } // end while

  // Check whether all necessary values are specified and whether combinations are correct:
  if (!specified_name)
    throw readerror(is.line(), "Name of model needs to be specified");

  if (specified_ncat && !(specified_shape || specified_cat_rates) )
  {
    cerr << specified_ncat << " " << specified_shape << " " << specified_cat_rates << endl;
    throw readerror(is.line(), "Specifying the number of rate categories (ncat) only makes sense if a shape paramter has been specified for gamma distributed site rates, or if cat_rates have been specified.");
  }

  if (specified_propinv && specified_cat_rates)
  {
    ratetype = ratetype_inv_cat_rates;
  }
  else if (specified_cat_rates)
  {
    ratetype = ratetype_cat_rates;
  }
  else if (specified_propinv && specified_distfunction)
  {
    ratetype = ratetype_distfunction_inv;
  }
  else if (specified_distfunction)
  {
    ratetype = ratetype_distfunction;
  }
  else if (specified_shape && specified_propinv)
  {
    ratetype = ratetype_invgamma;
  }
  else if (specified_shape)
  {
    ratetype = ratetype_gamma;
  }
  else if (specified_propinv)
  {
    ratetype = ratetype_propinv;
  }
  else
  {
    ratetype = ratetype_equal;
  }

  set_model_specific_parameters(&is,
                                input_modeltypename,
                                input_tstv,
                                specified_tstv,
                                specified_rrates,
                                specified_base_frequencies);

  complete_relRateMatrix();
  normalize_rrates();

  // Warn if sum of basefreq deviates from 1 only if basefrequences
  // have been specified by the user.
  // For all aa-models the frequencies are corrected slighly.
  if (!normalize_basefreq(specified_base_frequencies))
    throw readerror(is.line(), "Base frequencies do not add up to 1.");

  //     ///XXXXXXXXXXXXXXXXXXXXXXXxxxxxxxxxxxxxxxxxxxxxxxxx
  //   set_model(tmodelname, tmodeltype, tratetype, ttstv,
  // 	    trAC, trAG, trAT, trCG, trCT, trGT, tshape, tinv,
  // 	    tpia, tpic, tpig, tpit);
}


// Seems to be not fully implemented:
// -does not set ratetype
// -does not handle distribution function
template <int N>
void molecular_model<N>::set_model(faststring modelname_param,
                                   faststring modeltype_param,
                                   vector<double> *rrates_param, // Parameters are supplied as in an upper triangular matrix.
                                   vector<double> *base_param,
                                   double         shape_param,
                                   double         pinv_param,
                                   unsigned       ncat_param,
                                   double         *tstv_param

                                   )

{

  //temporary variables
  //  faststring        tmodelname;                                     // required parameter
  //  staticSquareMatrix<N>                    inputRates;
  //  staticVector<N>                          inputFrequencies;
  //  double    tshape=0, tinv=0;                               // default 0,0

  faststring        input_modeltypename;
  double        input_tstv=1;

  //   bool          specified_name              = false;
  //   bool          specified_modeltype         = false;
  bool          specified_rrates            = false;
  bool          specified_tstv              = false;
  //   bool          specified_shape             = false;
  //   bool          specified_ncat              = false;
  //   bool          specified_propinv           = false;
  bool          specified_base_frequencies  = false;

  relRates.assign(1);
  pi.assign(0.25);

  modelname           = modelname_param;
  input_modeltypename = modeltype_param;

  //XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


  int expected_rates = ((N-1)*N)/2;// nuc: N=4, (N-1)*N/2=6 rates, aa: N=20, (N-1)*N/2=190 rates
  int i, j, k;

  if (rrates_param != NULL)
  {
    short number_of_rates = rrates_param->size();

    if (number_of_rates < expected_rates ||       // < 6
        number_of_rates > expected_rates)         // > 6
      throw setmodelerror("Wrong number of rates");

    // Copy rrates
    k=0;
    for (i=0; i<N; ++i) {           // for all rows i  // We fill the upper triangualr matrix with: relRates_ij
      for (j=i+1; j<N; ++j, ++k) {  // for all cols j  //   -  r_01=a r_02=b r_03=c
        //	  if (k == number_of_tokens )                    //         -   r_12=d r_13=e
        //	    break;                                       //                -   r_23=f
        relRates(i,j) = (*rrates_param)[k];        //                       -
      }
    }
    lower_or_upper_matrix_in_input = 'u';
    specified_rrates = true;
  }

  shape = shape_param;
  if (shape <= 0)
    throw setmodelerror("Bad shape value: Shape must be >0.");

  ncat = ncat_param;
  if (ncat <= 0)
    throw setmodelerror("Bad ncat value: Ncat must be >0.");

  inv = pinv_param;
  if (inv < 0 || inv >= 1)
    throw setmodelerror("Bad inv value: Inv must be in the range 0 <= inv < 1.");

  if (base_param != NULL)
  {
    for (i=0; i < N; ++i)
    {
      pi[i]  = (*base_param)[i];
    }
    specified_base_frequencies = true;
  }

  if (tstv_param != NULL)
  {
    input_tstv = *tstv_param;
    specified_tstv = true;
  }

  set_model_specific_parameters(NULL,
                                modelname,
                                input_tstv,
                                specified_tstv,
                                specified_rrates,
                                specified_base_frequencies);

  complete_relRateMatrix();
  normalize_rrates();

  if (!normalize_basefreq(specified_base_frequencies))
    throw setmodelerror("Base frequencies do not add up to 1.");

  //     ///XXXXXXXXXXXXXXXXXXXXXXXxxxxxxxxxxxxxxxxxxxxxxxxx
  //   set_model(tmodelname, tmodeltype, tratetype, ttstv,
  // 	    trAC, trAG, trAT, trCG, trCT, trGT, tshape, tinv,
  // 	    tpia, tpic, tpig, tpit);
}



template <int N>
void molecular_model<N>::get_random_sequence(faststring& s, unsigned len) const
{
  double                 x;
  unsigned               pos, i;
  const unsigned char*   index_to_symbol = get_index_to_symbol();
  double                 probability_sum[N];

  probability_sum[0] = pi[0];
  for (i=1; i<N; ++i)
    probability_sum[i] = probability_sum[i-1] + pi[i];

  s="";
  s.reserve(len);

  for (pos=0; pos < len; ++pos)
  {
    x = (*random_lf_co)();          // must be a random number in interval [0,1)

    for (i=0; i<N; ++i)
    {
      if (x <= probability_sum[i])
        break;
    }
    s.push_back(index_to_symbol[i]);
  }
}

//***************************************************
//** nuc_model
//***************************************************

// For some models, some parameters need to be set before calling this function: tstv
void nuc_model::set_rates_nst12()
{
  // We follow Felsenstein, Inferring Phylogenies, First print, page 200-203
  // Starting from the Timura Nei model, we obtain everything we need:

  double rho;  // = alpha_R/alpha_Y
  double alpha_R;
  double alpha_Y;
  double beta;
  double piR = pi[nA] + pi[nG];
  double piY = pi[nC] + pi[nT];

  switch ((enummodeltype)modeltype) {
    case nuc_model::JC:
      pi[nA] = pi[nC] = pi[nG] = pi[nT] = 0.25;
      tstv = 0.5;
      rho = 1;
      break;

    case nuc_model::F81:
      tstv = (pi[nA]*pi[nG] + pi[nC]*pi[nT])/(piR*piY);
      rho = piR/piY;
      break;

    case nuc_model::K2P:    // take tstv as supplied by caller
      pi[nA] = pi[nC] = pi[nG] = pi[nT] = 0.25;
      rho = 1;
      break;

    case nuc_model::F84:    // take tstv as supplied by caller
      rho = 1;
      break;

    case nuc_model::HKY:  // take tstv as supplied by caller
      rho = piR/piY;
      break;

    default:
      throw 1;
  }

  alpha_Y = ( piR*piY*tstv-pi[nA]*pi[nG]-pi[nC]*pi[nT] )/
  ( 2.0*(1+tstv)*(piY*pi[nA]*pi[nG]*rho + piR*pi[nC]*pi[nT]) );
  alpha_R = rho * alpha_Y;
  beta    = 1.0/(2.0*piR*piY*(1+tstv));

  // we set the upper triangular matrix
  lower_or_upper_matrix_in_input = 'u';

  relRates(nA,nC) =               beta;
  relRates(nA,nG) = alpha_R/piR + beta;
  relRates(nA,nT) =               beta;
  relRates(nC,nG) =               beta;
  relRates(nC,nT) = alpha_Y/piY + beta;
  relRates(nG,nT) =               beta;
}

// Does only need non-normalized upper triangular matrix
double nuc_model::compute_tstv() const
{
  return (  pi[nA]*pi[nG]*relRates(nA,nG) + pi[nC]*pi[nT]*relRates(nC,nT))/
  (  pi[nA]*pi[nC]*relRates(nA,nC) + pi[nA]*pi[nT]*relRates(nA,nT)
   + pi[nC]*pi[nG]*relRates(nC,nG) + pi[nG]*pi[nT]*relRates(nG,nT));
}

void nuc_model::print(ostream& os, unsigned flag) {

  //  DEBUGCODE( debug_print_model_symbols(); );

  (void) flag;

  os << "begin nuc-model" << endl;
  os << "name: " << modelname << endl;
  os << "modeltype= " << modeltypenames[modeltype] << endl;
  if ( (enummodeltype)modeltype == K2P || (enummodeltype)modeltype == F84 || (enummodeltype)modeltype == HKY )
  {
    os << "tstv: ";
    os << tstv << endl;
  }
  else if ( (enummodeltype)modeltype == GTR)
  {
#define GTnormalized
#ifdef  GTnormalized
    if(relRates(nG,nT) == 0) {
      exit(1);
    }
    os << "rrates: ";
    os << relRates(nA,nC)/relRates(nG,nT) << " "
    << relRates(nA,nG)/relRates(nG,nT) << " "
    << relRates(nA,nT)/relRates(nG,nT) << " "
    << relRates(nC,nG)/relRates(nG,nT) << " "
    << relRates(nC,nT)/relRates(nG,nT) << " "
    << relRates(nG,nT)/relRates(nG,nT) << endl;
#else
    os << relRates(nA,nC) << " " << relRates(nA,nG) << " " << relRates(nA,nT) << " "
    << relRates(nC,nG) << " " << relRates(nC,nT) << " " << relRates(nG,nT) << endl;
#endif
  }
  os << "# normalized rrates: ";
#ifdef GTnormalized
  os << relRates(nA,nC) << " " << relRates(nA,nG) << " " << relRates(nA,nT) << " "
  << relRates(nC,nG) << " " << relRates(nC,nT) << " " << relRates(nG,nT) << endl;
#else
  os << relRates(nA,nC)/relRates(nG,nT) << " " << relRates(nA,nG)/relRates(nG,nT) << " " << relRates(nA,nT)/relRates(nG,nT) << " "
  << relRates(nC,nG)/relRates(nG,nT) << " " << relRates(nC,nT)/relRates(nG,nT) << " " << relRates(nG,nT)/relRates(nG,nT) << endl;
#endif

  os << "# recomputed tstv: " << compute_tstv() << endl;

  if (ratetype == ratetype_gamma || ratetype == ratetype_invgamma )
    os << "shape: " << shape << endl;
  if (ncat > 0)
    os << "ncat: " << ncat << endl;
  if (ratetype == ratetype_propinv || ratetype == ratetype_invgamma )
    os << "inv: " << inv << endl;

  os << "basefreq: ";
  int i;
  for (i = nA; i < nT; ++i)
    os << pi[i] << " ";
  os << pi[i] << endl;
  os << "end model" << endl << endl;
}


void nuc_model::print(FILE *os, unsigned flag) {

  //  DEBUGCODE( debug_print_model_symbols(); );

  (void) flag;

  myPrint(os, "begin nuc-model\n");
  myPrint(os, "name: ", modelname.c_str(), "\n");
  myPrint(os, "modeltype= ", modeltypenames[modeltype], "\n");
  if ( (enummodeltype)modeltype == K2P || (enummodeltype)modeltype == F84 || (enummodeltype)modeltype == HKY )
  {
    myPrint(os, "tstv: ", faststring(tstv).c_str(), "\n");
  }
  else if ( (enummodeltype)modeltype == GTR)
  {
#define GTnormalized
#ifdef  GTnormalized
    if(relRates(nG,nT) == 0) {
      exit(1);
    }
    myPrint(os, "rrates: ");
    myPrint(os, faststring(relRates(nA,nC)/relRates(nG,nT)).c_str(), " ");
    myPrint(os, faststring(relRates(nA,nG)/relRates(nG,nT)).c_str(), " ");
    myPrint(os, faststring(relRates(nA,nT)/relRates(nG,nT)).c_str(), " ");
    myPrint(os, faststring(relRates(nC,nG)/relRates(nG,nT)).c_str(), " ");
    myPrint(os, faststring(relRates(nC,nT)/relRates(nG,nT)).c_str(), " ");
    myPrint(os, faststring(relRates(nG,nT)/relRates(nG,nT)).c_str(), "\n");
#else
    myPrint(os, faststring(relRates(nA,nC)).c_str(), " ");
    myPrint(os, faststring(relRates(nA,nG)).c_str(), " ");
    myPrint(os, faststring(relRates(nA,nT)).c_str(), " ");
    myPrint(os, faststring(relRates(nC,nG)).c_str(), " ");
    myPrint(os, faststring(relRates(nC,nT)).c_str(), " ");
    myPrint(os, faststring(relRates(nG,nT)).c_str(), "\n");
#endif
  }
  myPrint(os, "# normalized rrates: ");
#ifdef GTnormalized
  myPrint(os, faststring(relRates(nA,nC)).c_str(), " ");
  myPrint(os, faststring(relRates(nA,nG)).c_str(), " ");
  myPrint(os, faststring(relRates(nA,nT)).c_str(), " ");
  myPrint(os, faststring(relRates(nC,nG)).c_str(), " ");
  myPrint(os, faststring(relRates(nC,nT)).c_str(), " ");
  myPrint(os, faststring(relRates(nG,nT)).c_str(), "\n");
#else
  myPrint(os, faststring(relRates(nA,nC)/relRates(nG,nT)).c_str(), " ");
  myPrint(os, faststring(relRates(nA,nG)/relRates(nG,nT)).c_str(), " ");
  myPrint(os, faststring(relRates(nA,nT)/relRates(nG,nT)).c_str(), " ");
  myPrint(os, faststring(relRates(nC,nG)/relRates(nG,nT)).c_str(), " ");
  myPrint(os, faststring(relRates(nC,nT)/relRates(nG,nT)).c_str(), " ");
  myPrint(os, faststring(relRates(nG,nT)/relRates(nG,nT)).c_str(), "\n");
#endif

  myPrint(os, "# recomputed tstv: ", faststring(compute_tstv()).c_str(), "\n");

  if (ratetype == ratetype_gamma || ratetype == ratetype_invgamma )
    myPrint(os, "shape: ", faststring(shape).c_str(), "\n");

  if (ncat > 0)
    myPrint(os, "ncat: ", faststring(ncat).c_str(), "\n");
  if (ratetype == ratetype_propinv || ratetype == ratetype_invgamma )
    myPrint(os, "inv: ", faststring(inv).c_str(), "\n");

  myPrint(os, "basefreq: ");
  int i;
  for (i = nA; i < nT; ++i)
    myPrint(os, faststring(pi[i]).c_str(), " ");
  myPrint(os, faststring(pi[i]).c_str(), "\n");
  myPrint(os, "end model\n\n");
}




void nuc_model::set_model_specific_parameters(CFile  *is,
                                              faststring& input_modeltypename,
                                              double  input_tstv,
                                              bool    specified_tstv,
                                              bool    specified_rrates,
                                              bool    specified_base_frequencies)
{
  int i;

  for (i=JC; i<=GTR; ++i)
  {
    if ( input_modeltypename == modeltypenames[i] )
      break;
  }
  if (i > GTR)
  {
    if (is == NULL)
    {
      throw setmodelerror("Unknown modeltype");
    }
    else
    {
      throw readerror(is->line(), "Unknown modeltype");
    }
  }

  modeltype = (int)i;

  if ( (modeltype == JC || modeltype == F81) && (specified_tstv || specified_rrates) )
  {
    if (is == NULL)
    {
      throw setmodelerror("In the model ending on this line: For a JC or F81 model, tstv and rrates cannot be specified");
    }
    else
    {
      throw readerror(is->line()-1, "In the model ending on this line: For a JC or F81 model, tstv and rrates cannot be specified");
    }
  }

  if ( (modeltype == JC || modeltype == K2P) && (specified_base_frequencies) )
  {
    if (is == NULL)
    {
      throw setmodelerror("In the model ending on this line: For a JC or K2P base frequencies cannot be specified");
    }
    else
    {
      throw readerror(is->line()-1, "In the model ending on this line: For a JC or K2P base frequencies cannot be specified");
    }
  }

  if ( (modeltype == K2P || modeltype == F84 || modeltype == HKY) && (specified_rrates) )
  {
    if (is == NULL)
    {
      throw setmodelerror("In the model ending on this line: For a K2P or F84 or HKY model, rrates cannot be specified");
    }
    else
    {
      throw readerror(is->line()-1, "In the model ending on this line: For a K2P or F84 or HKY model, rrates cannot be specified");
    }
  }

  if (modeltype == GTR && (specified_tstv) )
  {
    if (is == NULL)
    {
      throw setmodelerror("In the model ending on this line: For a GTR model, tstv cannot be specified");
    }
    else
    {
      throw readerror(is->line()-1, "In the model ending on this line: For a GTR model, tstv cannot be specified");
    }
  }

  if (modeltype == GTR )
  {
    if (!specified_rrates)
    {
      if (is == NULL)
      {
        throw setmodelerror("In the model ending on this line: For a GTR model rrates need to be specified.");
      }
      else
      {
        throw readerror(is->line()-1, "In the model ending on this line: For a GTR model rrates need to be specified.");
      }
    }
  }


  // For these models, the user is expected to specify base frequencies
  // What if he/she hasn't done so?
  // Should we thow an error or use default frequencies??
  //   if (modeltype == F81 || modeltype == F84 || modeltype == HKY || modeltype == GTR)
  //   {

  //   }

  if (specified_tstv)
  {
    tstv = input_tstv;
  }

  // the nst == 1 and 2 models are treated here:
  if ( modeltype == JC   || modeltype == F81 ||
      modeltype == K2P  || modeltype == F84 || modeltype == HKY )
  {
    set_rates_nst12();
  }
  else if (modeltype == GTR)
  {
    tstv = compute_tstv();
  }
  else
  {
    throw 15;
  }
}


int nuc_model::get_modeltype() const
{
  return modeltype;
}

faststring  nuc_model::get_modeltypename() const
{
  return modeltypenames[modeltype];
}

// const char[][6] nuc_model::get_modeltypenames() const
// {
//   return modeltypenames;
// }

double nuc_model::get_rAC() const {
  return relRates(nA,nC);
}

double nuc_model::get_rAG() const {
  return relRates(nA,nG);
}

double nuc_model::get_rAT() const {
  return relRates(nA,nT);
}

double nuc_model::get_rCG() const {
  return relRates(nC,nG);
}

double nuc_model::get_rCT() const {
  return relRates(nC,nT);
}

double nuc_model::get_rGT() const {
  return relRates(nG,nT);
}

double nuc_model::get_tstv() const {
  return tstv;
}

double nuc_model::get_PI_A() const {
  return pi[nA];
}

double nuc_model::get_PI_G() const {
  return pi[nG];
}

double nuc_model::get_PI_T() const {
  return pi[nT];
}

double nuc_model::get_PI_C() const {
  return pi[nC];
}



//***************************************************
//** a_model
//***************************************************
int  aa_model::get_modeltype() const
{
  return modeltype;
}

faststring  aa_model::get_modeltypename() const
{
  return modeltypenames[modeltype];
}

void aa_model::set_model_specific_parameters(CFile *is,
                                             faststring& input_modeltypename,
                                             double  input_tstv,  // This parameter has no meaning for aa models. Its used merely for compatibility in the class hierarchy.
                                             bool    specified_tstv,
                                             bool    specified_rrates,
                                             bool    specified_base_frequencies)
{
  int i,j,k;

  (void) input_tstv;

  if (is == NULL)
  {
    cerr << "Internal problem 33. Please report the occurence of this problem to the programmer. " << endl;
    exit(-33);
  }


  for (i=USER; i <= DAY; ++i)
  {
    if ( input_modeltypename == modeltypenames[i] )
      break;
  }

  if (i > DAY)
    throw readerror(is->line(), "Unknown modeltype");
  modeltype = i;

  if ( specified_tstv )
    throw readerror(is->line()-1, "In the model ending on this line: tstv cannot be specified in protein model.");

  if ( modeltype > USER && specified_rrates )
    throw readerror(is->line()-1, "In the model ending on this line: A relative rate matrix for protein models can only be specified in the USER model.");

  if ( modeltype > USER && specified_base_frequencies )
    throw readerror(is->line()-1, "In the model ending on this line: Base frequencies for protein models can only be specified in the USER model.");

  double *relR, *aaFreq;
  bool   is_lower_triangular_matrix = true;
  // There is a compiler warning that these can be used uninitialised. I think
  // this is due to complex if, the compiler does not understand. But anyway.
  // We want to be sure.

  bool   resR_and_aaFreq_initialised = false;

  if (modeltype == JTT)
  {
    relR   = jttRelativeRates;
    aaFreq = jttFrequencies;
    resR_and_aaFreq_initialised = true;
  }
  else if (modeltype == WAG_OLD)
  {
    relR   = wag_oldRelativeRates;
    aaFreq = wag_oldFrequencies;
    is_lower_triangular_matrix = false;
    resR_and_aaFreq_initialised = true;
  }
  else if (modeltype == WAG)
  {
    relR   = wagRelativeRates;
    aaFreq = wagFrequencies;
    resR_and_aaFreq_initialised = true;
  }
  else if (modeltype == WAG_STAR)
  {
    relR   = wagstarRelativeRates;
    aaFreq = wagstarFrequencies;
    resR_and_aaFreq_initialised = true;
  }
  else if (modeltype == LG)
  {
    relR   = lgRelativeRates;
    aaFreq = lgFrequencies;
    resR_and_aaFreq_initialised = true;
  }
  else if (modeltype == DAY)
  {
    relR   = dayRelativeRates;
    aaFreq = dayFrequencies;
    resR_and_aaFreq_initialised = true;
  }
  else if (modeltype != USER)
  {
    aaFreq = NULL;
    faststring err = "No data is available for this model: " + input_modeltypename;
    throw readerror(is->line()-1, err.c_str() );
  }

  if (modeltype == USER) // USER type
  {
    if (!specified_rrates)
    {
      throw readerror(is->line()-1, "In the model ending on this line: If the model is of type USER, rrates must be specified in the model.");
    }

    if (!specified_base_frequencies)
    {
      throw readerror(is->line()-1, "In the model ending on this line: If the model is of type USER, basefreq must be specified in the model.");
    }
  }

  if (modeltype != USER) // relR and aaFreq still need to be copied from the models specified.
  {
    if (!resR_and_aaFreq_initialised)
    {
      cerr << "ERROR: relR and aaFreq are used uninitialized. This has to be fixed.\n";
      exit(0);
    }

    if ( !is_lower_triangular_matrix ) // for the elements, row by row, of an upper triangular matrix
    {
      k=0;
      for (i=0; i < 20; ++i) {           // for all rows i
        for (j=i+1; j < 20; ++j, ++k) {  // for all cols j
          relRates(i,j) = relR[k];       // relRates_ij in order (0,1) (0,2) (0,3) ... (1,2) (1,3) ...
        }
      }
    }
    else  // For elements, row by row, from a lower triangular matrix, e.g. the PAML standard representation.
    {     // Since the rate matrix is symmetric, the lower diagonal can be inserted into the upper diagonal.
      k=0;
      for (i=0; i < 20; ++i) {         // For all rows i (which become cols in relRates)
        for (j=0; j < i; ++j, ++k) {   // For all cols j (which become rows in relRates)
          relRates(j,i) = relR[k];     // We read an (i,j) from the lower triangle and store it as (j,i)
        }                              // in the upper triangle.
      }
    }

    for (i=0; i < 20; ++i) {
      pi[i] = aaFreq[i];
    }
  }
}

void aa_model::print(ostream& os, unsigned flag) {

  //  DEBUGCODE( debug_print_model_symbols(); );

  //  cerr << "Print aa-model has been called." << endl;

  (void) flag;
  int N = 20;

  os << "begin aa-model" << endl;
  os << "name: " << modelname << endl;
  os << "modeltype= " << modeltypenames[modeltype] << endl;

  if (ratetype == ratetype_gamma || ratetype == ratetype_invgamma )
    os << "shape: " << shape << endl;
  if (ncat > 0)
    os << "ncat: " << ncat << endl;
  if (ratetype == ratetype_propinv || ratetype == ratetype_invgamma )
    os << "inv: " << inv << endl;

  int i, j;
  os << "rrates: ";
  for (i=0; i < N; ++i) {            // for all rows i
    for (j=i+1; j < N; ++j) {        // for all cols j in the upper triangle
      os << relRates(i,j) << " ";
    }
    os << " ";
  }
  os << endl;

  os << "# rrates_lower: ";
  for (i=0; i < N; ++i) {            // for all rows i
    for (j=0; j < i; ++j) {          // for all cols j in the lower triangle
      os << relRates(i,j) << " ";
    }
    os << " ";
  }
  os << endl;

  os << "basefreq: ";
  for (i = 0; i < N-1; ++i)
    os << pi[i] << " ";
  os << pi[i] << endl;
  os << "end model" << endl << endl;
}


void aa_model::print(FILE *os, unsigned flag) {

  //  DEBUGCODE( debug_print_model_symbols(); );

  //  cerr << "Print aa-model has been called." << endl;

  (void) flag;
  int N = 20;

  myPrint(os, "begin aa-model\n");
  myPrint(os , "name: ", modelname.c_str(), "\n");
  myPrint(os, "modeltype= ", modeltypenames[modeltype], "\n");

  if (ratetype == ratetype_gamma || ratetype == ratetype_invgamma )
    myPrint(os, "shape: ", faststring(shape).c_str(), "\n");
  if (ncat > 0)
    myPrint(os, "ncat: ", faststring(ncat).c_str(), "\n");
  if (ratetype == ratetype_propinv || ratetype == ratetype_invgamma )
    myPrint(os, "inv: ", faststring(inv).c_str(), "\n");

  int i, j;
  myPrint(os, "rrates: ");
  for (i=0; i < N; ++i) {            // for all rows i
    for (j=i+1; j < N; ++j) {        // for all cols j in the upper triangle
      myPrint(os, faststring(relRates(i,j)).c_str(), " ");
    }
    myPrint(os, " ");
  }
  myPrint(os, "\n");

  myPrint(os, "# rrates_lower: ");
  for (i=0; i < N; ++i) {            // for all rows i
    for (j=0; j < i; ++j) {          // for all cols j in the lower triangle
      myPrint(os, faststring(relRates(i,j)).c_str(), " ");
    }
    myPrint(os, " ");
  }
  myPrint(os, "\n");

  myPrint(os, "basefreq: ");
  for (i = 0; i < N-1; ++i)
    myPrint(os, faststring(pi[i]).c_str(), " ");
  myPrint(os, faststring(pi[i]).c_str(), "\n");
  myPrint(os, "end model\n\n");
}



//***************************************************
// End of file
//***************************************************

// //default constructor creates JC;
// nuc_model<N>::nuc_model():siterates(NULL), seq_length(0)
// {
//   set_model("JC", nuc_model::JC, ratetype_equal, 0.5,
// 	    1, 1, 1, 1, 1, 1, 0, 0, 0.25, 0.25, 0.25, 0.25);
// }

/*  function pointer is recieved when siterates have to be created
 //needs function pointer to function which creates random numbers
 nuc_model::nuc_model(double (*ran_op)(), double (*ran_hop)()) {
 ranf_op = ran_op;
 ranf_halfop = ran_hop;
 set_model("JC", 1, _equal_, 0.5, 1, 1, 1, 1, 1, 1, 0, 0, 0.25, 0.25, 0.25, 0.25);
 }*/


// //************************************************************************
// // This routine, and routine called form here are the only ones that
// // should be allowed to change tstv and the rates.
// //************************************************************************
// void nuc_model::set_model(const faststring &name,  enummodeltype            xmodeltype, enumratetype r, double xtstv,
// 			      double xrAC,   double xrAG,  double xrAT,
// 			      double xrCG,   double xrCT,  double xrGT,
// 			      double xshape, double xinv,
// 			double xPI_A,  double xPI_C, double xPI_G, double xPI_T);



// void molecular_model::check_model()
// {
//   if (name == "")
//     throw setmodelerror();
//   modelname = name;

//   if(r == ratetype_equal    || r == ratetype_gamma ||
//      r == ratetype_invgamma || r == ratetype_propinv) {
//     ratetype = r;
//   } else {
//     throw setmodelerror();
//   }

//   if (xmodeltype == JC    || xmodeltype == F81 ||
//       xmodeltype == K2P   || xmodeltype == F84 || xmodeltype == HKY ||
//       xmodeltype == GTR) {
//     modeltype = xmodeltype;
//   } else {
//     throw setmodelerror();
//   }


//   // the nst == 1 and 2 models are treated here:
//   if ( modeltype == JC   || modeltype == F81 ||
//        modeltype == K2P  || modeltype == F84 || modeltype == HKY )
//   {
//     set_rates_nst12(modeltype, tstv, PI_A, PI_C, PI_G, PI_T,
// 		    rAC, rAG, rAT, rCG, rCT, rGT);
//   }
//   else if (modeltype == GTR)
//   {
//     normalize_rrates( PI_A,  PI_C,  PI_G,  PI_T, rAC,  rAG,  rAT, rCG,  rCT,  rGT);
//     tstv = compute_tstv( PI_A,  PI_C,  PI_G,  PI_T, rAC,  rAG,  rAT, rCG,  rCT,  rGT);
//   }
//   else
//   {
//     throw 15;
//   }

//   double PI_all = PI_A + PI_C + PI_G + PI_T;

//   DEBUGCODE( cerr << "EPS: " << EPS << " all " << PI_all << endl);
//   if (PI_all < 1-EPS || PI_all > 1+EPS)
//     throw setmodelerror();
// }



template class molecular_model<4>;
template class molecular_model<20>;

