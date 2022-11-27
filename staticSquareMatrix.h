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

#ifndef staticSquareMatrixH
#define staticSquareMatrixH

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <cfloat>
//#include <String>

#include <cstring>

#include "staticVector.h"


template<int N>
class staticSquareMatrix
{
 private:
  double matrix[N][N];
  static void Jacobi_Cyclic_Method(double[], double *, double *);

 public:
  staticSquareMatrix(){};
  //  staticSquareMatrix(double x);
  staticSquareMatrix(const double m[N][N]) {matrix = m; };
  staticSquareMatrix(const double *m)      {std::memcpy(matrix, m, sizeof(matrix));};

  void    print(std::ostream& = std::cout, int=4);

  void    assign(double);
  void    assign(double *);
  void    assign(double [N][N]);

  void    assign_diagonal(double);
  void    assign_diagonal(double *);
  void    assign_diagonal(const staticVector<N>&);

  void    setArrayToThis(double [N][N]);
  void    setArrayToThis(double *);
  
  void    transpose();

  void    add(double);
  void    add(const staticSquareMatrix&);
  void    add_toDiagonal(double);

  void    sub(double);
  void    sub(const staticSquareMatrix&);
  void    sub_fromDiagonal(double);

  void    mult(double);
  void    setToProductOf(const staticSquareMatrix&, const staticSquareMatrix&);
  double  absmax();
  double  determinant();
  double  trace();
/*   void    invert(); */
/*   double  determinant2x2(double, double, double, double); */
/*   double  determinant3x3(double, double, double, double, double, double, double, double, double); */

  void    EigenVectorsValues_JCM(staticSquareMatrix<N> &, staticVector<N> &);
  void    setToExp_series(const staticSquareMatrix&, double, int, double);

  bool    orthogonal();
  bool    symmetric();
  bool    antisymmetric();

  staticSquareMatrix& operator += (double);
  staticSquareMatrix& operator -= (double);
  staticSquareMatrix& operator *= (double);

  staticSquareMatrix& operator += (const staticSquareMatrix &);
  staticSquareMatrix& operator -= (const staticSquareMatrix &);
  staticSquareMatrix& operator *= (const staticSquareMatrix &);

        double& operator()(int i, int j)       { return matrix[i][j]; }
  const double& operator()(int i, int j) const { return matrix[i][j]; }

  staticVector<N> mult_with_vector(const staticVector<N>&);
};



/* template<int N> */
/* inline staticSquareMatrix<N>::staticSquareMatrix(double x) { // creates diagonal matrix with one value */
/*   int i, j; */

/*   std::memset(matrix, 0, sizeof(matrix)); */
/*   for(i=0; i < N; ++i) { */
/*     matrix[i][i] = x; */
/*   } */
/* } */

template<int N>
inline void staticSquareMatrix<N>::print(std::ostream &os, int pres) {
  int i, j;
  os.setf(std::ios::fixed);
  os.precision(pres);
  for (i = 0; i < N; ++i) {
    for (j = 0; j < N-1; ++j) {
      std::cout << matrix[i][j] << " ";
    }
    std::cout << matrix[i][j] << std::endl;
  }
}

template <int N>
inline void staticSquareMatrix<N>::assign(double x) {
  std::memset(matrix, x, sizeof(matrix));
}

template <int N>
inline void staticSquareMatrix<N>::assign(double *m) {
  std::memcpy(matrix, m, sizeof(matrix));
}

template <int N>
inline void staticSquareMatrix<N>::assign(double m[N][N]) {
  std::memcpy(matrix, (double*)m, sizeof(matrix));
}

template <int N>
inline void staticSquareMatrix<N>::assign_diagonal(double x) {
  std::memset(matrix, 0, sizeof(matrix));

  for (int i=0; i<N; ++i)
    matrix[i][i] = x;
}

template <int N>
inline void staticSquareMatrix<N>::assign_diagonal(double *m) {
  std::memset(matrix, 0, sizeof(matrix));

  for (int i=0; i<N; ++i)
    matrix[i][i] = m[i];
}

template <int N>
inline void staticSquareMatrix<N>::assign_diagonal(const staticVector<N>& v) {
  std::memset(matrix, 0, sizeof(matrix));

  for (int i=0; i<N; ++i)
    matrix[i][i] = v.bvector[i];
}

template <int N>
inline void staticSquareMatrix<N>::setArrayToThis(double m[N][N]) {
  std::memcpy((double*)m, matrix, sizeof(matrix));
}

template <int N>
inline void staticSquareMatrix<N>::setArrayToThis(double *m) {
  std::memcpy((double*)m, matrix, sizeof(matrix));
}

template <int N>
inline void staticSquareMatrix<N>::transpose() {
  int i, j;
  double tmp;

  for(i=0; i < N; ++i) {
    for (j=i+1; j < N; ++j) {
      tmp          = matrix[i][j];
      matrix[i][j] = matrix[j][i];
      matrix[j][i] = tmp;
    }
  }
}

template <int N>
inline void staticSquareMatrix<N>::add(double x) {
  int i, j;

  for(i=0; i < N; ++i) {
    for (j=0; j < N; ++j) {
      matrix[i][j] += x;
    }
  }
}

template <int N>
inline void staticSquareMatrix<N>::add(const staticSquareMatrix &x) {
  int i, j;

  for(i=0; i < N; ++i) {
    for (j=0; j < N; ++j) {
      matrix[i][j] += x.matrix[i][j];
    }
  }
}

template <int N>
inline void staticSquareMatrix<N>::add_toDiagonal(double x) {
  int i;

  for(i=0; i < N; ++i) {
      matrix[i][i] += x;
  }
}

template <int N>
inline void staticSquareMatrix<N>::sub(double x) {
  int i, j;

  for(i=0; i < N; ++i) {
    for (j=0; j < N; ++j) {
      matrix[i][j] -= x;
    }
  }
}

template <int N>
inline void staticSquareMatrix<N>::sub(const staticSquareMatrix &x) {
  int i, j;

  for(i=0; i < N; ++i) {
    for (j=0; j < N; ++j) {
      matrix[i][j] -= x.matrix[i][j];
    }
  }
}

template <int N>
inline void staticSquareMatrix<N>::sub_fromDiagonal(double x) {
  int i;

  for(i=0; i < N; ++i) {
      matrix[i][i] -= x;
  }
}

template <int N>
inline void staticSquareMatrix<N>::mult(double x) {
  int i, j;

  for(i=0; i < N; ++i) {
    for (j=0; j < N; ++j) {
      matrix[i][j] *= x;
    }
  }
}

template <int N>
void staticSquareMatrix<N>::setToProductOf(const staticSquareMatrix& m1,
					   const staticSquareMatrix& m2) {
  int    i, j, k;
  double tmp;

  for (i = 0; i < N; ++i) {
    for (j = 0; j < N; ++j) {
      tmp = 0; 
      for (k = 0; k < N; ++k) {
	tmp += m1.matrix[i][k] * m2.matrix[k][j];
      }
      //      if(tmp < 0.000000000001) {
      //	tmp = 0.0;
      //      }
      matrix[i][j] = tmp;
    }
  }
}

template <int N>
inline staticSquareMatrix<N>& staticSquareMatrix<N>::operator+=(double x) {
  int i, j;

  for(i=0; i < N; ++i) {
    for (j=0; j < N; ++j) {
      matrix[i][j] += x;
    }
  }
}

template <int N>
inline staticSquareMatrix<N>& staticSquareMatrix<N>::operator-=(double x) {
  int i, j;

  for(i=0; i < N; ++i) {
    for (j=0; j < N; ++j) {
      matrix[i][j] -= x;
    }
  }
}

template <int N>
inline staticSquareMatrix<N>& staticSquareMatrix<N>::operator*=(double x) {
  int i, j;

  for(i=0; i < N; ++i) {
    for (j=0; j < N; ++j) {
      matrix[i][j] *= x;
    }
  }
  return *this;
}

template <int N>
inline staticSquareMatrix<N>& staticSquareMatrix<N>::operator+=(const staticSquareMatrix &x) {
  int i, j;

  for(i=0; i < N; ++i) {
    for (j=0; j < N; ++j) {
      matrix[i][j] += x.matrix[i][j];
    }
  }
  return *this;
}

template <int N>
inline staticSquareMatrix<N>& staticSquareMatrix<N>::operator-=(const staticSquareMatrix &x) {
  int i, j;

  for(i=0; i < N; ++i) {
    for (j=0; j < N; ++j) {
      matrix[i][j] -= x.matrix[i][j];
    }
  }
  return *this;
}

template <int N>
inline staticSquareMatrix<N>& staticSquareMatrix<N>::operator*=(const staticSquareMatrix &x) {
  staticSquareMatrix<N> old = *this;
  setToProductOf(old,x);
  return *this;
}

template <int N>
inline double staticSquareMatrix<N>::absmax() {
  int i, j;

  if (N == 0)
  {
    return 0;
  }

  double m   = fabs(matrix[0][0]);
  double tmp;
 
  for(i=0; i < N; ++i) {
    for (j=0; j < N; ++j) {
      tmp = fabs(matrix[i][j]);
      if (tmp > m)
	m = tmp;
    }
  }
  return m;
}

template <int N>
void staticSquareMatrix<N>::EigenVectorsValues_JCM(staticSquareMatrix<N> &e, staticVector<N> &v)
{
  staticSquareMatrix<N> tmp = *this;
  Jacobi_Cyclic_Method(v.bvector, (double*)e.matrix, (double*)tmp.matrix);
}

template <int N>
void staticSquareMatrix<N>::setToExp_series(const staticSquareMatrix<N> &A,
					    double t, int k, double eps)
{
  // 1+At+1/2AAtt...
  staticSquareMatrix<N> tmp(A);
  int                   i;

  tmp *= t;
  assign_diagonal(1);
  add(tmp);
  for (i=2; i<k; ++i)
  {
    tmp *= A;
    tmp *= (t/i);
    add(tmp);
/*     std::cout << "for tmp" << i << std::endl; */
/*     tmp.print(std::cout, 15); */
  }

  while (tmp.absmax() > eps)
  {
    tmp *= A;
    tmp *= (t/i);
    add(tmp);
    //   std::cout << "while tmp" << i << std::endl;
    //   tmp.print(std::cout, 15);
    ++i;
  }
}


template<int N>
double staticSquareMatrix<N>::determinant()
{
  if (N < 1)
  {
    return 0;
  }
  if (N == 1)
  {
    return matrix[0][0];
  }
  if (N==2)
  {
    return matrix[0][0]*matrix[1][1]-matrix[0][1]*matrix[1][0];
  }
  if (N==3)
  {
    return  matrix[0][0]*matrix[1][1]*matrix[2][2]
           +matrix[0][1]*matrix[1][2]*matrix[2][0]
           +matrix[0][2]*matrix[1][0]*matrix[2][1]
           -matrix[0][0]*matrix[1][2]*matrix[2][1]
           -matrix[0][1]*matrix[1][0]*matrix[2][2]
           -matrix[0][2]*matrix[1][1]*matrix[2][0];
  }
  if (N==4)
  {

  }
}

template<int N>
double staticSquareMatrix<N>::trace()
{
  double t=0;
  int    i;

  for (i=0; i < N; ++i)
  {
    t += matrix[i][i];
  }
  return t;
}

template<int N>
staticVector<N> staticSquareMatrix<N>::mult_with_vector(const staticVector<N>&v)
{
  staticVector<N> a;
  int     i, j;
  double  tmp;

  for(i=0; i<N; ++i) { //rows
    tmp = 0;
    for(j=0; j<N; ++j) { //cols
      tmp += matrix[i][j]*v.bvector[j];
    }
    a.bvector[i] = tmp;
  }
  return a;
}

template<int N>
bool staticSquareMatrix<N>::orthogonal()
{
  staticSquareMatrix<N> trans(*this);
  trans.transpose();
  staticSquareMatrix<N> res;
  res.setToProductOf(*this, trans);
  res.sub_fromDiagonal(1);
  return (res.absmax() < 0.0000001);
}

template<int N>
bool staticSquareMatrix<N>::symmetric()
{
  int i, j;

  for(i=0; i < N; ++i) {
    for (j=i+1; j < N; ++j) {
      if (fabs(matrix[i][j] - matrix[j][i]) > 0.0000001)
	return false;
    }
  }
  return true;
}

template<int N>
bool staticSquareMatrix<N>::antisymmetric()
{
  int i, j;

  for(i=0; i < N; ++i) {
    for (j=i+1; j < N; ++j) {
      if (fabs(matrix[i][j] + matrix[j][i]) > 0.0000001)
	return false;
    }
  }
  return true;
}



////////////////////////////////////////////////////////////////////////////////
// File: jacobi_cyclic_method.c                                               //
// Routines:                                                                  //
//    Jacobi_Cyclic_Method                                                    //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  void Jacobi_Cyclic_Method                                                 //
//            (double eigenvalues[], double *eigenvectors, double *A)  //
//                                                                            //
//  Description:                                                              //
//     Find the eigenvalues and eigenvectors of a symmetric N x N matrix A    //
//     using the Jacobi method. Upon return, the input matrix A will have     //
//     been modified.                                                         //
//     The Jacobi procedure for finding the eigenvalues and eigenvectors of a //
//     symmetric matrix A is based on finding a similarity transformation     //
//     which diagonalizes A.  The similarity transformation is given by a     //
//     product of a sequence of orthogonal (rotation) matrices each of which  //
//     annihilates an off-diagonal element and its transpose.  The rotation   //
//     effects only the rows and columns containing the off-diagonal element  //
//     and its transpose, i.e. if a[i][j] is an off-diagonal element, then    //
//     the orthogonal transformation rotates rows a[i][] and a[j][], and      //
//     equivalently it rotates columns a[][i] and a[][j], so that a[i][j] = 0 //
//     and a[j][i] = 0.                                                       //
//     The cyclic Jacobi method considers the off-diagonal elements in the    //
//     following order: (0,1),(0,2),...,(0,N-1),(1,2),...,(N-2,N-1).  If the  //
//     the magnitude of the off-diagonal element is greater than a treshold,  //
//     then a rotation is performed to annihilate that off-diagnonal element. //
//     The process described above is called a sweep.  After a sweep has been //
//     completed, the threshold is lowered and another sweep is performed     //
//     with the new threshold. This process is completed until the final      //
//     sweep is performed with the final threshold.                           //
//     The orthogonal transformation which annihilates the matrix element     //
//     a[k][m], k != m, is Q = q[i][j], where q[i][j] = 0 if i != j, i,j != k //
//     i,j != m and q[i][j] = 1 if i = j, i,j != k, i,j != m, q[k][k] =       //
//     q[m][m] = cos(phi), q[k][m] = -sin(phi), and q[m][k] = sin(phi), where //
//     the angle phi is determined by requiring a[k][m] -> 0.  This condition //
//     on the angle phi is equivalent to                                      //
//               cot(2 phi) = 0.5 * (a[k][k] - a[m][m]) / a[k][m]             //
//     Since tan(2 phi) = 2 tan(phi) / (1.0 - tan(phi)^2),                    //
//               tan(phi)^2 + 2cot(2 phi) * tan(phi) - 1 = 0.                 //
//     Solving for tan(phi), choosing the solution with smallest magnitude,   //
//       tan(phi) = - cot(2 phi) + sgn(cot(2 phi)) sqrt(cot(2phi)^2 + 1).     //
//     Then cos(phi)^2 = 1 / (1 + tan(phi)^2) and sin(phi)^2 = 1 - cos(phi)^2 //
//     Finally by taking the sqrts and assigning the sign to the sin the same //
//     as that of the tan, the orthogonal transformation Q is determined.     //
//     Let A" be the matrix obtained from the matrix A by applying the        //
//     similarity transformation Q, since Q is orthogonal, A" = Q'AQ, where Q'//
//     is the transpose of Q (which is the same as the inverse of Q).  Then   //
//         a"[i][j] = Q'[i][p] a[p][q] Q[q][j] = Q[p][i] a[p][q] Q[q][j],     //
//     where repeated indices are summed over.                                //
//     If i is not equal to either k or m, then Q[i][j] is the Kronecker      //
//     delta.   So if both i and j are not equal to either k or m,            //
//                                a"[i][j] = a[i][j].                         //
//     If i = k, j = k,                                                       //
//        a"[k][k] =                                                          //
//           a[k][k]*cos(phi)^2 + a[k][m]*sin(2 phi) + a[m][m]*sin(phi)^2     //
//     If i = k, j = m,                                                       //
//        a"[k][m] = a"[m][k] = 0 =                                           //
//           a[k][m]*cos(2 phi) + 0.5 * (a[m][m] - a[k][k])*sin(2 phi)        //
//     If i = k, j != k or m,                                                 //
//        a"[k][j] = a"[j][k] = a[k][j] * cos(phi) + a[m][j] * sin(phi)       //
//     If i = m, j = k, a"[m][k] = 0                                          //
//     If i = m, j = m,                                                       //
//        a"[m][m] =                                                          //
//           a[m][m]*cos(phi)^2 - a[k][m]*sin(2 phi) + a[k][k]*sin(phi)^2     //
//     If i= m, j != k or m,                                                  //
//        a"[m][j] = a"[j][m] = a[m][j] * cos(phi) - a[k][j] * sin(phi)       //
//                                                                            //
//     If X is the matrix of normalized eigenvectors stored so that the ith   //
//     column corresponds to the ith eigenvalue, then AX = X Lamda, where     //
//     Lambda is the diagonal matrix with the ith eigenvalue stored at        //
//     Lambda[i][i], i.e. X'AX = Lambda and X is orthogonal, the eigenvectors //
//     are normalized and orthogonal.  So, X = Q1 Q2 ... Qs, where Qi is      //
//     the ith orthogonal matrix,  i.e. X can be recursively approximated by  //
//     the recursion relation X" = X Q, where Q is the orthogonal matrix and  //
//     the initial estimate for X is the identity matrix.                     //
//     If j = k, then x"[i][k] = x[i][k] * cos(phi) + x[i][m] * sin(phi),     //
//     if j = m, then x"[i][m] = x[i][m] * cos(phi) - x[i][k] * sin(phi), and //
//     if j != k and j != m, then x"[i][j] = x[i][j].                         //
//                                                                            //
//  Arguments:                                                                //
//     double  eigenvalues                                                    //
//        Array of dimension N, which upon return contains the eigenvalues of //
//        the matrix A.                                                       //
//     double* eigenvectors                                                   //
//        Matrix of eigenvectors, the ith column of which contains an         //
//        eigenvector corresponding to the ith eigenvalue in the array        //
//        eigenvalues.                                                        //
//     double* A                                                              //
//        Pointer to the first element of the symmetric N x N matrix A. The   //
//        input matrix A is modified during the process.                      //
//                                                                            //
//  Return Values:                                                            //
//     Function is of type void.                                              //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double A[N][N], double eigenvalues[N], double eigenvectors[N][N]       //
//                                                                            //
//     (your code to initialize the matrix A )                                //
//                                                                            //
//     Jacobi_Cyclic_Method(eigenvalues, (double*)eigenvectors,               //
//                                                          (double *) A, N); //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //

template <int N>
void staticSquareMatrix<N>::Jacobi_Cyclic_Method(double eigenvalues[], double *eigenvectors, double *A){
  int i, j, k, m;
  double *pAk, *pAm, *p_r, *p_e;
  double threshold_norm;
  double threshold;
  double tan_phi, sin_phi, cos_phi, tan2_phi, sin2_phi, cos2_phi;
  double sin_2phi; // , cos_2phi, 
  double cot_2phi;
  double dum1;
  double dum2;
  double dum3;
  //	double r;
  double max;

  // Take care of trivial cases

  if ( N < 1) return;
  if ( N == 1) {
    eigenvalues[0] = *A;
    *eigenvectors = 1.0;
    return;
  }

  // Initialize the eigenvalues to the identity matrix.

  for (p_e = eigenvectors, i = 0; i < N; i++)
    for (j = 0; j < N; p_e++, j++)
      if (i == j) *p_e = 1.0; else *p_e = 0.0;

  // Calculate the threshold and threshold_norm.

  for (threshold = 0.0, pAk = A, i = 0; i < ( N - 1 ); pAk += N, i++) 
    for (j = i + 1; j < N; j++) threshold += *(pAk + j) * *(pAk + j);
  threshold = sqrt(threshold + threshold);
  threshold_norm = threshold * DBL_EPSILON;
  max = threshold + 1.0;
  while (threshold > threshold_norm) {
    threshold /= 10.0;
    if (max < threshold) continue;
    max = 0.0;
    for (pAk = A, k = 0; k < (N-1); pAk += N, k++) {
      for (pAm = pAk + N, m = k + 1; m < N; pAm += N, m++) {
	if ( fabs(*(pAk + m)) < threshold ) continue;

	// Calculate the sin and cos of the rotation angle which
	// annihilates A[k][m].
	
	cot_2phi = 0.5 * ( *(pAk + k) - *(pAm + m) ) / *(pAk + m);
	dum1 = sqrt( cot_2phi * cot_2phi + 1.0);
	if (cot_2phi < 0.0) dum1 = -dum1;
	tan_phi = -cot_2phi + dum1;
	tan2_phi = tan_phi * tan_phi;
	sin2_phi = tan2_phi / (1.0 + tan2_phi);
	cos2_phi = 1.0 - sin2_phi;
	sin_phi = sqrt(sin2_phi);
	if (tan_phi < 0.0) sin_phi = - sin_phi;
	cos_phi = sqrt(cos2_phi); 
	sin_2phi = 2.0 * sin_phi * cos_phi;
	//	cos_2phi = cos2_phi - sin2_phi;
	
	// Rotate columns k and m for both the matrix A 
	//     and the matrix of eigenvectors.
	
	p_r = A;
	dum1 = *(pAk + k);
	dum2 = *(pAm + m);
	dum3 = *(pAk + m);
	*(pAk + k) = dum1 * cos2_phi + dum2 * sin2_phi + dum3 * sin_2phi;
	*(pAm + m) = dum1 * sin2_phi + dum2 * cos2_phi - dum3 * sin_2phi;
	*(pAk + m) = 0.0;
	*(pAm + k) = 0.0;
	for (i = 0; i < N; p_r += N, i++) {
	  if ( (i == k) || (i == m) ) continue;
	  if ( i < k ) dum1 = *(p_r + k); else dum1 = *(pAk + i);
	  if ( i < m ) dum2 = *(p_r + m); else dum2 = *(pAm + i);
	  dum3 = dum1 * cos_phi + dum2 * sin_phi;
	  if ( i < k ) *(p_r + k) = dum3; else *(pAk + i) = dum3;
	  dum3 = - dum1 * sin_phi + dum2 * cos_phi;
	  if ( i < m ) *(p_r + m) = dum3; else *(pAm + i) = dum3;
	}
	for (p_e = eigenvectors, i = 0; i < N; p_e += N, i++) {
	  dum1 = *(p_e + k);
	  dum2 = *(p_e + m);
	  *(p_e + k) = dum1 * cos_phi + dum2 * sin_phi;
	  *(p_e + m) = - dum1 * sin_phi + dum2 * cos_phi;
	}
      }
      for (i = 0; i < N; i++)
	if ( i == k ) continue;
	else if ( max < fabs(*(pAk + i))) max = fabs(*(pAk + i));
    }
  }
  for (pAk = A, k = 0; k < N; pAk += N, k++) eigenvalues[k] = *(pAk + k); 
}



//*************************************************************

#endif
