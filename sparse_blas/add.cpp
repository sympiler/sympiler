//
// Created by kazem on 2/9/20.
//

#include <def.h>
#include <sparse_utilities.h>
#include "sparse_blas_lib.h"

namespace sym_lib {

 CSC* add(CSC *A, CSC *B, double alpha, double beta, bool sort){
  if(A->stype != B->stype)
   return NULLPNTR;
  int p, j, nz = 0, anz, *Cp, *Ci, *Bp, m, n, bnz, *w ;
  bool values = A->is_pattern || B->is_pattern; //true if either is pattern
  double *x, *Bx, *Cx ;
  m = A->m ; anz = A->p [A->n] ;
  n = B->n ; Bp = B->p ; Bx = B->x ; bnz = Bp [n] ;
  w = new int[m]();
  x = !values ? new double[m]() : NULLPNTR;
  CSC *C = new CSC(m,n,anz+bnz,values);
  C->stype = A->stype;
  C->is_pattern = A->is_pattern;
  Cp = C->p ; Ci = C->i ; Cx = C->x ;
  for (j = 0 ; j < n ; j++){
   Cp [j] = nz ;                   /* column j of C starts here */
   nz = scatter (A, j, alpha, w, x, j+1, C, nz) ;   /* alpha*A(:,j)*/
   nz = scatter (B, j, beta, w, x, j+1, C, nz) ;    /* beta*B(:,j) */
   if(sort) std::sort(&C->i[Cp[j]], &C->i[nz]);
   if (!values) for (p = Cp [j] ; p < nz ; p++) Cx [p] = x [Ci [p]] ;
  }
  Cp [n] = nz ;
  delete []x;
  delete []w;
  return C;
 }

 void csc_add_diag(int n, const int *Ap, const int *Ai, const double *Ax,
   const double* d, int *Bp, int *Bi, double *Bx){
  for (int i=0; i<n; i++){
   Bp[i] = Ap[i];
   for(int j=Ap[i]; j<=Ap[i+1]; j++){
    Bi[j] = Ai[j];
    if( i == Ai[j] )
     Bx[j] = (d[i]+Ax[j]);
    else
     Bx[j] = Ax[j];
   }}}


 void csr_add_diag(int n, const int *Ap, const int *Ai, const double *Ax,
                   const double* d, int *Bp, int *Bi, double *Bx){
  for (int i=0; i<n; i++){
   Bp[i] = Ap[i];
   for(int j=Ap[i]; j<=Ap[i+1]; j++){
    Bi[j] = Ai[j];
    if(i == Ai[j])
     Bx[j] = (d[i]+Ax[j]);
    else
     Bx[j] = Ax[j];
   }}}

}