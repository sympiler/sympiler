//
// Created by kazem on 2020-05-20.
//
#include <omp.h>
#include "def.h"

namespace sym_lib{

 void spilu0_csr(int n, int nnz, const int *Ap, const int *Ai,
   const double *Ax, const int *A_diag, double *l){
  int *iw = new int[n];
  for (int j = 0; j < n; j++ ){
   iw[j] = EMPTY;
  }
  int jrow=0, jw;
  double tl;
  // copying A
  for (int k = 0; k < nnz; k++ ){
   l[k] = Ax[k];
  }
  for (int i = 0; i < n; i++ ){
  // iw points to the nonzero entries in row i.
   for (int k = Ap[i]; k < Ap[i+1]; k++ ){
    iw[Ai[k]] = k;
   }
   int j;
   for (j = Ap[i]; j < Ap[i+1]; ++j) {
    jrow = Ai[j];
    if ( i <= jrow ){
     break;
    }
    tl = l[j] / l[A_diag[jrow]];
    l[j] = tl;
    for (int jj = A_diag[jrow] + 1; jj < Ap[jrow+1]; jj++ ){
     jw = iw[Ai[jj]];
     if ( jw != EMPTY ){
      l[jw] = l[jw] - tl * l[jj];
     }
    }
   }
   assert(A_diag[i] == j);
   assert(jrow == i);//if not, it mean it does not have a diagonal value
   assert(l[j] != 0.0); //zero diagonal check
   for (int k = Ap[i]; k < Ap[i+1]; k++ ){
    iw[Ai[k]] = EMPTY;
   }
  }
  delete []iw;
 }


 void spilu0_csr_lbc(int n, int nnz, const int *Ap, const int *Ai,
                     const double *Ax, const int *A_diag, double *l,
                     int level_no, const int *level_ptr,
                     const int *par_ptr, const int *partition,
                     int *tempvecs) {
  // copying A
  for (int k = 0; k < nnz; k++ ){
   l[k] = Ax[k];
  }
  for (int i1 = 0; i1 < level_no; ++i1) {
#pragma omp parallel
   {
    int id = omp_get_thread_num();
    int *iw = tempvecs+(id*n);
#pragma omp  for schedule(auto)
    for (int j1 = level_ptr[i1]; j1 < level_ptr[i1 + 1]; ++j1) {
     for (int k1 = par_ptr[j1]; k1 < par_ptr[j1 + 1]; ++k1) {
      int i = partition[k1];

      for (int k = Ap[i]; k < Ap[i+1]; k++ ){
       iw[Ai[k]] = k;
      }
      int j, jrow, jw;
      for (j = Ap[i]; j < Ap[i+1]; ++j) {
       jrow = Ai[j];
       if ( i <= jrow ){
        break;
       }
       double tl = l[j] / l[A_diag[jrow]];
       l[j] = tl;
       for (int jj = A_diag[jrow] + 1; jj < Ap[jrow+1]; jj++ ){
        jw = iw[Ai[jj]];
        if ( jw != EMPTY ){
         l[jw] = l[jw] - tl * l[jj];
        }
       }
      }
      assert(A_diag[i] == j);
      assert(jrow == i);//if not, it mean it does not have a diagonal value
      assert(l[j] != 0.0); //zero diagonal check
      for (int k = Ap[i]; k < Ap[i+1]; k++ ){
       iw[Ai[k]] = EMPTY;
      }
     }
    }
   //delete []iw;
  }
  }
 }


}