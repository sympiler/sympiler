//
// Created by george on 2019-10-09.
//
#include <omp.h>

#include "sparse_blas_lib.h"

namespace sym_lib {

 void sptrsv_csr(int n, int *Lp, int *Li, double *Lx, double *x) {
  int i, j;
  for (i = 0; i < n; i++) {
   for (j = Lp[i]; j < Lp[i + 1] - 1; j++) {
    x[i] -= Lx[j] * x[Li[j]];
   }
   x[i] /= Lx[Lp[i + 1] - 1];
  }
 }



 void sptrsv_csr_levelset(int n, const int *Lp, const int *Li, const double *Lx,
                          double *x,
                          int levels, const int *levelPtr,
                          const int *levelSet) {
  for (int l = 0; l < levels; l++) {
#pragma omp parallel for default(shared) schedule(auto)
   for (int k = levelPtr[l]; k < levelPtr[l + 1]; ++k) {
    int i = levelSet[k];
    for (int j = Lp[i]; j < Lp[i + 1] - 1; j++) {
     x[i] -= Lx[j] * x[Li[j]];
    }
    x[i] /= Lx[Lp[i + 1] - 1];
   }
  }
 }


 void sptrsv_csr_lbc(int n, int *Lp, int *Li, double *Lx, double *x,
                     int level_no, int *level_ptr,
                     int *par_ptr, int *partition) {
  for (int i1 = 0; i1 < level_no; ++i1) {
#pragma omp parallel
   {
#pragma omp  for schedule(auto)
    for (int j1 = level_ptr[i1]; j1 < level_ptr[i1 + 1]; ++j1) {
     for (int k1 = par_ptr[j1]; k1 < par_ptr[j1 + 1]; ++k1) {
      int i = partition[k1];
      for (int j = Lp[i]; j < Lp[i + 1] - 1; j++) {
       x[i] -= Lx[j] * x[Li[j]];
      }
      x[i] /= Lx[Lp[i + 1] - 1];
     }
    }
   }
  }
 }




}