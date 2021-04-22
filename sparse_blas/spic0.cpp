//
// Created by kazem on 2020-04-23.
//
#include <cmath>
#include <cassert>
#include <iostream>

namespace sym_lib {


 bool spic0_csc_left(int n, const int* Ap, const int* Ai,
   const double* Ax, const int *prunePtr, const int *pruneSet,
   double* Lx) {
  double *f = new double[n]();
  for (int colNo = 0; colNo < n; ++colNo) {
   //Uncompressing a col into a 1D array
   for (int nzNo = Ap[colNo]; nzNo < Ap[colNo + 1]; ++nzNo) {
    f[Ai[nzNo]] = Ax[nzNo];//Copying nonzero of the col
   }
   for (int i = prunePtr[colNo]; i < prunePtr[colNo + 1]; ++i) {
    double tmp=0;
    int spCol = pruneSet[i];
    bool sw=false;
    for (int l = Ap[spCol]; l < Ap[spCol + 1]; ++l) {
     if (Ai[l] == colNo && !sw) {
      tmp = Lx[l];
      sw=true;
     }
     if(sw){
      f[Ai[l]] -= Lx[l] * tmp;
     }
    }
   }
   if (f[colNo] <= 0)
    return false; //The matrix is not SPD
   double tmpSqrt = sqrt(f[colNo]);
   f[colNo] = 0;
   Lx[Ap[colNo]] = tmpSqrt;
   for (int j = Ap[colNo] + 1; j < Ap[colNo + 1]; ++j) {
    Lx[j] = f[Ai[j]] / tmpSqrt;
    f[Ai[j]] = 0;
   }
  }
  delete[]f;
  return true;
 }

 bool spic0_csc_left_LBC(int n, const int* Ap, const int* Ai,
                     const double* Ax, const int *prunePtr,
                     const int *pruneSet,const int nLevels,
                     const int *levelPtr,
                     const int *parPtr, const int *partition,
                     double* Lx) {
  bool ret = true;
  for (int i1 = 0; i1 < nLevels; ++i1) {
#pragma omp parallel
   {
#pragma omp  for schedule(static)
    for (int j1 = levelPtr[i1]; j1 < levelPtr[i1 + 1]; ++j1) {
     auto *f = new double[n]();
     for (int k1 = parPtr[j1]; k1 < parPtr[j1 + 1]; ++k1) {
      int colNo = partition[k1];

      //double *f = new double[n]();
      //for (int colNo = 0; colNo < n; ++colNo) {
      //Uncompressing a col into a 1D array
      for (int nzNo = Ap[colNo]; nzNo < Ap[colNo + 1]; ++nzNo) {
       f[Ai[nzNo]] = Ax[nzNo];//Copying nonzero of the col
      }
      for (int i = prunePtr[colNo]; i < prunePtr[colNo + 1]; ++i) {
       double tmp = 0;
       int spCol = pruneSet[i];
       bool sw = false;
       for (int l = Ap[spCol]; l < Ap[spCol + 1]; ++l) {
        if (Ai[l] == colNo && !sw) {
         tmp = Lx[l];
         sw = true;
        }
        if (sw) {
         f[Ai[l]] -= Lx[l] * tmp;
        }
       }
      }
      if (f[colNo] <= 0)
       ret = false; //The matrix is not SPD
      double tmpSqrt = sqrt(f[colNo]);
      f[colNo] = 0;
      Lx[Ap[colNo]] = tmpSqrt;
      for (int j = Ap[colNo] + 1; j < Ap[colNo + 1]; ++j) {
       Lx[j] = f[Ai[j]] / tmpSqrt;
       f[Ai[j]] = 0;
      }
     }
     delete []f;
    }
   }
  }
  return ret;
 }


 /*
  * Accepts upper triangular of A in CSR
  * or lower triangular of A in CSC
  */
 void spic0_csr(int n, double *Lx, int *Lp, int *Li)  {
  for (int i = 0; i < n; i++) {
   assert(Lx[Lp[i]] >0);
   Lx[Lp[i]] = sqrt(Lx[Lp[i]]);//S1
   double d = Lx[Lp[i]];

   for (int m = Lp[i] + 1; m < Lp[i + 1]; m++) {
    Lx[m] = Lx[m] / d;//S2
   }

   for (int m = Lp[i]+1 ; m < Lp[i + 1]; m++) {
    for (int k = Lp[Li[m]]; k < Lp[Li[m] + 1]; k++) {
     for (int l = m; l < Lp[i + 1] && Li[l+1] <= Li[k]; l++) {
      //if (Li[l] == Li[k] && Li[l + 1] <= Li[k]) {
      if (Li[l] == Li[k] ) {
       Lx[k] -= Lx[m] * Lx[l]; //S3
      }
     }
    }
   }
  }
 }


 void spic0_csr_lbc(int n, double *Lx, int *Lp, int *Li,
                    int level_no, int *level_ptr,
                    int *par_ptr, int *partition) {
  for (int i1 = 0; i1 < level_no; ++i1) {
#pragma omp parallel
   {
#pragma omp  for schedule(auto)
    for (int j1 = level_ptr[i1]; j1 < level_ptr[i1 + 1]; ++j1) {
     for (int k1 = par_ptr[j1]; k1 < par_ptr[j1 + 1]; ++k1) {
      int i = partition[k1];
      Lx[Lp[i]] = sqrt(Lx[Lp[i]]);//S1
      double d = Lx[Lp[i]];

      for (int m = Lp[i] + 1; m < Lp[i + 1]; m++) {
       Lx[m] = Lx[m] / d;//S2
      }

      for (int m = Lp[i]+1 ; m < Lp[i + 1]; m++) {
       for (int k = Lp[Li[m]]; k < Lp[Li[m] + 1]; k++) {
        for (int l = m; l < Lp[i + 1] && Li[l+1] <= Li[k]; l++) {
         //if (Li[l] == Li[k] && Li[l + 1] <= Li[k]) {
         if (Li[l] == Li[k] ) {
#pragma omp atomic
          Lx[k] -= Lx[m] * Lx[l]; //S3
         }

        }
       }
      }
     }}}}//LBC outermost
 }

}