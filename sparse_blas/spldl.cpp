//
// Created by kazem on 6/12/20.
//

#include <algorithm>
#include <cassert>
#include <sparse_inspector.h>
#include <omp.h>

namespace sym_lib{

 int ldl_left_simplicial_02(int n, int *c, int *r, double *values,
                            int *cT, int *rT,
                            int *lC, int *lR, double *&lValues,
                            double *d,
#ifdef PPRUNE
   int *prunePtr, int *pruneSet,
#endif
                            int *eTree,
                            double *ws, int *ws_int) {
  int top = 0;
  double *f = ws;
  int *finger = ws_int;
  int *xi = ws_int + n;
  std::fill_n(f, n, 0);
  std::fill_n(xi, 2 * n, 0);
  //Determining the column pointer
  for (int k = 0; k < n; ++k) {
   finger[k] = lC[k];
  }
  int spCol = 0;
  for (int colNo = 0; colNo < n; ++colNo) {
   //emptying col colNo in L
   //Uncompress a col into a 1D array
   for (int nzNo = c[colNo]; nzNo < c[colNo + 1]; ++nzNo) {
    f[r[nzNo]] = values[nzNo];//Copying nonzero of the col
   }
#ifdef PPRUNE
   for (int i = prunePtr[colNo]; i < prunePtr[colNo + 1]-1; ++i) {
    spCol = pruneSet[i];
#endif
   top = ereach(n, cT, rT, colNo, eTree, xi, xi + n);
   //std::cout<<n-top<<";\n";
   for (int i = top; i < n; ++i) {
    spCol = xi[i];
    bool sw = false;
    double facing_val = 0, tmp = 0;
    int facing_idx = -1;
    for (int l = lC[spCol]; l < lC[spCol + 1]; ++l) {
     if (lR[l] == colNo) {
      facing_val = lValues[l];
      tmp = facing_val * d[spCol];
      facing_idx = l;
      break;
     }
    }
    assert(facing_idx >= 0);
    for (int l = facing_idx + 1; l < lC[spCol + 1]; ++l) {
     f[lR[l]] -= lValues[l] * tmp;
    }
    d[colNo] += facing_val * tmp;
   }
   d[colNo] = f[colNo] - d[colNo];
   double diag = d[colNo];
   if(!diag)
    return 0;
   f[colNo] = 0;
   lValues[lC[colNo]] = 1;
   for (int j = lC[colNo] + 1; j < lC[colNo + 1]; ++j) {
    lValues[j] = f[lR[j]] / diag;
    f[lR[j]] = 0;
   }
  }
  return 1;
 } //end LDL


 int ldl_parallel_left_simplicial_01 (int n, int* c, int* r, double* values,
                                      int* cT, int* rT,
                                      int* lC, int* lR, double* &lValues,
                                      double *d,
#if 0
   int *prunePtr, int *pruneSet,
#endif
                                      int *eTree,
                                      int nLevels, int *levelPtr,
                                      int *parPtr, int *partition) {
//omp_set_num_threads(1);

  omp_set_nested(1);
  for (int i1 = 0; i1 < nLevels; ++i1) {
#pragma omp parallel
   {
#pragma omp  for schedule(static)
    for (int j1 = levelPtr[i1]; j1 < levelPtr[i1 + 1]; ++j1) {
     auto *f = new double[n]();
     auto *xi = new int[2 * n]();
     //int pls = levelSet[j1];
     for (int k1 = parPtr[j1]; k1 < parPtr[j1 + 1]; ++k1) {
      int colNo = partition[k1];
      int spCol = 0;
      //Uncompress a col into a 1D array
      for (int nzNo = c[colNo]; nzNo < c[colNo + 1]; ++nzNo) {
       f[r[nzNo]] = values[nzNo];//Copying nonzero of the col
      }
#if 0
      for (int i = prunePtr[colNo]; i < prunePtr[colNo + 1]-1; ++i) {
   spCol = pruneSet[i];
#endif
      int top = ereach(n, cT, rT, colNo, eTree, xi, xi + n);
      //std::cout<<n-top<<";\n";
      for (int i = top; i < n; ++i) {
       spCol = xi[i];
       bool sw = false;
       double facing_val = 0, tmp = 0;
       int facing_idx = -1;
       for (int l = lC[spCol]; l < lC[spCol + 1]; ++l) {
        if (lR[l] == colNo) {
         facing_val = lValues[l];
         tmp = facing_val * d[spCol];
         facing_idx = l;
         break;
        }
       }
       assert(facing_idx >= 0);
       for (int l = facing_idx + 1; l < lC[spCol + 1]; ++l) {
        f[lR[l]] -= lValues[l] * tmp;
       }
       d[colNo] += facing_val * tmp;
      }
      d[colNo] = f[colNo] - d[colNo];
      double diag = d[colNo];
      f[colNo] = 0;
      lValues[lC[colNo]] = 1;
      for (int j = lC[colNo] + 1; j < lC[colNo + 1]; ++j) {
       lValues[j] = f[lR[j]] / diag;
       f[lR[j]] = 0;
      }
     }
     delete[]f;
     delete[]xi;
    }
   }
  }
  return 1;
 }


}