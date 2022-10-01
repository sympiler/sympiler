//
// Created by kazem on 10/17/21.
//

#include <aggregation/sparse_inspector.h>
#include <valarray>

namespace sym_lib {
 /*
 * Left looking simplicial cholesky with reach function coupled
 */
 bool cholesky_left_serial(const int n, const int* c, const int* r, const
 double* values, const int* cT, const int* rT, const int* lC, const int*
 lR, double* &lValues, const int *eTree) {
  /*
   * Performs a Cholesky decomposition on a given matrix (c, r, values), i.e.
   * stored in compressed column format, and produces L, which are
   * stored in column compressed format. The column sparsity pattern of L
   * is provided in lC
   * The row sparsity pattern of L is provided in lR
   * No ordering is performed.
   */
  double *f = new double[n]();
  int top=n;
  int *xi = new int[2*n]();
  int *rowIdx = new int[n];
  //Determining the column pointer
  for (int k = 0; k < n; ++k) {
   rowIdx[k]=lC[k];
  }
  int spCol=0;
  for (int colNo = 0; colNo < n; ++colNo) {
   //Uncompressing a col into a 1D array
   for (int nzNo = c[colNo]; nzNo < c[colNo + 1]; ++nzNo) {
    // if(r[nzNo]>=colNo)
    f[r[nzNo]] = values[nzNo];//Copying nonzero of the col
   }
   top = ereach(n,cT,rT,colNo,eTree,xi,xi+n);
   for (int i = top; i < n; ++i) {
    spCol=xi[i];
    for (int l = lC[spCol]; l < lC[spCol+1] ; ++l) {
     if(lR[l]>=colNo){
      f[lR[l]] -= lValues[l]*lValues[rowIdx[spCol]];
     }
    }
    rowIdx[spCol]++;
   }
   rowIdx[colNo]++;//Skip diagonal row
   if(f[colNo]<=0)
    return false; //The matrix is not SPD
   double tmpSqrt=sqrt(f[colNo]);
   f[colNo]=0;
   lValues[lC[colNo]] = tmpSqrt;
   for (int j = lC[colNo]+1; j < lC[colNo+1]; ++j) {
    lValues[j] = f[lR[j]]/tmpSqrt;
    f[lR[j]]=0;
   }
  }
  delete []rowIdx;
  delete []f;
  delete []xi;
  return true;
 }

}
