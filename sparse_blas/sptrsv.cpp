//
// Created by george on 2019-10-09.
//
#include <omp.h>

#include "sympiler/sparse_blas_lib.h"
#include <sympiler/dense_blas/BLAS.h>

#ifdef MKL
#include "mkl.h"
#endif

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

/// for unit diagonal L-factor from ILU0
 void sptrsv_csr_ilu(int n, const int *Lp, const int *Li, const double *Lx,
   const int*diags, double *x) {
  int i, j;
  for (i = 0; i < n; i++) {
   int dia = diags[i];
   for (j = Lp[i]; j < dia; j++) {
    x[i] -= Lx[j] * x[Li[j]];
   }
   //x[i] /= Lx[dia];
  }
 }

 void sptrsv_csr_levelset(int n, const int *Lp, const int *Li, const double *Lx,
                          double *x,
                          int levels, const int *levelPtr,
                          const int *levelSet) {
  for (int l = 0; l < levels; l++) {
#pragma omp  parallel for default(shared) schedule(auto)
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
#pragma omp  parallel
   {
#pragma omp for schedule(auto)
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


 void sptrsv_csr_lbc_redundant(int n, int *Lp, int *Li, double *Lx,
   double *x,
   double *b,
   int level_no, int *level_ptr,
   int *par_ptr, int *partition, double *tmp) {
  for (int i1 = 0; i1 < level_no; ++i1) {
#pragma omp  parallel
   {
#pragma omp for schedule(auto)
    for (int j1 = level_ptr[i1]; j1 < level_ptr[i1 + 1]; ++j1) {
#ifdef ENABLE_OPENMP
     int tid = omp_get_thread_num();
#else
     int tid=0;
#endif
     double *tt = tmp + tid*n;
     std::fill_n(tt,n,0); //TODO not sure if needed
     for (int k1 = par_ptr[j1]; k1 < par_ptr[j1 + 1]; ++k1) {
      int i = partition[k1];
      for (int j = Lp[i]; j < Lp[i + 1] - 1; j++) {
       tt[i] -= Lx[j] * x[Li[j]];
      }
      tt[i] = (tt[i] + b[i]) / Lx[Lp[i + 1] - 1];
      x[i] = tt[i];
     }
    }
   }
  }
 }


 void sptrsv_csr_ilu_lbc(int n, const int *Lp,
                         const int *Li, const double *Lx, const int*diags,
                         double *x, const int level_no, const int *level_ptr,
                         const int *par_ptr, const int *partition) {
  for (int i1 = 0; i1 < level_no; ++i1) {
#pragma omp parallel
   {
#pragma omp  for schedule(auto)
    for (int j1 = level_ptr[i1]; j1 < level_ptr[i1 + 1]; ++j1) {
     for (int k1 = par_ptr[j1]; k1 < par_ptr[j1 + 1]; ++k1) {
      int i = partition[k1];
      int dia = diags[i];
      for (int j = Lp[i]; j < dia; j++) {
       x[i] -= Lx[j] * x[Li[j]];
      }
      x[i] /= Lx[dia];
     }
    }}}}


 void sptrsv_csc(int n, int *Lp, int *Li, double *Lx, double *x) {
  int i, j;
  for (i = 0; i < n; i++) {
   x[i] /= Lx[Lp[i]];
   for (j = Lp[i] + 1; j < Lp[i + 1]; j++) {
    x[Li[j]] -= Lx[j] * x[i];
   }
  }
 }

 void
 sptrsv_csc_levelset(int n, int *Lp, int *Li, double *Lx, double *x, int levels,
                     const int *levelptr, const int *levelset) {
  for (int i = 0; i < levels; i++) {
#pragma omp parallel
   {
#pragma omp for schedule(auto)
    for (int p = levelptr[i]; p < levelptr[i + 1]; p++) {
     int node = levelset[p];
     x[node] /= Lx[Lp[node]];
     for (int j = Lp[node] + 1; j < Lp[node + 1]; j++) {
#pragma omp atomic
      x[Li[j]] -= Lx[j] * x[node];
     }
    }
   }
  }
 }

 void
 sptrsv_csc_lbc(int n, int *Lp, int *Li, double *Lx, double *x, int level_no, int *level_ptr,
                int *par_ptr, int *partition) {
  for (int i1 = 0; i1 < level_no; ++i1) {
#pragma omp parallel
   {
#pragma omp  for schedule(auto)
    for (int j1 = level_ptr[i1]; j1 < level_ptr[i1 + 1]; ++j1) {
     for (int k1 = par_ptr[j1]; k1 < par_ptr[j1 + 1]; ++k1) {
      int node = partition[k1];
     x[node] /= Lx[Lp[node]];
     for (int j = Lp[node] + 1; j < Lp[node + 1]; j++) {
      int tmp = Li[j];
#pragma omp atomic
      x[tmp] -= Lx[j] * x[node];
     }
    }}}}}

 void sptrsv_csc_block(int n, const int *Lp, const int *Li, const int *nrows,
                         double *Lx, double *x, int num_nodes,
                         const int *supernodes) {
  int i, p, k;
  auto *tempvec = new double[n]();

  for (i = 0; i < num_nodes; i++) {
   int super = supernodes[i];

   int width = supernodes[i + 1] - super;
   int nrow = nrows[i];

   lsolve_BLAS(nrow, width, &Lx[Lp[i]], &x[super]);
   matvec_BLAS(nrow, nrow - width, width, &Lx[Lp[i] + width], &x[super], tempvec);

   for (p = Lp[i] + width, k = 0; p < Lp[i] + nrow; p++, k++) {
    int idx = Li[p];
    x[idx] -= tempvec[k];
    tempvec[k] = 0;
   }
  }

  delete[]tempvec;
 }

#ifdef MKL
 void sptrsv_csc_block_mkl(int n, const int *Lp, const int *Li, const int *nrows,
                       double *Lx, double *x, int num_nodes,
                       const int *supernodes) {
  int i, p, k;
  double one[2], zero[2];
  one[0] = 1.0;
  one[1] = 0.;
  zero[0] = 0.;
  zero[1] = 0.;
  int ione = 1;

  auto *tempvec = new double[n]();

  for (i = 0; i < num_nodes; i++) {
   int super = supernodes[i];

   int width = supernodes[i + 1] - super;
   int nrow = nrows[i];

   dtrsm("L", "L", "N", "N", &width,&ione,one,&Lx[Lp[i]],
               &nrow,&x[super],&n);

   int tmpRow=nrow - width;
     dgemv("N",&tmpRow,&width,one,&Lx[Lp[i] + width],&nrow,&x[super],&ione,
               zero,tempvec,&ione);

   for (p = Lp[i] + width, k = 0; p < Lp[i] + nrow; p++, k++) {
    int idx = Li[p];
    x[idx] -= tempvec[k];
    tempvec[k] = 0;
   }
  }

  delete[]tempvec;
 }
#endif

 void sptrsv_csc_levelset_block(int n, const int *Lp, const int *Li,
                                  const int *nrows, double *Lx, double *x,
                                  const int *levelPtr,
                                  const int *levels, int n_lev,
                                  const int *supernodes, const int *sup2node,
                                  double **tempvecs) {
  int i, j, p, k;

  for (i = 0; i < n_lev; i++) {
#pragma omp parallel default(shared) private(j, k, p)
   {
#ifdef ENABLE_OPENMP
    int id = omp_get_thread_num();
#else
    int id = 0;
#endif
    double *tempvec = tempvecs[id];

#pragma omp for schedule(auto)
    for (j = levelPtr[i]; j < levelPtr[i + 1]; j++) {
     int super = levels[j];

     // sup2node[j] stores the index of current supernode in the block set
     int index = sup2node[j];
     int width = supernodes[index + 1] - super;
     int nrow = nrows[index];

     lsolve_BLAS(nrow, width, &Lx[Lp[index]], &x[super]);
     matvec_BLAS(nrow, nrow - width, width, &Lx[Lp[index] + width], &x[super],
                 tempvec);

     for (p = Lp[index] + width, k = 0; p < Lp[index] + nrow; p++, k++) {
      int idx = Li[p];
#pragma omp atomic
      x[idx] -= tempvec[k];
      tempvec[k] = 0;
     }
    }
   }
  }
 }

 void sptrsv_csc_bn_multi(int n, const int *Lp, const int *Li, const int *nrows,
   double *Lx, double *x, const int *levelPtr, const int *levels, int n_lev,
   const int *supernodes, const int *sup2node, int num_nodes, double **tempvecs)
 {
  int i, j, p, k, base;

  for (i = 0; i < n_lev; i++) {
#pragma omp parallel default(shared) private(j, k, p, base)
   {
#ifdef ENABLE_OPENMP
    int id = omp_get_thread_num();
#else
    int id = 0;
#endif
    double *tempvec = tempvecs[id];
#pragma omp for schedule(auto)
    for (j = levelPtr[i]; j < levelPtr[i + 1]; j++) {
     int super = levels[j];
     int index = sup2node[j];

     // check which layer
     // FIXME: this is not too correct...
     int kern_num = index / num_nodes;
     base = kern_num * n;
     index -= kern_num * num_nodes;

     int width = supernodes[index + 1] - super + base;
     int nrow = nrows[index];

     // copy entire block data
     if (base != 0)
      for (int m = 0; m < width; m++)
       x[super + m] += x[super - n + m];

     lsolve_BLAS(nrow, width, &Lx[Lp[index]], &x[super]);
     matvec_BLAS(nrow, nrow - width, width, &Lx[Lp[index] + width], &x[super], tempvec);

     for (p = Lp[index] + width, k = 0; p < Lp[index] + nrow; p++, k++) {
      int idx = Li[p] + base;

#pragma omp atomic
      x[idx] -= tempvec[k];
      tempvec[k] = 0;
     }
    }
   }
  }
 }


/// triangular solve a diagonal matrix with a vector
/// \param L input two dim array with size of rowNoxdNum
/// \param v input two dim array with size ofrowNo x 1
/// \param dNum input number of diagonals
/// \param rowNo input number of rows
/// \param result output two dim array with size of rowNo x 1
/// \return
    int BandedMatrixTriangularSolve(const double **L, const double **v, const int dNum,
                                    const int rowNo, double **&result)
    {
        int k;
        double s = 0;
        int l=0;
        int n = dNum/2 + 1;
        double **rowIdx;
        for(int i=0; i<rowNo; i++){
            result[i][0] = v[i][0];
        }
        for(int i=0; i<n; i++){
            k=i;
            for(int j=0; j< rowNo-l; j++){
                rowIdx[i][j] = k;
                k++;
            }
            l++;
        }
        result[0][0] = v[0][0] / L[0][0];
        for(int i=1; i<rowNo; i++){
            l=1;
            for(int j=1; j< n; j++){
                for(int k=0; k<rowNo-l; k++){
                    if (rowIdx[j][k] == i){
                        s += L[j][k] * result[k][0];
                    }
                }
                l++;
            }
            result[i][0] = (v[i][0] - s) / L[0][i];
            s = 0;
        }
        return 1;
    }

}
