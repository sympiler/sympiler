//
// Created by kazem on 10/9/19.
//

#include <omp.h>
#include "sparse_blas_lib.h"

namespace sym_lib {

 void spmv_csr(int n, const int *Ap, const int *Ai, const double *Ax,
               const double *x, double *y) {
  int i, j;
  for (i = 0; i < n; i++) {
   for (j = Ap[i]; j < Ap[i + 1]; j++) {
    y[i] += Ax[j] * x[Ai[j]];
   }
  }
 }

 void spmv_csr(int n, const int *Ap, const int *Ai, const double *Ax,
               const double *x, const double *y, double *z, double a,
               double b) {
  int i, j;
#pragma omp parallel for default(shared) schedule(auto)
  for (i = 0; i < n; i++) {
   for (j = Ap[i]; j < Ap[i + 1]; j++) {
    z[i] += Ax[j] * x[Ai[j]];
   }
   z[i] = a * z[i] + b * y[i];
  }
 }

#define Prof11
 void spmv_csr_parallel(int n, const int *Ap, const int *Ai, const double *Ax,
                        const double *x, double *y) {
#pragma omp parallel
  {
#ifdef Prof1
   double wtime = omp_get_wtime();
#endif
#pragma omp for schedule(static) nowait
  for (int i = 0; i < n; i++) {
   for (int j = Ap[i]; j < Ap[i + 1]; j++) {
    y[i] += Ax[j] * x[Ai[j]];
   }
  }
#ifdef Prof1
   wtime = omp_get_wtime() - wtime;
      printf( "Time taken by thread %d is %f\n", omp_get_thread_num(), wtime );
#endif
  }
 // std::cout<<"\n";
 }

 void spmv_csr_block(int n, const int *Ap, const int *Ai, const double *Ax,
                     int nodes, const int *supernodes, const double *x,
                     double *y) {
#pragma omp parallel for default(shared) schedule(auto)
  for (int i = 0; i < nodes; i++) {
   for (int j = supernodes[i]; j < supernodes[i + 1]; j++) {
    y[j] = 0;
    for (int k = Ap[j]; k < Ap[j + 1]; k++) {
     y[j] += Ax[k] * x[Ai[k]];
    }
   }
  }
 }

 void spmv_csc(int n, const int *Ap, const int *Ai, const double *Ax,
               const double *x, double *y) {
  int i, j;
  for (i = 0; i < n; i++) {
   for (j = Ap[i]; j < Ap[i + 1]; j++) {
    y[Ai[j]] += Ax[j] * x[i];
   }
  }
 }



 void spmv_csc(int n, const int *Ap, const int *Ai, const double *Ax,
               const double *x, const double *y, double *z, double a,
               double b) {
  for (int i = 0; i < n; i++) {
   for (int j = Ap[i]; j < Ap[i + 1]; j++) {
    z[Ai[j]] += Ax[j] * x[i];
   }
  }
  for (int i = 0; i < n; i++)
   z[i] = a * z[i] + b * y[i];
 }

 void spmv_csc(int n, const int *Ap, const int *Ai, const double *Ax,
               const double *x, const double *y, double *z, double *a,
               double *b) {
  for (int i = 0; i < n; i++) {
   for (int j = Ap[i]; j < Ap[i + 1]; j++) {
    z[Ai[j]] += Ax[j] * x[i];
   }
  }

#pragma omp parallel for schedule(auto)
  for (int i = 0; i < n; i++)
   z[i] = a[i] * z[i] + b[i] * y[i];
 }

 void spmv_csc_parallel(int n, const int *Ap, const int *Ai, const double *Ax,
                        const double *x, double *y) {
#pragma omp parallel for default(shared) schedule(auto)
  for (int i = 0; i < n; i++) {
   for (int j = Ap[i]; j < Ap[i + 1]; j++) {
#pragma omp atomic
    y[Ai[j]] += Ax[j] * x[i];
   }
  }
 }

 void spmv_csc_lbc(int n_parts, const int *sets, const int *ids,
                   int n, const int *Ap,
                   const int *Ai,
                   const double *Ax, const double *x, double *y) {
  double *red_tmp = new double[n_parts*n]();
#pragma omp parallel for default(shared) schedule(auto)
  for (int l = 0; l < n_parts; l++) {
   for (int k = sets[l]; k < sets[l + 1]; ++k) {
    int i = ids[k];
    for (int j = Ap[i]; j < Ap[i + 1]; j++) {
     red_tmp[l*n + Ai[j]] += Ax[j] * x[i];
    }
   }
  }
  for (int m = 0; m < n_parts; ++m) {
   for (int i = 0; i < n; ++i) {
    y[i] += red_tmp[m*n + i];
   }
  }

  delete []red_tmp;
 }


 void spmv_csc_lbc(const int *set_ptr, const int *sets, const int *ids,
                   int n, const int *Ap,
                   const int *Ai,
                   const double *Ax, const double *x, double *y) {
  int n_parts = set_ptr[1];
  double *red_tmp = new double[n_parts*n]();
#pragma omp parallel for default(shared) schedule(auto)
  for (int l = 0; l < n_parts; l++) {
   for (int k = sets[l]; k < sets[l + 1]; ++k) {
    int i = ids[k];
    for (int j = Ap[i]; j < Ap[i + 1]; j++) {
     red_tmp[l*n + Ai[j]] += Ax[j] * x[i];
    }
   }
  }
  for (int m = 0; m < n_parts; ++m) {
   for (int i = 0; i < n; ++i) {
    y[i] += red_tmp[m*n + i];
   }
  }

  delete []red_tmp;
 }


 void
 spmv_csc_block(int n, const int *Ap, const int *Ai, const double *Ax,
                const double *x, double *y, int nodes,
                const int *supernodes) {

//#pragma omp parallel for default(shared) schedule(auto)
  for (int i = 0; i < nodes; i++) {
   for (int j = supernodes[i]; j < supernodes[i + 1]; j++) {
    for (int k = Ap[j]; k < Ap[j + 1]; k++) {

//#pragma omp atomic
     y[Ai[k]] += Ax[k] * x[j];
    }
   }
  }
 }
}
