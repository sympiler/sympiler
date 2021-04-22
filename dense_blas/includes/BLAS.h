#ifndef BLAS_H
#define BLAS_H


void lsolve_BLAS(int dim, int num_col, double *M, double *b);

void ltsolve_BLAS(int dim, int num_col, double *M, double *b);

void matvec_BLAS(int dim, int nrow, int ncol, double *M, double *x, double *b);

void mattvec_BLAS(int dim, int nrow, int ncol, double *M, double *x, double *b);

double DDOT(int n, double *v0, double *v1);

void DGEMM(int N, int M, int K, const double *A, const double *B, double *C);

void GEMV(int N, int M, const double *A, const double *x, double *b);

void DAXPY(int n, double a, double *x, double *y, double *z);

void AtA(int N, int M, const double *A, double *C);

#endif