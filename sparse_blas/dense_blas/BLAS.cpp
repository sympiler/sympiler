#include <dense_blas/BLAS.h>

void lsolve_BLAS(int dim, int num_col, double *M, double *b) {
 int k;
 double x0, x1, x2, x3, x4, x5, x6, x7;
 double *M0;
 register double *Mi0, *Mi1, *Mi2, *Mi3, *Mi4, *Mi5, *Mi6, *Mi7;

 register int col = 0;

 M0 = M;
 while (col < num_col - 7) {
  Mi0 = M0;
  Mi1 = Mi0 + dim + 1;
  Mi2 = Mi1 + dim + 1;
  Mi3 = Mi2 + dim + 1;
  Mi4 = Mi3 + dim + 1;
  Mi5 = Mi4 + dim + 1;
  Mi6 = Mi5 + dim + 1;
  Mi7 = Mi6 + dim + 1;

  x0 = b[col] / *Mi0++;
  x1 = (b[col + 1] - x0 * *Mi0++) / *Mi1++;
  x2 = (b[col + 2] - x0 * *Mi0++ - x1 * *Mi1++) / *Mi2++;
  x3 = (b[col + 3] - x0 * *Mi0++ - x1 * *Mi1++ - x2 * *Mi2++) / *Mi3++;
  x4 = (b[col + 4] - x0 * *Mi0++ - x1 * *Mi1++ - x2 * *Mi2++ - x3 * *Mi3++) /
       *Mi4++;
  x5 = (b[col + 5] - x0 * *Mi0++ - x1 * *Mi1++ - x2 * *Mi2++ - x3 * *Mi3++ -
        x4 * *Mi4++) / *Mi5++;
  x6 = (b[col + 6] - x0 * *Mi0++ - x1 * *Mi1++ - x2 * *Mi2++ - x3 * *Mi3++ -
        x4 * *Mi4++ - x5 * *Mi5++) / *Mi6++;
  x7 = (b[col + 7] - x0 * *Mi0++ - x1 * *Mi1++ - x2 * *Mi2++ - x3 * *Mi3++ -
        x4 * *Mi4++ - x5 * *Mi5++ -
        x6 * *Mi6++) / *Mi7++;

  b[col++] = x0;
  b[col++] = x1;
  b[col++] = x2;
  b[col++] = x3;
  b[col++] = x4;
  b[col++] = x5;
  b[col++] = x6;
  b[col++] = x7;

  for (k = col; k < num_col; k++)
   b[k] = b[k] - x0 * *Mi0++ - x1 * *Mi1++ - x2 * *Mi2++ - x3 * *Mi3++
          - x4 * *Mi4++ - x5 * *Mi5++ - x6 * *Mi6++ - x7 * *Mi7++;

  M0 += 8 * dim + 8;
 }

 while (col < num_col - 3) {
  Mi0 = M0;
  Mi1 = Mi0 + dim + 1;
  Mi2 = Mi1 + dim + 1;
  Mi3 = Mi2 + dim + 1;

  x0 = b[col] / *Mi0++;
  x1 = (b[col + 1] - x0 * *Mi0++) / *Mi1++;
  x2 = (b[col + 2] - x0 * *Mi0++ - x1 * *Mi1++) / *Mi2++;
  x3 = (b[col + 3] - x0 * *Mi0++ - x1 * *Mi1++ - x2 * *Mi2++) / *Mi3++;

  b[col++] = x0;
  b[col++] = x1;
  b[col++] = x2;
  b[col++] = x3;

  for (k = col; k < num_col; k++)
   b[k] = b[k] - x0 * *Mi0++ - x1 * *Mi1++ - x2 * *Mi2++ - x3 * *Mi3++;

  M0 += 4 * dim + 4;
 }

 if (col < num_col - 1) {
  Mi0 = M0;
  Mi1 = Mi0 + dim + 1;

  x0 = b[col] / *Mi0++;
  x1 = (b[col + 1] - x0 * *Mi0++) / *Mi1++;

  b[col++] = x0;
  b[col++] = x1;

  for (k = col; k < num_col; k++) {
   b[k] = b[k] - x0 * *Mi0++ - x1 * *Mi1++;
  }

  M0 += 2 * dim + 2;
 }

 if (col == num_col - 1) {
  Mi0 = M0;
  x0 = b[col] / *Mi0;
  b[col] = x0;
 }
}

void ltsolve_BLAS(int dim, int num_col, double *M, double *b) {
 int k;
 double x0, x1, x2, x3, x4, x5, x6, x7;
 double *M0, *M1, *M2, *M3, *M4, *M5, *M6, *M7;
 register double *Mi0, *Mi1, *Mi2, *Mi3, *Mi4, *Mi5, *Mi6, *Mi7;
 register int col = num_col - 1;

 M0 = M;
 while (col >= 7) {
  M1 = M0 - dim - 1;
  M2 = M1 - dim - 1;
  M3 = M2 - dim - 1;
  M4 = M3 - dim - 1;
  M5 = M4 - dim - 1;
  M6 = M5 - dim - 1;
  M7 = M6 - dim - 1;


  Mi0 = M0 + 1;
  Mi1 = M1 + 2;
  Mi2 = M2 + 3;
  Mi3 = M3 + 4;
  Mi4 = M4 + 5;
  Mi5 = M5 + 6;
  Mi6 = M6 + 7;
  Mi7 = M7 + 8;

  for (k = num_col - 1; k >= col + 1; k--) {
   b[col] -= *Mi0++ * b[k];
   b[col - 1] -= *Mi1++ * b[k];
   b[col - 2] -= *Mi2++ * b[k];
   b[col - 3] -= *Mi3++ * b[k];
   b[col - 4] -= *Mi4++ * b[k];
   b[col - 5] -= *Mi5++ * b[k];
   b[col - 6] -= *Mi6++ * b[k];
   b[col - 7] -= *Mi7++ * b[k];
  }

  x0 = b[col] / *M0;
  x1 = (b[col - 1] - x0 * *(M1 + 1)) / *M1;
  x2 = (b[col - 2] - x0 * *(M2 + 1) - x1 * *(M2 + 2)) / *M2;
  x3 = (b[col - 3] - x0 * *(M3 + 1) - x1 * *(M3 + 2) - x2 * *(M3 + 3)) / *M3;
  x4 = (b[col - 4] - x0 * *(M4 + 1) - x1 * *(M4 + 2) - x2 * *(M4 + 3) -
        x3 * *(M4 + 4)) / *M4;
  x5 = (b[col - 5] - x0 * *(M5 + 1) - x1 * *(M5 + 2) - x2 * *(M5 + 3) -
        x3 * *(M5 + 4) - x4 * *(M5 + 5)) / *M5;
  x6 = (b[col - 6] - x0 * *(M6 + 1) - x1 * *(M6 + 2) - x2 * *(M6 + 3) -
        x3 * *(M6 + 4) - x4 * *(M6 + 5) - x5 * *(M6 + 6)) / *M6;
  x7 = (b[col - 7] - x0 * *(M7 + 1) - x1 * *(M7 + 2) - x2 * *(M7 + 3) -
        x3 * *(M7 + 4) - x4 * *(M7 + 5) - x5 * *(M7 + 6) - x6 * *(M7 + 7)) /
       *M7;

  b[col--] = x0;
  b[col--] = x1;
  b[col--] = x2;
  b[col--] = x3;
  b[col--] = x4;
  b[col--] = x5;
  b[col--] = x6;
  b[col--] = x7;

  M0 = M0 - 8 * dim - 8;
 }

 while (col >= 3) {
  M1 = M0 - dim - 1;
  M2 = M1 - dim - 1;
  M3 = M2 - dim - 1;

  Mi0 = M0 + 1;
  Mi1 = M1 + 2;
  Mi2 = M2 + 3;
  Mi3 = M3 + 4;

  for (k = num_col - 1; k >= col + 1; k--) {
   b[col] -= *Mi0++ * b[k];
   b[col - 1] -= *Mi1++ * b[k];
   b[col - 2] -= *Mi2++ * b[k];
   b[col - 3] -= *Mi3++ * b[k];
  }

  x0 = b[col] / *M0;
  x1 = (b[col - 1] - x0 * *(M1 + 1)) / *M1;
  x2 = (b[col - 2] - x0 * *(M2 + 1) - x1 * *(M2 + 2)) / *M2;
  x3 = (b[col - 3] - x0 * *(M3 + 1) - x1 * *(M3 + 2) - x2 * *(M3 + 3)) / *M3;

  b[col--] = x0;
  b[col--] = x1;
  b[col--] = x2;
  b[col--] = x3;

  M0 = M0 - 4 * dim - 4;
 }

 while (col >= 1) {
  M1 = M0 - dim - 1;
  Mi0 = M0 + 1;
  Mi1 = M1 + 2;

  for (k = num_col - 1; k >= col + 1; k--) {
   b[col] -= *Mi0++ * b[k];
   b[col - 1] -= *Mi1++ * b[k];
  }

  x0 = b[col] / *M0;
  x1 = (b[col - 1] - x0 * *(M1 + 1)) / *M1;

  b[col--] = x0;
  b[col--] = x1;

  M0 = M0 - 2 * dim - 2;
 }

 if (col == 0) {
  Mi0 = M0 + 1;
  for (k = num_col - 1; k >= 1; k--) {
   b[col] -= *Mi0++ * b[k];
  }

  x0 = b[col] / *M0;
  b[col--] = x0;

  M0 = M0 - dim - 1;
 }
}

void
mattvec_BLAS(int dim, int nrow, int ncol, double *M, double *x, double *b) {
 double *M0 = M;
 register double *Mi0, *Mi1, *Mi2, *Mi3, *Mi4, *Mi5, *Mi6, *Mi7;
 register int col = 0;

 int k;

 while (col < ncol - 7) {
  Mi0 = M0;
  Mi1 = Mi0 + dim;
  Mi2 = Mi1 + dim;
  Mi3 = Mi2 + dim;
  Mi4 = Mi3 + dim;
  Mi5 = Mi4 + dim;
  Mi6 = Mi5 + dim;
  Mi7 = Mi6 + dim;

  for (k = 0; k < nrow; k++) {
   b[col] += *Mi0++ * x[k];
   b[col + 1] += *Mi1++ * x[k];
   b[col + 2] += *Mi2++ * x[k];
   b[col + 3] += *Mi3++ * x[k];
   b[col + 4] += *Mi4++ * x[k];
   b[col + 5] += *Mi5++ * x[k];
   b[col + 6] += *Mi6++ * x[k];
   b[col + 7] += *Mi7++ * x[k];
  }

  M0 += 8 * dim;
  col += 8;
 }

 while (col < ncol - 3) {
  Mi0 = M0;
  Mi1 = Mi0 + dim;
  Mi2 = Mi1 + dim;
  Mi3 = Mi2 + dim;

  for (k = 0; k < nrow; k++) {
   b[col] += *Mi0++ * x[k];
   b[col + 1] += *Mi1++ * x[k];
   b[col + 2] += *Mi2++ * x[k];
   b[col + 3] += *Mi3++ * x[k];
  }

  M0 += 4 * dim;
  col += 4;
 }

 while (col < ncol) {
  Mi0 = M0;
  for (k = 0; k < nrow; k++) {
   double tmp = *Mi0++ * x[k];

   b[col] += tmp;
  }
  M0 += dim;
  col++;
 }
}

void matvec_BLAS(int dim, int nrow, int ncol, double *M, double *x, double *b) {
 double x0, x1, x2, x3, x4, x5, x6, x7;
 double *M0 = M;
 register double *Mi0, *Mi1, *Mi2, *Mi3, *Mi4, *Mi5, *Mi6, *Mi7;
 register int col = 0;

 int k;

 while (col < ncol - 7) {
  Mi0 = M0;
  Mi1 = Mi0 + dim;
  Mi2 = Mi1 + dim;
  Mi3 = Mi2 + dim;
  Mi4 = Mi3 + dim;
  Mi5 = Mi4 + dim;
  Mi6 = Mi5 + dim;
  Mi7 = Mi6 + dim;

  x0 = x[col++];
  x1 = x[col++];
  x2 = x[col++];
  x3 = x[col++];
  x4 = x[col++];
  x5 = x[col++];
  x6 = x[col++];
  x7 = x[col++];

  for (k = 0; k < nrow; k++) {
   b[k] += x0 * *Mi0++ + x1 * *Mi1++ + x2 * *Mi2++ + x3 * *Mi3++
           + x4 * *Mi4++ + x5 * *Mi5++ + x6 * *Mi6++ + x7 * *Mi7++;
  }

  M0 += 8 * dim;
 }

 while (col < ncol - 3) {
  Mi0 = M0;
  Mi1 = Mi0 + dim;
  Mi2 = Mi1 + dim;
  Mi3 = Mi2 + dim;

  x0 = x[col++];
  x1 = x[col++];
  x2 = x[col++];
  x3 = x[col++];

  for (k = 0; k < nrow; k++) {
   b[k] += x0 * *Mi0++ + x1 * *Mi1++ + x2 * *Mi2++ + x3 * *Mi3++;
  }

  M0 += 4 * dim;
 }

 while (col < ncol) {
  Mi0 = M0;
  x0 = x[col++];

  for (k = 0; k < nrow; k++)
   b[k] += x0 * *Mi0++;

  M0 += dim;
 }
}


double DDOT(int n, double *v0, double *v1) {
 register int i = 0;
 register double sum = 0.0;
 register double *M0 = v0;
 register double *M1 = v1;

 while(i < n - 4) {
  sum += *M0++ * *M1++ + *M0++ * *M1++ + *M0++ * *M1++ + *M0++ * *M1++;
  i += 4;
 }

 while(i < n) {
  sum += *M0++ * *M1++;
  i++;
 }
 return sum;
}


// TODO: make this fast currently column majored
void DGEMM(int N, int M, int K, const double *A, const double *B, double *C) {
 // A is N * M
 // B is M * K
 for(int i = 0; i < N; i++) {
  for(int k = 0; k < K; k++) {
   for(int j = 0; j < M; j++) {
    C[i+k*N] += A[i+j*N] * B[j+k*M];
   }
  }
 }
}


void GEMV(int N, int M, const double *A, const double *x, double *b) {
 for(int i = 0; i < M; i++) {
  for(int j = 0; j < N; j++) {
   b[j] += A[j+i*N] * x[i];
  }
 }
}


// TODO: verify
void DAXPY(int n, double a, double *x, double *y, double *z) {
 register int i = 0;
 register double *M0 = x;
 register double *M1 = y;
 register double *M2 = z;

 while(i < n - 8) {
  M2[0] = a * *M0++ + *M1++;
  M2[1] = a * *M0++ + *M1++;
  M2[2] = a * *M0++ + *M1++;
  M2[3] = a * *M0++ + *M1++;
  M2[4] = a * *M0++ + *M1++;
  M2[5] = a * *M0++ + *M1++;
  M2[6] = a * *M0++ + *M1++;
  M2[7] = a * *M0++ + *M1++;
  M2 += 8;
  i += 8;
 }

 while(i < n - 4) {
  M2[0] = a * *M0++ + *M1++;
  M2[1] = a * *M0++ + *M1++;
  M2[2] = a * *M0++ + *M1++;
  M2[3] = a * *M0++ + *M1++;
  M2 += 4;
  i += 4;
 }

 while(i < n - 2) {
  M2[0] = a * *M0++ + *M1++;
  M2[1] = a * *M0++ + *M1++;
  M2 += 2;
  i += 2;
 }

 while(i < n) {
  M2[0] = a * *M0++ + *M1++;
  M2++;
  i++;
 }
}


void AtA(int N, int M, const double *A, double *C) {
 // resulting matrix is MxM
 for(int i = 0; i < M; i++) {
  for(int j = 0; j < M; j++) {
   for(int k = 0; k < N; k++) {
    C[i+j*M] += A[k+i*N] * A[k+j*N];
   }
  }
 }
}
