//
// Created by George Huang on 2019-11-21.
//

#include <immintrin.h>

double dotprod01(int n, double *v0, double *v1) {
 register int i = 0;
 register double sum = 0.0;
 register double *M0 = v0;
 register double *M1 = v1;

 __m256d a0, a1, b0;

 while(i < n-4) {
  a0 = _mm256_loadu_pd(M0);
  a1 = _mm256_loadu_pd(M1);
  b0 = _mm256_mul_pd(a0, a1);

  sum += b0[0] + b0[1] + b0[2] + b0[3];
  M0 += 4;
  M1 += 4;
  i += 4;
 }

 while(i < n) {
  sum += *M0++ * *M1++;
  i++;
 }

 return sum;
}

double dotprod02(int n, double *v0, double *v1) {
 register int i = 0;
 register double sum = 0.0;
 register double *M0 = v0;
 register double *M1 = v1;

 while(i < n - 4) {
  sum += *M0++ * *M1++ + *M0++ * *M1++ + *M0++ * *M1++ + *M0++ * *M1++;
       //+ *M0++ * *M1++ + *M0++ * *M1++ + *M0++ * *M1++ + *M0++ * *M1++;
  i += 4;
 }

 while(i < n) {
  sum += *M0++ * *M1++;
  i++;
 }
 return sum;
}