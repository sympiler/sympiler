//
// Created by kazem on 6/16/20.
//

#include <def.h>

namespace sym_lib{

 void csr_add_diag(int n, const int *Ap, const int *Ai, const double *Ax,
                   double* d){
  for (int i=0; i<n; i++){
//   Bp[i] = Ap[i];
   for(int j=Ap[i]; j<=Ap[i+1]; j++){
    d[i] += Ax[j];
   }}}




void mat_premult_diag(CSC *A, const double *d) {
 int j, i;

 for (j = 0; j < A->n; j++) {                // Cycle over columns
  for (i = A->p[j]; i < A->p[j + 1]; i++) { // Cycle every row in the column
   A->x[i] *= d[A->i[i]];                  // Scale by corresponding element
   // of d for row i
  }
 }
}


void mat_premult_diag(CSR *A, const double *d) {
 int j, i;

 for (j = 0; j < A->n; j++) {                // Cycle over columns
  for (i = A->p[j]; i < A->p[j + 1]; i++) { // Cycle every row in the column
   A->x[i] *= d[A->i[i]];                  // Scale by corresponding element
   // of d for row i
  }
 }
}

 void mat_premult_diag_parallel(CSR *A, const double *d) {
#pragma omp parallel for
  for (int j = 0; j < A->n; j++) {                // Cycle over columns
   for (int i = A->p[j]; i < A->p[j + 1]; i++) { // Cycle every row in the column
    A->x[i] *= d[A->i[i]];                  // Scale by corresponding element
    // of d for row i
   }
  }
 }

 void mat_premult_diag_parallel(CSC *A, const double *d) {
#pragma omp parallel for
  for (int j = 0; j < A->n; j++) {                // Cycle over columns
   for (int i = A->p[j]; i < A->p[j + 1]; i++) { // Cycle every row in the column
    A->x[i] *= d[A->i[i]];                  // Scale by corresponding element
    // of d for row i
   }
  }
 }

void mat_postmult_diag(CSC *A, const double *d) {
 int j, i;

 for (j = 0; j < A->n; j++) {                // Cycle over columns j
  for (i = A->p[j]; i < A->p[j + 1]; i++) { // Cycle every row i in column j
   A->x[i] *= d[j];                        // Scale by corresponding element
   // of d for column j
  }
 }
}


 void mat_postmult_diag(CSR *A, const double *d) {
  int j, i;

  for (j = 0; j < A->n; j++) {                // Cycle over columns j
   for (i = A->p[j]; i < A->p[j + 1]; i++) { // Cycle every row i in column j
    A->x[i] *= d[j];                        // Scale by corresponding element
    // of d for column j
   }
  }
 }


 void mat_postmult_diag_parallel(CSR *A, const double *d) {
#pragma omp parallel for
  for (int j = 0; j < A->n; j++) {                // Cycle over columns j
   for (int i = A->p[j]; i < A->p[j + 1]; i++) { // Cycle every row i in column j
    A->x[i] *= d[j];                        // Scale by corresponding element
    // of d for column j
   }
  }
 }
void mat_postmult_diag_parallel(CSC *A, const double *d) {
#pragma omp parallel for
 for (int j = 0; j < A->n; j++) {                // Cycle over columns j
  for (int i = A->p[j]; i < A->p[j + 1]; i++) { // Cycle every row i in column j
   A->x[i] *= d[j];                        // Scale by corresponding element
   // of d for column j
  }
 }
}

}

