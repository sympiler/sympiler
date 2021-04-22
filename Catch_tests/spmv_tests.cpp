//
// Created by kazem on 10/9/19.
//

//Please do not copy following line if
// you use this file to create a new test case
#define CATCH_CONFIG_MAIN

#include <omp.h>
#include <test_utils.h>
#include <utility>
#include "catch.hpp"
#include "sparse_blas_lib.h"
#include "sparse_utilities.h"

using namespace sym_lib;
using namespace std;

TEST_CASE("Check small general cases", "[spmvCorrectnessChecks]") {

 SECTION("Small SPMV1") {
  // CSR version
  int n = 4;
  int Ap[5] = {0, 2, 3, 4, 6};
  int Ai[6] = {0, 1, 3, 3, 1, 3};
  double Ax[6] = {1.0, 2.0, 1.0, 1.0, 1.0, 1.0};
  double x[4] = {1, 1, 1, 1};
  double y[4] = {};
  double z[4] = {3, 1, 1, 2};
  spmv_csr(n, Ap, Ai, Ax, x, y);
  for (int i = 0; i < n; ++i) {
   CHECK(y[i] == z[i]);
  }

  // CSC version
  int *nAp, *nAi;
  double *nAx;
  csr_to_csc(n, n, Ap, Ai, Ax, nAp, nAi, nAx);
  y[0] = y[1] = y[2] = y[3] = 0.0;
  spmv_csc(n, nAp, nAi, nAx, x, y);
  for (int i = 0; i < n; ++i) {
   CHECK(y[i] == z[i]);
  }
  delete[]nAp;
  delete[]nAi;
  delete[]nAx;
 }

 SECTION("Small SPMV2") {
  // Symmetric matrix

  // CSR version
  int n = 10;
  int Ap[11] = {0, 5, 8, 12, 16, 19, 24, 29, 34, 39, 44};
  int Ai[44] = {0, 1, 5, 7, 8, 0, 1, 7, 2, 6, 7, 9, 3, 4, 5, 8, 3, 4, 6, 0, 3,
                5, 6, 8, 2, 4, 5, 6, 9, 0, 1, 2, 7, 9, 0, 3, 5, 8, 9, 2, 6, 7,
                8, 9};
  double Ax[44];
  for(double & i : Ax)
   i = 1.0;

  double x[10] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  double y0[10] = {};
  double y1[10] = {};
  double z[10] = {5, 3, 4, 4, 3, 5, 5, 5, 5, 5};
  spmv_csr(n, Ap, Ai, Ax, x, y0);
  spmv_csc(n, Ap, Ai, Ax, x, y1);
  for (int i = 0; i < n; ++i) {
   CHECK(y0[i] == z[i]);
   CHECK(y1[i] == z[i]);
  }
 }

 SECTION("SpMV CSR, parallel vs. serial") {
  vector<pair<int,double>> configs = {make_pair(10,0.2),
                                      make_pair(100,0.03),
                                      make_pair(500,0.004),
                                      make_pair(1000,0.0005),
                                      make_pair(10000,.0001)};
  for (auto i : configs) {
   size_t n = i.first;
   double density = i.second;
   CSC *A = random_square_sparse(n,density);
   CSR *A_csr = csc_to_csr(A);

   auto *x = new double[A->n]();
   auto *y_serial = new double[A->n]();
   auto *y_parallel = new double[A->n]();
   std::fill_n(x,A->n,1.0);
   spmv_csr(A_csr->n,A_csr->p,A_csr->i,A_csr->x,x,y_serial);

   std::fill_n(x,A->n,1.0);
   spmv_csr_parallel(A_csr->n,A_csr->p,A_csr->i,A_csr->x,x,y_parallel);

   CHECK(is_equal(0,A->n,y_serial,y_parallel));

   delete []x;
   delete []y_serial;
   delete []y_parallel;
   delete A;
   delete A_csr;
  }
 }

}
