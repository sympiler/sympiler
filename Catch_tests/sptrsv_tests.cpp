//
// Created by george on 2019-10-09.
//

#include <test_utils.h>
#include <omp.h>
#include <sparse_io.h>
#include <sparse_inspector.h>
#include "catch.hpp"
#include "sparse_blas_lib.h"
#include "sparse_utilities.h"


using namespace sym_lib;
using namespace std;

TEST_CASE("Check small lower triangular cases", "[sptrsvCorrectnessChecks]") {

 SECTION("Small SPTRSV1") {
  // CSC version
  int n = 10;
  int Lp[11] = {0, 3, 6, 9, 13, 16, 18, 20, 23, 25, 26};
  int Li[26] = {0, 6, 9, 1, 2, 5, 2, 5, 8, 3, 4, 5, 8, 4, 5, 8, 5, 7, 6, 7, 7,
                8, 9, 8, 9, 9};
  double Lx[26];
  for (auto &i : Lx)
   i = 1.0;
  double b0[10] = {1, 1, 2, 1, 2, 5, 2, 3, 5, 4};
  double b1[10] = {1, 1, 2, 1, 2, 5, 2, 3, 5, 4};
  sptrsv_csc(n, Lp, Li, Lx, b0);
  for (int i = 0; i < n; ++i) {
   CHECK(b0[i] == 1.0);
  }

  // CSR version
  int *nLp, *nLi;
  double *nLx;
  csc_to_csr(n, n, Lp, Li, Lx, nLp, nLi, nLx);
  sptrsv_csr(n, nLp, nLi, nLx, b1);
  for (int i = 0; i < n; ++i) {
   CHECK(b1[i] == 1.0);
  }
  delete[]nLp;
  delete[]nLi;
  delete[]nLx;
 }

 SECTION("Sptrsv CSR, parallel vs. serial") {
  vector<pair<int,double>> configs = {make_pair(10,0.2),
                                      make_pair(100,0.03),
                                      make_pair(500,0.004),
                                      make_pair(1000,0.0005),
                                      make_pair(10000,.0001)};
  //omp_set_num_threads(1);
  for (auto i : configs) {
   size_t n = i.first;
   double density = i.second;
   int *level_set, *level_ptr, level_no;
   CSC *A = random_square_sparse(n,density);
   CSC *L = make_half(A->n,A->p,A->i,A->x);
   level_no = build_levelSet_CSC(L->n, L->p, L->i, level_ptr, level_set);
   CSR *L_csr = csc_to_csr(L);

   auto *y_serial = new double[A->n]();
   auto *y_parallel = new double[A->n]();
   std::fill_n(y_serial,A->n,1.0);
   sptrsv_csr(L_csr->n,L_csr->p,L_csr->i,L_csr->x,y_serial);

   std::fill_n(y_parallel,A->n,1.0);
   sptrsv_csr_levelset(L_csr->n,L_csr->p,L_csr->i,L_csr->x,y_parallel,
     level_no,level_ptr,level_set);

   CHECK(is_equal(0,A->n,y_serial,y_parallel));

   delete []y_serial;
   delete []y_parallel;
   delete []level_ptr;
   delete []level_set;
   delete A;
   delete L;
   delete L_csr;
  }
 }
}
