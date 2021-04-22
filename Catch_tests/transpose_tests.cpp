//
// Created by Kazem on 10/12/19.
//

#include <def.h>
#include <test_utils.h>
#include <sparse_utilities.h>
#include <sparse_blas_lib.h>
#include "catch.hpp"

using namespace sym_lib;
using namespace std;

TEST_CASE("Transpose random cases", "[transposeRandomChecks]") {

 SECTION("A^T^T = A case") {
  vector<pair<int,double>> configs = {make_pair(10,0.2),
                                      make_pair(100,0.03),
                                      make_pair(200,0.04),
                                      make_pair(400,0.06),
                                      make_pair(500,0.004),
                                      make_pair(1000,0.0005),
                                      make_pair(10000,.0001)};
  for (auto i : configs) {
   size_t n = i.first;
   double density = i.second;
   CSC *A = random_square_sparse(n,density);
   CSC *L = make_half(A->n,A->p,A->i,A->x);
   CSC *B = transpose_symmetric(L, NULLPNTR);
   CSC *C = transpose_symmetric(B, NULLPNTR);
   CHECK(is_equal(L,C));
   delete A;
   delete L;
  }
 }

 SECTION("PAP^T = P^-1BP^-1 case") {
  vector<pair<int,double>> configs = {make_pair(10,0.2),
                                      make_pair(100,0.03),
                                      make_pair(200,0.04),
                                      make_pair(400,0.06),
                                      make_pair(500,0.004),
                                      make_pair(1000,0.0005),
                                      make_pair(10000,.0001)};
  for (auto i : configs) {
   int n = i.first;
   double density = i.second;
   std::vector<int> perm_vec;
   generate_uniq_rand_vector(n,n,perm_vec);
   auto *perm = new int[n];
   auto *perm_inv = new int[n];
   for (int j = 0; j < n; ++j) {
    perm[j] = perm_vec[j];
   }
   compute_inv_perm(n,perm,perm_inv);
   CSC *A = random_square_sparse(n,density);
   CSC *L = make_half(A->n,A->p,A->i,A->x);
   CSC *BPT = transpose_symmetric(L, perm);
   CSC *BP = transpose_symmetric(BPT, NULLPNTR);
   CSC *BT = transpose_symmetric(BP, perm_inv);
   CSC* B = transpose_symmetric(BT, NULLPNTR);
   CHECK(is_equal(L,B));
   delete A;
   delete L;
   delete BPT;
   delete BP;
   delete BT;
   delete B;
   delete []perm;
   delete []perm_inv;
  }
 }

}