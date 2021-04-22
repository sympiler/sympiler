//
// Created by Kazem on 10/11/19.
//
#include <random>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include "includes/test_utils.h"
namespace sym_lib{



 double double_rand(double d_min, double d_max) {
  double f = (double)rand() / RAND_MAX;
  return d_min + f * (d_max - d_min);
 }


 void generate_uniq_rand_vector(int m, int n, std::vector<int> &rand_vec,
                                unsigned seed){
  std::vector<int> numbers;
  for(int i=0; i<n; i++)
   numbers.push_back(i);
  if(seed == 0)
   seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::shuffle(numbers.begin(), numbers.end(), std::default_random_engine(seed));
  for (int j = 0; j < m; ++j) {
   rand_vec.push_back(numbers[j]);
  }
 }



 CSC *random_square_sparse(size_t n, double density, double max_val,
                           unsigned seed){
  // Warning: for large n, e.g., n>5e5, this is very slow
  int nnz_per_col = density * n + 1;
  nnz_per_col = nnz_per_col> n ? n : nnz_per_col;
  size_t  nnz = n*nnz_per_col;
  if(nnz > INT32_MAX)
   return NULLPNTR;
  CSC *A = new CSC(n,n,nnz);
  A->p[0]=0;
  for (int i = 0; i < n; ++i) {
   A->p[i+1] = A->p[i] + nnz_per_col;
   std::vector<int> rand_vec;
   generate_uniq_rand_vector(nnz_per_col,n,rand_vec,seed);
   if( std::find(std::begin(rand_vec),std::end(rand_vec),i) ==
   std::end(rand_vec) ){//if diagonal idx is not there,
    rand_vec[0] = i; //Put the diagonal in the first location
   }
   std::sort(rand_vec.begin(),rand_vec.end());
   for (int j = A->p[i], k=0; j < A->p[i + 1]; ++j, ++k) {
    A->i[j] = rand_vec[k];
    A->x[j] = double_rand(0,max_val);
   }
  }
  return A;
 }

 CSC *random_symmetric_sparse(size_t n, double density, double max_val, unsigned seed) {
  // Warning: for large n, e.g., n>5e5, this is very slow
  int nnz_per_col = density * n + 1;
  nnz_per_col = nnz_per_col> n ? n : nnz_per_col;
  size_t  nnz = n*nnz_per_col;
  if(nnz > INT32_MAX)
   return NULLPNTR;
  CSC *A = new CSC(n,n,nnz);
  // TODO
  return A;
 }

 bool is_equal(CSC *A, CSC *B, double eps){
  if(A->n != B->n || A->nnz != B->nnz){
   return false;
  }
  for (int j = 0; j < A->n+1; ++j) {
   if(A->p[j] != B->p[j])
    return false;
  }
  for (int i = 0; i < A->nnz; ++i) {
   if(A->i[i] != B->i[i])
    return false;
  }
  for (int i = 0; i < A->nnz; ++i) {
   if(std::abs(A->x[i] -B->x[i]) > eps)
    return false;
  }
  return true;
 }


 void rhs_init(int n, int *Ap, int *Ai, double *Ax, double *b){
  /*generating a rhs that produces a result of all 1 vector*/
  for (int j = 0; j < n; ++j) {
   b[j]=0;
  }
  for (int c = 0; c < n ; ++c) {
   for (int cc = Ap[c]; cc < Ap[c + 1]; ++cc) {
    b[Ai[cc]]+=Ax[cc];
   }
  }
 }


 void rhs_init_blocked(size_t n, size_t nBlocks, size_t *Ap, int *Ai,
   size_t *AiP, double *Ax, double *b){
  /*generating a rhs that produces a result of all 1 vector*/
  for (int j = 0; j < n; ++j) {
   b[j]=0;
  }
  for (int c = 0; c < n ; ++c) {
   for (int cc = Ap[c], j=0; cc < Ap[c + 1]; ++cc, ++j) {
    size_t curRow = Ai[AiP[c]+j];
    b[curRow]+=Ax[cc];
   }
  }
 }

 bool test_unique(int n, int *vec){
  auto *tmp = new bool[n]();
  for (int i = 0; i < n; ++i) {
   assert(vec[i] < n && vec[i] >= 0);
   tmp[vec[i]] =  true;
  }
  for (int j = 0; j < n; ++j) {
   if(!tmp[j]){
    delete []tmp;
    return false;
   }
  }
  delete []tmp;
  return true;
 }

}
