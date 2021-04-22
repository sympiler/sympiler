//
// Created by Kazem on 10/10/19.
//

#ifndef PROJECT_DEF_H
#define PROJECT_DEF_H

#include <chrono>
#include <vector>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <queue>
#include <cstring>


#define NULLPNTR nullptr
#define EMPTY -1

#ifdef DBG_LOG
#define PRINT_LOG(x) std::cout << (x)
#else
#define PRINT_LOG(x)
#endif

#if 1 //def CSV_LOG
#define PRINT_CSV(x) std::cout <<(x)<<","
#else
#define PRINT_CSV(x)
#endif

#define CS_FLIP(i) (-(i)-2)
#define CS_MARKED(w, j) (w [j] < 0)
#define CS_MARK(w, j) { w [j] = CS_FLIP (w [j]) ; }

/// Measuring time of a kernel
///
struct timing_measurement {
 double elapsed_time;
 std::chrono::time_point<std::chrono::system_clock> start_time, end_time;
 std::chrono::duration<double> elapsed_seconds;
 std::vector<std::chrono::time_point<std::chrono::system_clock> > t_array;

 timing_measurement() : elapsed_time(0) {}
 ~timing_measurement()= default;

 void start_timer() {
  start_time = std::chrono::system_clock::now();
  t_array.push_back(start_time);
 }

 double measure_elapsed_time() {
  assert(t_array.size() > 0);
  end_time = std::chrono::system_clock::now();
  elapsed_seconds = end_time - t_array[0];
  elapsed_time = elapsed_seconds.count();
  t_array.push_back(end_time);
  return elapsed_time;
 }

 void reset (){
  elapsed_time=0; t_array.clear();
 }

 void print_t_array() {
  for (unsigned long i = 1; i < t_array.size(); ++i) {
   std::chrono::duration<double> et = t_array[i] - t_array[i - 1];
   PRINT_CSV(et.count());
  }
 }

};

/// The struct for storing CSC or CSR format
///
struct CSC {
 size_t m; // rows
 size_t n; // columns
 size_t nnz; // nonzeros
 int stype;
 bool is_pattern;
 bool pre_alloc; //if memory is allocated somewhere other than const.
 int *p; // Column pointer array
 int *i; // Row index array
 double *x;

 CSC(size_t M, size_t N, size_t NNZ) : m(M), n(N), nnz(NNZ) {
  is_pattern = false;
  pre_alloc = false;
  if (N > 0)
   p = new int[N + 1]();
  else
   p = NULLPNTR;
  if (NNZ > 0) {
   i = new int[NNZ]();
   x = new double[NNZ]();
  } else {
   i = NULLPNTR;
   x = NULLPNTR;
  }
  stype = 0;
 };

 CSC(size_t M, size_t N, size_t NNZ, bool ip) :
   m(M), n(N), nnz(NNZ), is_pattern(ip) {
  is_pattern = ip;
  pre_alloc = false;
  if (N > 0)
   p = new int[N + 1]();
  else
   p = NULLPNTR;
  if (NNZ > 0) {
   i = new int[NNZ]();
   if (!is_pattern)
    x = new double[NNZ]();
   else
    x = NULLPNTR;
  } else {
   i = NULLPNTR;
   x = NULLPNTR;
  }
  stype = 0;
 };

 CSC(size_t M, size_t N, size_t NNZ, int *Ap, int *Ai, double *Ax) {
  is_pattern = false;
  pre_alloc = true;
  m = M;
  n = N;
  nnz = NNZ;
  p = Ap;
  i = Ai;
  x = Ax;
 }

 CSC(size_t M, size_t N, size_t NNZ, int *Ap, int *Ai, int st) {
  is_pattern = true;
  pre_alloc = true;
  m = M;
  n = N;
  nnz = NNZ;
  p = Ap;
  i = Ai;
  x = NULLPNTR;
  stype = st;
 }

 ~CSC() {
  if(!pre_alloc){
   if (n > 0)
    delete[]p;
   if (nnz > 0) {
    delete[]i;
    if (!is_pattern)
     delete[]x;
   }
  }
 }

};

struct CSR {
 size_t m; // rows
 size_t n; // columns
 size_t nnz; // nonzeros
 int stype;
 bool is_pattern;
 bool pre_alloc; //if memory is allocated somewhere other than const.
 int *p; // Row pointer array
 int *i; // Column index array
 double *x;

 CSR(size_t M, size_t N, size_t NNZ) : m(M), n(N), nnz(NNZ) {
  is_pattern = false;
  pre_alloc=false;
  if (M > 0)
   p = new int[M + 1]();
  else
   p = NULLPNTR;
  if (NNZ > 0) {
   i = new int[NNZ]();
   x = new double[NNZ]();
  } else {
   i = NULLPNTR;
   x = NULLPNTR;
  }
  stype = 0;
 }

 CSR(size_t M, size_t N, size_t NNZ, bool ip) :
   m(M), n(N), nnz(NNZ), is_pattern(ip) {
  is_pattern = ip;
  pre_alloc=false;
  if (M > 0)
   p = new int[M + 1]();
  else
   p = NULLPNTR;
  if (NNZ > 0) {
   i = new int[NNZ]();
   if (!is_pattern)
    x = new double[NNZ]();
   else
    x = NULLPNTR;
  } else {
   i = NULLPNTR;
   x = NULLPNTR;
  }
  stype = 0;
 };

 CSR(size_t M, size_t N, size_t NNZ, int *Ap, int *Ai, int st) {
  is_pattern = true;
  pre_alloc = true;
  m = M;
  n = N;
  nnz = NNZ;
  p = Ap;
  i = Ai;
  x = NULLPNTR;
  stype = st;
 }

 ~CSR() {
  if(!pre_alloc){
   if (m > 0)
    delete[]p;
   if (nnz > 0) {
    delete[]i;
    delete[]x;
   }
  }
 }

};


struct BCSC {
 size_t m; // rows
 size_t n; // cols
 size_t nnz;
 int *p;
 int *i;
 double *x;

 // block format specifics
 size_t nodes;
 int *supernodes;
 int *nrows;

 BCSC (CSC *A) : m(A->m), n(A->n) {
  assert(A != nullptr);
 }

 ~BCSC() {
  delete[]p;
  delete[]i;
  delete[]x;
  delete[]supernodes;
  delete[]nrows;
 }
};

struct Dense {
 size_t row;
 size_t col;
 size_t lda;
 double *a;

 Dense(size_t M, size_t N, size_t LDA) : row(M), col(N), lda(LDA) {
  assert(lda == 1 || lda == N);
  if (row > 0 && col > 0)
   a = new double[row * col]();
 }

 ~Dense() {
  if (row > 0 && col > 0)
   delete[]a;
 }
};

#endif //PROJECT_DEF_H
