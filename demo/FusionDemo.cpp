//
// Created by Kazem on 11/10/19.
//
#define DBG_LOG
#define CSV_LOG
#include <def.h>
#include <test_utils.h>
#include <utils.h>
#include "FusionDemo.h"
#include <sparse_utilities.h>
#ifdef METIS
#include <metis_interface.h>
#endif

namespace sym_lib{

 FusionDemo::FusionDemo():L1_csr_(NULLPNTR), L1_csc_(NULLPNTR),
                          L2_csr_(NULLPNTR), L2_csc_(NULLPNTR),
                          A_csr_(NULLPNTR),A_csc_(NULLPNTR),
                          x_(NULLPNTR),
                          x_in_(NULLPNTR), correct_x_(NULLPNTR){
  num_test_=5;
  redundant_nodes_=0;
#ifdef PROFILE
  pw_ = NULLPNTR;
#endif
  }

  FusionDemo::FusionDemo(int n, std::string name):FusionDemo() {
  n_ = n;
  name_ = name;
  x_in_ = new double[n]();
  x_ = new double[n]();
 }



#ifdef PROFILE
  FusionDemo::FusionDemo(int n, std::string name, PAPIWrapper *pw):FusionDemo(n, name){
   pw_ = pw;
 }

 FusionDemo::FusionDemo(CSR *L, CSC* L_csc, CSR *A, CSC *A_csc,
                        double *correct_x, std::string name, PAPIWrapper *pw):
   FusionDemo(L->n,name, pw){
  L1_csr_ = L;
  L1_csc_ = L_csc;
  A_csr_ = A;
  A_csc_ = A_csc;
  correct_x_ = correct_x;
 }
#endif

  FusionDemo::~FusionDemo() {
  delete []x_in_;
  delete []x_;
 }

 void FusionDemo::setting_up() {
  std::fill_n(x_in_,n_,1);
  std::fill_n(x_,n_,0.0);
 }

 void FusionDemo::testing() {
  if(correct_x_)
   if (!is_equal(0, n_, correct_x_, x_,1e-6))
    PRINT_LOG(name_ + " code != reference solution.\n");
 }


 timing_measurement FusionDemo::evaluate() {
  timing_measurement median_t;
  std::vector<timing_measurement> time_array;
  analysis_time_.start_timer();
  build_set();
  analysis_time_.measure_elapsed_time();
  for (int i = 0; i < num_test_; ++i) {
   setting_up();
#ifdef PROFILE
   if(pw_)
    pw_->begin_profiling();
#endif
   timing_measurement t1 = fused_code();
#ifdef PROFILE
   if(pw_)
    pw_->finish_profiling();
#endif
   time_array.emplace_back(t1);
  }
  testing();

  median_t = time_median(time_array);
  return median_t;
 }

 timing_measurement FusionDemo::analysisTime() {
  return analysis_time_;
 }



 void generate_matrices_from_mtx(CSC *L1_csc_in,
   CSR *&L2_csr, CSC *&B, CSR *&B_csr){
  CSC *L1_csc = NULLPNTR;
#ifdef METIS
  //We only reorder L since dependency matters more in l-solve.
  auto *perm = new int[L1_csc_in->n]();
  CSC *L1_csc_full = make_full(L1_csc_in);
  metis_perm_general(L1_csc_full, perm);
  L1_csc = make_half(L1_csc_full->n, L1_csc_full->p, L1_csc_full->i,
                     L1_csc_full->x);
  CSC *Lt = transpose_symmetric(L1_csc, perm);
  CSC *L1_ord = transpose_symmetric(Lt, NULLPNTR);
  delete L1_csc;
  L1_csc = L1_ord;
  delete Lt;
  delete L1_csc_full;
  delete[]perm;
#endif
  L2_csr = csc_to_csr(L1_csc);
  if(true){
   B = make_full(L1_csc);
   B_csr = csc_to_csr(B);
  } else {
   B = diagonal(L1_csc->n,1.0);
   B_csr = csc_to_csr(B);
  }
 }

 void print_common_header(){
 PRINT_CSV("Matrix Name,A Dimension,A Nonzero,L Nonzero,Code Type,Data Type,"
           "Metis Enabled,Number of Threads");
 }

 void print_common(std::string matrix_name, std::string variant, std::string strategy,
   CSC *B, CSC *L, int num_threads){
  PRINT_CSV(matrix_name);
  PRINT_CSV(B->m);
  PRINT_CSV(B->nnz);
  if(L)
   PRINT_CSV(L->nnz);
  PRINT_CSV(variant);
  PRINT_CSV(strategy);
#ifdef METIS
  PRINT_CSV("Metis");
#else
  PRINT_CSV("No Metis");
#endif
  PRINT_CSV(num_threads);
 }

}