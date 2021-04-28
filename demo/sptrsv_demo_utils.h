//
// Created by kazem on 3/16/20.
//

#ifndef FUSION_SPTRSV_DEMO_UTILS_H
#define FUSION_SPTRSV_DEMO_UTILS_H

#include <sparse_inspector.h>
#include <lbc.h>

#include "FusionDemo.h"
#include "sparse_blas_lib.h"

namespace sym_lib {

 class SptrsvSerial : public FusionDemo {
 protected:
  timing_measurement fused_code() override {
   timing_measurement t1;
   t1.start_timer();
   sptrsv_csr(n_, L1_csr_->p, L1_csr_->i, L1_csr_->x, x_in_);
   t1.measure_elapsed_time();
   copy_vector(0,n_,x_in_,x_);
   return t1;
  }

 public:
  SptrsvSerial(CSR *L, CSC *L_csc,
                           double *correct_x, std::string name) :
    FusionDemo(L->n, name) {
   L1_csr_ = L;
   L1_csc_ = L_csc;
   correct_x_ = correct_x;
  };

  ~SptrsvSerial() override {};
 };

 class SptrsvLevelSet : public SptrsvSerial {
 protected:
  int *level_set, *level_ptr, level_no;
  void build_set() override {

   level_no = build_levelSet_CSC(L1_csc_->n, L1_csc_->p, L1_csc_->i,
                                 level_ptr, level_set);
  }

  timing_measurement fused_code() override {
   timing_measurement t1;
   t1.start_timer();
   sptrsv_csr_levelset(n_, L1_csr_->p, L1_csr_->i, L1_csr_->x, x_in_,
                       level_no, level_ptr, level_set);
   t1.measure_elapsed_time();
   copy_vector(0,n_,x_in_,x_);
   return t1;
  }

 public:
  SptrsvLevelSet (CSR *L, CSC *L_csc,
                  double *correct_x, std::string name) :
    SptrsvSerial(L, L_csc, correct_x, name) {
   L1_csr_ = L;
   L1_csc_ = L_csc;
   correct_x_ = correct_x;
  };

  ~SptrsvLevelSet () override {
   delete []level_ptr;
   delete []level_set;
  };
 };

 class SptrsvLBC : public SptrsvSerial {
 protected:
  int final_level_no, *fina_level_ptr, *final_part_ptr, *final_node_ptr;
  int part_no;
  int lp_, cp_, ic_;
  void build_set() override {
   auto *cost = new double[n_]();
   for (int i = 0; i < n_; ++i) {
    cost[i] = L1_csr_->p[i+1] - L1_csr_->p[i];
   }
   get_coarse_levelSet_DAG_CSC_tree(n_, L1_csr_->p, L1_csr_->i,
     L1_csr_->stype,
     final_level_no,
     fina_level_ptr,part_no,
     final_part_ptr,final_node_ptr,
     lp_,cp_, ic_, cost);
   delete []cost;
  }

  timing_measurement fused_code() override {
   timing_measurement t1;
   t1.start_timer();
   sptrsv_csr_lbc(n_, L1_csr_->p, L1_csr_->i, L1_csr_->x, x_in_,
                  final_level_no, fina_level_ptr,
                  final_part_ptr, final_node_ptr);
   t1.measure_elapsed_time();
   copy_vector(0,n_,x_in_,x_);
   return t1;
  }

 public:
  SptrsvLBC (CSR *L, CSC *L_csc,
               double *correct_x, std::string name,
               int lp, int cp, int ic) :
    SptrsvSerial(L, L_csc, correct_x, name) {
   L1_csr_ = L;
   L1_csc_ = L_csc;
   correct_x_ = correct_x;
   lp_=lp; cp_=cp; ic_=ic;
  };

  ~SptrsvLBC () override {
   delete []fina_level_ptr;
   delete []final_part_ptr;
   delete []final_node_ptr;
  };
 };



}
#endif //FUSION_SPTRSV_DEMO_UTILS_H
