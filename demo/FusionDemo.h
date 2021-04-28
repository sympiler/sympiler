//
// Created by Kazem on 11/10/19.
//

#ifndef PROJECT_FUSIONDEMO_H
#define PROJECT_FUSIONDEMO_H


#include <string>

#ifdef PAPI
#include "PAPIWrapper.h"
#else
#undef PROFILE
#endif


namespace sym_lib{
 class FusionDemo {

 protected:
  int n_{};
  double *correct_x_{};
  double *x_{}, *x_in_{};
  std::string name_{};

  int num_test_{};
  int num_threads_{};
  CSR *L1_csr_, *L2_csr_, *A_csr_;
  CSC *L1_csc_, *L2_csc_, *A_csc_;

  /// profiling info
  int redundant_nodes_{};
  double redundant_nnz_{};
  double cost_nnz_{}; double cost_unit_{};
  double avg_parallelism{}; double avg_iter_parallelism{};
  double critical_path_{}; double max_diff_{};
  timing_measurement analysis_time_{};

  std::vector<CSC*> dag_array_;

  virtual void build_set(){};
  virtual void setting_up();
  virtual timing_measurement fused_code() = 0;
  virtual void testing();
 public:
  FusionDemo();
  explicit  FusionDemo(int, std::string);
  virtual ~FusionDemo();

#ifdef PROFILE
  PAPIWrapper *pw_;
  explicit  FusionDemo(int, std::string, PAPIWrapper *pw);
  FusionDemo(CSR *L, CSC* L_csc, CSR *A, CSC *A_csc,
             double *correct_x, std::string name, PAPIWrapper *pw);
  void set_pw(PAPIWrapper *pw){pw_=pw;}
#endif

  int redundantNodes(){ return redundant_nodes_;}
  double redundantNNZ(){ return redundant_nnz_;}
  double *solution(){ return x_;}
  timing_measurement evaluate();

  void set_num_test(int nt){num_test_=nt;};
  void set_num_threads(int nt){num_threads_=nt;};
  std::string Name(){ return name_;}
  timing_measurement analysisTime();
  double CostUnit(){ return cost_unit_;}
  double CostNNZ(){ return cost_nnz_;}
  double avgParallleism(){ return avg_parallelism;}
  double avgIterParallelism(){ return avg_iter_parallelism;}
  double criticalPath(){ return critical_path_; }
  double maxDiff(){return max_diff_; }
 };

 void generate_matrices_from_mtx(CSC *L1_csc,
                                 CSR *&L2_csr, CSC *&B, CSR *&B_csr);

 void print_common_header();
 void print_common(std::string matrix_name, std::string variant, std::string strategy,
                   CSC *B, CSC *L, int num_thread);


}



#endif //PROJECT_FUSIONDEMO_H
