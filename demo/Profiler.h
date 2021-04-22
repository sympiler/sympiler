//
// Created by kazem on 2020-05-27.
//

#ifndef FUSION_PROFILER_H
#define FUSION_PROFILER_H

#include "FusionDemo.h"
#include "PAPIWrapper.h"

namespace sym_lib{

 class Profiler{
  std::vector<int> event_codes_;
  int bundle_size_;
  int num_instances_;
  std::string name_;
  std::vector<EventCounterBundle> counter_bundles_;
  double cost_unit_{};
  double cost_nnz_{};
  int redundant_nodes_;
  double redundant_nnz_;
  double avg_parallelism_, avg_iter_parallelism_;
  double critical_path_; double max_diff_;
 protected:
  virtual FusionDemo *construct_demo(CSR *L, CSC* L_csc, CSR *A, CSC *A_csc,
                             double *correct_x, std::string name,
                             PAPIWrapper *pw, int p1, int p2, int p3)=0;/*{
   auto *fd = new ProfilingSptrsvSpmvSerialNonFused(L2_csr, L1_csc, B_csr, B,
                                              NULLPNTR,name,pw);
   return fd;
  }*/

 public:
  Profiler(std::vector<int>  event_codes, int bnd_size, int num_inst):
    event_codes_(std::move(event_codes)),bundle_size_(bnd_size),
    num_instances_(num_inst),name_(""){}

  void profile(CSR *L, CSC* L_csc, CSR *A, CSC *A_csc,
               double *correct_x, std::string name, int p1=2, int p2=2, int p3=2){
   int num_iter = ceil(event_codes_.size() / (double)bundle_size_) ;
   for (int i = 0; i < num_iter; ++i) {
    std::vector<int> tmp;
    for (int j = i*bundle_size_;
    j < std::min((i + 1) * bundle_size_,(int)event_codes_.size()); ++j) {
     tmp.emplace_back(event_codes_[j]);
    }
    auto *pw = new PAPIWrapper(tmp, num_instances_);
    auto *fd = construct_demo(L, L_csc, A, A_csc, correct_x, name, pw, p1, p2, p3);
    fd->set_num_threads(p1);
    fd->evaluate();
    counter_bundles_.insert(counter_bundles_.end(),pw->counter_bundles_.begin(),
      pw->counter_bundles_.end());
    name_ = fd->Name();
    cost_unit_ = fd->CostUnit();
    cost_nnz_ = fd->CostNNZ();
    redundant_nodes_ = fd->redundantNodes();
    redundant_nnz_ = fd->redundantNNZ();
    avg_parallelism_ = fd->avgParallleism();
    avg_iter_parallelism_ = fd->avgIterParallelism();
    critical_path_ = fd->criticalPath();
    max_diff_ = fd->maxDiff();
    delete pw;
    delete fd;
   }
  }


  void print_headers(){
   for (int j = 0; j < counter_bundles_.size(); ++j) {
    PRINT_CSV(counter_bundles_[j].event_name);
   }
  }

  void print_counters(){
   for (int i = 0; i < num_instances_; ++i) {
    for (int j = 0; j < counter_bundles_.size(); ++j) {
     PRINT_CSV(counter_bundles_[j].counters[i]);
    }
   }
  }

  void print_counters(int i){
    for (int j = 0; j < counter_bundles_.size(); ++j) {
     PRINT_CSV(counter_bundles_[j].counters[i]);
    }
  }

  std::string Name(){ return name_;}
  double CostNNZ(){ return cost_nnz_;}
  double CostUnit(){ return cost_unit_;}
  int redundantNodes(){ return redundant_nodes_;}
  double redundantNNZ(){ return redundant_nnz_;}
  double avgParallelism(){ return avg_parallelism_;}
  double avgIterParallelism(){ return avg_iter_parallelism_;}
  double criticalPath(){ return critical_path_; }
  double maxDiff(){ return max_diff_; }
 };

}
#endif //FUSION_PROFILER_H
