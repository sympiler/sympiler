//
// Created by Kazem on 10/11/19.
//

#include <sparse_io.h>
#include <cmath>
#include "includes/utils.h"
namespace sym_lib{

 bool time_cmp(timing_measurement a, timing_measurement b){
  return a.elapsed_time<b.elapsed_time;}

 timing_measurement
 time_median(std::vector<timing_measurement> time_array){
  size_t n = time_array.size();
  if(n==0){
   timing_measurement t;
   t.elapsed_time=-1;
   return t;
  }
  std::sort(time_array.begin(),time_array.end(),time_cmp);
  if(n==1)
   return time_array[0];
  return time_array[n/2];
 }



/* Safely compute a*k, where k should be small, and check for integer overflow.
 * If overflow occurs, return 0 and set OK to FALSE.  Also return 0 if OK is
 * FALSE on input. */
 size_t mult_size_t(size_t a, size_t k, int *ok) {
  size_t p = 0, s;
  while (*ok) {
   if (k % 2) {
    p = p + a;
    (*ok) = (*ok) && (p >= a);
   }
   k = k / 2;
   if (!k) return (p);
   s = a + a;
   (*ok) = (*ok) && (s >= a);
   a = s;
  }
  return (0);
 }


 size_t add_size_t (size_t a, size_t b, int *ok){
  size_t s = a + b ;
  (*ok) = (*ok) && (s >= a) ;
  return ((*ok) ? s : 0) ;
 }

 int *partition_by_weight(int n,const int *set, const double *weight,
   int n_parts,
   double *target_weight){
  double *even_weight;
  if(target_weight)
   even_weight = target_weight;
  else{
   even_weight = new double[n_parts];
   double even_w = std::ceil(sum_vector(n, weight) / n_parts);
   std::fill_n(even_weight, n_parts, even_w);
  }
  int *indices = new int[n_parts+1]();
  int j = 0;
  for (int i = 0; i < n_parts; ++i) {
  double c_wgt = 0;
   while (c_wgt < even_weight[i] && j < n){
    int c_n = set[j];
    c_wgt += weight[c_n];
    j++;
   }
   indices[i+1] = j;
  }
  if(!target_weight)
   delete []even_weight;
  return indices;
 }


 void partition_by_bfs(int n, CSC *df, const double *weight,
                       int &final_level_no, // will be one in this case
                       int* &fina_level_ptr, int part_no,
                       int *&final_part_ptr, int *&final_node_ptr,
                          double *target_weight){
  double *even_weight;
  final_level_no = 1;
  fina_level_ptr = new int[final_level_no+1]();
  fina_level_ptr[final_level_no] = part_no;
  final_part_ptr = new int[part_no+1]();
  final_node_ptr = new int[n]();
  if(target_weight)
   even_weight = target_weight;
  else{
   even_weight = new double[part_no];
   double even_w = std::ceil(sum_vector(n, weight) / part_no);
   std::fill_n(even_weight, part_no, even_w);
  }
  auto *is_visited = new bool[n](); int visited_nodes=0;
  int nxt_start = 0; std::vector<std::vector<int>> set_temp;
  set_temp.resize(part_no);
  std::vector<int> stack;
  for (int i = 0; i < part_no && visited_nodes<n; ++i) {
   assert(stack.size() < n);
   double c_wgt = 0;
   while (c_wgt < even_weight[i] && visited_nodes < n){
    if(stack.empty()){
     nxt_start=-1;
     for (int k = 0; k < n; ++k) {
      if(!is_visited[k]){
       nxt_start = k;
       break;
      }
     }
     if(nxt_start < 0) break;
     stack.push_back(nxt_start);
     is_visited[nxt_start] = true;
    }
    // do bfs per node nxt_start
    while(!stack.empty()){
     int cn = stack[0]; stack.erase(stack.begin());
     c_wgt += weight[cn]; set_temp[i].push_back(cn);visited_nodes++;
     for (int k = df->p[cn]; k < df->p[cn + 1]; ++k) {
      auto cnn = df->i[k];
      if(!is_visited[cnn]){
       stack.push_back(cnn);
       is_visited[cnn] = true;
      }
     }
     if(c_wgt >= even_weight[i])
      break;
    }
   }
  }
  // if anything is remained add it to last part
  for (int m = 0; m < stack.size(); ++m) {
   set_temp[part_no-1].push_back(stack[m]);
  }
  for (int k = 0; k < n; ++k) {
   if(!is_visited[k]){
    set_temp[part_no-1].push_back(k);
   }
  }
  // puting into the set
  for (int l = 0; l < set_temp.size(); ++l) {
   int ss = set_temp[l].size();
   int ip = final_part_ptr[l];
   final_part_ptr[l+1] = ip + ss;
   for (int i = 0; i < ss; ++i) {
    assert(ip < n);
    final_node_ptr[ip] = set_temp[l][i];
    ip++;
   }
   //assert(ip <= n);
  }
  if(!target_weight)
   delete []even_weight;
  delete []is_visited;
 }



}
