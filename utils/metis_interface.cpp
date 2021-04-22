//
// Created by kazem on 10/13/19.
//

#include <metis.h>
#include <cassert>
#include <functional>
#include <sparse_io.h>
#include "includes/def.h"
#include "includes/sparse_utilities.h"
#include "includes/metis_interface.h"

namespace sym_lib {
 // FIXME
 CSC *metis_partitioning(CSC *A, int npar, int *&partition) {
  std::function<int(CSC*, int*&)> metis_perm_func;
  std::function<int(CSC*, int*&, int)> metis_part_func;

  if(A->is_pattern) {
   metis_perm_func = metis_perm_symmetric;
   metis_part_func = metis_partition_symmetric;
  } else {
   metis_perm_func = metis_perm_general;
   metis_part_func = metis_partition_general;
  }

  CSC *A1 = nullptr, *A2 = nullptr;
  int *perm = nullptr;
  metis_perm_func(A, perm);
  if(A->is_pattern) {
   A1 = transpose_symmetric(A, perm);
   A2 = make_full(A1);
  } else {
   A2 = copy_sparse(A);
  }
  metis_part_func(A2, partition, npar);

  delete[]perm;
  delete (A1);

  return A2;
 }


 int metis_perm_symmetric(CSC *A, int *&perm) {
  if (A->stype == 0) {//asymmetric case is not supported!
   return -1;
  }

  perm = new int[A->n]();
  int nnzFull = 2 * A->nnz - A->n;
  auto Afullp = new idx_t[A->n + 1]();
  auto Afulli = new idx_t[nnzFull]();
  auto ind = new int[A->n]();
  for (size_t i = 0; i < A->n; i++) {
   for (size_t p = A->p[i]; p < A->p[i + 1]; p++) {
    int row = A->i[p];
    ind[i]++;
    if (row != i)
     ind[row]++;
   }
  }
  Afullp[0] = 0;
  for (size_t i = 0; i < A->n; i++)
   Afullp[i + 1] = Afullp[i] + ind[i];

  for (size_t i = 0; i < A->n; i++)
   ind[i] = 0;
  for (size_t i = 0; i < A->n; i++) {
   for (size_t p = A->p[i]; p < A->p[i + 1]; p++) {
    int row = A->i[p];
    if (row != i) {
     int index = Afullp[i] + ind[i];
     Afulli[index] = row;
     ind[i]++;
     if (row != i) {
      index = Afullp[row] + ind[row];
      Afulli[index] = i;
      ind[row]++;
     }
    }
   }
  }

  //Making the graph for passing it to metis, it should have
  //both upper and lower parts
  idx_t options1[METIS_NOPTIONS];
  METIS_SetDefaultOptions(options1);
  int ncol = A->n;
  idx_t ncolIDXT = ncol;
  auto weigt = new idx_t[ncol];
  auto permIDX = new idx_t[ncol];
  auto IpermIDX = new idx_t[ncol];
  for (int i = 0; i < ncol; ++i) {
   permIDX[i] = 0;
   IpermIDX[i] = 0;
   weigt[i] = 1;
  }

#if 0
  CSC *tmp_mat = new CSC(ncol,ncol,nnzFull,AFullp,AFulli,NULLPNTR)
  std::cout<<number_empty_col()
#endif
  int retMet = METIS_NodeND(&ncolIDXT, Afullp, Afulli, nullptr, options1,
                            permIDX, IpermIDX);
  assert(retMet == METIS_OK);
  for (int i = 0; i < ncol; ++i)
   perm[i] = permIDX[i];
  delete[]Afullp;
  delete[]Afulli;
  delete[]ind;
  delete[]weigt;
  delete[]permIDX;
  delete[]IpermIDX;
  return retMet;
 }


 int metis_perm_general(CSC *A, int *&perm) {
  perm = new int[A->n]();

  // assumes A is full matrix
  idx_t n = A->n;
  CSR *At = csc_to_csr(A);

  // now remove the diagonal elements
  auto Ap = new idx_t[At->n+1]();
  auto Ai = new idx_t[At->nnz ]();

  Ap[0] = 0;
  int cnt = 0;
  for (int i = 0; i < n; i++) {
   int row_cnt = 0;
   for (int p = At->p[i]; p < At->p[i + 1]; p++) {
    if (i == At->i[p]) // diagonal, ignore
     continue;
    row_cnt++;
    Ai[cnt] = At->i[p];
    cnt++;
   }
   Ap[i + 1] = Ap[i] + row_cnt;
  }

  idx_t options1[METIS_NOPTIONS];
  METIS_SetDefaultOptions(options1);
  int ncol = A->n;
  idx_t ncolIDXT = ncol;
  auto weigt = new idx_t[ncol];
  auto permIDX = new idx_t[ncol];
  auto IpermIDX = new idx_t[ncol];
  for (int i = 0; i < ncol; ++i) {
   permIDX[i] = 0;
   IpermIDX[i] = 0;
   weigt[i] = 1;
  }

  int retMet = METIS_NodeND(&ncolIDXT, Ap, Ai, nullptr, options1,
                            permIDX, IpermIDX);
  assert(retMet == METIS_OK);
  for (int i = 0; i < ncol; ++i)
   perm[i] = permIDX[i];
  delete At;
  delete[]Ap;
  delete[]Ai;
  delete[]weigt;
  delete[]permIDX;
  delete[]IpermIDX;

  return retMet;
 }

 int metis_partition_symmetric(CSC *A, int *&part, int k) {
  part = new int[A->n]();

  int n = A->n;
  int nnzFull = 2 * A->nnz - A->n;
  auto Afullp = new idx_t[n + 1]();
  auto Afulli = new idx_t[nnzFull]();
  auto ind = new int[n]();
  for (size_t i = 0; i < n; i++) {
   for (size_t p = A->p[i]; p < A->p[i + 1]; p++) {
    int row = A->i[p];
    ind[i]++;
    if (row != i)
     ind[row]++;
   }
  }
  Afullp[0] = 0;
  for (size_t i = 0; i < n; i++)
   Afullp[i + 1] = Afullp[i] + ind[i];

  for (size_t i = 0; i < n; i++)
   ind[i] = 0;
  for (size_t i = 0; i < n; i++) {
   for (size_t p = A->p[i]; p < A->p[i + 1]; p++) {
    int row = A->i[p];
    if (row != i) {
     int index = Afullp[i] + ind[i];
     Afulli[index] = row;
     ind[i]++;
     if (row != i) {
      index = Afullp[row] + ind[row];
      Afulli[index] = i;
      ind[row]++;
     }
    }
   }
  }

  idx_t ncolIDXT = n;
  idx_t nconIDX = 1;
  idx_t nparsIDX = k;
  idx_t objvalIDX = 0;

  idx_t options1[METIS_NOPTIONS];
  METIS_SetDefaultOptions(options1);

  auto *weigt = new idx_t[n];
  auto *partIDX = new idx_t[n];
  for (int i = 0; i < n; ++i) {
   partIDX[i] = 0;
   weigt[i] = 1;
  }

  idx_t retMet = METIS_PartGraphKway(&ncolIDXT, &nconIDX, Afullp, Afulli, nullptr,
                                     nullptr, nullptr, &nparsIDX, nullptr,
                                     nullptr, options1, &objvalIDX, partIDX);
  assert(retMet == METIS_OK);
  for (int i = 0; i < n; ++i)
   part[i] = partIDX[i];
  delete[]Afullp;
  delete[]Afulli;
  delete[]weigt;
  delete[]partIDX;

  return 0;
 }


 int metis_partition_general(CSC *A, int *&part, int k) {
  part = new int[A->n]();

  // assumes A is full matrix
  idx_t n = A->n;
  CSR *At = csc_to_csr(A);

  // now remove the diagonal elements

  auto Ap = new idx_t[At->n+1]();
  auto Ai = new idx_t[At->nnz ]();

  Ap[0] = 0;
  int cnt1 = 0;
  int cnt = 0;
  for (int i = 0; i < n; i++) {
   int row_cnt = 0;
   for (int p = At->p[i]; p < At->p[i + 1]; p++) {
    if (i != At->i[p]) {
     row_cnt++;
     Ai[cnt] = At->i[p];
     cnt++;
    }
    cnt1++;
   }
   Ap[i + 1] = Ap[i] + row_cnt;
  }

  idx_t ncolIDXT = A->n;
  idx_t nconIDX = 1;
  idx_t nparsIDX = k;
  idx_t objvalIDX = 0;

  idx_t options1[METIS_NOPTIONS];
  METIS_SetDefaultOptions(options1);

  auto *weigt = new idx_t[n];
  //auto *unit_wgt = new idx_t[n];
  auto *partIDX = new idx_t[n];
  for (int i = 0; i < n; ++i) {
   partIDX[i] = 0;
   weigt[i] = 1 + A->p[i+1]-A->p[i];
//   for (int j = A->p[i]; j < A->p[i+1]; ++j) {
  //  weigt[i] += (A->p[j+1] - A->p[j]);
    //weigt[i] += 1;
//   }
   //unit_wgt[i] = 1;
  }
/*
  for (int j = 0; j < ncolIDXT; ++j) {
   for (int i = Ap[j]; i < Ap[j + 1]; ++i) {
    std::cout<<j <<", "<< i << " "<<Ai[i]<<"\n";
   }
  }
*/
  idx_t retMet = METIS_PartGraphKway(&ncolIDXT, &nconIDX, Ap, Ai, weigt,
                                     nullptr, nullptr, &nparsIDX, nullptr,
                                     nullptr, options1, &objvalIDX, partIDX);
  assert(retMet == METIS_OK);
  int min_p_no = n; // it happens to start from a value other than zero, maybe a bug.
  for (int i = 0; i < n; ++i) {
   part[i] = partIDX[i];
   min_p_no = std::min(min_p_no,part[i]);
  }
  //normalizing it
  if(min_p_no != 0){
   for (int i = 0; i < n; ++i) {
    part[i] -= min_p_no;
   }
  }
  delete (At);
  delete[]Ap;
  delete[]Ai;
  delete[]weigt;
  delete[]partIDX;
  //delete []unit_wgt;

  return 0;
 }

 int metis_partition_coarsened(CSC *A, int *&part, int k, int coarsen){
//  timing_measurement t; t.reset();
//  t.start_timer();
  CSC *tmp = coarsen_k_times(A->n, A->nnz, A->p, A->i, A->stype, coarsen);
  int *part_cors;
 // t.start_timer();
  if(A->stype < 0){
   metis_partition_symmetric(tmp, part_cors, k);
  } else{
   metis_partition_general(tmp, part_cors, k);
  }
  //t.start_timer();
  part = new int[A->n]();
  for (int i = 0; i < tmp->n; ++i) {
   int i2p = part_cors[i];
   int bnd = i+coarsen; // std::max<int>(i + coarsen, A->n);
   for (int j = i; j < bnd; ++j) {
    part[j] = i2p;
   }
  }
  //t.measure_elapsed_time();
  //t.print_t_array();
  delete part_cors;
  delete tmp;
  return 0;
 }

}
