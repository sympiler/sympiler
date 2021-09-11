//
// Created by kazem on 04/30/21.
//


#include "cholesky_solver.h"

#include <algorithm>
#include <iostream>
#include <queue>
#include <common/Sym_BLAS.h>

#include "common/Norm.h"
#include "common/Transpose.h"
#include "common/TreeUtils.h"
#include "common/Util.h"
#include "gmres/mgmres.hpp"
#include "matrixVector/spmv_CSC.h"
#include "cholesky/parallel_PB_Cholesky_05.h"
#include "cholesky/sequential_PB_Cholesky.h"
#include "linear_solver/solve_phase.h"
#include "symbolic/symbolic_phase.h"



namespace sym_lib {
 namespace parsy {

  profiling_solver_info::profiling_solver_info(int nt) : fact_time(0), analysis_time(0),
                                                         solve_time(0), iter_time(0),
                                                         ordering_time(0), update_time(0),
                                                         piv_reord(0) {
   timing_chol = new double[4 + nt]();
  }

  profiling_solver_info::profiling_solver_info(int nt, double ft, double at, double st, double it, double ot, double ut,
                                               double pr) :
    fact_time(ft), analysis_time(at),
    solve_time(st), iter_time(it),
    ordering_time(ot), update_time(ut),
    piv_reord(pr) {
   timing_chol = new double[4 + nt]();
  }

  profiling_solver_info::~profiling_solver_info() {
   delete[]timing_chol;
  }

  std::chrono::time_point<std::chrono::system_clock> profiling_solver_info::tic() {
   return std::chrono::system_clock::now();
  }

  std::chrono::time_point<std::chrono::system_clock> profiling_solver_info::toc() {
   return std::chrono::system_clock::now();
  }

  double profiling_solver_info::elapsed_time(std::chrono::time_point<std::chrono::system_clock> beg,
                                             std::chrono::time_point<std::chrono::system_clock> lst) {
   double ret = 0;
   elapsed_seconds = lst - beg;
   ret = elapsed_seconds.count();
   return ret;
  }

  void profiling_solver_info::print_profiling() {
   std::cout << "analysis time: " << analysis_time << ";";
   std::cout << "fact time: " << fact_time << ";";
   std::cout << "update time: " << update_time << ";";
   std::cout << "reordering pivot time: " << piv_reord << ";";
   std::cout << "solve time: " << solve_time << ";";
  }

  SolverSettings::SolverSettings(CSC *Amat) {
   default_setting();
   A = Amat;
   rhs = NULL;
   A_ord = NULL;//new CSC;
   AT_ord = NULL; //new CSC;
   SM = new CSC;
   psi = new profiling_solver_info(num_thread);

   solver_mode = 0;//basic mode
   remove_trans = 0;
   x = NULL;
   B = NULL;
   C = NULL;
   a_consistent = 1;
   to_del = 0;
   visible_sn = NULL;
   level_ptr = NULL;
   par_ptr = NULL;
   par_set = NULL;
  }

  SolverSettings::SolverSettings(CSC *Amat, double *rhs_in) {
   default_setting();
   A = Amat;
   rhs = rhs_in;
   A_ord = NULL;//new CSC;
   AT_ord = NULL; //new CSC;
   SM = new CSC;
   psi = new profiling_solver_info(num_thread);

   solver_mode = 0;//basic mode
   remove_trans = 0;
   x = NULL;
   B = NULL;
   C = NULL;
   a_consistent = 1;
   to_del = 0;
   visible_sn = NULL;
   level_ptr = NULL;
   par_ptr = NULL;
   par_set = NULL;
  }

  SolverSettings::SolverSettings(CSC *Amat, double *rhs_in, CSC *Bmat) :
    SolverSettings(Amat, rhs_in) { //This an expert routine
   B = Bmat;
   BT = new CSC;
   transpose_unsym(B->nrow, B->ncol, B->p, B->i, B->x,
                   BT->nrow, BT->ncol, BT->p, BT->i, BT->x);
   BT->stype = B->stype;
   BT->xtype = B->xtype;
   BT->packed = B->packed;
   BT->nz = B->nz;
   BT->sorted = B->sorted;
   remove_trans = 1;
   solver_mode = 1;
  }

  SolverSettings::SolverSettings(CSC *Amat, double *rhs_in, CSC *Bmat, double *b) :
    SolverSettings(Amat, rhs_in, Bmat) { //This an expert routine
   extra_rhs = b;
  }

  SolverSettings::SolverSettings(CSC *Amat, double *rhs_in, CSC *Bmat, CSC *BTmat) :
    SolverSettings(Amat, rhs_in) { //This an expert routine
   B = Bmat;
   BT = BTmat;
   solver_mode = 1;
  }

  SolverSettings::SolverSettings(CSC *Amat, double *rhs_in, CSC *Bmat, CSC *BTmat, CSC *Cmat, CSC *CTmat) :
    SolverSettings(Amat, rhs_in, Bmat, BTmat) { //This an expert routine
   C = Cmat;
   CT = CTmat;
   solver_mode = 1;
  }

  SolverSettings::~SolverSettings() {
   allocateAC(A_ord, 0, 0, 0, FALSE);
   allocateAC(AT_ord, 0, 0, 0, FALSE);
   allocateLC(L, FALSE);
   if (simplicial_alloc) {
    delete[]L->p_s;
    delete[]L->x_s;
   }
   delete psi;
   delete[]ws;
 //  delete[]ws_int;
//   delete[]ws_zeroed;
   //delete []pinv;
   delete[]perm_piv;
   delete[]marked;
   delete[]valL;
   delete[]d_val;
   delete[]visible_cnt;
   delete []x;
   if (remove_trans) {
    allocateAC(BT, 0, 0, 0, FALSE);
   }

   if (solver_mode == 1 || solver_mode == 2) {
#ifdef SYM_REMOV
    delete []child_sn_ptr;
     delete []child_sn_no;
     delete []num_sn_child;
#endif
    delete[]atree;
    delete[]visible_sn;
#if 0
    delete []child_ptr;
     delete []child_no;
     delete []num_child;
     delete []etree_mod;
#endif
   }

   delete[]level_ptr;
   delete[]par_ptr;
   delete[]par_set;
   if (simplicial_alloc) {
    //delete []level_ptr_s;
    delete[]par_ptr_s;
    delete[]par_set_s;
   }

   if (solver_mode == 2) {
    delete[]sm_rhs;
   }
   if (solver_mode == 1) {
    delete[]sm_rhs;
    delete[]sm_solution;
   }
   delete[]extra_cols;

   if (num_thread > thread_thresh && solver_mode == 1) {
    delete[]s_level_ptr;
    delete[]s_level_set;
   }

   // deleting matrices
   delete[]L->ColCount;
   delete L;
   if (solver_mode != 0)
    allocateAC(SM, 0, 0, 0, FALSE);
   else
    delete SM;
  }

  void SolverSettings::default_setting() {
   sym_order = SYM_ORDER::S_AMD;
#ifdef OPENMP
   ldl_variant = 4;
#else
   ldl_variant = 1;
#endif
   ldl_update_variant = 2;
#ifdef OPENBLAS
   num_thread = openblas_get_num_procs();
   openblas_set_num_threads(1);
#else
   num_thread = mkl_get_max_threads();
   MKL_Domain_Set_Num_Threads(1, MKL_DOMAIN_BLAS);
#endif
   chunk = 1;
   cost_param = num_thread;
   level_param = -3;
   final_seq_node = 50;// should be greater than 1
   n_relax[0] = 4;
   n_relax[1] = 16;
   n_relax[2] = 48;
   z_relax[0] = 0.8;
   z_relax[1] = 0.1;
   z_relax[2] = 0.05;
   //z_relax[0]=0.9; z_relax[1]=0.5; z_relax[2]=0.05;
   //refinemen
   req_ref_iter = 2;
   max_iter = req_ref_iter;
   max_inner_iter = 2;
   tol_abs = 1e-15;
   tol_rel = 1e-15;
   //
   reg_diag = 1e-8;
   is_super = 1;
   regularization = 1;
   thread_thresh = 25;
   simplicial_alloc = 0;
   is_solved = is_factorized = false;
   x = NULL; in_perm = NULL;
  }

  int SolverSettings::build_super_matrix() {
   int status = 0;
   size_t SM_nz;
   int *SMp;
   int *SMi;
   double *SMx;
   size_t SM_size = 0;
   if (C == NULL && B->nrow == 0) {
    SM_size = A->ncol;
    SMp = new int[SM_size + 1]();
    SMp[0] = 0;
    for (int i = 1; i < A->ncol + 1; ++i) {
     SMp[i] = SMp[i - 1] + (A->p[i] - A->p[i - 1]);
    }
   } else if (C == NULL && B->nrow > 0) {
    SM_size = A->ncol + B->nrow;
    SMp = new int[SM_size + 1]();
    SMp[0] = 0;
    for (int i = 1; i < A->ncol + 1; ++i) {
     SMp[i] = SMp[i - 1] + (A->p[i] - A->p[i - 1]) +
              (B->p[i] - B->p[i - 1]);
    }
   } else if (C->nrow > 0 && B->nrow == 0) {
    SM_size = A->ncol + C->nrow;
    SMp = new int[SM_size + 1]();
    SMp[0] = 0;
    for (int i = 1; i < A->ncol + 1; ++i) {
     SMp[i] = SMp[i - 1] + (A->p[i] - A->p[i - 1]) +
              (C->p[i] - C->p[i - 1]);
    }
   } else {// both are not null
    SM_size = A->ncol + B->nrow + C->nrow;
    SMp = new int[SM_size + 1]();
    SMp[0] = 0;
    for (int i = 1; i < A->ncol + 1; ++i) {
     SMp[i] = SMp[i - 1] + (A->p[i] - A->p[i - 1]) +
              (B->p[i] - B->p[i - 1]) +
              (C->p[i] - C->p[i - 1]);
    }
   }
   //Adding diagonal for columns with zero values.
   for (int k = A->ncol + 1; k < SM_size + 1; ++k) {
    SMp[k] = SMp[k - 1] + 1;
   }
   SM_nz = SMp[SM_size];
   SMi = new int[SM_nz]();
   SMx = new double[SM_nz]();

   int base1 = A->ncol;
   size_t stp = 0;
   for (int j = 0; j < A->ncol; ++j) {
    stp = SMp[j];
    //Adding Hessian
    for (int i = A->p[j]; i < A->p[j + 1]; ++i) {
     SMi[stp] = A->i[i];
     SMx[stp] = A->x[i];
     stp++;
    }
    base1 = A->ncol;
    //Adding equalities
    if (C != NULL) {
     if (C->nrow > 0) {
      for (int i = C->p[j]; i < C->p[j + 1]; ++i) {
       SMi[stp] = base1 + C->i[i];
       SMx[stp] = C->x[i];
       //std::cout<<"Eq: "<< base1 + C->i[i]<<"; "<<SMx[stp]<<"\n";
       stp++;
      }
     }
     base1 = A->ncol + C->nrow;
    }
    //Adding inequalities
    if (B->nrow > 0) {
     for (int i = B->p[j]; i < B->p[j + 1]; ++i) {
      SMi[stp] = base1 + B->i[i];
      //SMx[stp] = B->x[i];
      stp++;
     }
    }
    assert(stp == SMp[j + 1]);
   }
   //Putting a small value in diagonals
   base1 = A->ncol;
   for (int l = SMp[A->ncol], j = 0; l < SM_nz; ++l, ++j) {
    //SMx[l] = 1e-6;
    SMi[l] = base1 + j;
   }

   SM->ncol = SM->nrow = SM_size;
   SM->p = SMp;
   SM->i = SMi;
   SM->x = SMx;
   SM->stype = -1;
   SM->xtype = SYMPILER_REAL;
   SM->packed = TRUE;
   SM->sorted = TRUE;
   SM->nzmax = SM_nz;
   sm_rhs = new double[SM->ncol]();
   sm_solution = new double[SM->ncol]();
//  CSC *sKKTt = ptranspose(SM,2,NULL,NULL,0,status);
//  print_csc("skkt\n",SM->ncol,SM->p,SM->i,SM->x);
   //print_csc("\nskkt Trans\n",sKKTt->ncol,sKKTt->p,sKKTt->i,sKKTt->x);
   //std::cout<<"\n";
   return status;
  }

  void SolverSettings::find_perturbation(double tol) {
   for (int i = 0; i < AorSM->ncol; ++i) {
    perturbed_value pv;
    if (std::abs(AorSM->x[AorSM->p[i]]) < tol) {
     pv.col_idx = i;
     pv.per_value = reg_diag;
     AorSM->x[AorSM->p[i]] = reg_diag;
     perturbed_diags.push_back(pv);
    }
   }
   // print_csc("skkt\n",AorSM->ncol,AorSM->p,AorSM->i,AorSM->x);
  }

  void SolverSettings::apply_perturbation(double tol) {
   for (int i = 0; i < A->ncol; ++i) {
    AorSM->x[AorSM->p[i]] += tol;
   }
   for (int i = A->ncol; i < AorSM->ncol; ++i) {
    AorSM->x[AorSM->p[i]] -= tol;
   }
  }

  void SolverSettings::add_perturbation(double tol) {
   for (int i = 0; i < A_ord->ncol; ++i) {
    int ord_col = L->IPerm[i];
    if (ord_col < A->ncol) {
     A_ord->x[A_ord->p[ord_col]] += tol;
    } else {
     A_ord->x[A_ord->p[ord_col]] -= tol;
    }
   }
  }

  void SolverSettings::remove_perturbation(double tol) {
   for (int i = 0; i < A_ord->ncol; ++i) {
    int ord_col = L->IPerm[i];
    if (ord_col < A->ncol) {
     A_ord->x[A_ord->p[ord_col]] -= tol;
    } else {
     A_ord->x[A_ord->p[ord_col]] += tol;
    }
   }
  }

  int SolverSettings::symbolic_analysis() {
   psi->start = psi->tic();

   if (solver_mode == 0) {
    AorSM = A;
    base = 0;
    //x = new double[AorSM->ncol];
    if (ldl_variant == 6 || ldl_variant == 7) {
     is_super = 0;
     simplicial_alloc = 1;
    }
   } else { // Mode 1
    is_super = 1; //FIXME: do it in analysis part
    simplicial_alloc = 0;//TODO: this is for low-rank updates
    if (ldl_variant == 6 || ldl_variant == 7
        || ldl_update_variant == 6 || ldl_update_variant == 7) {
     is_super = 0;
     simplicial_alloc = 1;
    }
    build_super_matrix();

    for (int i = 0; i < A->ncol; ++i) {
     sm_rhs[i] = rhs[i];
    }
    rhs = sm_rhs; //FIXME rhs should be read only!
    x = sm_solution;
    AorSM = SM;
    if (C == NULL)
     base = A->ncol;
    else
     base = A->ncol + C->nrow;

   }
   extra_cols = new int[AorSM->ncol]();
#ifdef SYM_REMOV
   if(solver_mode == 1){
     for (int j = base; j < AorSM->ncol; ++j) {
      extra_cols[j] = 1;
     }
    }
#endif
//  for (int j = 147; j < AorSM->ncol; ++j) {
//   extra_cols[j] = 1;
//  }
   //print_vec("Extra: ",0,AorSM->ncol,extra_cols);
   //Fill diagonals with perturbed value
   if(ldl_variant == 1)
    cost_param = 1;
   L = symbolic_analysis_lin_solve(1, AorSM, NULL, NULL, n_relax, z_relax,
                                   AorSM->ncol, prune_ptr, prune_set,
                                   n_level, level_ptr, level_set,
                                   n_par, par_ptr, par_set,
                                   n_level_s, level_ptr_s,
                                   n_par_s, par_ptr_s, par_set_s,
                                   cost_param, level_param, final_seq_node,
                                   status, max_sup_wid, max_col, psi->ordering_time,
                                   simplicial_alloc, extra_cols, in_perm, sym_order);
   psi->end = psi->toc();
   psi->analysis_time = psi->elapsed_time(psi->start, psi->end);
   if (L == NULL)
    return 0;
   //
   // print_vec("ordering: ",0,AorSM->ncol,L->Perm);
   //print_vec("inv ordering: ",0,AorSM->ncol,L->IPerm);
   //print_vec("supernode: ",0,L->nsuper,L->super);
   valL = new double[L->xsize]();
   d_val = new double[2 * AorSM->ncol]();
   visible_cnt = new int[L->nsuper]();
   allocate_workspace();
   //ws = new double[2*AorSM->ncol]();
   //ws_int = new int[3*AorSM->ncol]();
   //pinv = new int[AorSM->ncol];
   perm_piv = new int[AorSM->ncol]();
   marked = new bool[AorSM->ncol]();
   if (simplicial_alloc) {
    L->x_s = new double[L->xsize_s]();
   }
   AT_ord = ptranspose(AorSM, 2, L->Perm, NULL, 0, status);
   A_ord = ptranspose(AT_ord, 2, NULL, NULL, 0, status);

   etree = L->Parent;
   if (solver_mode == 0) {
    atree = L->sParent;
    if (simplicial_alloc) {
     etree_mod = new int[AorSM->ncol]();
     for (int k = 0; k < AorSM->ncol; ++k) {
      etree_mod[k] = L->Parent[k];
     }
     //print_vec("after: \n", 0, AorSM->ncol,etree_mod);
     l_pb = L->p_s; //new int [AorSM->ncol];
     l_pe = NULL; //new int[AorSM->ncol];
     l_i = L->i; //new int[L->xsize]();
     l_x = L->x_s; //new double[L->xsize]();
    }
    //print_vec("\nordering: \n", 0, AorSM->ncol,L->Perm);
   } else { // Mode 1
    //allocate for simplicial factor for low rank update
    l_pb = L->p_s; //new int [AorSM->ncol];
    l_pe = NULL; //new int[AorSM->ncol];
    l_i = L->i; //new int[L->xsize]();
    l_x = L->x_s; //new double[L->xsize]();
    atree = new int[L->nsuper]();
    visible_sn = new bool[L->nsuper]();
    for (int i = 0; i < L->nsuper; ++i) {
     atree[i] = L->sParent[i];
    }
    for (int l = 0; l < AorSM->ncol; ++l) {
     marked[l] = true;
    }
#ifdef SYM_REMOV
    //Compute the second representation of the tree
     child_sn_ptr = new int[L->nsuper+1];
     child_sn_no = new int[L->nsuper];
     num_sn_child = new int[L->nsuper]();
     populateChildren(L->nsuper,atree,child_sn_ptr,child_sn_no,num_sn_child);
     compressed_set_to_vector(L->nsuper,child_sn_ptr,child_sn_no,children_vec);
     // deleting extra rows symbolically, numerics are already set to zero
     //print_vec("before: \n", 0, L->nsuper,atree);
     //print_vec("\nordering: \n", 0, AorSM->ncol,L->Perm);
     //print_vec("\nsupernode bnds: \n", 0, L->nsuper+1,L->super);
     //print_vec("etree before: \n", 0, AorSM->ncol,L->Parent);
     //int hh = getTreeHeightBruteForce(L->nsuper,atree);
     //print_csc("\noriginal \n",A_ord->ncol,A_ord->p,A_ord->i,A_ord->x);
     //std::cout<<" : "<<hh<<"\n";
     for (int j = base; j < AorSM->ncol; ++j) {//hiding extra cols/rows
      delete_node_tree_simple(j);
     }
    // print_vec("after: \n", 0, L->nsuper,atree);
#endif
    for (int m = 0; m < L->nsuper; ++m) {
     if (marked[m]) {
      modified_sns.push_back(m);
     }
     visible_sn[m] = marked[m];
    }
    architecture_related_params();
   }
   return 1;
  }

  void SolverSettings::reset_symbolic_factor() {
   std::fill(valL, valL + L->xsize, 0);
   std::fill(d_val, d_val + 2 * AorSM->ncol, 0);
   std::fill(ws, ws + ws_dbl_size, 0);
  // std::fill(ws_int, ws_int + ws_int_size, 0);
 //  std::fill(ws_zeroed, ws_zeroed + (num_thread * AorSM->ncol), 0);

   //TODO reseting other symbolic info

  }

  void SolverSettings::architecture_related_params() {
   if (num_thread > thread_thresh) {
    s_level_ptr = new int[L->nsuper + 1]();
    s_level_set = new int[L->nsuper]();
    s_level_no = getLevelSet(L->nsuper, L->sParent, s_level_ptr, s_level_set);
   } else {
    s_level_no = -1;
   }
  }

  void SolverSettings::allocate_workspace() {
   //size_t ws_int_size;
   //size_t ws_dbl_size;
   // for update func
   //int workspace: 4*super_max + 2*n + 3*supNo
   // double workspace: 2 * super_max*col_max
   size_t upd_int = 4 * (max_sup_wid + 1) + 2 * AorSM->ncol + 3 * L->nsuper;
   size_t upd_dbl = 2 * (max_sup_wid + 1) * (max_col + 1);

   // temporary used
   size_t tmp_int = 3 * AorSM->ncol;
   size_t tmp_dbl = 2 * AorSM->ncol;

   ws_int_size = std::max(upd_int, tmp_int);
   ws_dbl_size = std::max(upd_dbl, tmp_dbl);

   // ws_size for solve_only= 2*A_ord->ncol
   // ws_size for triangular solves = num_thread*A_ord->ncol
   // ws_size for iter_ref = 4*(max_inner_iter+1) +
   // max_inner_iter*(max_inner_iter+1) + n + max_inner_iter*(n-1) +
   // solve phase
   size_t slve_dbl = 2 * AorSM->ncol;
   if (max_iter > 0) {
    slve_dbl += 4 * (max_inner_iter + 1) +
                (max_inner_iter + 1) * (max_inner_iter + 2) + 1
                + AorSM->ncol +
                (max_inner_iter + 1) * (AorSM->ncol) + 1;
   }
   ws_dbl_size = std::max(slve_dbl, ws_dbl_size);


   //allocating
   ws = new double[ws_dbl_size]();
   //ws_int = new int[ws_int_size]();
   //ws_zeroed = new double[num_thread * AorSM->ncol]();


  }

  int SolverSettings::numerical_factorization(sym_lib::parsy::CSC *A_in) {
   A = A_in;
   AorSM = A;
   base = 0;
   //x = new double[AorSM->ncol];
   is_factorized = false; is_solved = false;
   reset_symbolic_factor();
   //delete AT_ord;
   //delete A_ord;
   allocateAC(A_ord, 0, 0, 0, FALSE);
   allocateAC(AT_ord, 0, 0, 0, FALSE);
   AT_ord = ptranspose(AorSM, 2, L->Perm, NULL, 0, status);
   A_ord = ptranspose(AT_ord, 2, NULL, NULL, 0, status);
   etree = L->Parent;
   atree = L->sParent;
   auto ret = numerical_factorization();
   return ret;
  }

  int SolverSettings::numerical_factorization() {
   int ret_val = 0;
   if(is_factorized)
    return 1;
   //print_csc("\nORdered: ",A_ord->ncol,A_ord->p,A_ord->i,A_ord->x);
   switch (ldl_variant) {
    case 1:
//    MKL_Domain_Set_Num_Threads(num_thread, MKL_DOMAIN_BLAS);
     SET_BLAS_THREAD(num_thread);
     psi->start = psi->tic();
     ret_val = cholesky_left_sn_07(A_ord->ncol, A_ord->p, A_ord->i, A_ord->x,
                              L->p, L->s, L->i_ptr, valL,
                              L->super, L->nsuper, psi->timing_chol,
#ifndef PRUNE
                              atree, AT_ord->p, AT_ord->i, L->col2Sup,
#else
       prune_ptr,prune_set,
#endif
                              max_sup_wid + 1, max_col + 1);
     psi->end = psi->toc();
     psi->fact_time += psi->elapsed_time(psi->start, psi->end);
     //MKL_Domain_Set_Num_Threads(1, MKL_DOMAIN_BLAS);
     SET_BLAS_THREAD(1);
     is_factorized = true;
     break;
    case 2:
     //MKL_Domain_Set_Num_Threads(num_thread, MKL_DOMAIN_BLAS);
     SET_BLAS_THREAD(num_thread);
     psi->start = psi->tic();
/*     ret_val = ldl_left_sn_02_v2(A_ord->ncol, A_ord->p, A_ord->i, A_ord->x,
                                 L->p, L->s, L->i_ptr, valL,
                                 d_val,
                                 L->super, L->nsuper, psi->timing_chol,
#ifndef PRUNE
                                 atree, AT_ord->p, AT_ord->i, L->col2Sup,
#else
       prune_ptr,prune_set,
#endif
                                 max_sup_wid + 1, max_col + 1, num_pivot, perm_piv,
                                 L->sParent);*/
     psi->end = psi->toc();
     psi->fact_time += psi->elapsed_time(psi->start, psi->end);
     //MKL_Domain_Set_Num_Threads(1, MKL_DOMAIN_BLAS);
     SET_BLAS_THREAD(1);
     break;
    case 3://parallel static
     psi->start = psi->tic();
#ifdef OPENMP
/*    ret_val = ldl_left_sn_parallel_01(A_ord->ncol, A_ord->p, A_ord->i,
                                      A_ord->x, L->p, L->s, L->i_ptr, valL,
                                      d_val,
                                      L->super, L->nsuper, psi->timing_chol,
#ifndef PRUNE
                                      atree, AT_ord->p, AT_ord->i,
                                      L->col2Sup,
#else
      prune_ptr,prune_set,
#endif
                                      n_level, level_ptr, level_set,
                                      n_par, par_ptr, par_set,
                                      chunk, num_thread,
                                      max_sup_wid + 1, max_col + 1, num_pivot,
                                      reg_diag);*/
#endif

     psi->end = psi->toc();
     psi->fact_time += psi->elapsed_time(psi->start, psi->end);
     is_factorized = true;
     break;
    case 4://Parallel SBK
     psi->start = psi->tic();
#ifdef OPENMP
    ret_val = cholesky_left_par_05(A_ord->ncol, A_ord->p, A_ord->i, A_ord->x,
                                      L->p, L->s, L->i_ptr, valL,
                                      L->super, L->nsuper, psi->timing_chol,
#ifndef PRUNE
                                      atree, AT_ord->p, AT_ord->i, L->col2Sup,
#else
      prune_ptr,prune_set,
#endif
                                      n_level, level_ptr, level_set,
                                      n_par, par_ptr, par_set,
                                      chunk, num_thread,
                                      max_sup_wid + 1, max_col + 1);
#endif
     psi->end = psi->toc();
     psi->fact_time += psi->elapsed_time(psi->start, psi->end);
     //print_vec("simpl-o: ",0,A_ord->ncol,d_val);
     is_factorized = true;
     break;
    default:
     std::cout << " Wrong algorithm type! \n";
     return -1;
   }
   return ret_val;
  }


  double *SolverSettings::solve_only() {
   // workspace needed for solve:
   // ws_size for solve_only= 2*A_ord->ncol
   // ws_size for triangular solves = num_thread*A_ord->ncol
   // ws_size for iter_ref = 4*(max_inner_iter+1) +
   // max_inner_iter*(max_inner_iter+1) + n + max_inner_iter*(n-1) +
   //
   if(is_solved)
    return x;
   if (!rhs)
    return NULL;//rhs is not given
   psi->start = psi->toc();
   x = new double[A_ord->ncol]();
   double *x_ord = ws;
   double *rhs_ord = ws + A_ord->ncol;
   for (int i = 0; i < A_ord->ncol; ++i) {
    assert(L->Perm[i] < A_ord->ncol);
    x_ord[i] = rhs[L->Perm[i]];
   }
   //print_csc("\nA reo: \n",A_ord->ncol,A_ord->p,A_ord->i,A_ord->x);
   //print_vec<double >("\n\nrrhhssss\n: ",0,A_ord->ncol,rhs);
   //print_vec<double >("x reordered: ",0,A->ncol,x_ord);
   //print_vec<int >("ordering: ",0,A_ord->ncol,L->Perm);
   //print_vec<int >("inverse ordering: ",0,A->ncol,L->IPerm);

   if (ldl_variant == 1) {
    solve_phase_ll_blocked(A_ord->ncol, x_ord,
                                    L->col2Sup, L->super,
                                    L->p, L->s, valL, L->i_ptr,
                                    L->nsuper, L->nzmax);

   }  else if(ldl_variant == 4){
    solve_phase_ll_blocked_parallel(A_ord->ncol, x_ord,
                                              L->col2Sup, L->super,
                                              L->p, L->s, valL, L->i_ptr,
                                              L->nsuper, L->nzmax,
                                              n_level, level_ptr, level_set,
                                              n_par, par_ptr, par_set, chunk);
   } else{
    assert(false); // Wrong input
   }

   for (int i = 0; i < A_ord->ncol; ++i) {
    x[i] = x_ord[L->IPerm[i]];
   }
   //print_vec<double >("x orig order: ",0,A->ncol,x);
   //print_vec("dval \n",0,2*AorSM->ncol,d_val);
   psi->end = psi->toc();
   psi->solve_time += psi->elapsed_time(psi->start, psi->end);
   is_solved=true;
   return x;
  }


  double *SolverSettings::solve_only(const double *rhs_in, const int n_rhs) {
   // workspace needed for solve:
   // ws_size for solve_only= 2*A_ord->ncol
   // ws_size for triangular solves = num_thread*A_ord->ncol
   // ws_size for iter_ref = 4*(max_inner_iter+1) +
   // max_inner_iter*(max_inner_iter+1) + n + max_inner_iter*(n-1) +
   //
   const double *rhs = rhs_in;
   delete []x;
   x = new double[n_rhs*A_ord->ncol]();
   psi->start = psi->toc();
   double *x_ord = new double[n_rhs*A_ord->ncol]();
   double *rhs_ord = ws + A_ord->ncol;
#pragma omp parallel for
   for (int i = 0; i < A_ord->ncol; ++i) {
    assert(L->Perm[i] < A_ord->ncol);
    int tmp = L->Perm[i];
    for (int j = 0; j < n_rhs; ++j) {
     int idx = j *A_ord->ncol;
     x_ord[i + idx] = rhs[tmp + idx];
    }
   }
   //print_csc("\nA reo: \n",A_ord->ncol,A_ord->p,A_ord->i,A_ord->x);
   //print_vec<double >("\n\nrrhhssss\n: ",0,A_ord->ncol,rhs);
   //print_vec<double >("x reordered: ",0,A->ncol,x_ord);
   //print_vec<int >("ordering: ",0,A_ord->ncol,L->Perm);
   //print_vec<int >("inverse ordering: ",0,A->ncol,L->IPerm);

   if (ldl_variant == 1) {
    if(n_rhs == 1)
    solve_phase_ll_blocked(A_ord->ncol, x_ord,
                           L->col2Sup, L->super,
                           L->p, L->s, valL, L->i_ptr,
                           L->nsuper, L->nzmax);
    else
     solve_phase_ll_blocked_nrhs(A_ord->ncol, x_ord, n_rhs,
                            L->col2Sup, L->super,
                            L->p, L->s, valL, L->i_ptr,
                            L->nsuper, L->nzmax, max_col);


   }  else if(ldl_variant == 4){
    if(n_rhs == 1)
     solve_phase_ll_blocked_parallel(A_ord->ncol, x_ord,
                                     L->col2Sup, L->super,
                                     L->p, L->s, valL, L->i_ptr,
                                     L->nsuper, L->nzmax,
                                     n_level, level_ptr, level_set,
                                     n_par, par_ptr, par_set, chunk);
    else
     // TODO: replace the following with the parallel version
     solve_phase_ll_blocked_nrhs(A_ord->ncol, x_ord, n_rhs,
                                 L->col2Sup, L->super,
                                 L->p, L->s, valL, L->i_ptr,
                                 L->nsuper, L->nzmax, max_col);
/*     solve_phase_ll_blocked_parallel_nrhs(A_ord->ncol, x_ord, n_rhs,
                                     L->col2Sup, L->super,
                                     L->p, L->s, valL, L->i_ptr,
                                     L->nsuper, L->nzmax,
                                     n_level, level_ptr, level_set,
                                     n_par, par_ptr, par_set, chunk, max_col);*/
   } else{
    assert(false); // Wrong input
   }
#pragma omp parallel for
   for (int i = 0; i < A_ord->ncol; ++i) {
    int tmp = L->IPerm[i];
    for (int j = 0; j < n_rhs; ++j) {
     auto idx = j *A_ord->ncol;
     x[i + idx] = x_ord[tmp + idx];
    }
   }
   //print_vec<double >("x orig order: ",0,A->ncol,x);
   //print_vec("dval \n",0,2*AorSM->ncol,d_val);
   psi->end = psi->toc();
   psi->solve_time += psi->elapsed_time(psi->start, psi->end);
   delete []x_ord;
   return x;
  }


  void SolverSettings::compute_norms(double *rhs_in) {
   rhs = rhs_in != NULL ? rhs_in : rhs;
   double alp[2] = {1.0, 0};
   double bet[2] = {0.0, 0};
   int norm_type = 0;
   CSC *TMP = ptranspose(A_ord, 2, L->IPerm, NULL, 0, status);
   CSC *TMP2 = ptranspose(TMP, 2, NULL, NULL, 0, status);
   double *res = new double[TMP2->ncol]();
   x_l1 = norm_dense(1, TMP2->ncol, x, norm_type);
   rhs_l1 = norm_dense(1, TMP2->ncol, rhs, norm_type);
   spmv_csc_sym_one_int(TMP2->ncol, TMP2->p, TMP2->i, TMP2->x, -1, alp, bet,
                        1, x, res);
   //print_vec("res mult: ",0,TMP2->ncol,res);
   for (int i = 0; i < TMP2->ncol; ++i) {
    res[i] = rhs[i] - res[i];
   }
   //print_vec("res: ",0,TMP2->ncol,res);
   res_l1 = norm_dense(1, TMP2->ncol, res, norm_type);
   A_l1 = norm_sparse_int(TMP2->ncol, TMP2->p, TMP2->i, TMP2->x, -1, norm_type);
   delete[]res;
   allocateAC(TMP, 0, 0, 0, FALSE);
   allocateAC(TMP2, 0, 0, 0, FALSE);
  }

  double SolverSettings::backward_error() {
   compute_norms();
   if (A_l1 * x_l1 + rhs_l1 > 0)
    bwd_err = res_l1 / (A_l1 * x_l1 + rhs_l1);
   else
    bwd_err = 0;
   std::cout << "d: " << res_l1 << "\n";
   std::cout << "A l1: " << A_l1 << "\n";
   std::cout << "x l1: " << x_l1 << "\n";
   std::cout << "rhs l1: " << rhs_l1 << "\n";
   return bwd_err;
  }

  void SolverSettings::convert_supernode_to_simplicial() {
   size_t actualNNZ = 0;
   for (int i = 0; i < L->nsuper; ++i) {
    int curCol = L->super[i];
    int nxtCol = L->super[i + 1];
    int supWdt = nxtCol - curCol;
    assert(supWdt > 0);
    for (int j = curCol; j < nxtCol; ++j) {
     l_pb[j] = L->p[j] + (j - curCol);
     l_pe[j] = L->p[j + 1];
     for (int k = l_pb[j],
            kk = L->i_ptr[curCol] + (j - curCol);
          k < l_pe[j]; ++k, ++kk) { // copy row indices
      l_i[k] = L->i[kk];
     }
    }
   }
  }

  void SolverSettings::copy_lsuper_to_l() {

  }

  void SolverSettings::set_diags(double d) {
   for (int i = 0; i < perturbed_diags.size(); ++i) {
    int cc = perturbed_diags[i].col_idx;
    int ord_col = L->IPerm[cc];
    double pv = d == 0 ? 0 : perturbed_diags[i].per_value;
    A_ord->x[A_ord->p[ord_col]] = pv;
   }
/*  for (int i = A->ncol; i < A_ord->ncol; ++i) {
   int ord_col = L->IPerm[i];
   //std::cout<<"Diag: "<<A_ord->x[A_ord->p[ord_col]]<<","<<A_ord->i[A_ord->p[ord_col]]<<"\n";
   A_ord->x[A_ord->p[ord_col]] = d;
  }*/
  }

  int SolverSettings::check_ldlt_factor() { // TODO move it to unit test later

   size_t *ia = new size_t[A_ord->ncol + 1]();
   int *ja = new int[L->xsize]();
   double *a = new double[L->xsize]();
   size_t *A2p = new size_t[A_ord->ncol + 1]();
   for (int m = 0; m <= A_ord->ncol; ++m) {
    A2p[m] = static_cast<size_t >(A_ord->p[m]);
   }
   //Converting BCSC to CSC
   bcsc2csc_aggressive(A_ord->ncol, L->nsuper, L->p, L->s, L->i_ptr,
                       L->super, valL, ia, ja, a);
   int *tmp = new int[A_ord->ncol + 1]();
   for (int m = 0; m <= A_ord->ncol; ++m) {
    tmp[m] = ia[m];
   }
   print_csc("\n%%MatrixMarket matrix coordinate real general\n",
             A_ord->ncol, tmp, ja, a);
   //print_vec("\nD: ",0,2*A_ord->ncol,d_val);
   print_csc("\n%%MatrixMarket matrix coordinate real symmetric\n",
             A_ord->ncol, A_ord->p, A_ord->i, A_ord->x);
   delete[]tmp;
   double nrm = norm_sparse(A_ord->ncol, ia, ja, a, -1, 1);
   bool check = true; //ldlt_check(A_ord->ncol, ia, ja, a, d_val, A2p,
//    A_ord->i,A_ord->x);
   delete[]ia;
   delete[]ja;
   delete[]a;
   delete[]A2p;

   return check;
  }
 }
}