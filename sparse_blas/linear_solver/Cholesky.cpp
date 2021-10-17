//
// Created by kazem on 10/17/21.
//

#include <sparse_blas_lib.h>
#include "Cholesky.h"

namespace sym_lib{

 Cholesky::Cholesky(CSC *A): LinearSolver(A){}

 Cholesky::~Cholesky(){}

 bool Cholesky::symbolic_analysis(){
  assert(_in_A->m == _in_A->n);
  bool sym_correct = true;//TODO: check for memory alloc's success
  _l_factor_info = build_symbolic_simplicial_lfactor(_in_A, _in_A_ordered,
                                                     _in_At_ordered);
  return sym_correct;
 }

 bool Cholesky::factorize(){
  bool is_spd;
  CSC *L = _l_factor_info->L;
  is_spd = cholesky_left_serial(_in_A->m, _in_A_ordered->p, _in_A_ordered->i,
                                _in_A_ordered->x,
                                _in_At_ordered->p, _in_At_ordered->i, L->p,
                                L->i, L->x, _l_factor_info->parent);
  return is_spd;
 }

 bool Cholesky::solve(Dense *b){
  bool is_solved = false;
  CSC *L = _l_factor_info->L;
  sptrsv_csc(_dim, L->p, L->i, L->x, _solution->a);
  ltsolve(_dim, L->p, L->i, L->x, _solution->a);
  return is_solved;
 }

}