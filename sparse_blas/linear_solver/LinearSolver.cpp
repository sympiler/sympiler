//
// Created by kazem on 10/17/21.
//

#include "LinearSolver.h"

namespace sym_lib{

 LinearSolver::LinearSolver(CSC *A){
  _in_A = A; _rhs=NULLPNTR;
  _in_A_ordered = NULLPNTR; _in_At_ordered=NULLPNTR;
  _dim = _in_A->m;
 }

 LinearSolver::~LinearSolver(){
  delete _l_factor_info;
  delete _solution;
  delete _in_At_ordered;
  delete _in_A_ordered;
 }

}