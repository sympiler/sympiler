//
// Created by kazem on 10/17/21.
//

#ifndef SYMPILER_PROJECT_LINEARSOLVER_H
#define SYMPILER_PROJECT_LINEARSOLVER_H

#include <def.h>
#include <lfactor_creation.h>

namespace sym_lib {

 class LinearSolver {

 protected:
  int _dim;
  CSC *_in_A;
  CSC *_in_A_ordered, *_in_At_ordered;
  LFactorSymbolicInfo *_l_factor_info;
  Dense *_rhs, *_solution;

 public:
  LinearSolver(CSC *A);

  ~LinearSolver();

  virtual bool symbolic_analysis() = 0;

  virtual bool factorize() = 0;

  virtual bool solve(sym_lib::Dense *b) = 0;


 };

}
#endif //SYMPILER_PROJECT_LINEARSOLVER_H
