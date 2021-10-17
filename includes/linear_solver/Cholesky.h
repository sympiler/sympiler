//
// Created by kazem on 10/17/21.
//

#ifndef SYMPILER_PROJECT_SIMPLICIALLEFTCHOLESKY_H
#define SYMPILER_PROJECT_SIMPLICIALLEFTCHOLESKY_H

#include "LinearSolver.h"

namespace sym_lib {
 class Cholesky : public LinearSolver{

 protected:

 public:
  Cholesky(CSC *A);

  ~Cholesky();

  bool symbolic_analysis() override;

  bool factorize() override;

  bool solve(Dense *b) override;

 };
}

#endif //SYMPILER_PROJECT_SIMPLICIALLEFTCHOLESKY_H
