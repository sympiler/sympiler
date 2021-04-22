//
// Created by kazem on 17/02/16.
//

#ifndef DSLPROJECT_OPTIMIZATION_H
#define DSLPROJECT_OPTIMIZATION_H

#include "Expr.h"
#include "Func.h"
namespace Sympiler {
namespace Internal {
Stmt Optimization(Func &output);
}
}

#endif //DSLPROJECT_OPTIMIZATION_H
