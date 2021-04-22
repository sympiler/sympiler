//
// Created by kazem on 5/2/17.
//

#ifndef SYMPILER_PROJ_LOWER_H
#define SYMPILER_PROJ_LOWER_H

#include "Expr.h"
#include "Kernel.h"

namespace Sympiler {
namespace Internal {
Stmt lower(Kernel *ker, SymbolicObject *sym);

}
}
#endif //SYMPILER_PROJ_LOWER_H
