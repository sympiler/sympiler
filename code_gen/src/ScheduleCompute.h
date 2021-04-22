//
// Created by kazem on 3/16/16.
//

#ifndef CPSL_HALIDE_SCHEDULECOMPUTE_H
#define CPSL_HALIDE_SCHEDULECOMPUTE_H

#include "IRMutator.h"
#include "Func.h"

namespace Sympiler {
namespace Internal {

class ScheduleCompute : public IRMutator {
    Func func;


public:
    ScheduleCompute(Stmt s, Func f){};
    //ScheduleCompute(Func f, Stmt s);
    virtual void visit(const Add *);

};

}
}


#endif //CPSL_HALIDE_SCHEDULECOMPUTE_H
