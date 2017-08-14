//
// Created by kazem on 5/8/17.
//

#ifndef SYMPILER_PROJ_VIPRUNE_H
#define SYMPILER_PROJ_VIPRUNE_H

#include "IRMutator.h"
#include "Substitute.h"

namespace Sympiler{
namespace Internal {
class VIPrune : public IRMutator {
    std::string setVal;
    Expr lb, ub;
    void visit(const For *for_loop){
        if(for_loop->for_type == ForType::Pruned){
            Stmt body = mutate(for_loop->body);
            std::string loopStr=for_loop->name+"P";
            Expr newIdx = Pointer::make(halide_type_t(halide_type_int,32),setVal,
                                        Variable::make(halide_type_t(halide_type_int,32),
                                                       loopStr));
            Stmt newbody = substitute(for_loop->name,newIdx,body); //substitute current for loop name with the new index
            Stmt prunedFor = For::make(for_loop->name+"P", lb, ub,
                                       ForType::Serial, for_loop->device_api, newbody);

            stmt = prunedFor;
        } else{
            IRMutator::visit(for_loop);
        }

        /*Expr min = mutate(for_loop->min);
        Expr extent = mutate(for_loop->extent);
        Stmt body = mutate(for_loop->body);
        if (min.same_as(for_loop->min) &&
            extent.same_as(for_loop->extent) &&
            body.same_as(for_loop->body)) {
            stmt = for_loop;
        } else {
            stmt = For::make(for_loop->name, min, extent, for_loop->for_type, for_loop->device_api, body);
        }*/

    }

public:
    VIPrune(Expr Lb, std::string val, Expr Ub):
            lb(Lb),setVal(val), ub(Ub){}
};

}
}


#endif //SYMPILER_PROJ_VIPRUNE_H
