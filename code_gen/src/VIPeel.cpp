//
// Created by george on 2020-02-16.
//

#include "VIPeel.h"

namespace Sympiler {
namespace Internal {

void VIPeel::visit(const For *for_loop) {
    if(for_loop->for_type == ForType::Peeled && num_peel > 0) {
        Stmt body = mutate(for_loop->body);
        std::string loopStr = for_loop->name;

        Stmt s = substitute(for_loop->name, for_loop->min, body);
        for(int i = 1; i < num_peel; i++) {
            Expr imm0 = IntImm::make(halide_type_t(halide_type_int, 32), i);
            Expr idx = Add::make(for_loop->min, imm0);
            Expr cond = LT::make(idx, for_loop->extent);
            Stmt si = substitute(for_loop->name, idx, body);

            Stmt ite = IfThenElse::make(cond, si);
            s = Block::make(s, ite);
        }

        Expr imm0 = IntImm::make(halide_type_t(halide_type_int, 32), num_peel);
        Expr new_min = Add::make(for_loop->min, imm0);
        Stmt peeledFor = For::make(for_loop->name, new_min, for_loop->extent, ForType::Serial, for_loop->device_api, body);

        stmt = Block::make(s, peeledFor);
    } else {
        IRMutator::visit(for_loop);
    }
}

}
}