//
// Created by Shujian Qian on 2021-05-16.
//

#include "PeelIters.h"

#include "Substitute.h"

namespace Sympiler {
namespace Internal {

void
PeelIters::visit(const Sympiler::Internal::For *for_loop) {
    if (for_loop->for_type != ForType::Peeled) {
        IRMutator::visit(for_loop);
        return;
    }

    if (m_num_iters_to_peel > 0) {
        // If there are iterations to peel, peel them

        // The body of the for loop is not mutated since this visitor only peels
        // the outermost for loop
        Stmt body = for_loop->body;

        // Create for loop before the first iter to peel
        std::string name_pre = for_loop->name + "_pre";
        Expr var_pre = Variable::make(halide_type_t(halide_type_int, 32),
                                      name_pre);
        Stmt body_pre = substitute(for_loop->name, var_pre, body);
        Expr imm_ext_pre = IntImm::make(halide_type_t(halide_type_int, 32),
                                        m_iters_to_peel[0]);
        Stmt s = For::make(name_pre, for_loop->min, imm_ext_pre,
                           ForType::Serial, for_loop->device_api, body_pre);

        // Create peeled iterations and loops in between
        int peel_cnt;
        for (peel_cnt = 0; peel_cnt < m_num_iters_to_peel; ++peel_cnt) {
            Expr imm_iter = IntImm::make(halide_type_t(halide_type_int, 32),
                                         m_iters_to_peel[peel_cnt]);
            Stmt peeled_body = substitute(for_loop->name, imm_iter, body);
            s = Block::make(s, peeled_body);

            if (peel_cnt < m_num_iters_to_peel - 1 &&
                m_iters_to_peel[peel_cnt + 1] > m_iters_to_peel[peel_cnt] + 1) {
                std::string name_loop =
                        for_loop->name + "_" + std::to_string(peel_cnt);
                Expr var_loop = Variable::make(
                        halide_type_t(halide_type_int, 32), name_loop);
                Expr imm1 = IntImm::make(halide_type_t(halide_type_int, 32), 1);
                Expr imm_min_loop = Add::make(imm_iter, imm1);
                Expr imm_ext_loop = IntImm::make(
                        halide_type_t(halide_type_int, 32),
                        m_iters_to_peel[peel_cnt + 1]);

                Stmt body_loop = substitute(for_loop->name, var_loop, body);
                Stmt btw_loop = For::make(name_loop, imm_min_loop, imm_ext_loop,
                                          ForType::Serial, for_loop->device_api,
                                          body_loop);

                s = Block::make(s, btw_loop);
            }
        }

        // Create for loop after the last iteration to peel
        std::string name_aft = for_loop->name + "_aft";
        Expr var_aft = Variable::make(halide_type_t(halide_type_int, 32),
                                      name_aft);
        Stmt body_aft = substitute(for_loop->name, var_aft, body);
        Expr imm_min_aft = Add::make(
                IntImm::make(halide_type_t(halide_type_int, 32),
                             m_iters_to_peel[m_num_iters_to_peel - 1]),
                IntImm::make(halide_type_t(halide_type_int, 32), 1));
        Stmt for_aft = For::make(name_aft, imm_min_aft, for_loop->extent,
                                 ForType::Serial, for_loop->device_api,
                                 body_aft);

        stmt = Block::make(s, for_aft);
    } else {
        // If no iterations to peel, change the for_type to Serial and skip
        stmt = For::make(for_loop->name, for_loop->min, for_loop->extent,
                         ForType::Serial, for_loop->device_api, for_loop->body);
    }
}

}
}

