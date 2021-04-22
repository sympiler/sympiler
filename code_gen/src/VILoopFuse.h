//
// Created by george on 2020-02-19.
//

#ifndef FUSION_VILOOPFUSE_H
#define FUSION_VILOOPFUSE_H

#include "IREquality.h"
#include "IRMutator.h"
#include "Substitute.h"

namespace Sympiler {
 namespace Internal {
  class VILoopFuse : IRMutator {

   void visit(const For *loop1, const For *loop2) {
    // Check for loop bounds first
    Stmt body1 = mutate(loop1->body);
    Stmt body2 = mutate(loop2->body);

    Expr min1 = loop1->min;
    Expr min2 = loop2->min;
    Expr ext1 = loop1->extent;
    Expr ext2 = loop2->extent;

    if(!equal(min1, min2) || !equal(ext1, ext2)) {
     Stmt sep_loop1 = For::make(loop1->name, min1, ext1, ForType::Serial, loop1->device_api, body1);
     Stmt sep_loop2 = For::make(loop2->name, min2, ext2, ForType::Serial, loop2->device_api, body2);
     stmt = Block::make(sep_loop1, sep_loop1);
     return;
    }

    // TODO: check dependency

    // Merge naive bodies
    /* for (i = 0; i < n; i++)
     *  y[i] = f(x)
     * for (i = 0; i < n; i++)
     *  z[i] = g(x)
     */
    std::string loop_str = loop1->name + "_" + loop2->name + "f";
    Expr idx = Variable::make(halide_type_t(halide_type_int, 32), loop_str);
    Stmt newbody1 = substitute(loop1->name, idx, body1);
    Stmt newbody2 = substitute(loop2->name, idx, body2);
    Stmt fused_bdy = Block::make(newbody1, newbody2);

    Stmt fused_loop = For::make(loop_str, min1, ext1, ForType::Serial, loop1->device_api, fused_bdy);
    stmt = fused_loop;

    // TODO: expand this to
    /* for (i = 0; i < n; i++)
     *  y[i] = f(x)
     * z = y
     * for (i = 0; i < n; i++)
     *  w[i] = g(z)
     *
     * To fuse this, need w[i] depend only on z[0:i]
     * Maybe to expand this with loop peeling if
     *  w[i] depends on z[0:i+k]
     */
   }
  };
 }
}


#endif //FUSION_VILOOPFUSE_H
