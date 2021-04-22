//
// Created by kazem on 07/03/16.
//

#ifndef CPSL_HALIDE_IRMUTATOR_H
#define CPSL_HALIDE_IRMUTATOR_H

#include "IRVisitor.h"

namespace Sympiler{
namespace Internal{
/** A base class for passes over the IR which modify it
 * (e.g. replacing a variable with a value (Substitute.h), or
 * constant-folding).
 *
 * Your mutate should override the visit methods you care about. Return
 * the new expression by assigning to expr or stmt. The default ones
 * recursively mutate their children. To mutate sub-expressions and
 * sub-statements you should the mutate method, which will dispatch to
 * the appropriate visit method and then return the value of expr or
 * stmt after the call to visit.
 */
class IRMutator: public IRVisitor {
public:
    virtual Expr mutate(Expr expr);
    virtual Stmt mutate(Stmt stmt);
protected:
    Expr expr;
    Stmt stmt;

    virtual void visit(const Variable *v);
    virtual void visit(const Pointer *v);
    virtual void visit(const IntImm *i);
    virtual void visit(const UIntImm *);
    virtual void visit(const FloatImm *);
    virtual void visit(const StringImm *);
    virtual void visit(const Add *op);
    virtual void visit(const Sub *);
    virtual void visit(const Mul *);
    virtual void visit(const Div *);
    virtual void visit(const Mod *);
    virtual void visit(const Min *);
    virtual void visit(const Max *);
    virtual void visit(const EQ *);
    virtual void visit(const NE *);
    virtual void visit(const LT *);
    virtual void visit(const LE *);
    virtual void visit(const GT *);
    virtual void visit(const GE *);
    virtual void visit(const And *);
    virtual void visit(const Or *);
    virtual void visit(const Not *);
    virtual void visit(const Select *);
    virtual void visit(const For *op);
    virtual void visit(const Let *op);
    virtual void visit(const LetStmt *op);
    virtual void visit(const LetStmtPtr *op);
    virtual void visit(const Evaluate *op);
    virtual void visit(const Block *op);
    virtual void visit(const Load *op);
    virtual void visit(const Store *op);
    virtual void visit(const Call *op);
    virtual void visit(const CallX *op);
    virtual void visit(const Allocate *);
    virtual void visit(const Free *);
    virtual void visit(const IfThenElse *);
};

}
}


#endif //CPSL_HALIDE_IRMUTATOR_H
