//
// Created by kazem on 07/03/16.
//

#include <vector>
#include "IRMutator.h"
#include "IR.h"


namespace Sympiler{
namespace Internal{

Expr IRMutator::mutate(Expr e) {
    if (e.defined()) {
        e.accept(this);
    } else {
        expr = Expr();
    }
    //stmt = Stmt();
    return expr;
}

Stmt IRMutator::mutate(Stmt s) {
    if (s.defined()) {
        s.accept(this);
    } else {
        stmt = Stmt();
    }
    expr = Expr();
    return stmt;
}


namespace {
template<typename T>
void mutate_binary_operator(IRMutator *mutator, const T *op, Expr *expr, Stmt *stmt) {
    Expr a = mutator->mutate(op->a);
    Expr b = mutator->mutate(op->b);
    if (a.same_as(op->a) &&
        b.same_as(op->b)) {
        *expr = op;
    } else {
        *expr = T::make(a, b);
    }
    *stmt = NULL;
}
}
void IRMutator::visit(const IntImm *op)   {expr = op;}
void IRMutator::visit(const UIntImm *op)   {expr = op;}
void IRMutator::visit(const FloatImm *op) {expr = op;}
void IRMutator::visit(const StringImm *op) {expr = op;}

void IRMutator::visit(const Pointer *p) {
//    if(p->indirect)
//        p->idx->accept(this);
    Expr value;
    if(p->idx.defined())
        value = mutate(p->idx);
    if (value.same_as(p->idx)) {
        expr = p;
    } else {
        expr = Pointer::make(p->type, p->name, value);
    }

}


void IRMutator::visit(const Add *op)     {mutate_binary_operator(this, op, &expr, &stmt);}
void IRMutator::visit(const Sub *op)     {mutate_binary_operator(this, op, &expr, &stmt);}
void IRMutator::visit(const Mul *op)     {mutate_binary_operator(this, op, &expr, &stmt);}
void IRMutator::visit(const Div *op)     {mutate_binary_operator(this, op, &expr, &stmt);}
void IRMutator::visit(const Mod *op)     {mutate_binary_operator(this, op, &expr, &stmt);}
void IRMutator::visit(const Min *op)     {mutate_binary_operator(this, op, &expr, &stmt);}
void IRMutator::visit(const Max *op)     {mutate_binary_operator(this, op, &expr, &stmt);}
void IRMutator::visit(const EQ *op)      {mutate_binary_operator(this, op, &expr, &stmt);}
void IRMutator::visit(const NE *op)      {mutate_binary_operator(this, op, &expr, &stmt);}
void IRMutator::visit(const LT *op)      {mutate_binary_operator(this, op, &expr, &stmt);}
void IRMutator::visit(const LE *op)      {mutate_binary_operator(this, op, &expr, &stmt);}
void IRMutator::visit(const GT *op)      {mutate_binary_operator(this, op, &expr, &stmt);}
void IRMutator::visit(const GE *op)      {mutate_binary_operator(this, op, &expr, &stmt);}
void IRMutator::visit(const And *op)     {mutate_binary_operator(this, op, &expr, &stmt);}
void IRMutator::visit(const Or *op)      {mutate_binary_operator(this, op, &expr, &stmt);}

void IRMutator::visit(const Not *op) {
    Expr a = mutate(op->a);
    if (a.same_as(op->a)) expr = op;
    else expr = Not::make(a);
}

void IRMutator::visit(const Select *op)  {
    Expr cond = mutate(op->condition);
    Expr t = mutate(op->true_value);
    Expr f = mutate(op->false_value);
    if (cond.same_as(op->condition) &&
        t.same_as(op->true_value) &&
        f.same_as(op->false_value)) {
        expr = op;
    } else {
        expr = Select::make(cond, t, f);
    }
}

void IRMutator::visit(const Variable *op) {expr = op;}
void IRMutator::visit(const Load *op) {
    Expr index = mutate(op->ptr);
    if (index.same_as(op->ptr)) {
        expr = op;
    } else {
        expr = Load::make(op->type, op->name, op->ptr);
    }
}
void IRMutator::visit(const For *op) {
    Expr min = mutate(op->min);
    Expr extent = mutate(op->extent);
    Stmt body = mutate(op->body);
    if (min.same_as(op->min) &&
        extent.same_as(op->extent) &&
        body.same_as(op->body)) {
        stmt = op;
    } else {
        stmt = For::make(op->name, min, extent, op->for_type, op->device_api, body);
    }
}

void IRMutator::visit(const Call *op) {
    std::vector<Expr > new_args(op->args.size());
    bool changed = false;

    // Mutate the args
    for (size_t i = 0; i < op->args.size(); i++) {
        Expr old_arg = op->args[i];
        Expr new_arg = mutate(old_arg);
        if (!new_arg.same_as(old_arg)) changed = true;
        new_args[i] = new_arg;
    }

    if (!changed) {
        expr = op;
    } else {
        expr = Call::make(op->type,op->name, new_args);
    }
}

void IRMutator::visit(const CallX *op) {
    std::vector<Expr > new_args(op->args.size());
    std::vector<Argument > new_argums(op->argums.size());
    for (int j = 0; j < op->argums.size(); ++j) {
        new_argums[j] = op->argums[j];
    }

    bool changed = false;

    // Mutate the args
    for (size_t i = 0; i < op->args.size(); i++) {
        Expr old_arg = op->args[i];
        Expr new_arg = mutate(old_arg);
        if (!new_arg.same_as(old_arg)) changed = true;
        new_args[i] = new_arg;
    }

    if (!changed) {
        stmt = op;
    } else {
        stmt = CallX::make(op->name, new_args, new_argums);
    }
}

void IRMutator::visit(const Let *op) {
    Expr value = mutate(op->value);
    Expr body = mutate(op->body);
    if (value.same_as(op->value) &&
        body.same_as(op->body)) {
        expr = op;
    } else {
        expr = Let::make(op->name, value, body);
    }
}

void IRMutator::visit(const LetStmt *op) {
    Expr value = mutate(op->value);
    Stmt body = mutate(op->body);
    if (value.same_as(op->value) &&
        body.same_as(op->body)) {
        stmt = op;
    } else {
        stmt = LetStmt::make(op->name, value, body);
    }
}

void IRMutator::visit(const LetStmtPtr *op) {
    Expr lhs = mutate(op->lhs);
    Expr value = mutate(op->value);
    Stmt body = mutate(op->body);
    if (lhs.same_as(op->lhs) &&
        value.same_as(op->value) &&
        body.same_as(op->body)) {
        stmt = op;
    } else {
        stmt = LetStmtPtr::make(lhs, value, body);
    }
}

void IRMutator::visit(const Store *op) {
    Expr value = mutate(op->value);
    Expr index = mutate(op->ptr);
    if (value.same_as(op->value) && index.same_as(op->ptr)) {
        stmt = op;
    } else {
        stmt = Store::make(index,op->name, value);
    }
}

void IRMutator::visit(const Allocate *op) {
    std::vector<Expr> new_extents;
    bool all_extents_unmodified = true;
    for (size_t i = 0; i < op->extents.size(); i++) {
        new_extents.push_back(mutate(op->extents[i]));
        all_extents_unmodified &= new_extents[i].same_as(op->extents[i]);
    }
    Stmt body = mutate(op->body);
    Expr condition = mutate(op->condition);
    Expr new_expr;
    if (op->new_expr.defined()) {
        new_expr = mutate(op->new_expr);
    }
    if (all_extents_unmodified &&
        body.same_as(op->body) &&
        condition.same_as(op->condition) &&
        new_expr.same_as(op->new_expr)) {
        stmt = op;
    } else {
        stmt = Allocate::make(op->name, op->type, new_extents, condition, body, new_expr, op->free_function);
    }
}

void IRMutator::visit(const Free *op) {
    stmt = op;
}

void IRMutator::visit(const Evaluate *op) {
    Expr v = mutate(op->value);
    if (v.same_as(op->value)) {
        stmt = op;
    } else {
        stmt = Evaluate::make(v);
    }
}

void IRMutator::visit(const Block *op) {
    Stmt first = mutate(op->first);
    Stmt rest = mutate(op->rest);
    if (first.same_as(op->first) &&
        rest.same_as(op->rest)) {
        stmt = op;
    } else {
        stmt = Block::make(first, rest);
    }
}

void IRMutator::visit(const IfThenElse *op) {
    Expr condition = mutate(op->condition);
    Stmt then_case = mutate(op->then_case);
    Stmt else_case = mutate(op->else_case);
    if (condition.same_as(op->condition) &&
        then_case.same_as(op->then_case) &&
        else_case.same_as(op->else_case)) {
        stmt = op;
    } else {
        stmt = IfThenElse::make(condition, then_case, else_case);
    }
}
}
}