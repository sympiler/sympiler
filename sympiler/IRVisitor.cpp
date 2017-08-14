//
// Created by kazem on 12/02/16.
//

#include "IRVisitor.h"

namespace Sympiler {
namespace Internal {
IRVisitor::~IRVisitor() {
}

void IRVisitor::visit(const Variable *v) {

}

void IRVisitor::visit(const Pointer *p) {
//    if(p->indirect)
//        p->idx->accept(this);
}
void IRVisitor::visit(const IntImm *i) {

}

void IRVisitor::visit(const UIntImm *) {
}

void IRVisitor::visit(const FloatImm *) {
}

void IRVisitor::visit(const StringImm *) {
}

void IRVisitor::visit(const Add *op) {
    op->a.accept(this);
    op->b.accept(this);
}

void IRVisitor::visit(const Div *op) {
    op->a.accept(this);
    op->b.accept(this);
}

void IRVisitor::visit(const Sub *op) {
    op->a.accept(this);
    op->b.accept(this);
}

void IRVisitor::visit(const Mul *op) {
    op->a.accept(this);
    op->b.accept(this);
}


void IRVisitor::visit(const Mod *op) {
    op->a.accept(this);
    op->b.accept(this);
}

void IRVisitor::visit(const Min *op) {
    op->a.accept(this);
    op->b.accept(this);
}

void IRVisitor::visit(const Max *op) {
    op->a.accept(this);
    op->b.accept(this);
}

void IRVisitor::visit(const EQ *op) {
    op->a.accept(this);
    op->b.accept(this);
}

void IRVisitor::visit(const NE *op) {
    op->a.accept(this);
    op->b.accept(this);
}

void IRVisitor::visit(const LT *op) {
    op->a.accept(this);
    op->b.accept(this);
}

void IRVisitor::visit(const LE *op) {
    op->a.accept(this);
    op->b.accept(this);
}

void IRVisitor::visit(const GT *op) {
    op->a.accept(this);
    op->b.accept(this);
}

void IRVisitor::visit(const GE *op) {
    op->a.accept(this);
    op->b.accept(this);
}

void IRVisitor::visit(const And *op) {
    op->a.accept(this);
    op->b.accept(this);
}

void IRVisitor::visit(const Or *op) {
    op->a.accept(this);
    op->b.accept(this);
}

void IRVisitor::visit(const Not *op) {
    op->a.accept(this);
}

void IRVisitor::visit(const Select *op) {
    op->condition.accept(this);
    op->true_value.accept(this);
    op->false_value.accept(this);
}

void IRVisitor::visit(const For *op) {
    op->min.accept(this);
    op->extent.accept(this);
    op->body.accept(this);
}

void IRVisitor::visit(const Let *op) {
    op->value.accept(this);
    op->body.accept(this);
}

void IRVisitor::visit(const LetStmt *op) {
    op->value.accept(this);
    op->body.accept(this);
}

void IRVisitor::visit(const LetStmtPtr *op) {
    op->value.accept(this);
    op->body.accept(this);
}

void IRVisitor::visit(const Evaluate *op) {
    op->value.accept(this);
}

void IRVisitor::visit(const Block *op){
    op->first.accept(this);
    if (op->rest.defined()) {
        op->rest.accept(this);
    }
}

void IRVisitor::visit(const Load *op){
    op->ptr.accept(this);
}

void IRVisitor::visit(const Store *op){
    op->value.accept(this);
    op->ptr.accept(this);
}

void IRVisitor::visit(const Call *op) {
    for (size_t i = 0; i < op->args.size(); i++) {
        op->args[i].accept(this);
    }

    // Consider extern call args
/*    Function f = op->func;
    if (op->call_type == Call::Halide && f.has_extern_definition()) {
        for (size_t i = 0; i < f.extern_arguments().size(); i++) {
            ExternFuncArgument arg = f.extern_arguments()[i];
            if (arg.is_expr()) {
                arg.expr.accept(this);
            }
        }
    }*/
}

void IRVisitor::visit(const CallX *op){
    for (size_t i = 0; i < op->args.size(); i++) {
        op->args[i].accept(this);
    }

}

void IRVisitor::visit(const Allocate *op) {
    for (size_t i = 0; i < op->extents.size(); i++) {
        op->extents[i].accept(this);
    }
    op->condition.accept(this);
    if (op->new_expr.defined()) {
        op->new_expr.accept(this);
    }
    op->body.accept(this);
}

void IRVisitor::visit(const Free *op) {
}

void IRVisitor::visit(const IfThenElse *op) {
    op->condition.accept(this);
    op->then_case.accept(this);
    if (op->else_case.defined()) {
        op->else_case.accept(this);
    }
}

}
}