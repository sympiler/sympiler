//
// Created by kazem on 07/02/16.
//

#include "IR.h"
#include "IRVisitor.h"
#include "Expr.h"

namespace Sympiler {
namespace Internal {

Expr Add::make(Expr a, Expr b) {

    Add *node = new Add;
    node->type = a.type();
    node->a = a;
    node->b = b;
    return node;
}

Expr Sub::make(Expr a, Expr b) {
    //internal_assert(a.defined()) << "Sub of undefined\n";
    //internal_assert(b.defined()) << "Sub of undefined\n";
    //internal_assert(a.type() == b.type()) << "Sub of mismatched types\n";

    Sub *node = new Sub;
    node->type = a.type();
    node->a = a;
    node->b = b;
    return node;
}

Expr Mul::make(Expr a, Expr b) {
  //  internal_assert(a.defined()) << "Mul of undefined\n";
  //  internal_assert(b.defined()) << "Mul of undefined\n";
  //  internal_assert(a.type() == b.type()) << "Mul of mismatched types\n";

    Mul *node = new Mul;
    node->type = a.type();
    node->a = a;
    node->b = b;
    return node;
}

Expr Div::make(Expr a, Expr b) {
   // internal_assert(a.defined()) << "Div of undefined\n";
   // internal_assert(b.defined()) << "Div of undefined\n";
   // internal_assert(a.type() == b.type()) << "Div of mismatched types\n";

    Div *node = new Div;
    node->type = a.type();
    node->a = a;
    node->b = b;
    return node;
}

Expr Mod::make(Expr a, Expr b) {
/*
    internal_assert(a.defined()) << "Mod of undefined\n";
    internal_assert(b.defined()) << "Mod of undefined\n";
    internal_assert(a.type() == b.type()) << "Mod of mismatched types\n";
*/

    Mod *node = new Mod;
    node->type = a.type();
    node->a = a;
    node->b = b;
    return node;
}

Expr Min::make(Expr a, Expr b) {
/*
    internal_assert(a.defined()) << "Min of undefined\n";
    internal_assert(b.defined()) << "Min of undefined\n";
    internal_assert(a.type() == b.type()) << "Min of mismatched types\n";
*/

    Min *node = new Min;
    node->type = a.type();
    node->a = a;
    node->b = b;
    return node;
}

Expr Max::make(Expr a, Expr b) {
/*
    internal_assert(a.defined()) << "Max of undefined\n";
    internal_assert(b.defined()) << "Max of undefined\n";
    internal_assert(a.type() == b.type()) << "Max of mismatched types\n";
*/

    Max *node = new Max;
    node->type = a.type();
    node->a = a;
    node->b = b;
    return node;
}

Expr EQ::make(Expr a, Expr b) {
/*
    internal_assert(a.defined()) << "EQ of undefined\n";
    internal_assert(b.defined()) << "EQ of undefined\n";
    internal_assert(a.type() == b.type()) << "EQ of mismatched types\n";
*/

    EQ *node = new EQ;
    node->type = Bool(a.type().lanes());
    node->a = a;
    node->b = b;
    return node;
}

Expr NE::make(Expr a, Expr b) {
/*
    internal_assert(a.defined()) << "NE of undefined\n";
    internal_assert(b.defined()) << "NE of undefined\n";
    internal_assert(a.type() == b.type()) << "NE of mismatched types\n";
*/

    NE *node = new NE;
    node->type = Bool(a.type().lanes());
    node->a = a;
    node->b = b;
    return node;
}

Expr LT::make(Expr a, Expr b) {
/*
    internal_assert(a.defined()) << "LT of undefined\n";
    internal_assert(b.defined()) << "LT of undefined\n";
    internal_assert(a.type() == b.type()) << "LT of mismatched types\n";
*/

    LT *node = new LT;
    node->type = Bool(a.type().lanes());
    node->a = a;
    node->b = b;
    return node;
}


Expr LE::make(Expr a, Expr b) {
/*
    internal_assert(a.defined()) << "LE of undefined\n";
    internal_assert(b.defined()) << "LE of undefined\n";
    internal_assert(a.type() == b.type()) << "LE of mismatched types\n";
*/

    LE *node = new LE;
    node->type = Bool(a.type().lanes());
    node->a = a;
    node->b = b;
    return node;
}

Expr GT::make(Expr a, Expr b) {
/*
    internal_assert(a.defined()) << "GT of undefined\n";
    internal_assert(b.defined()) << "GT of undefined\n";
    internal_assert(a.type() == b.type()) << "GT of mismatched types\n";
*/

    GT *node = new GT;
    node->type = Bool(a.type().lanes());
    node->a = a;
    node->b = b;
    return node;
}


Expr GE::make(Expr a, Expr b) {
/*
    internal_assert(a.defined()) << "GE of undefined\n";
    internal_assert(b.defined()) << "GE of undefined\n";
    internal_assert(a.type() == b.type()) << "GE of mismatched types\n";
*/

    GE *node = new GE;
    node->type = Bool(a.type().lanes());
    node->a = a;
    node->b = b;
    return node;
}

Expr And::make(Expr a, Expr b) {
/*
    internal_assert(a.defined()) << "And of undefined\n";
    internal_assert(b.defined()) << "And of undefined\n";
    internal_assert(a.type().is_bool()) << "lhs of And is not a bool\n";
    internal_assert(b.type().is_bool()) << "rhsDense of And is not a bool\n";
    internal_assert(a.type() == b.type()) << "And of mismatched types\n";
*/

    And *node = new And;
    node->type = Bool(a.type().lanes());
    node->a = a;
    node->b = b;
    return node;
}

Expr Or::make(Expr a, Expr b) {
/*      internal_assert(a.defined()) << "Or of undefined\n";
    internal_assert(b.defined()) << "Or of undefined\n";
    internal_assert(a.type().is_bool()) << "lhs of Or is not a bool\n";
    internal_assert(b.type().is_bool()) << "rhsDense of Or is not a bool\n";
    internal_assert(a.type() == b.type()) << "Or of mismatched types\n";
*/
    Or *node = new Or;
    node->type = Bool(a.type().lanes());
    node->a = a;
    node->b = b;
    return node;
}

Expr Not::make(Expr a) {
/*        internal_assert(a.defined()) << "Not of undefined\n";
    internal_assert(a.type().is_bool()) << "argument of Not is not a bool\n";*/

    Not *node = new Not;
    node->type = Bool(a.type().lanes());
    node->a = a;
    return node;
}

Expr Select::make(Expr condition, Expr true_value, Expr false_value) {
 /*   internal_assert(condition.defined()) << "Select of undefined\n";
    internal_assert(true_value.defined()) << "Select of undefined\n";
    internal_assert(false_value.defined()) << "Select of undefined\n";
    internal_assert(condition.type().is_bool()) << "First argument to Select is not a bool: " << condition.type() << "\n";
    internal_assert(false_value.type() == true_value.type()) << "Select of mismatched types\n";
    internal_assert(condition.type().is_scalar() ||
                    condition.type().lanes() == true_value.type().lanes())
    << "In Select, vector lanes of condition must either be 1, or equal to vector lanes of arguments\n";
*/
    Select *node = new Select;
    node->type = true_value.type();
    node->condition = condition;
    node->true_value = true_value;
    node->false_value = false_value;
    return node;
}

Expr Variable::make(Type type, std::string name) {
    Variable *node = new Variable;
    node->type = type;
    node->name = name;
    return node;
}

Expr Pointer::make(Type type, std::string name, Expr idx) {
    Pointer *node = new Pointer;
    node->name = name;
    node->type = type;
    if(idx.defined())
        node->indirect=true;//indirect memory access
    else
        node->indirect= false;
    node->idx = idx;
    return node;
}
//bool Pointer::defined() { return name!=""; }

Stmt For::make(std::string name, Expr min, Expr extent, ForType for_type, DeviceAPI device_api, Stmt body) {
    For *node = new For;
    node->name = name;
    node->min = min;
    node->extent = extent;
    node->for_type = for_type;
    node->device_api = device_api;
    node->body = body;
    return node;
}

Expr Let::make(std::string name, Expr value, Expr body) {

    Let *node = new Let;
//    node->type = body.type();
    node->name = name;
    node->value = value;
    node->body = body;
    return node;
}

Stmt LetStmt::make(std::string name, Expr value, Stmt body) {
    LetStmt *node = new LetStmt;
    node->name = name;
    node->value = value;
    node->body = body;
    return node;
}

Stmt LetStmtPtr::make(Expr lhs, Expr value, Stmt body) {
    LetStmtPtr *node = new LetStmtPtr;
    node->lhs = lhs;
    node->value = value;
    node->body = body;
    return node;
}

Stmt Evaluate::make(Expr v) {
    Evaluate *node = new Evaluate;
    node->value = v;
    return node;
}

Stmt Block::make(const Stmt &first, const Stmt &rest) {
  //  internal_assert(first.defined()) << "Block of undefined\n";
  //  internal_assert(rest.defined()) << "Block of undefined\n";

    Block *node = new Block;
    node->first = first;
    node->rest = rest;
    return node;
}

Stmt Block::make(const std::vector<Stmt> &stmts) {
    if (stmts.empty()) {
        return Stmt();
    }
    Stmt result = stmts.back();
    for (size_t i = stmts.size()-1; i > 0; i--) {
        result = Block::make(stmts[i-1], result);
    }
    return result;
}

Expr Load::make(Type type, std::string name, Expr ptr) {
    //internal_assert(index.defined()) << "Load of undefined\n";
    //internal_assert(type.lanes() == index.type().lanes()) << "Vector lanes of Load must match vector lanes of index\n";

    Load *node = new Load;
    node->type = type;
    node->ptr = ptr;
    node->name = name;
    //node->image = image;
    //node->param = param;
    return node;
}

Stmt Store::make(Expr ptr, std::string name, Expr value) {
//    internal_assert(value.defined()) << "Store of undefined\n";
//    internal_assert(index.defined()) << "Store of undefined\n";

    Store *node = new Store;
    node->ptr = ptr;
    node->value = value;
    node->name = name;
    return node;
}

Expr Call::make(Type t, std::string name, const std::vector<Expr> &args){
    Call *node = new Call;
    node->type = t;
    node->call_type = CallType::Intrinsic;//FIXME
    node->name=name;
    node->args = args;
    return node;
}

Stmt CallX::make(std::string name, const std::vector<Expr> &args){
    CallX *node = new CallX;
    node->name=name;
    node->args = args;
    node->call_type=CallType::Sympiler;//FIXME
    return node;
}

Stmt CallX::make(std::string name, const std::vector<Expr> &args,
                 const std::vector<Argument> &argums){
    CallX *node = new CallX;
    node->name=name;
    node->args = args;
    node->argums=argums;
    node->call_type=CallType::Sympiler;//FIXME
    return node;
}

Stmt Allocate::make(std::string name, Type type, const std::vector<Expr> &extents,
                    Expr condition, Stmt body,
                    Expr new_expr, std::string free_function) {
    for (size_t i = 0; i < extents.size(); i++) {
/*                internal_assert(extents[i].defined()) << "Allocate of undefined extent\n";
        internal_assert(extents[i].type().is_scalar() == 1) << "Allocate of vector extent\n";
    }
    internal_assert(body.defined()) << "Allocate of undefined\n";
    internal_assert(condition.defined()) << "Allocate with undefined condition\n";
    internal_assert(condition.type().is_bool()) << "Allocate condition is not boolean\n";*/

        Allocate *node = new Allocate;
        node->name = name;
        node->type = type;
        node->extents = extents;
        node->new_expr = new_expr;
        node->free_function = free_function;
        node->condition = condition;
        node->body = body;
        return node;
    }
}

Stmt Free::make(std::string name) {
    Free *node = new Free;
    node->name = name;
    return node;
}

Stmt IfThenElse::make(Expr condition, Stmt then_case, Stmt else_case) {
//    internal_assert(condition.defined() && then_case.defined()) << "IfThenElse of undefined\n";
    // else_case may be null.

    IfThenElse *node = new IfThenElse;
    node->condition = condition;
    node->then_case = then_case;
    node->else_case = else_case;
    return node;
}

template<> void ExprNode<IntImm>::accept(IRVisitor *v) const { v->visit((const IntImm *)this); }
template<> void ExprNode<UIntImm>::accept(IRVisitor *v) const { v->visit((const UIntImm *)this); }
template<> void ExprNode<FloatImm>::accept(IRVisitor *v) const { v->visit((const FloatImm *)this); }
template<> void ExprNode<StringImm>::accept(IRVisitor *v) const { v->visit((const StringImm *)this); }
template<> void ExprNode<Variable>::accept(IRVisitor *v) const { v->visit((const Variable *) this); }
template<> void ExprNode<Pointer>::accept(IRVisitor *v) const { v->visit((const Pointer *) this); }
template<> void ExprNode<Add>::accept(IRVisitor *v) const { v->visit((const Add *) this); }
template<> void ExprNode<Sub>::accept(IRVisitor *v) const { v->visit((const Sub *) this); }
template<> void ExprNode<Mul>::accept(IRVisitor *v) const { v->visit((const Mul *) this); }
template<> void ExprNode<Div>::accept(IRVisitor *v) const { v->visit((const Div *)this); }
template<> void ExprNode<Mod>::accept(IRVisitor *v) const { v->visit((const Mod *)this); }
template<> void ExprNode<Min>::accept(IRVisitor *v) const { v->visit((const Min *)this); }
template<> void ExprNode<Max>::accept(IRVisitor *v) const { v->visit((const Max *)this); }
template<> void ExprNode<EQ>::accept(IRVisitor *v) const { v->visit((const EQ *)this); }
template<> void ExprNode<NE>::accept(IRVisitor *v) const { v->visit((const NE *)this); }
template<> void ExprNode<LT>::accept(IRVisitor *v) const { v->visit((const LT *)this); }
template<> void ExprNode<LE>::accept(IRVisitor *v) const { v->visit((const LE *)this); }
template<> void ExprNode<GT>::accept(IRVisitor *v) const { v->visit((const GT *)this); }
template<> void ExprNode<GE>::accept(IRVisitor *v) const { v->visit((const GE *)this); }
template<> void ExprNode<And>::accept(IRVisitor *v) const { v->visit((const And *)this); }
template<> void ExprNode<Or>::accept(IRVisitor *v) const { v->visit((const Or *)this); }
template<> void ExprNode<Not>::accept(IRVisitor *v) const { v->visit((const Not *)this); }
template<> void ExprNode<Select>::accept(IRVisitor *v) const { v->visit((const Select *)this); }
template<> void ExprNode<Let>::accept(IRVisitor *v) const { v->visit((const Let *)this); }

        template<> void StmtNode<For>::accept(IRVisitor *v) const { v->visit((const For *)this); }
template<> void StmtNode<LetStmt>::accept(IRVisitor *v) const { v->visit((const LetStmt *)this); }
template<> void StmtNode<LetStmtPtr>::accept(IRVisitor *v) const { v->visit((const LetStmtPtr *)this); }
template<> void StmtNode<Evaluate>::accept(IRVisitor *v) const { v->visit((const Evaluate *)this); }
template<> void StmtNode<Block>::accept(IRVisitor *v) const { v->visit((const Block *)this); }
template<> void ExprNode<Load>::accept(IRVisitor *v) const { v->visit((const Load *)this); }
template<> void StmtNode<Store>::accept(IRVisitor *v) const { v->visit((const Store *)this); }
template<> void ExprNode<Call>::accept(IRVisitor *v) const { v->visit((const Call *)this); }
template<> void StmtNode<IfThenElse>::accept(IRVisitor *v) const { v->visit((const IfThenElse *)this); }
template<> void StmtNode<Allocate>::accept(IRVisitor *v) const { v->visit((const Allocate *)this); }
template<> void StmtNode<Free>::accept(IRVisitor *v) const { v->visit((const Free *)this); }
template<> void StmtNode<CallX>::accept(IRVisitor *v) const { v->visit((const CallX *)this); }

        template<> IRNodeType ExprNode<IntImm>::_type_info = {};
        template<> IRNodeType ExprNode<UIntImm>::_type_info = {};
        template<> IRNodeType ExprNode<FloatImm>::_type_info = {};
        template<> IRNodeType ExprNode<StringImm>::_type_info = {};
        //template<> IRNodeType ExprNode<Cast>::_type_info = {};
        template<> IRNodeType ExprNode<Variable>::_type_info = {};
        template<> IRNodeType ExprNode<Pointer>::_type_info = {};
        template<> IRNodeType ExprNode<Add>::_type_info = {};
        template<> IRNodeType ExprNode<Sub>::_type_info = {};
        template<> IRNodeType ExprNode<Mul>::_type_info = {};
        template<> IRNodeType ExprNode<Div>::_type_info = {};
        template<> IRNodeType ExprNode<Mod>::_type_info = {};
        template<> IRNodeType ExprNode<Min>::_type_info = {};
        template<> IRNodeType ExprNode<Max>::_type_info = {};
        template<> IRNodeType ExprNode<EQ>::_type_info = {};
        template<> IRNodeType ExprNode<NE>::_type_info = {};
        template<> IRNodeType ExprNode<LT>::_type_info = {};
        template<> IRNodeType ExprNode<LE>::_type_info = {};
        template<> IRNodeType ExprNode<GT>::_type_info = {};
        template<> IRNodeType ExprNode<GE>::_type_info = {};
        template<> IRNodeType ExprNode<And>::_type_info = {};
        template<> IRNodeType ExprNode<Or>::_type_info = {};
        template<> IRNodeType ExprNode<Not>::_type_info = {};
        template<> IRNodeType ExprNode<Select>::_type_info = {};
        template<> IRNodeType ExprNode<Load>::_type_info = {};
        //template<> IRNodeType ExprNode<Ramp>::_type_info = {};
        //template<> IRNodeType ExprNode<Broadcast>::_type_info = {};
        template<> IRNodeType ExprNode<Call>::_type_info = {};
        template<> IRNodeType StmtNode<CallX>::_type_info = {};
        template<> IRNodeType ExprNode<Let>::_type_info = {};
        template<> IRNodeType StmtNode<LetStmt>::_type_info = {};
        template<> IRNodeType StmtNode<LetStmtPtr>::_type_info = {};
        //template<> IRNodeType StmtNode<AssertStmt>::_type_info = {};
        //template<> IRNodeType StmtNode<ProducerConsumer>::_type_info = {};
        template<> IRNodeType StmtNode<For>::_type_info = {};
        template<> IRNodeType StmtNode<Store>::_type_info = {};
        //template<> IRNodeType StmtNode<Provide>::_type_info = {};
        template<> IRNodeType StmtNode<Allocate>::_type_info = {};
        template<> IRNodeType StmtNode<Free>::_type_info = {};
        //template<> IRNodeType StmtNode<Realize>::_type_info = {};
        template<> IRNodeType StmtNode<Block>::_type_info = {};
        template<> IRNodeType StmtNode<IfThenElse>::_type_info = {};
        template<> IRNodeType StmtNode<Evaluate>::_type_info = {};
    }
}