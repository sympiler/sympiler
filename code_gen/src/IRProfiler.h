//
// Created by george on 2020-02-22.
//

#ifndef FUSION_IRPROFILER_H
#define FUSION_IRPROFILER_H

#include <map>
#include <ostream>
#include <fstream>
#include <unordered_map>
#include <sstream>
#include <stack>
#include <unordered_set>
#include "IRVisitor.h"
#include "Module.h"
#include "IRPrinter.h"

#define PRINT_AST

namespace Sympiler {
namespace Internal {


typedef struct {
    bool array;
    int def;
    int dead;
    std::vector<int> uses;
} var;

typedef struct {
    const For *l;
    int start;
    int end;
    int scope; // outmost is 0
    std::vector<const For *> innerloops;
    std::unordered_set<std::string> reads;
    std::unordered_set<std::string> writes;
} loop;


class IRProfiler : public IRVisitor {

protected:
    // total line numbers
    std::string id, id_ptr; // most recently used Expr
    std::string new_v; // most recently created _x variable
    std::vector<std::string> scope_var;

    int cur_lines;
    std::map<std::string, var*> variables;
    std::unordered_map<std::string, std::string> expr_map;

    // for loops
    bool parseloop;
    int loopscope;
    std::stack<loop *> loopstack;

    std::stringstream ss;
    IRPrinter *printer;

    std::ofstream stream;
    bool print_ast;
    int op_count;

public:
    std::vector<const loop*> loops;

    IRProfiler();

    void start_analysis(Module module);
    void start_analysis(Kernel& ker);

    void print_trace();

    void print_AST(const std::string &fname, Kernel &ker);

    void clear_map();

    std::string visit_expr(Expr e);
    void visit_stmt(Stmt s);

    /// generate new variable and place in table
    void use_var(int line, const std::string& name);
    std::string store_name(const std::string& s);
    var *new_var(const std::string& s);

    void visit_binop(Type t, Expr a, Expr b, const char *op);

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

#endif //FUSION_IRPROFILER_H
