/**
 * This file generates C code for the given program
 */
#ifndef DSLPROJECT_CODEGEN_C_H
#define DSLPROJECT_CODEGEN_C_H


#include <iosfwd>
#include <fstream>
#include <map>
#include "Module.h"
#include "IRPrinter.h"
#include "IR.h"
#include "Util.h"
#include "Scope.h"

namespace Sympiler{
namespace Internal{

class CodeGen_C : public IRPrinter {
    //std::ostream &outFile;

protected:
    /** An ID for the most recently generated ssa variable */
    std::string id, ptr_id;
    bool loopParsing;
    bool isHeader;
    /** A cache of generated values in scope */
    std::map<std::string, std::string> cache;
    std::vector<std::string> tmpPtrId;
    int openBracketNo;
    enum AppendSpaceIfNeeded {
        DoNotAppendSpace,
        AppendSpace,
    };

    /** Open a new C scope (i.e. throw in a brace, increase the indent) */
    void open_scope();

    /** Close a C scope (i.e. throw in an end brace, decrease the indent) */
    void close_scope(const std::string &comment);

    struct Allocation {
        Type type;
        std::string free_function;
    };

    /** Track the types of allocations to avoid unnecessary casts. */
    Scope<Allocation> allocations;

    /** Track which allocations actually went on the heap. */
    Scope<int> heap_allocations;

    /** True if there is a void * __user_context parameter in the arguments. */
    bool have_user_context;

public:
    //CodeGen_C();
    CodeGen_C(std::ostream& file, bool is_header=false, bool is_IG=true, const std::string &include_guard = "");
    ~CodeGen_C();
    void compile(Module module);
    void compile(Kernel& ker);

    void visit(const Variable *v );
    void visit(const Pointer *p);
    void visit(const IntImm *);
    void visit(const UIntImm *);
    void visit(const StringImm *);
    void visit(const FloatImm *);
    void visit(const Add *op);
    void visit(const Sub *);
    void visit(const Mul *);
    void visit(const Div *);
    void visit(const Mod *);
    void visit(const Min *);
    void visit(const Max *);
    void visit(const EQ *);
    void visit(const NE *);
    void visit(const LT *);
    void visit(const LE *);
    void visit(const GT *);
    void visit(const GE *);
    void visit(const And *);
    void visit(const Or *);
    void visit(const Not *);
    void visit(const Select *);
    void visit(const For *op);
    void visit(const Evaluate *op);
    void visit(const Let *);
    void visit(const LetStmt *);
    void visit(const LetStmtPtr *);
    void visit(const Load *);
    void visit(const Store *);
    void visit_binop(Type t, Expr a, Expr b, const char * op);
    void visit(const Call *op);
    void visit(const CallX *op);
    void visit(const Allocate *);
    void visit(const Free *);
    void visit(const IfThenElse *);
    std::string print_expr(Expr e);
    std::string print_ptr(Pointer e);
    void print_stmt(Stmt e);
    std::string print_type(Type );
    std::string print_assignment(Type t, const std::string &rhs);
    //std::string print_type(Type type, AppendSpaceIfNeeded space_option = DoNotAppendSpace);
    std::string print_name(const std::string &name);

};
}
}



#endif //DSLPROJECT_CODEGEN_C_H
