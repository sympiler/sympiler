//
// Created by kazem on 13/02/16.
//

#ifndef DSLPROJECT_IRPRINTER_H
#define DSLPROJECT_IRPRINTER_H

#include <ostream>
#include "IRVisitor.h"
namespace Sympiler {
namespace Internal {

/** Emit an expression on an output stream (such as std::cout) in a
 * human-readable form */
        std::ostream &operator<<(std::ostream &stream, const Expr &);

/** Emit a halide type on an output stream (such as std::cout) in a
 * human-readable form */
        std::ostream &operator<<(std::ostream &stream, const Type &);

/** Emit a halide Module on an output stream (such as std::cout) in a
 * human-readable form */
//        EXPORT std::ostream &operator<<(std::ostream &stream, const Module &);

/** Emit a halide device api type in a human readable form */
        std::ostream &operator<<(std::ostream &stream, const DeviceAPI &);

        namespace Internal {

/** Emit a halide statement on an output stream (such as std::cout) in
 * a human-readable form */
            std::ostream &operator<<(std::ostream &stream, const Stmt &);

/** Emit a halide for loop type (vectorized, serial, etc) in a human
 * readable form */
            std::ostream &operator<<(std::ostream &stream, const ForType &);

        }

class IRPrinter : public IRVisitor {
public:
    /** Construct an IRPrinter pointed at a given output stream
     * (e.g. std::cout, or a std::ofstream) */

    IRPrinter(std::ostream &s);


    /** emit an expression on the output stream */
    void print(Expr);
    /** emit a statement on the output stream */
    void print(Stmt);

protected:
    /** The stream we're outputting on */
    std::ostream &stream;

    /** The current indentation level, useful for pretty-printing
     * statements */
    int indent;

    /** Emit spaces according to the current indentation level */
    void do_indent();
    void visit(const Variable *v);
    void visit(const Pointer *p);
    void visit(const IntImm *);
    void visit(const UIntImm *);
    void visit(const FloatImm *);
    void visit(const StringImm *);
    void visit(const Add *);
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
    void visit(const Evaluate *);
    void visit(const Block *op);
    void visit(const Let *);
    void visit(const LetStmt *);
    void visit(const Load *);
    void visit(const Store *);
    void visit(const Call *op);
    void visit(const CallX *op);
    void visit(const Allocate *);
    void visit(const Free *);
    void visit(const IfThenElse *);
};
}
}
#endif //DSLPROJECT_IRPRINTER_H
