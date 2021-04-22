//
// Created by kazem on 07/02/16.
//

#ifndef DSLPROJECT_EXPR_H
#define DSLPROJECT_EXPR_H

#include "IntrusivePtr.h"
#include "Type.h"
#include "Util.h"
//#include "Error.h"

namespace Sympiler {
namespace Internal {

class IRVisitor;
/** A class representing a type of IR node (e.g. Add, or Mul, or
 * For). We use it for rtti (without having to compile with rtti). */
    struct IRNodeType {};
  /*  enum class IRNodeType {
        IntImm,
        UIntImm,
        FloatImm,
        StringImm,
        Cast,
        Variable,
        Add,
        Sub,
        Mul,
        Div,
        Mod,
        Min,
        Max,
        EQ,
        NE,
        LT,
        LE,
        GT,
        GE,
        And,
        Or,
        Not,
        Select,
        Load,
        Ramp,
        Broadcast,
        Call,
        Let,
        LetStmt,
        AssertStmt,
        ProducerConsumer,
        For,
        Store,
        Provide,
        Allocate,
        Free,
        Realize,
        Block,
        IfThenElse,
        Evaluate,
        Shuffle,
        Prefetch,
    };*/

    enum sympiler_argument_kind_t {
        sympiler_argument_kind_input_scalar = 0,
        sympiler_argument_kind_input_buffer = 1,
        sympiler_argument_kind_output_buffer = 2
    };

/** The abstract base classes for a node in the Halide IR. */
struct IRNode {
    /** We use the visitor pattern to traverse IR nodes throughout the
     * compiler, so we have a virtual accept method which accepts
     * visitors.
     */
    IRNode() { };

    virtual void accept(IRVisitor *v) const = 0;

    virtual ~IRNode() { }
    /** These classes are all managed with intrusive reference
       counting, so we also track a reference count. It's mutable
       so that we can do reference counting even through const
       references to IR nodes. */
    mutable RefCount ref_count;

    /** Each IR node subclass should return some unique pointer. We
 * can compare these pointers to do runtime type
 * identification. We don't compile with rtti because that
 * injects run-time type identification stuff everywhere (and
 * often breaks when linking external libraries compiled
 * without it), and we only want it for IR nodes. */
    virtual const IRNodeType *type_info() const = 0;

};

template<>
inline RefCount &ref_count<IRNode>(const IRNode *n) {
    return n->ref_count; }

template<>
inline void destroy<IRNode>(const IRNode *n) {
    delete n; }

/** A base class for statement nodes. They have no properties or
methods beyond base IR nodes for now */
struct BaseStmtNode : public IRNode {
};

struct BaseExprNode : public IRNode {
      Type type;
};

template<typename T>
struct ExprNode : public BaseExprNode {

    void accept(IRVisitor *v) const;
    virtual IRNodeType *type_info() const override{return &_type_info;}
    static IRNodeType _type_info;
};

template<typename T>
struct StmtNode : public BaseStmtNode {
    void accept(IRVisitor *v) const;
    virtual IRNodeType *type_info() const override{return &_type_info;}
    static IRNodeType _type_info;
};


struct IRHandle : public IntrusivePtr<const IRNode> {
    IRHandle() : IntrusivePtr<const IRNode>() { }

    IRHandle(const IRNode *p) : IntrusivePtr<const IRNode>(p) { }

    /** Dispatch to the correct visitor method for this node. E.g. if
     * this node is actually an Add node, then this will call
     * IRVisitor::visit(const Add *) */
    void accept(IRVisitor *v) const {
        ptr->accept(v);
    }

    /** Downcast this ir node to its actual type (e.g. Add, or
     * Select). This returns NULL if the node is not of the requested
     * type. Example usage:
     *
     * if (const Add *add = node->as<Add>()) {
     *   // This is an add node
     * }
     */
    template<typename T> const T *as() const {
        if (ptr->type_info() == &T::_type_info) {
            return (const T *)ptr;
        }
        return NULL;
    }
};

/** Integer constants */
struct IntImm : public ExprNode<IntImm> {
    int64_t value;

    static IntImm *make(Type t, int64_t value) {
      //  internal_assert(t.is_int() && t.is_scalar()) << "IntImm must be a scalar Int\n";
      //  internal_assert(t.bits() == 8 || t.bits() == 16 || t.bits() == 32 || t.bits() == 64)
      //  << "IntImm must be 8, 16, 32, or 64-bit\n";

     /*   if (t.bits() == 32 && value >= -8 && value <= 8 &&
            !small_int_cache[(int)value + 8].ref_count.is_zero()) {
            return &small_int_cache[(int)value + 8];
        }*/

        IntImm *node = new IntImm;
        node->type = t;
        // Normalize the value by dropping the high bits
        value <<= (64 - t.bits());
        // Then sign-extending to get them back
        value >>= (64 - t.bits());
        node->value = value;
        return node;
    }

private:
    static IntImm small_int_cache[17];
    //static const IRNodeType _type_info = IRNodeType::IntImm;//for RTTI
};
    struct UIntImm : public ExprNode<UIntImm> {
        uint64_t value;

        static UIntImm *make(Type t, uint64_t value) {
/*            internal_assert(t.is_uint() && t.is_scalar())
            << "UIntImm must be a scalar UInt\n";
            internal_assert(t.bits() == 1 || t.bits() == 8 || t.bits() == 16 || t.bits() == 32 || t.bits() == 64)
            << "UIntImm must be 1, 8, 16, 32, or 64-bit\n";
*/
            UIntImm *node = new UIntImm;
            node->type = t;
            // Normalize the value by dropping the high bits
            value <<= (64 - t.bits());
            value >>= (64 - t.bits());
            node->value = value;
            return node;
        }
    };

/** Floating point constants */
    struct FloatImm : public ExprNode<FloatImm> {
        double value;

        static FloatImm *make(Type t, double value) {
//            internal_assert(t.is_float() && t.is_scalar()) << "FloatImm must be a scalar Float\n";
            FloatImm *node = new FloatImm;
            node->type = t;
            switch (t.bits()) {
  /*              case 16:
                    node->value = (double)((float16_t)value);
                    break;*/ //TODO: not sure if we need this
                case 32:
                    node->value = (float)value;
                    break;
                case 64:
                    node->value = value;
                    break;
//                default: //FIXME
//                    internal_error << "FloatImm must be 16, 32, or 64-bit\n";
            }

            return node;
        }
    };

/** String constants */
struct StringImm : public ExprNode<StringImm> {
    std::string value;

    static StringImm *make(const std::string &val) {
        StringImm *node = new StringImm;
        node->type = Handle();
        node->value = val;
        return node;
    }
};

struct Expr : public IRHandle {
    Expr() : IRHandle() { };

    Expr(const BaseExprNode *n) : IRHandle(n) { }
    //Expr(int v){};
    /** Make an expression representing numeric constants of various types. */
    // @{
    EXPORT explicit Expr(int8_t x)    : IRHandle(Internal::IntImm::make(Int(8), x)) {}
    EXPORT explicit Expr(int16_t x)   : IRHandle(Internal::IntImm::make(Int(16), x)) {}
    EXPORT          Expr(int32_t x)   : IRHandle(Internal::IntImm::make(Int(32), x)) {}
    EXPORT explicit Expr(int64_t x)   : IRHandle(Internal::IntImm::make(Int(64), x)) {}
    EXPORT explicit Expr(uint8_t x)   : IRHandle(Internal::UIntImm::make(UInt(8), x)) {}
    EXPORT explicit Expr(uint16_t x)  : IRHandle(Internal::UIntImm::make(UInt(16), x)) {}
    EXPORT explicit Expr(uint32_t x)  : IRHandle(Internal::UIntImm::make(UInt(32), x)) {}
    EXPORT explicit Expr(uint64_t x)  : IRHandle(Internal::UIntImm::make(UInt(64), x)) {}
  //  EXPORT          Expr(float16_t x) : IRHandle(Internal::FloatImm::make(Float(16), (double)x)) {}
    EXPORT          Expr(float x)     : IRHandle(Internal::FloatImm::make(Float(32), x)) {}
    EXPORT explicit Expr(double x)    : IRHandle(Internal::FloatImm::make(Float(64), x)) {}
    // @}

    /** Make an expression representing a const string (i.e. a StringImm) */
//    EXPORT          Expr(const std::string &s) : IRHandle(Internal::StringImm::make(s)) {}

    /** Get the type of this expression node */
    Type type() const {
        return ((const Internal::BaseExprNode *)ptr)->type;
    }
};

/** An enum describing a type of device API. Used by schedules, and in
* the For loop IR node. */
enum class DeviceAPI {
    Parent, /// Used to denote for loops that inherit their device from where they are used, generally the default
    Host,
    Default_GPU,
    CUDA,
    OpenCL,
    GLSL,
    Renderscript,
    OpenGLCompute,
    Metal
};

/** An enum describing a type of loop traversal. Used in schedules,
* and in the For loop IR node. */
enum class ForType {
    Serial,
    Parallel,
    Vectorized,
    Unrolled,
    Pruned,
    Blocked,
    Peeled
};

/** A reference-counted handle to a statement node. */
struct Stmt : public IRHandle {
    Stmt() : IRHandle() { }

    Stmt(const BaseStmtNode *n) : IRHandle(n) { }

    /** This lets you use a Stmt as a key in a map of the form
     * map<Stmt, Foo, Stmt::Compare> */
    struct Compare {
        bool operator()(const Stmt &a, const Stmt &b) const {
            return a.ptr < b.ptr;
        }
    };
};
}
}
#endif //DSLPROJECT_EXPR_H
