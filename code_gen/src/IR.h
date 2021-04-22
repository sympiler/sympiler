//
// Created by kazem on 07/02/16.
//

#ifndef DSLPROJECT_IR_H
#define DSLPROJECT_IR_H

#include <string>
#include "Expr.h"
#include "Util.h"
#include "Argument.h"

namespace Sympiler {
namespace Internal {

struct Add : public ExprNode<Add> {
    Expr a, b;

    static Expr make(Expr a, Expr b);

};

/** The ratio of two expressions */
struct Div : public ExprNode<Div> {
    Expr a, b;

    static Expr make(Expr a, Expr b);
};

/** The difference of two expressions */
struct Sub : public ExprNode<Sub> {
    Expr a, b;

    static Expr make(Expr a, Expr b);
};

/** The product of two expressions */
struct Mul : public ExprNode<Mul> {
    Expr a, b;

    static Expr make(Expr a, Expr b);
};

/** The remainder of a / b. Mostly equivalent to '%' in C, except that
* the result here is always positive. For floats, this is equivalent
* to calling fmod. */
struct Mod : public ExprNode<Mod> {
    Expr a, b;

    static Expr make(Expr a, Expr b);
};

/** The lesser of two values. */
struct Min : public ExprNode<Min> {
    Expr a, b;

    static Expr make(Expr a, Expr b);
};

/** The greater of two values */
struct Max : public ExprNode<Max> {
    Expr a, b;

    static Expr make(Expr a, Expr b);
};

/** Is the first expression equal to the second */
struct EQ : public ExprNode<EQ> {
    Expr a, b;

    static Expr make(Expr a, Expr b);
};

/** Is the first expression not equal to the second */
struct NE : public ExprNode<NE> {
    Expr a, b;

    static Expr make(Expr a, Expr b);
};

/** Is the first expression less than the second. */
struct LT : public ExprNode<LT> {
    Expr a, b;

    static Expr make(Expr a, Expr b);
};

/** Is the first expression less than or equal to the second. */
struct LE : public ExprNode<LE> {
    Expr a, b;

    static Expr make(Expr a, Expr b);
};

/** Is the first expression greater than the second. */
struct GT : public ExprNode<GT> {
    Expr a, b;

    static Expr make(Expr a, Expr b);
};

/** Is the first expression greater than or equal to the second. */
struct GE : public ExprNode<GE> {
    Expr a, b;

    static Expr make(Expr a, Expr b);
};

/** Logical and - are both expressions true */
struct And : public ExprNode<And> {
    Expr a, b;

    static Expr make(Expr a, Expr b);
};

/** Logical or - is at least one of the expression true */
struct Or : public ExprNode<Or> {
    Expr a, b;

    static Expr make(Expr a, Expr b);
};

/** Logical not - true if the expression false */
struct Not : public ExprNode<Not> {
    Expr a;

    static Expr make(Expr a);
};

/** A ternary operator. Evalutes 'true_value' and 'false_value',
* then selects between them based on 'condition'. Equivalent to
* the ternary operator in C. */
struct Select : public ExprNode<Select> {
    Expr condition, true_value, false_value;

    static Expr make(Expr condition, Expr true_value, Expr false_value);
};

struct Variable : public ExprNode<Variable> {
    std::string name;

    //static Expr make(std::string name);
    static Expr make(Type type, std::string name); //{
        //return make(type, name);
    //}

};

struct Pointer : public ExprNode<Pointer> {
    std::string name;
    bool indirect;
    Expr idx;

    //bool defined();
    static Expr make(Type type, std::string name, Expr idx);

};

struct For : public StmtNode<For> {
    std::string name;
    Expr min, extent;
    ForType for_type;
    DeviceAPI device_api;
    Stmt body;

    static Stmt make(std::string name, Expr min, Expr extent, ForType for_type, DeviceAPI device_api,
                     Stmt body);
};

/** A let expression, like you might find in a functional
* language. Within the expression \ref Let::body, instances of the Var
* node \ref Let::name refer to \ref Let::value. */
struct Let : public ExprNode<Let> {
    std::string name;
    Expr value, body;

    static Expr make(std::string name, Expr value, Expr body);
};

struct LetStmt : public StmtNode<LetStmt> {
    std::string name;
    Expr value;
    Stmt body;

    static Stmt make(std::string name, Expr value, Stmt body);
};

    /** LetStmtPtr allows using an array in the left hand side. We cannot
     * make the arrays and indirect memory accesses in SSA form
     */
struct LetStmtPtr : public StmtNode<LetStmtPtr> {
    Expr lhs;
    Expr value;
    Stmt body;

    static Stmt make(Expr lhs, Expr value, Stmt body);
};
/** Evaluate and discard an expression, presumably because it has some side-effect. */
struct Evaluate : public StmtNode<Evaluate> {
    Expr value;

    static Stmt make(Expr v);
};


/** A sequence of statements to be executed in-order. 'rest' may be
* undefined. Used rest.defined() to find out. */
struct Block : public StmtNode<Block> {
    Stmt first, rest;

    static Stmt make(const Stmt &first, const Stmt &rest);
    /** Construct zero or more Blocks to invoke a list of statements in order.
     * This method may not return a Block statement if stmts.size() <= 1. */
    static Stmt make(const std::vector<Stmt> &stmts);

   // static const IRNodeType _type_info = IRNodeType::Block;
};

/** Load a value from a named buffer. The buffer is treated as an
* array of the 'type' of this Load node. That is, the buffer has
* no inherent type. */
struct Load : public ExprNode<Load> {
    std::string name;
    Expr ptr;

    // If it's a load from an image argument or compiled-in constant
    // image, this will point to that
    //Buffer image; //TODO change Buffer to Matrix?

    // If it's a load from an image parameter, this points to that
    //Parameter param; //TODO using this for Matrix?

    //EXPORT static Expr make(Type type, std::string name, Expr index, Buffer image, Parameter param);
    static Expr make(Type type, std::string name, Expr ptr);
};

/** Store a 'value' to the buffer called 'name' at a given
* 'index'. The buffer is interpreted as an array of the same type as
* 'value'. */
struct Store : public StmtNode<Store> {
    std::string name;
    Expr ptr;
    Expr value;

    static Stmt make(Expr ptr, std::string name, Expr value);
};

/** A function call. This can represent a call to some extern
* function (like sin), but it's also our multi-dimensional
* version of a Load, so it can be a load from an input image, or
* a call to another Sympiler kernel. The latter two types of call
* nodes don't survive all the way down to code generation - the
* lowering process converts them to Load nodes. */
struct Call : public ExprNode<Call> {
    std::string name;
    std::vector<Expr> args;
    typedef enum {Matrix, Extern, Sympiler, Intrinsic} CallType;
    CallType call_type;

    // Sympiler uses calls internally to represent certain operations
    // (instead of IR nodes). These are matched by name. Note that
    // these are deliberately char* (rather than std::string) so that
    // they can be referenced at static-initialization time without
    // risking ambiguous initalization order; we use a typedef to simplify
    // declaration.
    typedef const char* const ConstString;
    EXPORT static ConstString debug_to_file,
            shuffle_vector,
            interleave_vectors,
            reinterpret,
            bitwise_and,
            bitwise_not,
            bitwise_xor,
            bitwise_or,
            shift_left,
            shift_right,
            abs,
            absd,
            rewrite_buffer,
            random,
            lerp,
            create_buffer_t,
            copy_buffer_t,
            extract_buffer_min,
            extract_buffer_max,
            extract_buffer_host,
            set_host_dirty,
            set_dev_dirty,
            popcount,
            count_leading_zeros,
            count_trailing_zeros,
            undef,
            null_handle,
            address_of,
            return_second,
            if_then_else,
            trace,
            trace_expr,
            glsl_texture_load,
            glsl_texture_store,
            glsl_varying,
            image_load,
            image_store,
            make_struct,
            stringify,
            memoize_expr,
            copy_memory,
            likely,
            make_int64,
            make_float64,
            register_destructor;

    // If it's a call to another Sympiler function, this call node
    // holds onto a pointer to that function.
    //Kernel func;

    // If that function has multiple values, which value does this
    // call node refer to?
    int value_index;

    // If it's a call to an image, this call nodes hold a
    // pointer to that image's buffer
    //Buffer image; //TODO maybe a matrix

    // If it's a call to an image parameter, this call node holds a
    // pointer to that
    //Parameter param;
/*
    EXPORT static Expr make(Type type, std::string name, const std::vector<Expr> &args, CallType call_type,
                            Function func = Function(), int value_index = 0,
                            Buffer image = Buffer(), Parameter param = Parameter());

    /** Convenience constructor for calls to other Sympiler kernels */
/*        static Expr make(Function func, const std::vector<Expr> &args, int idx = 0) {
     /*   internal_assert(idx >= 0 &&
                        idx < func.outputs())
                << "Value index out of range in call to Sympiler kernel\n";
        internal_assert(func.has_pure_definition() || func.has_extern_definition())
                << "Call to undefined Sympiler kernel\n";*/
/*            return make(func.output_types()[(size_t)idx], func.name(), args, Sympiler, func, idx, Buffer(), Parameter());
    }

    /** Convenience constructor for loads from concrete images */
/*        static Expr make(Buffer image, const std::vector<Expr> &args) {
        return make(image.type(), image.name(), args, Image, Function(), 0, image, Parameter());
    }

    /** Convenience constructor for loads from images parameters */
/*        static Expr make(Parameter param, const std::vector<Expr> &args) {
        return make(param.type(), param.name(), args, Image, Function(), 0, Buffer(), param);
    }*/
    static Expr make(Type t, std::string name, const std::vector<Expr> &args);
};

struct CallX : public StmtNode<CallX> {
    std::string name;
    std::vector<Expr> args;
    std::vector<Argument> argums;
    typedef enum {BLAS, Sympiler} CallType;
    CallType call_type;

    //
    typedef const char* const ConstString;
    static ConstString debug_to_file,
            shuffle_vector,
            register_destructor;

    int value_index;
    static Stmt make(std::string name, const std::vector<Expr> &args);

    static Stmt make(std::string name, const std::vector<Expr> &args,
                     const std::vector<Argument> &argums );
};

/** An if-then-else block. 'else' may be undefined. */
struct IfThenElse : public StmtNode<IfThenElse> {
    Expr condition;
    Stmt then_case, else_case;

    static Stmt make(Expr condition, Stmt then_case, Stmt else_case = Stmt());
};

/** Allocate a scratch area called with the given name, type, and
* size. The buffer lives for at most the duration of the body
* statement, within which it is freed. It is an error for an allocate
* node not to contain a free node of the same buffer. Allocation only
* occurs if the condition evaluates to true. */
struct Allocate : public StmtNode<Allocate> {
    std::string name;
    Type type;
    std::vector<Expr> extents;
    Expr condition;

    // These override the code generator dependent malloc and free
    // equivalents if provided. If the new_expr succeeds, that is it
    // returns non-NULL, the function named be free_function is
    // guaranteed to be called. The free function signature must match
    // that of the code generator dependent free (typically
    // sympiler_free). If free_function is left empty, code generator
    // default will be called.
    Expr new_expr;
    std::string free_function;
    Stmt body;

    static Stmt make(std::string name, Type type, const std::vector<Expr> &extents,
                            Expr condition, Stmt body,
                            Expr new_expr = Expr(), std::string free_function = std::string());
};

/** Free the resources associated with the given buffer. */
struct Free : public StmtNode<Free> {
    std::string name;

    static Stmt make(std::string name);
};

}
}
#endif //DSLPROJECT_IR_H
