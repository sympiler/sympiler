//
// Created by kazem on 07/02/16.
//

#ifndef DSLPROJECT_IROPERATOR_H
#define DSLPROJECT_IROPERATOR_H

#include "Expr.h"
#include "IR.h"
#include "Var.h"
namespace Sympiler {
namespace Internal {

    /** Is the expression either an IntImm, a FloatImm, a StringImm, or a
 * Cast of the same, or a Ramp or Broadcast of the same. Doesn't do
 * any constant folding. */
    bool is_const(Expr e);

/** Is the expression an IntImm, FloatImm of a particular value, or a
 * Cast, or Broadcast of the same. */
        bool is_const(Expr e, int64_t v);

/** If an expression is an IntImm or a Broadcast of an IntImm, return
 * a pointer to its value. Otherwise returns NULL. */
        const int64_t *as_const_int(Expr e);

/** If an expression is a UIntImm or a Broadcast of a UIntImm, return
 * a pointer to its value. Otherwise returns NULL. */
        const uint64_t *as_const_uint(Expr e);

/** If an expression is a FloatImm or a Broadcast of a FloatImm,
 * return a pointer to its value. Otherwise returns NULL. */
        const double *as_const_float(Expr e);

/** Is the expression a constant integer power of two. Also returns
 * log base two of the expression if it is. Only returns true for
 * integer types. */
        bool is_const_power_of_two_integer(Expr e, int *bits);

/** Is the expression a const (as defined by is_const), and also
 * strictly greater than zero (in all lanes, if a vector expression) */
        bool is_positive_const(Expr e);

/** Is the expression a const (as defined by is_const), and also
 * strictly less than zero (in all lanes, if a vector expression) */
        bool is_negative_const(Expr e);

/** Is the expression a const (as defined by is_const), and also
 * strictly less than zero (in all lanes, if a vector expression) and
 * is its negative value representable. (This excludes the most
 * negative value of the Expr's type from inclusion. Intended to be
 * used when the value will be negated as part of simplification.)
 */
        bool is_negative_negatable_const(Expr e);

/** Is the expression a const (as defined by is_const), and also equal
 * to zero (in all lanes, if a vector expression) */
        bool is_zero(Expr e);

/** Is the expression a const (as defined by is_const), and also equal
 * to one (in all lanes, if a vector expression) */
        bool is_one(Expr e);

/** Is the expression a const (as defined by is_const), and also equal
 * to two (in all lanes, if a vector expression) */
        bool is_two(Expr e);

/** Is the statement a no-op (which we represent as either an
 * undefined Stmt, or as an Evaluate node of a constant) */
        bool is_no_op(Stmt s);

/** Construct an immediate of the given type from any numeric C++ type. */
// @{
        Expr make_const(Type t, int64_t val);
        Expr make_const(Type t, uint64_t val);
        Expr make_const(Type t, double val);
        inline Expr make_const(Type t, int32_t val)   {return make_const(t, (int64_t)val);}
        inline Expr make_const(Type t, uint32_t val)  {return make_const(t, (uint64_t)val);}
        inline Expr make_const(Type t, int16_t val)   {return make_const(t, (int64_t)val);}
        inline Expr make_const(Type t, uint16_t val)  {return make_const(t, (uint64_t)val);}
        inline Expr make_const(Type t, int8_t val)    {return make_const(t, (int64_t)val);}
        inline Expr make_const(Type t, uint8_t val)   {return make_const(t, (uint64_t)val);}
        inline Expr make_const(Type t, bool val)      {return make_const(t, (uint64_t)val);}
        inline Expr make_const(Type t, float val)     {return make_const(t, (double)val);}
      //  inline Expr make_const(Type t, float16_t val) {return make_const(t, (double)val);}
/** Check if a constant value can be correctly represented as the given type. */
        void check_representable(Type t, int64_t val);

/** Construct a boolean constant from a C++ boolean value.
 * May also be a vector if width is given.
 * It is not possible to coerce a C++ boolean to Expr because
 * if we provide such a path then char objects can ambiguously
 * be converted to Halide Expr or to std::string.  The problem
 * is that C++ does not have a real bool type - it is in fact
 * close enough to char that C++ does not know how to distinguish them.
 * make_bool is the explicit coercion. */
        Expr make_bool(bool val, int lanes = 1);

/** Construct the representation of zero in the given type */
        Expr make_zero(Type t);

/** Construct the representation of one in the given type */
        Expr make_one(Type t);

/** Construct the representation of two in the given type */
        Expr make_two(Type t);

/** Construct the constant boolean true. May also be a vector of
 * trues, if a lanes argument is given. */
        Expr const_true(int lanes = 1);

/** Construct the constant boolean false. May also be a vector of
 * falses, if a lanes argument is given. */
        Expr const_false(int lanes = 1);

Expr operator+(Expr a, Expr b);
//bool is_const_power_of_two_integer(const Expr &e, int *bits);//TODO
}
}
#endif //DSLPROJECT_IROPERATOR_H
