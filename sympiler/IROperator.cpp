//
// Created by kazem on 07/02/16.
//

#include "IROperator.h"
namespace Sympiler {


// Evaluate a float polynomial efficiently, taking instruction latency
// into account. The high order terms come first. n is the number of
// terms, which is the degree plus one.
/*        namespace {
            Expr evaluate_polynomial(Expr x, float *coeff, int n) {
                internal_assert(n >= 2);

                Expr x2 = x * x;

                Expr even_terms = coeff[0];
                Expr odd_terms = coeff[1];

                for (int i = 2; i < n; i++) {
                    if ((i & 1) == 0) {
                        if (coeff[i] == 0.0f) {
                            even_terms *= x2;
                        } else {
                            even_terms = even_terms * x2 + coeff[i];
                        }
                    } else {
                        if (coeff[i] == 0.0f) {
                            odd_terms *= x2;
                        } else {
                            odd_terms = odd_terms * x2 + coeff[i];
                        }
                    }
                }

                if ((n & 1) == 0) {
                    return even_terms * x + odd_terms;
                } else {
                    return odd_terms * x + even_terms;
                }
            }
        }*/

namespace Internal {


    bool is_const(Expr e) {
        if (e.as<IntImm>() ||
            e.as<UIntImm>() ||
            e.as<FloatImm>() ||
            e.as<StringImm>()) {
            return true;
/*                } else if (const Cast *c = e.as<Cast>()) {
            return is_const(c->value);
        } else if (const Ramp *r = e.as<Ramp>()) {
            return is_const(r->base) && is_const(r->stride);
        } else if (const Broadcast *b = e.as<Broadcast>()) {
            return is_const(b->value);*/
        } else {
            return false;
        }
    }

    bool is_const(Expr e, int64_t value) {
        if (const IntImm *i = e.as<IntImm>()) {
            return i->value == value;
        } else if (const UIntImm *i = e.as<UIntImm>()) {
            return (value >= 0) && (i->value == (uint64_t)value);
        } else if (const FloatImm *i = e.as<FloatImm>()) {
            return i->value == value;
/*                } else if (const Cast *c = e.as<Cast>()) {
            return is_const(c->value, value);
        } else if (const Broadcast *b = e.as<Broadcast>()) {
            return is_const(b->value, value);*/
        } else {
            return false;
        }
    }

    bool is_no_op(Stmt s) {
        if (!s.defined()) return true;
        const Evaluate *e = s.as<Evaluate>();
        return e && is_const(e->value);
    }

    const int64_t *as_const_int(Expr e) {
        if (!e.defined()) {
            return NULL;
/*                } else if (const Broadcast *b = e.as<Broadcast>()) {
            return as_const_int(b->value);*/
        } else if (const IntImm *i = e.as<IntImm>()) {
            return &(i->value);
        } else {
            return NULL;
        }
    }

    const uint64_t *as_const_uint(Expr e) {
        if (!e.defined()) {
            return NULL;
/*                } else if (const Broadcast *b = e.as<Broadcast>()) {
            return as_const_uint(b->value);*/
        } else if (const UIntImm *i = e.as<UIntImm>()) {
            return &(i->value);
        } else {
            return NULL;
        }
    }

    const double *as_const_float(Expr e) {
        if (!e.defined()) {
            return NULL;
        /*} else if (const Broadcast *b = e.as<Broadcast>()) {
            return as_const_float(b->value);*/
        } else if (const FloatImm *f = e.as<FloatImm>()) {
            return &(f->value);
        } else {
            return NULL;
        }
    }

    bool is_const_power_of_two_integer(Expr e, int *bits) {
        if (!(e.type().is_int() || e.type().is_uint())) return false;

/*                const Broadcast *b = e.as<Broadcast>();
        if (b) return is_const_power_of_two_integer(b->value, bits);

        const Cast *c = e.as<Cast>();
        if (c) return is_const_power_of_two_integer(c->value, bits);*/

        uint64_t val = 0;

        if (const int64_t *i = as_const_int(e)) {
            if (*i < 0) return false;
            val = (uint64_t)(*i);
        } else if (const uint64_t *u = as_const_uint(e)) {
            val = *u;
        }

        if (val && ((val & (val - 1)) == 0)) {
            *bits = 0;
            for (; val; val >>= 1) {
                if (val == 1) return true;
                (*bits)++;
            }
        }

        return false;
    }

    bool is_positive_const(Expr e) {
        if (const IntImm *i = e.as<IntImm>()) return i->value > 0;
        if (const UIntImm *u = e.as<UIntImm>()) return u->value > 0;
        if (const FloatImm *f = e.as<FloatImm>()) return f->value > 0.0f;
/*                if (const Cast *c = e.as<Cast>()) {
            return is_positive_const(c->value);
        }
        if (const Ramp *r = e.as<Ramp>()) {
            // slightly conservative
            return is_positive_const(r->base) && is_positive_const(r->stride);
        }
        if (const Broadcast *b = e.as<Broadcast>()) {
            return is_positive_const(b->value);
        }*/
        return false;
    }

    bool is_negative_const(Expr e) {
        if (const IntImm *i = e.as<IntImm>()) return i->value < 0;
        if (const FloatImm *f = e.as<FloatImm>()) return f->value < 0.0f;
/*                if (const Cast *c = e.as<Cast>()) {
            return is_negative_const(c->value);
        }
        if (const Ramp *r = e.as<Ramp>()) {
            // slightly conservative
            return is_negative_const(r->base) && is_negative_const(r->stride);
        }
        if (const Broadcast *b = e.as<Broadcast>()) {
            return is_negative_const(b->value);
        }*/
        return false;
    }

    bool is_negative_negatable_const(Expr e, Type T) {
        if (const IntImm *i = e.as<IntImm>()) {
            return (i->value < 0 && !T.is_min(i->value));
        }
        if (const FloatImm *f = e.as<FloatImm>()) return f->value < 0.0f;
/*                if (const Cast *c = e.as<Cast>()) {
            return is_negative_negatable_const(c->value, c->type);
        }
        if (const Ramp *r = e.as<Ramp>()) {
            // slightly conservative
            return is_negative_negatable_const(r->base) && is_negative_const(r->stride);
        }
        if (const Broadcast *b = e.as<Broadcast>()) {
            return is_negative_negatable_const(b->value);
        }*/
        return false;
    }

    bool is_negative_negatable_const(Expr e) {
        return is_negative_negatable_const(e, e.type());
    }

    bool is_zero(Expr e) {
        if (const IntImm *int_imm = e.as<IntImm>()) return int_imm->value == 0;
        if (const UIntImm *uint_imm = e.as<UIntImm>()) return uint_imm->value == 0;
        if (const FloatImm *float_imm = e.as<FloatImm>()) return float_imm->value == 0.0;
     //   if (const Cast *c = e.as<Cast>()) return is_zero(c->value);
    //    if (const Broadcast *b = e.as<Broadcast>()) return is_zero(b->value);
        return false;
    }

    bool is_one(Expr e) {
        if (const IntImm *int_imm = e.as<IntImm>()) return int_imm->value == 1;
        if (const UIntImm *uint_imm = e.as<UIntImm>()) return uint_imm->value == 1;
        if (const FloatImm *float_imm = e.as<FloatImm>()) return float_imm->value == 1.0;
    //    if (const Cast *c = e.as<Cast>()) return is_one(c->value);
    //    if (const Broadcast *b = e.as<Broadcast>()) return is_one(b->value);
        return false;
    }

    bool is_two(Expr e) {
        if (const IntImm *int_imm = e.as<IntImm>()) return int_imm->value == 2;
        if (const UIntImm *uint_imm = e.as<UIntImm>()) return uint_imm->value == 2;
        if (const FloatImm *float_imm = e.as<FloatImm>()) return float_imm->value == 2.0;
    //    if (const Cast *c = e.as<Cast>()) return is_two(c->value);
    //    if (const Broadcast *b = e.as<Broadcast>()) return is_two(b->value);
        return false;
    }

    namespace {
        template<typename T>
        Expr make_const_helper(Type t, T val) {
            if (t.is_vector()) {
              //  return Broadcast::make(make_const(t.element_of(), val), t.lanes());
            } else if (t.is_int()) {
                return IntImm::make(t, (int64_t)val);
            } else if (t.is_uint()) {
                return UIntImm::make(t, (uint64_t)val);
            } else if (t.is_float()) {
                return FloatImm::make(t, (double)val);
            } else {
              //  internal_error << "Can't make a constant of type " << t << "\n";
                return Expr();
            }
        }
    }

    Expr make_const(Type t, int64_t val) {
        return make_const_helper(t, val);
    }

    Expr make_const(Type t, uint64_t val) {
        return make_const_helper(t, val);
    }

    Expr make_const(Type t, double val) {
        return make_const_helper(t, val);
    }

    void check_representable(Type dst, int64_t x) {
/*        if (dst.is_handle()) {
            user_assert(dst.can_represent(x))
            << "Integer constant " << x
            << " will be implicitly coerced to type " << dst
            << ", but Halide does not support pointer arithmetic.\n";
        } else {
            user_assert(dst.can_represent(x))
            << "Integer constant " << x
            << " will be implicitly coerced to type " << dst
            << ", which changes its value to " << make_const(dst, x)
            << ".\n";
        }*/
    }

    Expr make_bool(bool val, int w) {
        return make_const(UInt(1, w), val);
    }

    Expr make_zero(Type t) {
        /*if (t.is_handle()) {
            return Call::make(Handle(), Call::null_handle, std::vector<Expr>(), Call::Intrinsic);
        } else */{
            return make_const(t, 0);
        }
    }

    Expr make_one(Type t) {
        return make_const(t, 1);
    }

    Expr make_two(Type t) {
        return make_const(t, 2);
    }

    Expr const_true(int w) {
        return make_one(UInt(1, w));
    }

    Expr const_false(int w) {
        return make_zero(UInt(1, w));
    }

    Expr operator+(Expr a, Expr b) {
        return Add::make(a, b);
    }
//TODO
/*bool is_const_power_of_two_integer(const Expr &e, int *bits) {
if (!(e.type().is_int() || e.type().is_uint())) return false;

const Broadcast *b = e.as<Broadcast>();
if (b) return is_const_power_of_two_integer(b->value, bits);

const Cast *c = e.as<Cast>();
if (c) return is_const_power_of_two_integer(c->value, bits);

uint64_t val = 0;

if (const int64_t *i = as_const_int(e)) {
    if (*i < 0) return false;
    val = (uint64_t)(*i);
} else if (const uint64_t *u = as_const_uint(e)) {
    val = *u;
}

if (val && ((val & (val - 1)) == 0)) {
    *bits = 0;
    for (; val; val >>= 1) {
        if (val == 1) return true;
        (*bits)++;
    }
}

return false;
}
*/
}
}