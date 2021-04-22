#include <sstream>
#include <limits>
#include "CodeGen_C.h"
#include "IROperator.h"
#include "Substitute.h"

namespace Sympiler {
namespace Internal {

namespace {
const std::string buffer_t_definition =
    "#ifndef HALIDE_ATTRIBUTE_ALIGN\n"
            "  #ifdef _MSC_VER\n"
            "    #define HALIDE_ATTRIBUTE_ALIGN(x) __declspec(align(x))\n"
            "  #else\n"
            "    #define HALIDE_ATTRIBUTE_ALIGN(x) __attribute__((aligned(x)))\n"
            "  #endif\n"
            "#endif\n"
            "#ifndef BUFFER_T_DEFINED\n"
            "#define BUFFER_T_DEFINED\n"
            "#include <stdbool.h>\n"
            "#include <stdint.h>\n"
            "typedef struct buffer_t {\n"
            "    uint64_t dev;\n"
            "    uint8_t* host;\n"
            "    int32_t extent[4];\n"
            "    int32_t stride[4];\n"
            "    int32_t min[4];\n"
            "    int32_t elem_size;\n"
            "    HALIDE_ATTRIBUTE_ALIGN(1) bool host_dirty;\n"
            "    HALIDE_ATTRIBUTE_ALIGN(1) bool dev_dirty;\n"
            "    HALIDE_ATTRIBUTE_ALIGN(1) uint8_t _padding[10 - sizeof(void *)];\n"
            "} buffer_t;\n"
            "#endif\n"
            "#define __user_context_ NULL\n";

const std::string headers =
    "#include <iostream>\n"
            "#include <math.h>\n"
            "#include <float.h>\n"
            "#include <assert.h>\n"
            "#include <string.h>\n"
            "#include <stdio.h>\n"
            "#include <stdint.h>\n"
            "#include <cholUtils.h>\n";

const std::string globals =//TODO add more feastures in future
    "extern \"C\" {\n"
            "void *sympiler_malloc(void *ctx, size_t s){return(malloc(s));}\n"
            "void sympiler_free(void *ctx, void *ptr){free(ptr);};\n"
       /*     "void *sympiler_print(void *ctx, const void *str);\n"
            "void *sympiler_error(void *ctx, const void *str);\n"
            "int sympiler_debug_to_file(void *ctx, const char *filename, void *data, int, int, int, int, int, int);\n"
            "int sympiler_start_clock(void *ctx);\n"
            "int64_t sympiler_current_time_ns(void *ctx);\n"
            "void sympiler_profiler_pipeline_end(void *, void *);\n"*/
            "}\n"
            "\n"

            // TODO: this next chunk is copy-pasted from posix_math.cpp. A
            // better solution for the C runtime would be nice.
            "#ifdef _WIN32\n"
            "float roundf(float);\n"
            "double round(double);\n"
            "#else\n"
            "inline float asinh_f32(float x) {return asinhf(x);}\n"
            "inline float acosh_f32(float x) {return acoshf(x);}\n"
            "inline float atanh_f32(float x) {return atanhf(x);}\n"
            "inline double asinh_f64(double x) {return asinh(x);}\n"
            "inline double acosh_f64(double x) {return acosh(x);}\n"
            "inline double atanh_f64(double x) {return atanh(x);}\n"
            "#endif\n"
            "inline float sqrt_f32(float x) {return sqrtf(x);}\n"
            "inline float sin_f32(float x) {return sinf(x);}\n"
            "inline float asin_f32(float x) {return asinf(x);}\n"
            "inline float cos_f32(float x) {return cosf(x);}\n"
            "inline float acos_f32(float x) {return acosf(x);}\n"
            "inline float tan_f32(float x) {return tanf(x);}\n"
            "inline float atan_f32(float x) {return atanf(x);}\n"
            "inline float sinh_f32(float x) {return sinhf(x);}\n"
            "inline float cosh_f32(float x) {return coshf(x);}\n"
            "inline float tanh_f32(float x) {return tanhf(x);}\n"
            "inline float hypot_f32(float x, float y) {return hypotf(x, y);}\n"
            "inline float exp_f32(float x) {return expf(x);}\n"
            "inline float log_f32(float x) {return logf(x);}\n"
            "inline float pow_f32(float x, float y) {return powf(x, y);}\n"
            "inline float floor_f32(float x) {return floorf(x);}\n"
            "inline float ceil_f32(float x) {return ceilf(x);}\n"
            "inline float round_f32(float x) {return roundf(x);}\n"
            "\n"
            "inline double sqrt_f64(double x) {return sqrt(x);}\n"
            "inline double sin_f64(double x) {return sin(x);}\n"
            "inline double asin_f64(double x) {return asin(x);}\n"
            "inline double cos_f64(double x) {return cos(x);}\n"
            "inline double acos_f64(double x) {return acos(x);}\n"
            "inline double tan_f64(double x) {return tan(x);}\n"
            "inline double atan_f64(double x) {return atan(x);}\n"
            "inline double sinh_f64(double x) {return sinh(x);}\n"
            "inline double cosh_f64(double x) {return cosh(x);}\n"
            "inline double tanh_f64(double x) {return tanh(x);}\n"
            "inline double hypot_f64(double x, double y) {return hypot(x, y);}\n"
            "inline double exp_f64(double x) {return exp(x);}\n"
            "inline double log_f64(double x) {return log(x);}\n"
            "inline double pow_f64(double x, double y) {return pow(x, y);}\n"
            "inline double floor_f64(double x) {return floor(x);}\n"
            "inline double ceil_f64(double x) {return ceil(x);}\n"
            "inline double round_f64(double x) {return round(x);}\n"
            "\n"
            "inline float nan_f32() {return NAN;}\n"
            "inline float neg_inf_f32() {return -INFINITY;}\n"
            "inline float inf_f32() {return INFINITY;}\n"
            "inline bool is_nan_f32(float x) {return x != x;}\n"
            "inline bool is_nan_f64(double x) {return x != x;}\n"
            "inline float float_from_bits(uint32_t bits) {\n"
            " union {\n"
            "  uint32_t as_uint;\n"
            "  float as_float;\n"
            " } u;\n"
            " u.as_uint = bits;\n"
            " return u.as_float;\n"
            "}\n"
            "inline int64_t make_int64(int32_t hi, int32_t lo) {\n"
            "    return (((int64_t)hi) << 32) | (uint32_t)lo;\n"
            "}\n"
            "inline double make_float64(int32_t i0, int32_t i1) {\n"
            "    union {\n"
            "        int32_t as_int32[2];\n"
            "        double as_double;\n"
            "    } u;\n"
            "    u.as_int32[0] = i0;\n"
            "    u.as_int32[1] = i1;\n"
            "    return u.as_double;\n"
            "}\n"
            "\n"
            "template<typename T> T max(T a, T b) {if (a > b) return a; return b;}\n"
            "template<typename T> T min(T a, T b) {if (a < b) return a; return b;}\n"

            // This may look wasteful, but it's the right way to do
            // it. Compilers understand memcpy and will convert it to a no-op
            // when used in this way. See http://blog.regehr.org/archives/959
            // for a detailed comparison of type-punning methods.
            "template<typename A, typename B> A reinterpret(B b) {A a; memcpy(&a, &b, sizeof(a)); return a;}\n"
            "\n"
            //We will need zero and one arrays for BLAS  calls
            "double one [2]={1.0,0.}, zero [2]={0.,0.};\n"
            "int sw = false, lb1 = 0, ub1 = 0; \n"
            "double *cur; int info=0; \n"    ;
}

CodeGen_C::CodeGen_C(std::ostream &s, bool is_header, bool IGType,
                     const std::string &guard) : IRPrinter(s),

                              id("$$ BAD ID $$"),
                              isHeader(is_header),
                              loopParsing(false),
                              ptr_id(""),
                              openBracketNo(0) {
    if (is_header) {
        // If it's a header, emit an include guard.
        stream << "#ifndef HALIDE_" << print_name(guard) << '\n'
        << "#define HALIDE_" << print_name(guard) << '\n';
    }

    if (!is_header) {
        stream << headers;
    }

    // Throw in a definition of a buffer_t
    stream << buffer_t_definition;
    if(IGType){
        stream << "#define BLOCKED\n";
    }

    // halide_filter_metadata_t just gets a forward declaration
    // (include HalideRuntime.h for the full goodness)
    stream << "struct halide_filter_metadata_t;\n";

    if (!is_header) {
        stream << globals;
    }



    // Everything from here on out is extern "C".
    stream << "#ifdef __cplusplus\n";
    stream << "extern \"C\" {\n";
    stream << "#endif\n";
}

CodeGen_C::~CodeGen_C() {
   // stream <<"}\n";
    stream << "#ifdef __cplusplus\n";
    stream << "}  // extern \"C\"\n";
    stream << "#endif\n";

    if (isHeader) {
        stream << "#endif\n";
    }
}

namespace {
std::string type_to_c_type(Type type) {
    std::ostringstream oss;
    //user_assert(type.lanes() == 1) << "Can't use vector types when compiling to C (yet)\n";
    if (type.is_float()) {
        if (type.bits() == 32) {
            oss << "float";
        } else if (type.bits() == 64) {
            oss << "double";
        } else {
      //      user_error << "Can't represent a float with this many bits in C: " << type << "\n";
        }

    } else if (type.is_handle()) {
        oss << "void *";
    } else {
        switch (type.bits()) {
            case 1:
                oss << "bool";
                break;
            case 8: case 16: case 32: case 64:
                if (type.is_uint()) oss << 'u';
                oss << "int" << type.bits() << "_t";
                break;
            default:
                oss << "ERR";
                //user_error << "Can't represent an integer with this many bits in C: " << type << "\n";
        }
    }
    return oss.str();
}
}

std::string CodeGen_C::print_type(Type type) {
    return type_to_c_type(type);
}

void CodeGen_C::compile(Module module) {
    //inserting headers ...
    compile(module.kers[0]);

}

void CodeGen_C::compile(Kernel& ker) {

    //Emit the function prototype
    //return value
    stream << print_type(ker.retVal.type) <<" " <<ker.Name() <<"(";
    for (int i = 0; i < ker.argType.size(); ++i) {
        stream << print_type(ker.argType[i].type)<<" ";
        if(ker.argType[i].is_output())
            stream << "&";
        if(ker.argType[i].is_pointer())
            stream << "*";
        stream << ker.argType[i].name;
        if(i < ker.argType.size()-1)
            stream <<", ";
    }
    if(isHeader)
        stream << ")\n";
    else{
        stream<< ") {\n";
        indent+=1;
        //Emit body
        print_stmt(ker.LoweredKer());
        do_indent();
        stream << "return 0;\n";
        indent -= 1;
        stream << "}\n";
    }
}

void CodeGen_C::visit(const Variable *v) {
    id = v->name;
}

void CodeGen_C::visit(const Pointer *p) {
    ptr_id += p->name; //define group member ptr_id
    if(p->idx.defined() ){
        if(p->indirect){
            ptr_id += "[";
            openBracketNo++;//To make sure that the brackets are paired
            ptr_id += print_expr(p->idx);
            ptr_id += "]";
            openBracketNo--;
        }
        else{
            ptr_id += print_expr(p->idx);
        }
    }
    id=ptr_id;
    ptr_id="";
}

void CodeGen_C::visit(const IntImm *op) {
    if (op->type == Int(32)) {
        id = std::to_string(op->value);
    } else {
        print_assignment(op->type, "(" + print_type(op->type) + ")(" + std::to_string(op->value) + ")");
    }
}

void CodeGen_C::visit(const UIntImm *op) {
    print_assignment(op->type, "(" + print_type(op->type) + ")(" + std::to_string(op->value) + ")");
}

void CodeGen_C::visit(const StringImm *op) {
    std::ostringstream oss;
    oss << Expr(op);
    id = oss.str();
}

// NaN is the only float/double for which this is true... and
// surprisingly, there doesn't seem to be a portable isnan function
// (dsharlet).
template <typename T>
static bool isnan(T x) { return x != x; }

template <typename T>
static bool isinf(T x)
{
    return std::numeric_limits<T>::has_infinity && (
            x == std::numeric_limits<T>::infinity() ||
            x == -std::numeric_limits<T>::infinity());
}

void CodeGen_C::visit(const FloatImm *op) {
    if (isnan(op->value)) {
        id = "nan_f32()";
    } else if (isinf(op->value)) {
        if (op->value > 0) {
            id = "inf_f32()";
        } else {
            id = "neg_inf_f32()";
        }
    } else {
        // Write the constant as reinterpreted uint to avoid any bits lost in conversion.
        union {
            uint32_t as_uint;
            float as_float;
        } u;
        u.as_float = op->value;

        std::ostringstream oss;
        oss << "float_from_bits(" << u.as_uint << " /* " << u.as_float << " */)";
        id = oss.str();
    }
}

void CodeGen_C::visit(const Add *op) {
    visit_binop(op->type, op->a, op->b, "+");
}

void CodeGen_C::visit(const Sub *op) {
    visit_binop(op->type, op->a, op->b, "-");
}

void CodeGen_C::visit(const Mul *op) {
    visit_binop(op->type, op->a, op->b, "*");
}

void CodeGen_C::visit(const Div *op) {
    int bits;
    if (is_const_power_of_two_integer(op->b, &bits)) {
        std::ostringstream oss;
        oss << print_expr(op->a) << " >> " << bits;
        print_assignment(op->type, oss.str());
    } else if (op->type.is_int()) {
        std::string a = print_expr(op->a);
        std::string b = print_expr(op->b);
        // q = a / b
        std::string q = print_assignment(op->type, a + " / " + b);
        // r = a - q * b
        std::string r = print_assignment(op->type, a + " - " + q + " * " + b);
        // bs = b >> (8*sizeof(T) - 1)
        std::string bs = print_assignment(op->type, b + " >> (" + print_type(op->type.element_of()) + ")" + std::to_string(op->type.bits() - 1));
        // rs = r >> (8*sizeof(T) - 1)
        std::string rs = print_assignment(op->type, r + " >> (" + print_type(op->type.element_of()) + ")" + std::to_string(op->type.bits() - 1));
        // id = q - (rs & bs) + (rs & bs)
        print_assignment(op->type, q + " - (" + rs + " & " + bs + ") + (" + rs + " & ~" + bs + ")");
    } else  {
        visit_binop(op->type, op->a, op->b, "/");
    }
}

void CodeGen_C::visit(const Mod *op) {
    int bits;
    if (is_const_power_of_two_integer(op->b, &bits)) {
        std::ostringstream oss;
        oss << print_expr(op->a) << " & " << ((1 << bits)-1);
        print_assignment(op->type, oss.str());
    /*} else if (op->type.is_int()) {//FIXME
        std::string a = print_expr(op->a);
        std::string b = print_expr(op->b);
        // r = a % b
        std::string r = print_assignment(op->type, a + " % " + b);
        // rs = r >> (8*sizeof(T) - 1)
        std::string rs = print_assignment(op->type, r + " >> (" + print_type(op->type.element_of()) + ")" + std::to_string(op->type.bits() - 1));
        // abs_b = abs(b)
        std::string abs_b = print_expr(cast(op->type, abs(op->b)));
        // id = r + (abs_b & rs)
        print_assignment(op->type, r + " + (" + abs_b + " & " + rs + ")");*/
    } else {
        visit_binop(op->type, op->a, op->b, "%");
    }
}

void CodeGen_C::visit(const Max *op) {
    print_expr(Call::make(op->type, "max", {op->a, op->b}));
}

void CodeGen_C::visit(const Min *op) {
    print_expr(Call::make(op->type, "min", {op->a, op->b}));
}

void CodeGen_C::visit(const EQ *op) {
    visit_binop(op->type, op->a, op->b, "==");
}

void CodeGen_C::visit(const NE *op) {
    visit_binop(op->type, op->a, op->b, "!=");
}

void CodeGen_C::visit(const LT *op) {
    visit_binop(op->type, op->a, op->b, "<");
}

void CodeGen_C::visit(const LE *op) {
    visit_binop(op->type, op->a, op->b, "<=");
}

void CodeGen_C::visit(const GT *op) {
    visit_binop(op->type, op->a, op->b, ">");
}

void CodeGen_C::visit(const GE *op) {
    visit_binop(op->type, op->a, op->b, ">=");
}

void CodeGen_C::visit(const Or *op) {
    visit_binop(op->type, op->a, op->b, "||");
}

void CodeGen_C::visit(const And *op) {
    visit_binop(op->type, op->a, op->b, "&&");
}

void CodeGen_C::visit(const Not *op) {
    print_assignment(op->type, "!(" + print_expr(op->a) + ")");
}

void CodeGen_C::visit(const Select *op) {
    std::ostringstream rhs;
    std::string tmp;
    tmp = ptr_id;
    ptr_id="";//There might be a value from printing a pointer
    std::string true_val = print_expr(op->true_value);
    std::string false_val = print_expr(op->false_value);
    std::string cond = print_expr(op->condition);
    ptr_id = tmp;
    rhs << "(" << print_type(op->type) << ")"
    << "(" << cond
    << " ? " << true_val
    << " : " << false_val
    << ")";
    print_assignment(op->type, rhs.str());
}

void CodeGen_C::visit(const For *op) {
    if (op->for_type == ForType::Parallel) {
        do_indent();
        stream << "#pragma omp parallel for\n";
    }
    loopParsing=false;
    std::string id_min = print_expr(op->min);
    std::string id_extent = print_expr(op->extent);
    loopParsing= false;
    do_indent();
    stream << "for (int "
    << print_name(op->name)
    << " = " << id_min
    << "; "
    << print_name(op->name)
    << " < " << id_extent
    //<< " < " << id_min
    //<< " + " << id_extent
    << "; "
    << print_name(op->name)
    << "++)\n";

    open_scope();
    op->body.accept(this);
    close_scope("for " + print_name(op->name));
}

void CodeGen_C::visit(const Evaluate *op){
    print(op->value);
}

void CodeGen_C::visit(const Let *op) {
    std::string id_value = print_expr(op->value);
    Expr body = op->body;
    if (op->value.type().is_handle()) {
        // The body might contain a Load that references this directly
        // by name, so we can't rewrite the name.
        do_indent();
        stream << print_type(op->value.type())
        << " " << print_name(op->name)
        << " = " << id_value << ";\n";
    } else {
        Expr new_var = Variable::make(op->value.type(), id_value);
        body = substitute(op->name, new_var, body);
    }
    print_expr(body);
}

void CodeGen_C::visit(const LetStmt *op) {
    std::string id_value = print_expr(op->value);
    Stmt body = op->body;
    if (op->value.type().is_handle()) {
        // The body might contain a Load or Store that references this
        // directly by name, so we can't rewrite the name.
        do_indent();
        stream << print_type(op->value.type())
        << " " << print_name(op->name)
        << " = " << id_value << ";\n";
    } else {
        Expr new_var = Variable::make(op->value.type(), id_value);
        body = substitute(op->name, new_var, body);
    }
    body.accept(this);
}

void CodeGen_C::visit(const LetStmtPtr *op) { //FIXME probably don't need this

    std::string id_value = print_expr(op->value);
    Expr lhs = op->lhs;
    Stmt body = op->body;
    if (op->value.type().is_handle()) {
        // The body might contain a Load or Store that references this
        // directly by name, so we can't rewrite the name.
        do_indent();
        //stream << print_type(op->value.type())
        stream << print_expr(lhs)
        << " = " << id_value << ";\n";
    } else {
//        Expr new_var = Variable::make(op->value.type(), id_value);
//        body = substitute(op->name, new_var, body);
    }
    body.accept(this);
}

void CodeGen_C::visit(const Load *op){
    Type t = op->type;
    bool type_cast_needed =
            !allocations.contains(op->name) ||
            allocations.get(op->name).type != t;

    std::ostringstream rhs;
    if (type_cast_needed) {
        rhs << "(("
        << print_type(op->type)
        << " *)"
        << print_expr(op->ptr)
        << ")";
    } else {
        rhs << print_expr(op->ptr);
    }

    print_assignment(op->type, rhs.str());
}

bool constant_allocation_size(const std::vector<Expr> &extents, const std::string &name, int32_t &size) {
    int64_t result = 1;

    for (size_t i = 0; i < extents.size(); i++) {
        if (const IntImm *int_size = extents[i].as<IntImm>()) {
            // Check if the individual dimension is > 2^31 - 1. Not
            // currently necessary because it's an int32_t, which is
            // always smaller than 2^31 - 1. If we ever upgrade the
            // type of IntImm but not the maximum allocation size, we
            // should re-enable this.
            /*
            if ((int64_t)int_size->value > (((int64_t)(1)<<31) - 1)) {
                user_error
                    << "Dimension " << i << " for allocation " << name << " has size " <<
                    int_size->value << " which is greater than 2^31 - 1.";
            }
            */
            result *= int_size->value;
            if (result > (static_cast<int64_t>(1)<<31) - 1) {
                /*user_error
                << "Total size for allocation " << name
                << " is constant but exceeds 2^31 - 1.\n";*/
            }
        } else {
            return false;
        }
    }

    size = static_cast<int32_t>(result);
    return true;
}


void CodeGen_C::visit(const Store *op) {

    Type t = op->value.type();

    /*bool type_cast_needed =
            t.is_handle() ||
            !allocations.contains(op->name) ||
            allocations.get(op->name).type != t;*/

    std::string id_index = print_expr(op->ptr);
    std::string id_value = print_expr(op->value);
    do_indent();

    /*if (type_cast_needed) {
        stream << "((const "
        << print_type(t)
        << " *)"
        << print_name(op->name)
        << ")";
    } else {
        stream << print_name(op->name);
    }*/
    //stream << "["
    stream << id_index
    << " = "
    << id_value
    << ";\n";

    cache.clear();
}

void CodeGen_C::visit(const Allocate *op) {
    open_scope();

    // For sizes less than 8k, do a stack allocation
    bool on_stack = false;
    int32_t constant_size;
    std::string size_id;
    if (op->new_expr.defined()) {
        Allocation alloc;
        alloc.type = op->type;
        alloc.free_function = op->free_function;
        allocations.push(op->name, alloc);
        heap_allocations.push(op->name, 0);
        stream << print_type(op->type) << "*" << print_name(op->name) << " = (" << print_expr(op->new_expr) << ");\n";
    } else {
        if (constant_allocation_size(op->extents, op->name, constant_size)) {
            int64_t stack_bytes = constant_size * op->type.bytes();

            if (stack_bytes > ((int64_t(1) << 31) - 1)) {
             //   user_error << "Total size for allocation "
             //   << op->name << " is constant but exceeds 2^31 - 1.\n";
            } else {
                size_id = print_expr(Expr(static_cast<int32_t>(constant_size)));
                if (stack_bytes <= 1024 * 8) {
                    on_stack = true;
                }
            }
        } else {
            // Check that the allocation is not scalar (if it were scalar
            // it would have constant size).
            //internal_assert(op->extents.size() > 0);

            size_id = print_assignment(Int(64), print_expr(op->extents[0]));

            for (size_t i = 1; i < op->extents.size(); i++) {
                // Make the code a little less cluttered for two-dimensional case
                std::string new_size_id_rhs;
                std::string next_extent = print_expr(op->extents[i]);
                if (i > 1) {
                    new_size_id_rhs =  "(" + size_id + " > ((int64_t(1) << 31) - 1)) ? " + size_id + " : (" + size_id + " * " + next_extent + ")";
                } else {
                    new_size_id_rhs = size_id + " * " + next_extent;
                }
                size_id = print_assignment(Int(64), new_size_id_rhs);
            }
            do_indent();
            stream << "if ((" << size_id << " > ((int64_t(1) << 31) - 1)) || ((" << size_id <<
            " * sizeof(" << print_type(op->type) << ")) > ((int64_t(1) << 31) - 1)))\n";
            open_scope();
            do_indent();
            /*stream << "sympiler_error("
            << (have_user_context ? "__user_context_" : "NULL")
            << ", \"32-bit signed overflow computing size of allocation "
            << op->name << "\\n\");\n";*/
            do_indent();
            stream << "return -1;\n";
            close_scope("overflow test " + op->name);
        }

        // Check the condition to see if this allocation should actually be created.
        // If the allocation is on the stack, the only condition we can respect is
        // unconditional false (otherwise a non-constant-sized array declaration
        // will be generated).
        if (!on_stack || is_zero(op->condition)) {
            Expr conditional_size = Select::make(op->condition,
                                                 //Var(size_id),
                                                 Variable::make(Int(64),size_id),
                                                 Expr(static_cast<int32_t>(0)));
            //conditional_size = simplify(conditional_size);
            size_id = print_assignment(Int(64), print_expr(conditional_size));
        }

        Allocation alloc;
        alloc.type = op->type;
        allocations.push(op->name, alloc);

        do_indent();
        stream << print_type(op->type) << ' ';

        if (on_stack) {
            stream << print_name(op->name)
            << "[" << size_id << "];\n";
        } else {
            stream << "*"
            << print_name(op->name)
            << " = ("
            << print_type(op->type)
            << " *)sympiler_malloc("
            << (have_user_context ? "__user_context_" : "NULL")
            << ", sizeof("
            << print_type(op->type)
            << ")*" << size_id << ");\n";
            heap_allocations.push(op->name, 0);
        }
    }

    op->body.accept(this);

    // Should have been freed internally
    //internal_assert(!allocations.contains(op->name));

    close_scope("alloc " + print_name(op->name));
}

void CodeGen_C::visit(const Free *op) {
    if (heap_allocations.contains(op->name)) {
        std::string free_function = allocations.get(op->name).free_function;
        if (free_function.empty()) {
            free_function = "sympiler_free";
        }

        do_indent();
        stream << free_function << "("
        << (have_user_context ? "__user_context_, " : "NULL, ")
        << print_name(op->name)
        << ");\n";
        heap_allocations.pop(op->name);
    }
    allocations.pop(op->name);
}

void CodeGen_C::visit(const IfThenElse *op) {
    std::string cond_id = print_expr(op->condition);

    do_indent();
    stream << "if (" << cond_id << ")\n";
    open_scope();
    op->then_case.accept(this);
    close_scope("if " + cond_id);

    if (op->else_case.defined()) {
        do_indent();
        stream << "else\n";
        open_scope();
        op->else_case.accept(this);
        close_scope("if " + cond_id + " else");
    }
}

std::string CodeGen_C::print_expr(Expr e) {
    id = "$$ BAD ID $$";
    e.accept(this);
    return id;
}

void CodeGen_C::print_stmt(Stmt s){
    s.accept(this);
}

std::string CodeGen_C::print_ptr(Pointer e) {
    id = "$$ BAD ID $$";
    e.accept(this);
    return id;
}



void CodeGen_C::visit_binop(Type t, Expr a, Expr b, const char *op) {
    std::string tmp = ptr_id;
    ptr_id="";
    std::string sa = print_expr(a);
    std::string sb = print_expr(b);
    ptr_id=tmp;
    if(loopParsing){
        ptr_id+=(sa+op+sb);
        id = "";
    }
    else{
        print_assignment(t, sa + " " + op + " " + sb);
    }
}

std::string CodeGen_C::print_name(const std::string &name) {
    std::ostringstream oss;

    // Prefix an underscore to avoid reserved words (e.g. a variable named "while")
    /*if (isalpha(name[0])) {
        oss << '_';
    }*/

    for (size_t i = 0; i < name.size(); i++) {
        if (name[i] == '.') {
            oss << '_';
        } else if (name[i] == '$') {
            oss << "__";
        } else if (name[i] != '_' && !isalnum(name[i])) {
            oss << "___";
        }
        else oss << name[i];
    }
    return oss.str();
}




void CodeGen_C::visit(const Call *op){
    std::ostringstream rhs;
    if(op->call_type == Call::Intrinsic){//FIXME: needs to make it general
        std::vector<std::string> args;
        for (size_t i = 0; i < op->args.size(); i++) {
            args.push_back(print_expr(op->args[i]));
        }
        rhs<<op->name<<"(";
        for (size_t i = 0; i < args.size()-1; i++) {
            rhs << args[i]  << ", ";
        }
        rhs << args[args.size()-1]<<")";
    }
    print_assignment(op->type, rhs.str());
}

void CodeGen_C::visit(const CallX *op){
    std::ostringstream rhs;
    std::string tmp;
    if(op->call_type == CallX::CallType::Sympiler){//FIXME: needs to make it general
        std::vector<std::string> args;
        for (size_t i = 0; i < op->args.size(); i++) {
            tmp = print_expr(op->args[i]);
            args.push_back(tmp);
        }
        do_indent();
        bool sw = op->argums.size()==op->args.size();
        stream<<op->name<<"(";
        for (size_t i = 0; i < args.size()-1; i++) {
            if(sw){
                if(op->argums[i].ref)
                    stream <<"&"<<args[i]<< ", ";
                else
                    stream <<args[i]<< ", ";
            }else{
                stream <<args[i]<< ", ";
            }
        }
        if(sw)
            stream <<"&"<<args[args.size()-1]<< ");\n";
        else
            stream << args[args.size()-1]<<");\n";
    }
}

std::string CodeGen_C::print_assignment(Type t, const std::string &rhs) {

    std::map<std::string, std::string>::iterator cached = cache.find(rhs);

    if (cached == cache.end()) {
        id = unique_name('_');
        do_indent();
        stream << print_type(t) << " " << id << " = " << rhs << ";\n";
        cache[rhs] = id;
    } else {
        id = cached->second;
    }
    return id;
}

void CodeGen_C::open_scope() {
    cache.clear();
    do_indent();
    indent++;
    stream << "{\n";
}

void CodeGen_C::close_scope(const std::string &comment) {
    cache.clear();
    indent--;
    do_indent();
    if (!comment.empty()) {
        stream << "} // " << comment << "\n";
    } else {
        stream << "}\n";
    }
}


}
}