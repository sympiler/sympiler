//
// Created by kazem on 15/02/16.
//

#ifndef DSLPROJECT_TYPE_H
#define DSLPROJECT_TYPE_H
#include <string>
/** A set of types to represent a C++ function signature. This allows
 * two things.  First, proper prototypes can be provided for Halide
 * generated functions, giving better compile time type
 * checking. Second, C++ name mangling can be done to provide link
 * time type checking for both Halide generated functions and calls
 * from Halide to external functions.
 *
 * These are intended to be constexpr producable, but we don't depend
 * on C++11 yet. In C++14, it is possible these will be replaced with
 * introspection/reflection facilities.
 *
 * halide_handle_traits has to go outside the Halide namespace due to template
 * resolution rules. TODO(zalman): Do all types need to be in global namespace?
 */
//@{

/** A structure to represent the (unscoped) name of a C++ composite type for use
 * as a single argument (or return value) in a function signature.
 *
 * Currently does not support the restrict qualifier, references, or
 * r-value references.  These features cannot be used in extern
 * function calls from Halide or in the generated function from
 * Halide, but their applicability seems limited anyway.
 */
struct halide_cplusplus_type_name {
    /// An enum to indicate whether a C++ type is non-composite, a struct, class, or union
    enum CPPTypeType {
        Simple, ///< "int"
        Struct, ///< "struct Foo"
        Class,  ///< "class Foo"
        Union,  ///< "union Foo"
        Enum,   ///< "enum Foo"
    } cpp_type_type;  // Note: order is reflected in map_to_name table in CPlusPlusMangle.cpp

    std::string name;

    halide_cplusplus_type_name(CPPTypeType cpp_type_type, const std::string &name)
            : cpp_type_type(cpp_type_type), name(name) {
    }

    bool operator==(const halide_cplusplus_type_name &rhs) const {
        return cpp_type_type == rhs.cpp_type_type &&
               name == rhs.name;
    }

    bool operator!=(const halide_cplusplus_type_name &rhs) const {
        return !(*this == rhs);
    }

    bool operator<(const halide_cplusplus_type_name &rhs) const {
        return cpp_type_type < rhs.cpp_type_type ||
               (cpp_type_type == rhs.cpp_type_type &&
                name < rhs.name);
    }
};

/** Types in the halide type system. They can be ints, unsigned ints,
* or floats (of various bit-widths), or a handle (which is always 64-bits).
* Note that the int/uint/float values do not imply a specific bit width
* (the bit width is expected to be encoded in a separate value).
*/
typedef enum halide_type_code_t
#if __cplusplus >= 201103L
        : uint8_t
#endif
{
    halide_type_int = 0,   //!< signed integers
    halide_type_uint = 1,  //!< unsigned integers
    halide_type_float = 2, //!< floating point numbers
    halide_type_handle = 3 //!< opaque pointer type (void *)
} halide_type_code_t;

// Note that while __attribute__ can go before or after the declaration,
// __declspec apparently is only allowed before.
#ifndef HALIDE_ATTRIBUTE_ALIGN
#ifdef _MSC_VER
#define HALIDE_ATTRIBUTE_ALIGN(x) __declspec(align(x))
#else
#define HALIDE_ATTRIBUTE_ALIGN(x) __attribute__((aligned(x)))
#endif
#endif

/** A runtime tag for a type in the halide type system. Can be ints,
* unsigned ints, or floats of various bit-widths (the 'bits'
* field). Can also be vectors of the same (by setting the 'lanes'
* field to something larger than one). This struct should be
* exactly 32-bits in size. */
struct halide_type_t {
    /** The basic type code: signed integer, unsigned integer, or floating point. */
#if __cplusplus >= 201103L
    HALIDE_ATTRIBUTE_ALIGN(1) halide_type_code_t code; // halide_type_code_t
#else
    HALIDE_ATTRIBUTE_ALIGN(1) uint8_t code; // halide_type_code_t
#endif

    /** The number of bits of precision of a single scalar value of this type. */
    HALIDE_ATTRIBUTE_ALIGN(1) uint8_t bits;

    /** How many elements in a vector. This is 1 for scalar types. */
    HALIDE_ATTRIBUTE_ALIGN(2) uint16_t lanes;

#ifdef __cplusplus

    /** Construct a runtime representation of a Halide type from:
     * code: The fundamental type from an enum.
     * bits: The bit size of one element.
     * lanes: The number of vector elements in the type. */
    halide_type_t(halide_type_code_t code, uint8_t bits, uint16_t lanes = 1)
            : code(code), bits(bits), lanes(lanes) {
    }

    /** Size in bytes for a single element, even if width is not 1, of this type. */
    size_t bytes() { return (bits + 7) / 8; }

#endif
};

template<typename T>
struct halide_c_type_to_name {
    static const bool known_type = false;
};

#define HALIDE_DECLARE_EXTERN_TYPE(TypeType, Type)                      \
    template<> struct halide_c_type_to_name<Type> {                     \
        static const bool known_type = true;                            \
        static halide_cplusplus_type_name name() {                      \
            return { halide_cplusplus_type_name::TypeType, #Type};      \
        }                                                               \
    }

#define HALIDE_DECLARE_EXTERN_SIMPLE_TYPE(T)     HALIDE_DECLARE_EXTERN_TYPE(Simple, T)
#define HALIDE_DECLARE_EXTERN_STRUCT_TYPE(T)     HALIDE_DECLARE_EXTERN_TYPE(Struct, T)
#define HALIDE_DECLARE_EXTERN_CLASS_TYPE(T)      HALIDE_DECLARE_EXTERN_TYPE(Class, T)
#define HALIDE_DECLARE_EXTERN_UNION_TYPE(T)      HALIDE_DECLARE_EXTERN_TYPE(Union, T)

HALIDE_DECLARE_EXTERN_SIMPLE_TYPE(bool);
HALIDE_DECLARE_EXTERN_SIMPLE_TYPE(int8_t);
HALIDE_DECLARE_EXTERN_SIMPLE_TYPE(uint8_t);
HALIDE_DECLARE_EXTERN_SIMPLE_TYPE(int16_t);
HALIDE_DECLARE_EXTERN_SIMPLE_TYPE(uint16_t);
HALIDE_DECLARE_EXTERN_SIMPLE_TYPE(int32_t);
HALIDE_DECLARE_EXTERN_SIMPLE_TYPE(uint32_t);
HALIDE_DECLARE_EXTERN_SIMPLE_TYPE(int64_t);
HALIDE_DECLARE_EXTERN_SIMPLE_TYPE(uint64_t);
HALIDE_DECLARE_EXTERN_SIMPLE_TYPE(float);
HALIDE_DECLARE_EXTERN_SIMPLE_TYPE(double);
/*HALIDE_DECLARE_EXTERN_STRUCT_TYPE(buffer_t);
HALIDE_DECLARE_EXTERN_STRUCT_TYPE(halide_buffer_t);
HALIDE_DECLARE_EXTERN_STRUCT_TYPE(halide_dimension_t);
HALIDE_DECLARE_EXTERN_STRUCT_TYPE(halide_device_interface_t);
HALIDE_DECLARE_EXTERN_STRUCT_TYPE(halide_filter_metadata_t);
*/
namespace Sympiler {
struct Expr;


/** Types in the halide type system. They can be ints, unsigned ints,
* or floats of various bit-widths (the 'bits' field). They can also
* be vectors of the same (by setting the 'lanes' field to something
* larger than one). Front-end code shouldn't use vector
* types. Instead vectorize a function. */
struct Type {
private:
    halide_type_t type;

public:
    /** Aliases for halide_type_code_t values for legacy compatibility
     * and to match the Halide internal C++ style. */
    // @{
    static const halide_type_code_t Int = halide_type_int;
    static const halide_type_code_t UInt = halide_type_uint;
    static const halide_type_code_t Float = halide_type_float;
    static const halide_type_code_t Handle = halide_type_handle;
    // @}

    /** The number of bytes required to store a single scalar value of this type. Ignores vector lanes. */
    int bytes() const { return (bits() + 7) / 8; }

    // Default ctor initializes everything to predictable-but-unlikely values
    Type() : type(Handle, 0, 0) { }


    /** Construct a runtime representation of a Halide type from:
     * code: The fundamental type from an enum.
     * bits: The bit size of one element.
     * lanes: The number of vector elements in the type. */
    Type(halide_type_code_t code, uint8_t bits, int lanes)
            : type(code, (uint8_t) bits, (uint16_t) lanes) {
    }

    /** Trivial copy constructor. */
    Type(const Type &that) = default;

    /** Type is a wrapper around halide_type_t with more methods for use
     * inside the compiler. This simply constructs the wrapper around
     * the runtime value. */
    Type(const halide_type_t &that) : type(that) { }

    /** Unwrap the runtime halide_type_t for use in runtime calls, etc.
     * Representation is exactly equivalent. */
    operator halide_type_t() const { return type; }

    /** Return the underlying data type of an element as an enum value. */
    halide_type_code_t code() const { return (halide_type_code_t) type.code; }

    /** Return the bit size of a single element of this type. */
    int bits() const { return type.bits; }

    /** Return the number of vector elements in this type. */
    int lanes() const { return type.lanes; }

    /** Return Type with same number of bits and lanes, but new_code for a type code. */
    Type with_code(halide_type_code_t new_code) const {
        return Type(new_code, bits(), lanes());
    }

    /** Return Type with same type code and lanes, but new_bits for the number of bits. */
    Type with_bits(uint8_t new_bits) const {
        return Type(code(), new_bits, lanes());
    }

    /** Return Type with same type code and number of bits,
     * but new_lanes for the number of vector lanes. */
    Type with_lanes(uint16_t new_lanes) const {
        return Type(code(), bits(), new_lanes);
    }

    /** Is this type boolean (represented as UInt(1))? */
    bool is_bool() const { return code() == UInt && bits() == 1; }

    /** Is this type a vector type? (lanes() != 1).
     * TODO(abadams): Decide what to do for lanes() == 0. */
    bool is_vector() const { return lanes() != 1; }

    /** Is this type a scalar type? (lanes() == 1).
     * TODO(abadams): Decide what to do for lanes() == 0. */
    bool is_scalar() const { return lanes() == 1; }

    /** Is this type a floating point type (float or double). */
    bool is_float() const { return code() == Float; }

    /** Is this type a signed integer type? */
    bool is_int() const { return code() == Int; }

    /** Is this type an unsigned integer type? */
    bool is_uint() const { return code() == UInt; }

    /** Is this type an opaque handle type (void *) */
    bool is_handle() const { return code() == Handle; }

    /** Compare two types for equality */
    bool operator==(const Type &other) const {
        return code() == other.code() && bits() == other.bits() && lanes() == other.lanes();
    }

    /** Compare two types for inequality */
    bool operator!=(const Type &other) const {
        return code() != other.code() || bits() != other.bits() || lanes() != other.lanes();
    }

    /** Produce the scalar type (that of a single element) of this vector type */
    Type element_of() const {
        return Type(code(), bits(), 1);
    }

    /** Can this type represent all values of another type? */
    bool can_represent(Type other) const;

    /** Can this type represent a particular constant? */
    // @{
    bool can_represent(double x) const;

    bool can_represent(int64_t x) const;

    bool can_represent(uint64_t x) const;
    // @}

    /** Check if an integer constant value is the maximum or minimum
     * representable value for this type. */
    // @{
    bool is_max(uint64_t) const;

    bool is_max(int64_t) const;

    bool is_min(uint64_t) const;

    bool is_min(int64_t) const;
    // @}

    /** Return an expression which is the maximum value of this type */
    Expr max() const;

    /** Return an expression which is the minimum value of this type */
    Expr min() const;
};

/** Constructing a signed integer type */
inline Type Int(int bits, int lanes = 1) {
    return Type(Type::Int, bits, lanes);
}

/** Constructing an unsigned integer type */
inline Type UInt(int bits, int lanes = 1) {
    return Type(Type::UInt, bits, lanes);
}

/** Construct a floating-point type */
inline Type Float(int bits, int lanes = 1) {
    return Type(Type::Float, bits, lanes);
}

/** Construct a boolean type */
inline Type Bool(int lanes = 1) {
    return UInt(1, lanes);
}

/** Construct a handle type */
inline Type Handle(int lanes = 1) {
    return Type(Type::Handle, 64, lanes);
}

/** Construct the halide equivalent of a C type */
/*template<typename T> Type type_of() {
return Type(halide_type_of<T>());
}*/


}
#endif //DSLPROJECT_TYPE_H
