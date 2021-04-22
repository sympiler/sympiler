#include <iostream>
#include <math.h>
#include <float.h>
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <cholUtils.h>
#ifndef HALIDE_ATTRIBUTE_ALIGN
  #ifdef _MSC_VER
    #define HALIDE_ATTRIBUTE_ALIGN(x) __declspec(align(x))
  #else
    #define HALIDE_ATTRIBUTE_ALIGN(x) __attribute__((aligned(x)))
  #endif
#endif
#ifndef BUFFER_T_DEFINED
#define BUFFER_T_DEFINED
#include <stdbool.h>
#include <stdint.h>
typedef struct buffer_t {
    uint64_t dev;
    uint8_t* host;
    int32_t extent[4];
    int32_t stride[4];
    int32_t min[4];
    int32_t elem_size;
    HALIDE_ATTRIBUTE_ALIGN(1) bool host_dirty;
    HALIDE_ATTRIBUTE_ALIGN(1) bool dev_dirty;
    HALIDE_ATTRIBUTE_ALIGN(1) uint8_t _padding[10 - sizeof(void *)];
} buffer_t;
#endif
#define __user_context_ NULL
#define BLOCKED
struct halide_filter_metadata_t;
extern "C" {
void *sympiler_malloc(void *ctx, size_t s){return(malloc(s));}
void sympiler_free(void *ctx, void *ptr){free(ptr);};
}

#ifdef _WIN32
float roundf(float);
double round(double);
#else
inline float asinh_f32(float x) {return asinhf(x);}
inline float acosh_f32(float x) {return acoshf(x);}
inline float atanh_f32(float x) {return atanhf(x);}
inline double asinh_f64(double x) {return asinh(x);}
inline double acosh_f64(double x) {return acosh(x);}
inline double atanh_f64(double x) {return atanh(x);}
#endif
inline float sqrt_f32(float x) {return sqrtf(x);}
inline float sin_f32(float x) {return sinf(x);}
inline float asin_f32(float x) {return asinf(x);}
inline float cos_f32(float x) {return cosf(x);}
inline float acos_f32(float x) {return acosf(x);}
inline float tan_f32(float x) {return tanf(x);}
inline float atan_f32(float x) {return atanf(x);}
inline float sinh_f32(float x) {return sinhf(x);}
inline float cosh_f32(float x) {return coshf(x);}
inline float tanh_f32(float x) {return tanhf(x);}
inline float hypot_f32(float x, float y) {return hypotf(x, y);}
inline float exp_f32(float x) {return expf(x);}
inline float log_f32(float x) {return logf(x);}
inline float pow_f32(float x, float y) {return powf(x, y);}
inline float floor_f32(float x) {return floorf(x);}
inline float ceil_f32(float x) {return ceilf(x);}
inline float round_f32(float x) {return roundf(x);}

inline double sqrt_f64(double x) {return sqrt(x);}
inline double sin_f64(double x) {return sin(x);}
inline double asin_f64(double x) {return asin(x);}
inline double cos_f64(double x) {return cos(x);}
inline double acos_f64(double x) {return acos(x);}
inline double tan_f64(double x) {return tan(x);}
inline double atan_f64(double x) {return atan(x);}
inline double sinh_f64(double x) {return sinh(x);}
inline double cosh_f64(double x) {return cosh(x);}
inline double tanh_f64(double x) {return tanh(x);}
inline double hypot_f64(double x, double y) {return hypot(x, y);}
inline double exp_f64(double x) {return exp(x);}
inline double log_f64(double x) {return log(x);}
inline double pow_f64(double x, double y) {return pow(x, y);}
inline double floor_f64(double x) {return floor(x);}
inline double ceil_f64(double x) {return ceil(x);}
inline double round_f64(double x) {return round(x);}

inline float nan_f32() {return NAN;}
inline float neg_inf_f32() {return -INFINITY;}
inline float inf_f32() {return INFINITY;}
inline bool is_nan_f32(float x) {return x != x;}
inline bool is_nan_f64(double x) {return x != x;}
inline float float_from_bits(uint32_t bits) {
 union {
  uint32_t as_uint;
  float as_float;
 } u;
 u.as_uint = bits;
 return u.as_float;
}
inline int64_t make_int64(int32_t hi, int32_t lo) {
    return (((int64_t)hi) << 32) | (uint32_t)lo;
}
inline double make_float64(int32_t i0, int32_t i1) {
    union {
        int32_t as_int32[2];
        double as_double;
    } u;
    u.as_int32[0] = i0;
    u.as_int32[1] = i1;
    return u.as_double;
}

template<typename T> T max(T a, T b) {if (a > b) return a; return b;}
template<typename T> T min(T a, T b) {if (a < b) return a; return b;}
template<typename A, typename B> A reinterpret(B b) {A a; memcpy(&a, &b, sizeof(a)); return a;}

double one [2]={1.0,0.}, zero [2]={0.,0.};
int sw = false, lb1 = 0, ub1 = 0; 
double *cur; int info=0; 
#ifdef __cplusplus
extern "C" {
#endif
int32_t Chol(int32_t n1, int32_t *Mp1, int32_t *Mi1, double *Mx1, int32_t *Mip1, int32_t n2, int32_t *Mp2, int32_t *Mi2, double *Mx2, int32_t *Mip2, int32_t *tempVec, double *finger, int32_t *setPtr, int32_t *setVal, int32_t ub, int32_t *blk2Col, int32_t blkNo) {
 int32_t _0 = blkNo + 1;
 for (int Cholf0 = 1; Cholf0 < _0; Cholf0++)
 {
  int32_t _1 = Cholf0 - 1;
  bool _2 = Cholf0 != 0;
  int32_t _3 = (int32_t)(_2 ? blk2Col[_1] : 0);
  for (int c0 = Mip2[_3]; c0 < Mip2[blk2Col[Cholf0]]; c0++)
  {
   int32_t _4 = Cholf0 - 1;
   bool _5 = Cholf0 != 0;
   int32_t _6 = (int32_t)(_5 ? blk2Col[_4] : 0);
   int32_t _7 = c0 - Mip2[_6];
   tempVec[Mi2[c0]] = _7;
  } // for c0
  int32_t _8 = Cholf0 - 1;
  bool _9 = Cholf0 != 0;
  int32_t _10 = (int32_t)(_9 ? blk2Col[_8] : 0);
  for (int c1 = _10; c1 < blk2Col[Cholf0]; c1++)
  {
   int32_t _11 = c1 + 1;
   for (int c2 = Mp1[c1]; c2 < Mp1[_11]; c2++)
   {
    int32_t _12 = Mp2[c1] + tempVec[Mi1[c2]];
    Mx2[_12] = Mx1[c2];
   } // for c2
  } // for c1
  int32_t _13 = Cholf0 - 1;
  for (int u0P = setPtr[_13]; u0P < setPtr[Cholf0]; u0P++)
  {
   ub1 = 0;
   lb1 = 0;
   sw = 1;
   int32_t _14 = setVal[u0P] + 1;
   for (int u1 = Mip2[blk2Col[setVal[u0P]]]; u1 < Mip2[blk2Col[_14]]; u1++)
   {
    int32_t _15 = Cholf0 - 1;
    bool _16 = Cholf0 != 0;
    int32_t _17 = (int32_t)(_16 ? blk2Col[_15] : 0);
    bool _18 = Mi2[u1] >= _17;
    bool _19 = _18 && sw;
    if (_19)
    {
     sw = 0;
     int32_t _20 = u1 - Mip2[blk2Col[setVal[u0P]]];
     lb1 = _20;
    } // if _19
    int32_t _21 = Cholf0 - 1;
    bool _22 = Cholf0 != 0;
    int32_t _23 = (int32_t)(_22 ? blk2Col[_21] : 0);
    int32_t _24 = blk2Col[Cholf0] - _23;
    int32_t _25 = _23 + _24;
    bool _26 = Mi2[u1] < _25;
    bool _27 = !(sw);
    bool _28 = _26 && _27;
    if (_28)
    {
     int32_t _29 = u1 - Mip2[blk2Col[setVal[u0P]]];
     ub1 = _29;
    } // if _28
   } // for u1
   int32_t _30 = ub1 - lb1;
   int32_t _31 = _30 + 1;
   int32_t _32 = setVal[u0P] + 1;
   int32_t _33 = blk2Col[_32] - blk2Col[setVal[u0P]];
   int32_t _34 = Mp2[blk2Col[setVal[u0P]]] + lb1;
   int32_t _35 = Mip2[blk2Col[_32]] - Mip2[blk2Col[setVal[u0P]]];
   int32_t _36 = _35 - lb1;
   syrk(&_31, &_33, one, &Mx2[_34], &_35, zero, finger, &_36);
   bool _37 = _31 > 0;
   if (_37)
   {
    int32_t _38 = setVal[u0P] + 1;
    int32_t _39 = Mip2[blk2Col[_38]] - Mip2[blk2Col[setVal[u0P]]];
    int32_t _40 = _39 - lb1;
    int32_t _41 = ub1 - lb1;
    int32_t _42 = _41 + 1;
    int32_t _43 = _40 - _42;
    int32_t _44 = blk2Col[_38] - blk2Col[setVal[u0P]];
    int32_t _45 = ub1 + 1;
    int32_t _46 = Mp2[blk2Col[setVal[u0P]]] + _45;
    int32_t _47 = Mp2[blk2Col[setVal[u0P]]] + lb1;
    gemm(&_43, &_42, &_44, one, &Mx2[_46], &_39, &Mx2[_47], &_39, zero, &finger[_42], &_40);
   } // if _37
   int32_t _48 = ub1 - lb1;
   int32_t _49 = _48 + 1;
   for (int uc0 = 0; uc0 < _49; uc0++)
   {
    int32_t _50 = setVal[u0P] + 1;
    int32_t _51 = Mip2[blk2Col[_50]] - Mip2[blk2Col[setVal[u0P]]];
    int32_t _52 = _51 - lb1;
    for (int uc1 = uc0; uc1 < _52; uc1++)
    {
     int32_t _53 = Mip2[blk2Col[setVal[u0P]]] + lb1;
     int32_t _54 = _53 + uc0;
     int32_t _55 = Cholf0 - 1;
     bool _56 = Cholf0 != 0;
     int32_t _57 = (int32_t)(_56 ? blk2Col[_55] : 0);
     int32_t _58 = Mip2[blk2Col[Cholf0]] - Mip2[_57];
     int32_t _59 = tempVec[Mi2[_54]] * _58;
     int32_t _60 = _53 + uc1;
     int32_t _61 = _59 + tempVec[Mi2[_60]];
     int32_t _62 = _61 + Mp2[_57];
     int32_t _63 = setVal[u0P] + 1;
     int32_t _64 = Mip2[blk2Col[_63]] - Mip2[blk2Col[setVal[u0P]]];
     int32_t _65 = _64 - lb1;
     int32_t _66 = uc0 * _65;
     int32_t _67 = _66 + uc1;
     double _68 = Mx2[_62] - finger[_67];
     Mx2[_62] = _68;
    } // for uc1
   } // for uc0
  } // for u0P
  int32_t _69 = Cholf0 - 1;
  bool _70 = Cholf0 != 0;
  int32_t _71 = (int32_t)(_70 ? blk2Col[_69] : 0);
  int32_t _72 = blk2Col[Cholf0] - _71;
  int32_t _73 = Mip2[blk2Col[Cholf0]] - Mip2[_71];
  potrf(&_72, &Mx2[Mp2[_71]], &_73, &info);
  int32_t _74 = _73 - _72;
  int32_t _75 = Mp2[_71] + _72;
  trsm(&_74, &_72, one, &Mx2[Mp2[_71]], &_73, &Mx2[_75], &_73);
 } // for Cholf0
 return 0;
}
#ifdef __cplusplus
}  // extern "C"
#endif
