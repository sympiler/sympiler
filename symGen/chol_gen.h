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
int32_t Chol(int32_t n1, int32_t *Mp1, int32_t *Mi1, double *Mx1, int32_t *Mip1, int32_t n2, int32_t *Mp2, int32_t *Mi2, double *Mx2, int32_t *Mip2, int32_t *setPtr, int32_t *setVal, int32_t ub, int32_t *blk2Col, int32_t blkNo) {
 {
  int64_t _0 = n1;
  if ((_0 > ((int64_t(1) << 31) - 1)) || ((_0 * sizeof(double)) > ((int64_t(1) << 31) - 1)))
  {
      return -1;
  } // overflow test finger
  bool _1 = (bool)(1);
  int64_t _2 = (int64_t)(_1 ? _0 : 0);
  int64_t _3 = _2;
  double *finger = (double *)sympiler_malloc(__user_context_, sizeof(double)*_3);
  {
   int64_t _4 = n1;
   if ((_4 > ((int64_t(1) << 31) - 1)) || ((_4 * sizeof(double)) > ((int64_t(1) << 31) - 1)))
   {
        return -1;
   } // overflow test tempVec
   bool _5 = (bool)(1);
   int64_t _6 = (int64_t)(_5 ? _4 : 0);
   int64_t _7 = _6;
   double *tempVec = (double *)sympiler_malloc(__user_context_, sizeof(double)*_7);
   for (int Cholf0 = 1; Cholf0 < n1; Cholf0++)
   {
    int16_t _8 = Cholf0 - 1;
    bool _9 = Cholf0 != 0;
    int32_t _10 = (int32_t)(_9 ? blk2Col[_8] : 0);
    for (int c0 = Mip2[_10]; c0 < Mip2[blk2Col[Cholf0]]; c0++)
    {
     int16_t _11 = Cholf0 - 1;
     bool _12 = Cholf0 != 0;
     int32_t _13 = (int32_t)(_12 ? blk2Col[_11] : 0);
     int32_t _14 = c0 - Mip2[_13];
     tempVec[Mi2[c0]] = _14;
    } // for c0
    int16_t _15 = Cholf0 - 1;
    bool _16 = Cholf0 != 0;
    int32_t _17 = (int32_t)(_16 ? blk2Col[_15] : 0);
    for (int c1 = _17; c1 < blk2Col[Cholf0]; c1++)
    {
     int16_t _18 = (int16_t)(1);
     int32_t _19 = c1 + _18;
     for (int c2 = Mp1[c1]; c2 < Mp1[_19]; c2++)
     {
      int32_t _20 = Mp2[c1] + tempVec[Mi1[c2]];
      Mx2[_20] = Mx1[c2];
     } // for c2
    } // for c1
    int16_t _21 = Cholf0 - 1;
    bool _22 = Cholf0 != 0;
    int32_t _23 = (int32_t)(_22 ? blk2Col[_21] : 0);
    int32_t _24 = _23 + 1;
    for (int u0P = setPtr[_23]; u0P < setPtr[_24]; u0P++)
    {
     sw = 1;
     int32_t _25 = setVal[u0P] + 1;
     for (int u1 = Mip2[blk2Col[setVal[u0P]]]; u1 < Mip2[blk2Col[_25]]; u1++)
     {
      int16_t _26 = Cholf0 - 1;
      bool _27 = Cholf0 != 0;
      int32_t _28 = (int32_t)(_27 ? blk2Col[_26] : 0);
      bool _29 = Mi2[u1] >= _28;
      bool _30 = _29 && sw;
      if (_30)
      {
       int32_t _31 = u1 - Mip2[blk2Col[setVal[u0P]]];
       lb1 = _31;
      } // if _30
      int16_t _32 = Cholf0 - 1;
      bool _33 = Cholf0 != 0;
      int32_t _34 = (int32_t)(_33 ? blk2Col[_32] : 0);
      int32_t _35 = blk2Col[Cholf0] - _34;
      int32_t _36 = _34 + _35;
      bool _37 = Mi2[u1] < _36;
      bool _38 = !(sw);
      bool _39 = _37 && _38;
      if (_39)
      {
       int32_t _40 = u1 - Mip2[blk2Col[setVal[u0P]]];
       ub1 = _40;
      } // if _39
     } // for u1
     int32_t _41 = ub1 - lb1;
     int32_t _42 = _41 + 1;
     int16_t _43 = Cholf0 - 1;
     bool _44 = Cholf0 != 0;
     int32_t _45 = (int32_t)(_44 ? blk2Col[_43] : 0);
     int32_t _46 = blk2Col[Cholf0] - _45;
     int32_t _47 = Mp2[blk2Col[setVal[u0P]]] + lb1;
     int32_t _48 = setVal[u0P] + 1;
     int32_t _49 = Mip2[blk2Col[_48]] - Mip2[blk2Col[setVal[u0P]]];
     int32_t _50 = _49 - lb1;
     syrk(&_42, &_46, one, &Mx2[_47], &_49, zero, finger, &_50);
     bool _51 = _42 > 0;
     if (_51)
     {
      int32_t _52 = setVal[u0P] + 1;
      int32_t _53 = Mip2[blk2Col[_52]] - Mip2[blk2Col[setVal[u0P]]];
      int32_t _54 = _53 - lb1;
      int32_t _55 = ub1 - lb1;
      int32_t _56 = _55 + 1;
      int32_t _57 = _54 - _56;
      int16_t _58 = Cholf0 - 1;
      bool _59 = Cholf0 != 0;
      int32_t _60 = (int32_t)(_59 ? blk2Col[_58] : 0);
      int32_t _61 = blk2Col[Cholf0] - _60;
      int32_t _62 = ub1 + 1;
      int32_t _63 = Mp2[blk2Col[setVal[u0P]]] + _62;
      int32_t _64 = Mp2[blk2Col[setVal[u0P]]] + lb1;
      gemm(&_57, &_56, &_61, one, &Mx2[_63], &_53, &Mx2[_64], &_53, zero, &finger[_56], &_54);
     } // if _51
     int32_t _65 = ub1 - lb1;
     int32_t _66 = _65 + 1;
     for (int uc0 = 0; uc0 < _66; uc0++)
     {
      int32_t _67 = setVal[u0P] + 1;
      int32_t _68 = Mip2[blk2Col[_67]] - Mip2[blk2Col[setVal[u0P]]];
      int32_t _69 = _68 - lb1;
      for (int uc1 = uc0; uc1 < _69; uc1++)
      {
       int32_t _70 = Mip2[blk2Col[setVal[u0P]]] + lb1;
       int32_t _71 = _70 + uc0;
       int16_t _72 = Cholf0 - 1;
       bool _73 = Cholf0 != 0;
       int32_t _74 = (int32_t)(_73 ? blk2Col[_72] : 0);
       int32_t _75 = Mip2[blk2Col[Cholf0]] - Mip2[_74];
       int32_t _76 = tempVec[Mi1[_71]] * _75;
       int32_t _77 = _70 + uc1;
       int32_t _78 = _76 + tempVec[Mi1[_77]];
       int32_t _79 = setVal[u0P] + 1;
       int32_t _80 = Mip2[blk2Col[_79]] - Mip2[blk2Col[setVal[u0P]]];
       int32_t _81 = _80 - lb1;
       int32_t _82 = uc0 * _81;
       int32_t _83 = _82 + uc1;
       double _84 = cur[_78] - tempVec[_83];
       cur[_78] = _84;
      } // for uc1
     } // for uc0
    } // for u0P
    int16_t _85 = Cholf0 - 1;
    bool _86 = Cholf0 != 0;
    int32_t _87 = (int32_t)(_86 ? blk2Col[_85] : 0);
    int32_t _88 = blk2Col[Cholf0] - _87;
    int32_t _89 = Mip2[blk2Col[Cholf0]] - Mip2[_87];
    potrf(&_88, &Mx1[Mp1[_87]], &_89, &info);
    int32_t _90 = _89 - _88;
    double _91 = Mx1[Mp1[_87]] + _88;
    trsm(&_90, &_88, one, &Mx1[Mp1[_87]], &_89, &_91, &_89);
   } // for Cholf0
  } // alloc tempVec
 } // alloc finger
 return 0;
}
#ifdef __cplusplus
}  // extern "C"
#endif
