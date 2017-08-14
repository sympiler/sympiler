#include <iostream>
#include <math.h>
#include <float.h>
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <myBLAS.h>
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
#ifdef __cplusplus
extern "C" {
#endif
int32_t trns(int32_t n1, int32_t *Mp1, int32_t *Mi1, double *Mx1, int32_t *Mip1, int32_t n2, int32_t *Mp2, int32_t *Mi2, double *Mx2, int32_t *, double *D2, int32_t *setPtr, int32_t *setVal, int32_t ub, int32_t *blk2Col, int32_t blkNo) {
 {
  int16_t _0 = (int16_t)(164);
  double _1 = n1 * _0;
  int64_t _2 = _1;
  if ((_2 > ((int64_t(1) << 31) - 1)) || ((_2 * sizeof(double)) > ((int64_t(1) << 31) - 1)))
  {
      return -1;
  } // overflow test tempVec
  bool _3 = (bool)(1);
  int64_t _4 = (int64_t)(_3 ? _2 : 0);
  int64_t _5 = _4;
  double *tempVec = (double *)sympiler_malloc(NULL, sizeof(double)*_5);
  int32_t _6 = Mp2[0] + 1;
  for (int copy = Mp2[0]; copy < Mp2[_6]; copy++)
  {
   D2[Mi2[copy]] = Mx2[copy];
  } // for copy
  for (int trnsf0P = 0; trnsf0P < ub; trnsf0P++)
  {
   int32_t _7 = setVal[trnsf0P] - 1;
   bool _8 = setVal[trnsf0P] != 0;
   int32_t _9 = (int32_t)(_8 ? blk2Col[_7] : 0);
   int32_t _10 = Mip1[blk2Col[setVal[trnsf0P]]] - Mip1[_9];
   int32_t _11 = blk2Col[setVal[trnsf0P]] - _9;
   trsm_blas(_10, _11, &Mx1[Mp1[_9]], &D2[_9]);
   int32_t _12 = _10 - _11;
   int32_t _13 = Mp1[_9] + _11;
   matvec(_10, _12, _11, &Mx1[_13], &D2[_9], &tempVec[0]);
   int32_t _14 = Mip1[_9] + _11;
   for (int trnsf1 = _14; trnsf1 < Mip1[blk2Col[setVal[trnsf0P]]]; trnsf1++)
   {
    int32_t _15 = setVal[trnsf0P] - 1;
    bool _16 = setVal[trnsf0P] != 0;
    int32_t _17 = (int32_t)(_16 ? blk2Col[_15] : 0);
    int32_t _18 = blk2Col[setVal[trnsf0P]] - _17;
    int32_t _19 = Mip1[_17] + _18;
    int16_t _20 = trnsf1 - _19;
    double _21 = D2[Mi1[trnsf1]] - tempVec[_20];
    D2[Mi1[trnsf1]] = _21;
    int32_t _22 = setVal[trnsf0P] - 1;
    bool _23 = setVal[trnsf0P] != 0;
    int32_t _24 = (int32_t)(_23 ? blk2Col[_22] : 0);
    int32_t _25 = blk2Col[setVal[trnsf0P]] - _24;
    int32_t _26 = Mip1[_24] + _25;
    int16_t _27 = trnsf1 - _26;
    int16_t _28 = (int16_t)(0);
    tempVec[_27] = _28;
   } // for trnsf1
  } // for trnsf0P
 } // alloc tempVec
 return 0;
}
#ifdef __cplusplus
}  // extern "C"
#endif
