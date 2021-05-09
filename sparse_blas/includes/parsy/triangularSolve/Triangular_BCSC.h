//
// Created by kazem on 7/17/17.
//

#ifndef TRIANGOPENMP_TRIANGULAR_H
#define TRIANGOPENMP_TRIANGULAR_H

#include <cstddef>

namespace sym_lib {
 namespace parsy {

#define MKL_BLAS

/*
 * Forward solve blocked
 */
  int blockedLsolve(int n, size_t *Lp, int *Li, double *Lx, int NNZ,
                    size_t *Li_ptr, int *col2sup, int *sup2col, int supNo, double *x);
  int blockedLsolve_mrhs(int n, size_t *Lp, int *Li, double *Lx, int NNZ, size_t *Li_ptr,
    int *col2sup, int *sup2col,
    int supNo, double *x, int n_rhs, int max_col=-1);
/*
 * Backward solve blocked, u
 */
  int blockedLTsolve(int n, size_t *Lp, int *Li, double *Lx, int NNZ, size_t *Li_ptr,
                     int *col2sup, int *sup2col, int supNo, double *x);
  int
  blockedLTsolve_mrhs(int n, size_t *Lp, int *Li, double *Lx, int NNZ, size_t *Li_ptr, int *col2sup, int *sup2col, int supNo,
                      double *x, int m_rhs, int max_col=-1);
/*
 *
 */
  int LeveledBlockedLTsolve_update(int n, size_t *Lp, int *Li, double *Lx,
                                   int NNZ, size_t *Li_ptr,
                                   int *col2sup, int *sup2col,
                                   int supNo, double *x,
                                   int levels, int *levelPtr,
                                   int *levelSet,
                                   int chunk,
                                   bool *marked, double *ws_dbl = NULL);


/*
 * Blocked Serial code
 */
//#define BLAS2
  int blockedPrunedLSolve(int n, int *Lp, int *Li, double *Lx, int NNZ, int *Li_ptr,
                          int *BPSet, int PBSetSize, int *sup2col, int supNo, double *x);

/*
 * Parallel Blocked
 */

  int leveledBlockedLsolve(int n, size_t *Lp, int *Li, double *Lx, int NNZ, size_t *Li_ptr,
                           int *col2sup, int *sup2col, int supNo, double *x,
                           int levels, int *levelPtr, int *levelSet, int chunk);

/*
 * Parallel Blocked with masked cols
 */

  int leveledBlockedLsolve_update(int n, size_t *Lp, int *Li, double *Lx,
                                  int NNZ, size_t *Li_ptr,
                                  int *col2sup, int *sup2col,
                                  int supNo, double *x,
                                  int levels, int *levelPtr, int *levelSet,
                                  int chunk, bool *mask, double *ws_dbl = NULL);


/*
 * Parallel H2 Blocked
 */
//#define MKL_BLAS
  int H2LeveledBlockedLsolve(int n, size_t *Lp, int *Li, double *Lx, int NNZ, size_t
  *Li_ptr,
                             int *col2sup, int *sup2col, int supNo, double *x,
                             int levels, int *levelPtr, int *levelSet,
                             int parts, int *parPtr, int *partition, int chunk);
  int
  H2LeveledBlockedLsolve_mrhs(int n, size_t *Lp, int *Li, double *Lx, int NNZ, size_t *Li_ptr, int *col2sup, int *sup2col,
                              int supNo, double *x, int n_rhs, int levels, int *levelPtr, int *levelSet, int parts, int *parPtr,
                              int *partition, int chunk, int max_col=-1);

  int H2LeveledBlockedLsolve_update(int n, size_t *Lp, int *Li, double *Lx, int NNZ, size_t
  *Li_ptr,
                                    int *col2sup, int *sup2col, int supNo, double *x,
                                    int levels, int *levelPtr, int *levelSet,
                                    int parts, int *parPtr, int *partition,
                                    int chunk, bool *mask, double *ws_dbl = NULL);
/*
 * Backward solve blocked, unit triangular
 */
//#define MKL_BLAS
  int H2LeveledBlockedLTsolve(int n, size_t *Lp, int *Li, double *Lx, int NNZ, size_t *Li_ptr,
                              int *col2sup, int *sup2col, int supNo, double *x,
                              int levels, int *levelPtr, int *levelSet,
                              int parts, int *parPtr, int *partition, int chunk);
  int
  H2LeveledBlockedLTsolve_mrhs(int n, size_t *Lp, int *Li, double *Lx, int NNZ, size_t *Li_ptr, int *col2sup, int *sup2col,
                               int supNo, double *x, int m_rhs, int levels, int *levelPtr, int *levelSet, int parts, int *parPtr,
                               int *partition, int chunk, int max_col=-1);

  int H2LeveledBlockedLTsolve_update(int n, size_t *Lp, int *Li,
                                     double *Lx, int NNZ, size_t *Li_ptr,
                                     int *col2sup, int *sup2col, int supNo,
                                     double *x,
                                     int levels, int *levelPtr, int *levelSet,
                                     int parts, int *parPtr, int *partition,
                                     int chunk, bool *mask, double *ws_dbl = NULL);



/*
 * Parallel H2 Blocked with peeling
 */
#undef MKL_BLAS

  int H2LeveledBlockedLsolve_Peeled(int n, size_t *Lp, int *Li, double *Lx, int NNZ, size_t
  *Li_ptr,
                                    int *col2sup, int *sup2col, int supNo, double *x,
                                    int levels, int *levelPtr, int *levelSet,
                                    int parts, int *parPtr, int *partition, int chunk, int threads);
 }
}
#endif //TRIANGOPENMP_TRIANGULAR_H
