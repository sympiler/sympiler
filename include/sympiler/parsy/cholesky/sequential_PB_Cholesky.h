//
// Created by kazem on 4/30/21.
//

#ifndef SYMPILER_PUBLIC_SEQUENTIAL_PB_CHOLESKY_H
#define SYMPILER_PUBLIC_SEQUENTIAL_PB_CHOLESKY_H
namespace sym_lib {
 namespace parsy {

  bool cholesky_left_sn_07(int n, int *c, int *r, double *values,
                           size_t *lC, int *lR, size_t *Li_ptr, double *lValues,
                           int *blockSet, int supNo, double *timing,

#ifndef PRUNE
                           int *aTree, int *cT, int *rT, int *col2Sup,
#else
    int *prunePtr, int *pruneSet,
#endif
                           int super_max,
                           int col_max);


 }
}
#endif //SYMPILER_PUBLIC_SEQUENTIAL_PB_CHOLESKY_H
