//
// Created by kazem on 4/30/21.
//

#ifndef SYMPILER_PUBLIC_PARALLEL_PB_CHOLESKY_05_H
#define SYMPILER_PUBLIC_PARALLEL_PB_CHOLESKY_05_H

namespace sym_lib {
 namespace parsy {

  bool cholesky_left_par_05(int n, int *c, int *r, double *values,
                            size_t *lC, int *lR, size_t *Li_ptr, double *lValues,
                            int *blockSet, int supNo, double *timing,
#ifndef PRUNE
                            int *aTree, int *cT, int *rT, int *col2Sup,
#else
    int *prunePtr, int *pruneSet,
#endif
                            int nLevels, int *levelPtr, int *levelSet,
                            int nPar, int *parPtr, int *partition,
                            int chunk, int threads, int super_max,
                            int col_max, double *nodCost = NULL);
 }
}
#endif //SYMPILER_PUBLIC_PARALLEL_PB_CHOLESKY_05_H
