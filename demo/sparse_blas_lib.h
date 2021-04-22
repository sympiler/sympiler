//
// Created by kazem on 10/9/19.
//

#ifndef FUSION_SPARSEBLASLIB_H
#define FUSION_SPARSEBLASLIB_H

#include "def.h"
namespace sym_lib {

 ///////////////////////// SPTRSV

 ///
 /// Forward-substitution on lower triangular matrix L
 /// \param n rank of matrix L
 /// \param Lp row pointers of matrix L
 /// \param Li column indices of matrix L
 /// \param Lx nonzero values of matrix L
 /// \param x initial RHS and outputs as result
 void
 sptrsv_csr(int n, int *Lp, int *Li, double *Lx, double *x);


 void
 sptrsv_csr_levelset(int n, const int *Lp, const int *Li, const double *Lx,
                     double *x,
                     int levels, const int *levelPtr, const int *levelSet);

 void
 sptrsv_csr_lbc(int n, int *Lp, int *Li, double *Lx, double *x,
                int level_no, int *level_ptr,
                int *par_ptr, int *partition);


}


#endif //FUSION_SPARSEBLASLIB_H
