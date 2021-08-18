//
// Created by kazem on 10/9/19.
//

#ifndef FUSION_SPARSEBLASLIB_H
#define FUSION_SPARSEBLASLIB_H

#include "def.h"
namespace sym_lib {

 ///
 /// \param n number of rows in matrix A
 /// \param Ap row pointer of matrix A
 /// \param Ai column indices of matrix A
 /// \param Ax nonzero values of matrix A
 /// \param x  multiplicand vector of size equal to number of columns
 /// \param y product vector of size n
 void
 spmv_csr(int n, const int *Ap, const int *Ai, const double *Ax,
          const double *x, double *y);

 /// full SpMV a*(Ax) + by = z
 /// \param n
 /// \param Ap
 /// \param Ai
 /// \param Ax
 /// \param x
 /// \param y
 /// \param z resulting vector
 /// \param a
 /// \param b
 void spmv_csr(int n, const int *Ap, const int *Ai, const double *Ax,
               const double *x, const double *y, double *z, double a, double b);

 void
 spmv_csr_parallel(int n, const int *Ap, const int *Ai, const double *Ax,
                   const double *x, double *y);

 void spmv_csr_block(int n, const int *Ap, const int *Ai, const double *Ax,
                     int nodes, const int *supernodes, const double *x,
                     double *y);

 ///
 /// \param n number of columns in matrix A
 /// \param Ap column pointers of matrix A
 /// \param Ai row indices of matrix A
 /// \param Ax nonzero values of matrix A
 /// \param x multiplicand vector of size n
 /// \param y product vector of size equal to number of rows
 void
 spmv_csc(int n, const int *Ap, const int *Ai, const double *Ax,
          const double *x, double *y);

 /// full SpMV: a*(Ax) + by = z
 /// \param n
 /// \param Ap
 /// \param Ai
 /// \param Ax
 /// \param x
 /// \param y
 /// \param z
 /// \param a
 /// \param b
 void spmv_csc(int n, const int *Ap, const int *Ai, const double *Ax,
               const double *x, const double *y, double *z, double a, double b);

 void spmv_csc(int n, const int *Ap, const int *Ai, const double *Ax,
               const double *x, const double *y, double *z, double *a,
               double *b);


 /// Uses a simple LBC paritioning with reduction
 /// \param n_parts
 /// \param sets
 /// \param ids
 /// \param n
 /// \param Ap
 /// \param Ai
 /// \param Ax
 /// \param x
 /// \param y
 void spmv_csc_lbc(int n_parts, const int *sets, const int *ids,
                   int n, const int *Ap,
                   const int *Ai,
                   const double *Ax, const double *x, double *y);
 void spmv_csc_lbc(const int *set_ptr, const int *sets, const int *ids,
                   int n, const int *Ap,
                   const int *Ai,
                   const double *Ax, const double *x, double *y);

 ///
 /// \param n
 /// \param Ap
 /// \param Ai
 /// \param Ax
 /// \param x
 /// \param y
 void spmv_csc_parallel(int n, const int *Ap, const int *Ai, const double *Ax,
                        const double *x, double *y);


 ///
 /// \param n
 /// \param Ap
 /// \param Ai
 /// \param Ax
 /// \param x
 /// \param y
 /// \param nodes
 /// \param supernodes
 void
 spmv_csc_block(int n, const int *Ap, const int *Ai, const double *Ax,
                const double *x, double *y, int nodes,
                const int *supernodes);


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


 /// Does a lsolve give a general A matrix with diagonal localtions
 /// \param n
 /// \param Lp General A, pointers
 /// \param Li
 /// \param Lx
 /// \param diags location of diagonal indices
 /// \param x
 void sptrsv_csr_ilu(int n, const int *Lp, const int *Li, const double *Lx,
                     const int*diags, double *x);
 void sptrsv_csr_ilu_lbc(int n, const int *Lp,
                         const int *Li, const double *Lx, const int*diags,
                         double *x, const int level_no, const int *level_ptr,
                         const int *par_ptr, const int *partition);


 void
 sptrsv_csr_levelset(int n, const int *Lp, const int *Li, const double *Lx,
                     double *x,
                     int levels, const int *levelPtr, const int *levelSet);

 void
 sptrsv_csr_lbc(int n, int *Lp, int *Li, double *Lx, double *x,
                int level_no, int *level_ptr,
                int *par_ptr, int *partition);

 // tmp size is nume_thread times n
 void
 sptrsv_csr_lbc_redundant(int n, int *Lp, int *Li, double *Lx,
                               double *x, double *b,
                               int level_no, int *level_ptr,
                               int *par_ptr, int *partition, 
                               double *tmp);
 ///
 /// Forward-substitution on lower triangular matrix L
 /// \param n rank of matrix L
 /// \param Lp column pointers of matrix L
 /// \param Li row indices of matrix L
 /// \param Lx nonzero values of matrix L
 /// \param x initial RHS and outputs as result
 void
 sptrsv_csc(int n, int *Lp, int *Li, double *Lx, double *x);


 void
 sptrsv_csc_levelset(int n, int *Lp, int *Li, double *Lx, double *x, int levels,
                     const int *levelptr, const int *levelset);
 void
 sptrsv_csc_lbc(int n, int *Lp, int *Li, double *Lx, double *x, int level_no, int *level_ptr,
                int *par_ptr, int *partition);


 /// ========= Blocked versions ===========

 ///
 /// \param n
 /// \param Lp
 /// \param Li
 /// \param nrows
 /// \param Lx
 /// \param x
 /// \param levelPtr
 /// \param levels
 /// \param n_lev
 /// \param supernodes
 /// \param sup2node
 /// \param tempvecs
 void sptrsv_csc_levelset_block(int n, const int *Lp, const int *Li,
                                const int *nrows, double *Lx, double *x,
                                const int *levelPtr,
                                const int *levels, int n_lev,
                                const int *supernodes, const int *sup2node,
                                double **tempvecs);

 ///
 /// \param n
 /// \param Lp
 /// \param Li
 /// \param nrows
 /// \param Lx
 /// \param x
 /// \param num_nodes
 /// \param supernodes
 void sptrsv_csc_block(int n, const int *Lp, const int *Li, const int *nrows,
                         double *Lx, double *x, int num_nodes,
                         const int *supernodes);

#ifdef MKL
 void sptrsv_csc_block_mkl(int n, const int *Lp, const int *Li, const int *nrows,
                           double *Lx, double *x, int num_nodes,
                           const int *supernodes);
#endif

 /////////////////////// Add
///
/// \param A
/// \param B
/// \param alpha
/// \param beta
/// \param sort : if true, sorts index arrays
/// \return
 CSC* add(CSC *A, CSC *B, double alpha, double beta, bool sort=true);



 /////////////////////////// ICholesky0

 /// Takes the lower part of a SPD matrix in CSC
 /// or the upper part is CSR and factorize it with incomplete Cholesky0
 /// \param n
 /// \param val
 /// \param rowPtr
 /// \param rowIdx
 void spic0_csr(int n, double *val, int *rowPtr, int *rowIdx);
 void spic0_csr_lbc(int n, double *val, int *rowPtr, int *colIdx,
                  int level_no, int *level_ptr,
                  int *par_ptr, int *partition);


  /*
   * Performs an incomplete Cholesky decomposition
   * on a given matrix (Ap, Ai, Ax), i.e.
   * stored in compressed column format,
   * and produces an approximate for L, which are
   * stored in column compressed format.
   * (Lx) : OUT : the nonzero values of the approximate L factor
   */
 bool spic0_csc_left(int n, const int* Ap, const int* Ai,
                const double* Ax, const int *prunePtr,
                const int *pruneSet,
                double* Lx);
 bool spic0_csc_left_LBC(int n, const int* Ap, const int* Ai,
                         const double* Ax, const int *prunePtr,
                         const int *pruneSet,const int nLevels,
                         const int *levelPtr,
                         const int *parPtr, const int *partition,
                         double* Lx);
 /////////////////////////// ILU0
 /// spilu0_csr computes the incomplete LU factorization of a matrix.
 /// The matrix A is assumed to be stored in compressed row format.
 /// \param n the matrix dimension
 /// \param nnz number of nonzero values
 /// \param Ap row pointer with size of n
 /// \param Ai column indices with size of nz_num
 /// \param Ax values with size of nonzeros
 /// \param A_diag location of diagonal with size of n
 /// \param l The factorized matix with size of nnz
 void spilu0_csr(int n, int nnz, const int *Ap, const int *Ai,
                 const double *Ax, const int *A_diag, double *l);
 //tempvec with size of n*num_threads, initilized with -1 (EMPTY)
 void spilu0_csr_lbc(int n, int nnz, const int *Ap, const int *Ai,
                     const double *Ax, const int *A_diag, double *l,
                     int level_no, const int *level_ptr,
                     const int *par_ptr, const int *partition,
                     int *tempvecs);



 ////////////////////////// LDL Simplicial (non-supernodal)

  /*
    Performs a LDL decomposition on a given matrix (c, r, values), i.e.
   * stored in compressed column format, and produces L, which are
   * stored in column compressed format.
   * (n, c, r, values) : IN : input matrix
   * (lC, lR) : IN : The column and rwo sparsity patterns of L
   * (lValues) : OUT : the nonzero values of the L factor
   * (pruneSet, prunePtr) : IN : the row sparsity pattern of the factor L
   * ws in: n, ws_int size of 3*n
   */

 int ldl_left_simplicial_02(int n, int *c, int *r, double *values,
                            int *cT, int *rT,
                            int *lC, int *lR, double *&lValues,
                            double *d,
#ifdef PPRUNE
   int *prunePtr, int *pruneSet,
#endif
                            int *eTree,
                            double *ws, int *ws_int);

 int ldl_parallel_left_simplicial_01 (int n, int* c, int* r, double* values,
                                      int* cT, int* rT,
                                      int* lC, int* lR, double* &lValues,
                                      double *d,
#if 0
   int *prunePtr, int *pruneSet,
#endif
                                      int *eTree,
                                      int nLevels, int *levelPtr,
                                      int *parPtr, int *partition);

 ///////////////////////// LT SOLVE

 /// Solves L^T x = b for x
 /// \param n
 /// \param Lp
 /// \param Li
 /// \param Lx
 /// \param x in/out pre-initilized with right-hand-side
 /// \return
 int ltsolve(int n, int *Lp, int *Li, double *Lx, double *x);



 void scale_vec_vec(int n, double *vec1, double *vec2);


//////////////////////////// Diagonal scaling
/// D*A
/// \param A in/out csc
/// \param d input diagonal, size n
 void mat_premult_diag(CSC *A, const double *d);
 void mat_premult_diag(CSR *A, const double *d);
 void mat_premult_diag_parallel(CSR *A, const double *d);
 void mat_premult_diag_parallel(CSC *A, const double *d);

 /// A*D^T
 /// \param A in/out
 /// \param d in
 void mat_postmult_diag(CSC *A, const double *d);
 void mat_postmult_diag(CSR *A, const double *d);
 void mat_postmult_diag_parallel(CSR *A, const double *d);
 void mat_postmult_diag_parallel(CSC *A, const double *d);

}


#endif //FUSION_SPARSEBLASLIB_H
