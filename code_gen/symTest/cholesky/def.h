//
// Created by kazem on 7/25/17.
//

#ifndef CHOLOPENMP_DEF_H
#define CHOLOPENMP_DEF_H

#include <assert.h>
#include "stddef.h"
#include "sys/param.h"
#define CS_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define CS_MIN(a,b) (((a) < (b)) ? (a) : (b))
#define CS_FLIP(i) (-(i)-2)
#define CS_UNFLIP(i) (((i) < 0) ? CS_FLIP(i) : (i))
#define CS_MARKED(w,j) (w [j] < 0)
#define CS_MARK(w,j) { w [j] = CS_FLIP (w [j]) ; }
#define HEAD(k,j) (ata ? head [k] : j)
#define NEXT(J)   (ata ? next [J] : -1)
#define  CS_MIN(a,b) (((a) < (b)) ? (a) : (b))
#define EMPTY -1
#define TRUE 1
#define FALSE 0
#define ASSERT(expression) (assert (expression))


#define CHOLMOD_AMD 2		/* use minimum degree (AMD) */

#define SIGN(x) (((x) < 0) ? (-1) : (((x) > 0) ? 1 : 0))


#define MATRIX_PATTERN 0	/* pattern only, no numerical values */
#define MATRIX_REAL 1		/* a real matrix */
#define Int_max INT_MAX


typedef struct sympiler_CSC
{
    size_t nrow ;	/* the matrix is nrow-by-ncol */
    size_t ncol ;
    size_t nzmax ;	/* maximum number of entries in the matrix */

    /* pointers to int or SuiteSparse_long: */
    int *p ;		/* p [0..ncol], the column pointers */
    int *i ;		/* i [0..nzmax-1], the row indices */

    /* for unpacked matrices only: */
    int *nz ;
    double *x ;		/* size nzmax or 2*nzmax, if present */

    int stype ;

    int xtype ;		/* pattern, real, complex, or zomplex */
    int sorted ;	/* TRUE if columns are sorted, FALSE otherwise */
    int packed ;	/* TRUE if packed (nz ignored), FALSE if unpacked
			 * (nz is required) */

} CSC ;

typedef struct Sympiler_BCSC
{

    size_t n ;		/* L is n-by-n */


    int *Perm ;	/* size n, permutation used */
    int *ColCount ;	/* size n, column counts for simplicial L */


    int *p ;		/* p [0..ncol], the column pointers */
    int *i ;		/* i [0..nzmax-1], the row indices */
    double *x ;		/* x [0..nzmax-1], the numerical values */
    void *nz ;

    size_t nsuper ;	/* number of supernodes */
    size_t ssize ;	/* size of s, integer part of supernodes */
    size_t xsize ;	/* size of x, real part of supernodes */

    int *super ;	/* size nsuper+1, first col in each supernode */ //SUP2Col
    int *pi ;		/* size nsuper+1, pointers to integer patterns */ //row indices
    int *px ;		/* size nsuper+1, pointers to real parts */ //
    int *s ;		/* size ssize, integer part of supernodes */


    int xtype ; /* pattern, real, complex, or zomplex */


} BCSC ;

int allocateLC(BCSC *L, int sw){
    int sNo = L->nsuper;
    if(sw){
        L->super = new int[sNo+1]();
        L->p = new int[sNo+1]();
        L->pi = new int[sNo+1]();
        L->i = new int[sNo+1]();
        L->px = new int[L->xsize]();
        L->s = new int[L->ssize]();
        L->x = new double[L->xsize]();
    }else{
        delete []L->super;
        delete []L->p;
        delete []L->pi;
        delete []L->i;
        delete []L->px;
        delete []L->s;
        delete []L->x;
    }

    return 1;
}

int allocateAC(CSC *A, int nrow, int nnz, int sytpe, int sw){
    if(sw){
        A->nrow=A->ncol=nrow;
        A->nzmax = nnz;
        A->stype = sytpe;
        A->xtype = MATRIX_REAL;//TODO removed later
        A->packed = TRUE; // Always
        A->p = new int[nrow+1]();
        A->i = new int[nnz]();
        A->x = new double[nnz]();
        A->nz = NULL;

    }else{
        delete []A->p;
        delete []A->i;
        delete []A->x;
    }

    return 1;
}
#endif //CHOLOPENMP_DEF_H
