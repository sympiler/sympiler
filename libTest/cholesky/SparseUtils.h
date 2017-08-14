//
// Created by kazem on 7/27/17.
//

#ifndef CHOLOPENMP_SPARSEUTILS_H
#define CHOLOPENMP_SPARSEUTILS_H

/* ========================================================================== */
/* === cholmod_nnz ========================================================== */
/* ========================================================================== */

#include "def.h"

/* Return the number of entries in a sparse matrix.
 *
 * workspace: none
 * integer overflow cannot occur, since the matrix is already allocated.
 */

long int getNNZ
        (
                /* ---- input ---- */
                int ncol,
                int *Ap,
                int *Anz,
                int packed,
                /* --------------- */
                int &status
        )
{
    size_t nz ;
    int j ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    /*RETURN_IF_NULL_COMMON (EMPTY) ;
    RETURN_IF_NULL (A, EMPTY) ;
    RETURN_IF_XTYPE_INVALID (A, CHOLMOD_PATTERN, CHOLMOD_ZOMPLEX, EMPTY) ;
    Common->status = CHOLMOD_OK ;*/
    status = TRUE;
    /* ---------------------------------------------------------------------- */
    /* return nnz (A) */
    /* ---------------------------------------------------------------------- */

    if (packed)
    {
        //Ap = A->p ;
        //RETURN_IF_NULL (Ap, EMPTY) ;
        nz = Ap [ncol] ;
    }
    else
    {
        //Anz = A->nz ;
        //RETURN_IF_NULL (Anz, EMPTY) ;
        nz = 0 ;
        for (j = 0 ; j < ncol ; j++)
        {
            nz += MAX (0, Anz [j]) ;
        }
    }
    return (nz) ;
}


/* ========================================================================== */
/* === cholmod_add_size_t =================================================== */
/* ========================================================================== */

/* Safely compute a+b, and check for integer overflow.  If overflow occurs,
 * return 0 and set OK to FALSE.  Also return 0 if OK is FALSE on input. */

size_t add_size_t (size_t a, size_t b, int *ok)
{
    size_t s = a + b ;
    (*ok) = (*ok) && (s >= a) ;
    return ((*ok) ? s : 0) ;
}


/* ========================================================================== */
/* === cholmod_mult_size_t ================================================== */
/* ========================================================================== */

/* Safely compute a*k, where k should be small, and check for integer overflow.
 * If overflow occurs, return 0 and set OK to FALSE.  Also return 0 if OK is
 * FALSE on input. */

size_t mult_size_t(size_t a, size_t k, int *ok)
{
    size_t p = 0, s ;
    while (*ok)
    {
        if (k % 2)
        {
            p = p + a ;
            (*ok) = (*ok) && (p >= a) ;
        }
        k = k / 2 ;
        if (!k) return (p) ;
        s = a + a ;
        (*ok) = (*ok) && (s >= a) ;
        a = s ;
    }
    return (0) ;
}

#endif //CHOLOPENMP_SPARSEUTILS_H
