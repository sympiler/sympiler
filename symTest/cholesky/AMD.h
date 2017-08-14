//
// Created by kazem on 7/31/17.
//

#ifndef CHOLOPENMP_AMD_H
#define CHOLOPENMP_AMD_H

#include "def.h"

/* ========================================================================= */
/* === AMD_aat ============================================================= */
/* ========================================================================= */

/* ------------------------------------------------------------------------- */
/* AMD, Copyright (c) Timothy A. Davis,					     */
/* Patrick R. Amestoy, and Iain S. Duff.  See ../README.txt for License.     */
/* email: DrTimothyAldenDavis@gmail.com                                      */
/* ------------------------------------------------------------------------- */

/* AMD_aat:  compute the symmetry of the pattern of A, and count the number of
 * nonzeros each column of A+A' (excluding the diagonal).  Assumes the input
 * matrix has no errors, with sorted columns and no duplicates
 * (AMD_valid (n, n, Ap, Ai) must be AMD_OK, but this condition is not
 * checked).
 */

size_t AMD_aat	/* returns nz in A+A' */
(
        int n,
const int Ap [ ],
const int Ai [ ],
        int Len [ ],	/* Len [j]: length of column j of A+A', excl diagonal*/
int Tp [ ],		/* workspace of size n */
double Info [ ]
)
{
    int p1, p2, p, i, j, pj, pj2, k, nzdiag, nzboth, nz ;
    double sym ;
    size_t nzaat ;

    #ifndef NDEBUG
//    AMD_debug_init ("AMD AAT") ;
    for (k = 0 ; k < n ; k++) Tp [k] = EMPTY ;
//    ASSERT (AMD_valid (n, n, Ap, Ai) == AMD_OK) ;
    #endif

/*    if (Info != (double *) NULL)
    {
    *//* clear the Info array, if it exists *//*
    for (i = 0 ; i < AMD_INFO ; i++)
    {
    Info [i] = EMPTY ;
    }
    Info [AMD_STATUS] = AMD_OK ;
    }*/

    for (k = 0 ; k < n ; k++)
    {
    Len [k] = 0 ;
    }

    nzdiag = 0 ;
    nzboth = 0 ;
    nz = Ap [n] ;

    for (k = 0 ; k < n ; k++)
    {
    p1 = Ap [k] ;
    p2 = Ap [k+1] ;
//    AMD_DEBUG2 (("\nAAT Column: "ID" p1: "ID" p2: "ID"\n", k, p1, p2)) ;

    /* construct A+A' */
    for (p = p1 ; p < p2 ; )
    {
    /* scan the upper triangular part of A */
    j = Ai [p] ;
    if (j < k)
    {
    /* entry A (j,k) is in the strictly upper triangular part,
     * add both A (j,k) and A (k,j) to the matrix A+A' */
    Len [j]++ ;
    Len [k]++ ;
//    AMD_DEBUG3 (("    upper ("ID","ID") ("ID","ID")\n", j,k, k,j));
    p++ ;
    }
    else if (j == k)
    {
    /* skip the diagonal */
    p++ ;
    nzdiag++ ;
    break ;
    }
    else /* j > k */
    {
    /* first entry below the diagonal */
    break ;
    }
    /* scan lower triangular part of A, in column j until reaching
     * row k.  Start where last scan left off. */
    ASSERT (Tp [j] != EMPTY) ;
    ASSERT (Ap [j] <= Tp [j] && Tp [j] <= Ap [j+1]) ;
    pj2 = Ap [j+1] ;
    for (pj = Tp [j] ; pj < pj2 ; )
    {
    i = Ai [pj] ;
    if (i < k)
    {
    /* A (i,j) is only in the lower part, not in upper.
     * add both A (i,j) and A (j,i) to the matrix A+A' */
    Len [i]++ ;
    Len [j]++ ;
//    AMD_DEBUG3 (("    lower ("ID","ID") ("ID","ID")\n",
    i,j, j,i)) ;
    pj++ ;
    }
    else if (i == k)
    {
    /* entry A (k,j) in lower part and A (j,k) in upper */
    pj++ ;
    nzboth++ ;
    break ;
    }
    else /* i > k */
    {
    /* consider this entry later, when k advances to i */
    break ;
    }
    }
    Tp [j] = pj ;
    }
    /* Tp [k] points to the entry just below the diagonal in column k */
    Tp [k] = p ;
    }

    /* clean up, for remaining mismatched entries */
    for (j = 0 ; j < n ; j++)
    {
    for (pj = Tp [j] ; pj < Ap [j+1] ; pj++)
    {
    i = Ai [pj] ;
    /* A (i,j) is only in the lower part, not in upper.
     * add both A (i,j) and A (j,i) to the matrix A+A' */
    Len [i]++ ;
    Len [j]++ ;
//    AMD_DEBUG3 (("    lower cleanup ("ID","ID") ("ID","ID")\n",
    i,j, j,i)) ;
    }
    }

    /* --------------------------------------------------------------------- */
    /* compute the symmetry of the nonzero pattern of A */
    /* --------------------------------------------------------------------- */

    /* Given a matrix A, the symmetry of A is:
     *	B = tril (spones (A), -1) + triu (spones (A), 1) ;
     *  sym = nnz (B & B') / nnz (B) ;
     *  or 1 if nnz (B) is zero.
     */

    if (nz == nzdiag)
    {
    sym = 1 ;
    }
    else
    {
    sym = (2 * (double) nzboth) / ((double) (nz - nzdiag)) ;
    }

    nzaat = 0 ;
    for (k = 0 ; k < n ; k++)
    {
    nzaat += Len [k] ;
    }

    //AMD_DEBUG1 (("AMD nz in A+A', excluding diagonal (nzaat) = %g\n",
    //(double) nzaat)) ;
    //AMD_DEBUG1 (("   nzboth: "ID" nz: "ID" nzdiag: "ID" symmetry: %g\n",
    //nzboth, nz, nzdiag, sym)) ;

/*    if (Info != (double *) NULL)
    {
    Info [AMD_STATUS] = AMD_OK ;
    Info [AMD_N] = n ;
    Info [AMD_NZ] = nz ;
    Info [AMD_SYMMETRY] = sym ;	    *//* symmetry of pattern of A *//*
    Info [AMD_NZDIAG] = nzdiag ;	    *//* nonzeros on diagonal of A *//*
    Info [AMD_NZ_A_PLUS_AT] = nzaat ;   *//* nonzeros in A+A' *//*
    }*/

    return (nzaat) ;
}


/* ========================================================================= */
/* === AMD_preprocess ====================================================== */
/* ========================================================================= */

/* AMD_preprocess does not check its input for errors or allocate workspace.
 * On input, the condition (AMD_valid (n,n,Ap,Ai) != AMD_INVALID) must hold.
 */

void AMD_preprocess
        (
                int n,		/* input matrix: A is n-by-n */
                const int Ap [ ],	/* size n+1 */
                const int Ai [ ],	/* size nz = Ap [n] */

                /* output matrix R: */
                int Rp [ ],		/* size n+1 */
                int Ri [ ],		/* size nz (or less, if duplicates present) */

                int W [ ],		/* workspace of size n */
                int Flag [ ]	/* workspace of size n */
        )
{

    /* --------------------------------------------------------------------- */
    /* local variables */
    /* --------------------------------------------------------------------- */

    int i, j, p, p2 ;

//    ASSERT (AMD_valid (n, n, Ap, Ai) != AMD_INVALID) ;

    /* --------------------------------------------------------------------- */
    /* count the entries in each row of A (excluding duplicates) */
    /* --------------------------------------------------------------------- */

    for (i = 0 ; i < n ; i++)
    {
        W [i] = 0 ;		/* # of nonzeros in row i (excl duplicates) */
        Flag [i] = EMPTY ;	/* Flag [i] = j if i appears in column j */
    }
    for (j = 0 ; j < n ; j++)
    {
        p2 = Ap [j+1] ;
        for (p = Ap [j] ; p < p2 ; p++)
        {
            i = Ai [p] ;
            if (Flag [i] != j)
            {
                /* row index i has not yet appeared in column j */
                W [i]++ ;	    /* one more entry in row i */
                Flag [i] = j ;	    /* flag row index i as appearing in col j*/
            }
        }
    }

    /* --------------------------------------------------------------------- */
    /* compute the row pointers for R */
    /* --------------------------------------------------------------------- */

    Rp [0] = 0 ;
    for (i = 0 ; i < n ; i++)
    {
        Rp [i+1] = Rp [i] + W [i] ;
    }
    for (i = 0 ; i < n ; i++)
    {
        W [i] = Rp [i] ;
        Flag [i] = EMPTY ;
    }

    /* --------------------------------------------------------------------- */
    /* construct the row form matrix R */
    /* --------------------------------------------------------------------- */

    /* R = row form of pattern of A */
    for (j = 0 ; j < n ; j++)
    {
        p2 = Ap [j+1] ;
        for (p = Ap [j] ; p < p2 ; p++)
        {
            i = Ai [p] ;
            if (Flag [i] != j)
            {
                /* row index i has not yet appeared in column j */
                Ri [W [i]++] = j ;  /* put col j in row i */
                Flag [i] = j ;	    /* flag row index i as appearing in col j*/
            }
        }
    }

/*#ifndef NDEBUG
    ASSERT (AMD_valid (n, n, Rp, Ri) == AMD_OK) ;
    for (j = 0 ; j < n ; j++)
    {
        ASSERT (W [j] == Rp [j+1]) ;
    }
#endif*/
}

/* ========================================================================= */
/* === AMD_order =========================================================== */
/* ========================================================================= */


int AMD_order
(
        int n,
const int Ap [ ],
const int Ai [ ],
        int P [ ],
double Control [ ],
double Info [ ]
)
{
    int *Len, *S, nz, i, *Pinv, info, status, *Rp, *Ri, *Cp, *Ci, ok ;
    size_t nzaat, slen ;
    double mem = 0 ;

    #ifndef NDEBUG
//    AMD_debug_init ("amd") ;
    #endif

    /* clear the Info array, if it exists */
    info = Info != (double *) NULL ;
    if (info)
    {
    for (i = 0 ; i < AMD_INFO ; i++)
    {
    Info [i] = EMPTY ;
    }
    Info [AMD_N] = n ;
    Info [AMD_STATUS] = AMD_OK ;
    }

    /* make sure inputs exist and n is >= 0 */
    if (Ai == (int *) NULL || Ap == (int *) NULL || P == (int *) NULL || n < 0)
    {
    if (info) Info [AMD_STATUS] = AMD_INVALID ;
    return (AMD_INVALID) ;	    /* arguments are invalid */
    }

    if (n == 0)
    {
    return (AMD_OK) ;	    /* n is 0 so there's nothing to do */
    }

    nz = Ap [n] ;
    if (info)
    {
    Info [AMD_NZ] = nz ;
    }
    if (nz < 0)
    {
    if (info) Info [AMD_STATUS] = AMD_INVALID ;
    return (AMD_INVALID) ;
    }

    /* check if n or nz will cause size_t overflow */
    if (((size_t) n) >= SIZE_T_MAX / sizeof (int)
    || ((size_t) nz) >= SIZE_T_MAX / sizeof (int))
    {
    if (info) Info [AMD_STATUS] = AMD_OUT_OF_MEMORY ;
    return (AMD_OUT_OF_MEMORY) ;	    /* problem too large */
    }

    /* check the input matrix:	AMD_OK, AMD_INVALID, or AMD_OK_BUT_JUMBLED */
    status = AMD_valid (n, n, Ap, Ai) ;

    if (status == AMD_INVALID)
    {
    if (info) Info [AMD_STATUS] = AMD_INVALID ;
    return (AMD_INVALID) ;	    /* matrix is invalid */
    }

    /* allocate two size-n integer workspaces */
//    Len  = SuiteSparse_malloc (n, sizeof (int)) ;
//    Pinv = SuiteSparse_malloc (n, sizeof (int)) ;
    mem += n ;
    mem += n ;
    if (!Len || !Pinv)
    {
    /* :: out of memory :: */
//    SuiteSparse_free (Len) ;
//    SuiteSparse_free (Pinv) ;
    if (info) Info [AMD_STATUS] = AMD_OUT_OF_MEMORY ;
    return (AMD_OUT_OF_MEMORY) ;
    }

    if (status == AMD_OK_BUT_JUMBLED)
    {
    /* sort the input matrix and remove duplicate entries */
//    AMD_DEBUG1 (("Matrix is jumbled\n")) ;
//    Rp = SuiteSparse_malloc (n+1, sizeof (int)) ;
//    Ri = SuiteSparse_malloc (nz,  sizeof (int)) ;
    mem += (n+1) ;
    mem += MAX (nz,1) ;
    if (!Rp || !Ri)
    {
    /* :: out of memory :: */
/*    SuiteSparse_free (Rp) ;
    SuiteSparse_free (Ri) ;
    SuiteSparse_free (Len) ;
    SuiteSparse_free (Pinv) ;*/
    if (info) Info [AMD_STATUS] = AMD_OUT_OF_MEMORY ;
    return (AMD_OUT_OF_MEMORY) ;
    }
    /* use Len and Pinv as workspace to create R = A' */
    AMD_preprocess (n, Ap, Ai, Rp, Ri, Len, Pinv) ;
    Cp = Rp ;
    Ci = Ri ;
    }
    else
    {
    /* order the input matrix as-is.  No need to compute R = A' first */
    Rp = NULL ;
    Ri = NULL ;
    Cp = (int *) Ap ;
    Ci = (int *) Ai ;
    }

    /* --------------------------------------------------------------------- */
    /* determine the symmetry and count off-diagonal nonzeros in A+A' */
    /* --------------------------------------------------------------------- */

    nzaat = AMD_aat (n, Cp, Ci, Len, P, Info) ;
 //   AMD_DEBUG1 (("nzaat: %g\n", (double) nzaat)) ;
 //   ASSERT ((MAX (nz-n, 0) <= nzaat) && (nzaat <= 2 * (size_t) nz)) ;

    /* --------------------------------------------------------------------- */
    /* allocate workspace for matrix, elbow room, and 6 size-n vectors */
    /* --------------------------------------------------------------------- */

    S = NULL ;
    slen = nzaat ;			/* space for matrix */
    ok = ((slen + nzaat/5) >= slen) ;	/* check for size_t overflow */
    slen += nzaat/5 ;			/* add elbow room */
    for (i = 0 ; ok && i < 7 ; i++)
    {
    ok = ((slen + n) > slen) ;	/* check for size_t overflow */
    slen += n ;			/* size-n elbow room, 6 size-n work */
    }
    mem += slen ;
    ok = ok && (slen < SIZE_T_MAX / sizeof (int)) ; /* check for overflow */
    ok = ok && (slen < int_MAX) ;	/* S[i] for int i must be OK */
    if (ok)
    {
    S = SuiteSparse_malloc (slen, sizeof (int)) ;
    }
    //AMD_DEBUG1 (("slen %g\n", (double) slen)) ;
    /*if (!S)
    {
    *//* :: out of memory :: (or problem too large) *//*
    SuiteSparse_free (Rp) ;
    SuiteSparse_free (Ri) ;
    SuiteSparse_free (Len) ;
    SuiteSparse_free (Pinv) ;
    if (info) Info [AMD_STATUS] = AMD_OUT_OF_MEMORY ;
    return (AMD_OUT_OF_MEMORY) ;
    }
    if (info)
    {
    *//* memory usage, in bytes. *//*
    Info [AMD_MEMORY] = mem * sizeof (int) ;
    }*/

    /* --------------------------------------------------------------------- */
    /* order the matrix */
    /* --------------------------------------------------------------------- */

    AMD_1 (n, Cp, Ci, P, Pinv, Len, slen, S, Control, Info) ;

    /* --------------------------------------------------------------------- */
    /* free the workspace */
    /* --------------------------------------------------------------------- */

    /*SuiteSparse_free (Rp) ;
    SuiteSparse_free (Ri) ;
    SuiteSparse_free (Len) ;
    SuiteSparse_free (Pinv) ;
    SuiteSparse_free (S) ;*/
    if (info) Info [AMD_STATUS] = status ;
    return (status) ;	    /* successful ordering */
}

#endif //CHOLOPENMP_AMD_H
