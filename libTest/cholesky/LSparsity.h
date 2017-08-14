//
// Created by kazem on 7/26/17.
//

#ifndef CHOLOPENMP_LSPARSITY_H
#define CHOLOPENMP_LSPARSITY_H

#include "def.h"
#include "PostOrder.h"
#include "Etree.h"
#include "ColumnCount.h"
#include "Inspection_BlockC.h"

//From CHOLMOD
#include "amd.h"
int permute_matrices
        (
                /* ---- input ---- */
                CSC *A,	/* matrix to permute */
                int ordering,	/* ordering method used */
                int *Perm,		/* fill-reducing permutation */
                int *fset,		/* subset of 0:(A->ncol)-1 */
                size_t fsize,	/* size of fset */
                int do_rowcolcounts,/* if TRUE, compute both S and F.  If FALSE, only
			 * S is needed for the symmetric case, and only F for
			 * the unsymmetric case */

                /* ---- output --- */
                CSC **A1_handle,	    /* see comments below for A1, A2, S, F */
                CSC **A2_handle,
                CSC **S_handle,
                CSC **F_handle,
                /* --------------- */
                //cholmod_common *Common
                int &status
        )
{
    CSC *A1, *A2, *S, *F ;

    *A1_handle = NULL ;
    *A2_handle = NULL ;
    *S_handle = NULL ;
    *F_handle = NULL ;
    A1 = NULL ;
    A2 = NULL ;

    if (ordering == CHOLMOD_NATURAL)
    {

        /* ------------------------------------------------------------------ */
        /* natural ordering of A */
        /* ------------------------------------------------------------------ */

        if (A->stype < 0)
        {
            /* symmetric lower case: A already in lower form, so S=A' */
            /* workspace: Iwork (nrow) */
            A2 = ptranspose(A, 0, NULL, NULL, 0, status) ;
            F = A ;
            S = A2 ;
        }


    }
    else
    {

        /* ------------------------------------------------------------------ */
        /* A is permuted */
        /* ------------------------------------------------------------------ */

        if (A->stype < 0)
        {
            /* symmetric lower case: S = tril (A (p,p))' and F = S' */
            /* workspace: Iwork (2*nrow) */
            A2 = ptranspose(A, 0, Perm, NULL, 0, status) ;
            S = A2 ;
            /* workspace: Iwork (nrow) */
            if (do_rowcolcounts)
            {
                /* F not needed for symmetric case if do_rowcolcounts FALSE */
                A1 = ptranspose(A2, 0, NULL, NULL, 0, status) ;
            }
            F = A1 ;
        }
    }

    /* If any cholmod_*transpose fails, one or more matrices will be NULL */
    *A1_handle = A1 ;
    *A2_handle = A2 ;
    *S_handle = S ;
    *F_handle = F ;
    return (TRUE) ;
}



int buildingSymbolicStructures
        (
                /* ---- input ---- */
                CSC *A,    /* matrix to analyze */
                int ordering,    /* ordering method used */
                int *Perm,        /* size n, fill-reducing permutation to analyze */
                int *fset,        /* subset of 0:(A->ncol)-1 */
                size_t fsize,    /* size of fset */
                /* ---- output --- */
                int *Parent,    /* size n, elimination tree */
                int *Post,        /* size n, postordering of elimination tree */
                int *ColCount,    /* size n, nnz in each column of L */

                /* ---- workspace  */
                int *First,        /* size n workspace for cholmod_postorder */
                int *Level,        /* size n workspace for cholmod_postorder */
                /* --------------- */
                // cholmod_common *Common
                int &status
        )
{
    CSC *A1, *A2, *S, *F ;
    int n, ok, do_rowcolcounts ;
    int fl, aatfl, lnz;

    /* check inputs */
 //   RETURN_IF_NULL_COMMON (FALSE) ;
 //   RETURN_IF_NULL (A, FALSE) ;

    n = A->nrow ;

    do_rowcolcounts = (ColCount != NULL) ;

    /* permute A according to Perm and fset */
    ok = permute_matrices (A, ordering, Perm, fset, fsize, do_rowcolcounts,
                           &A1, &A2, &S, &F, status) ;

    /* find etree of S (symmetric upper/lower case) or F (unsym case) */
    /* workspace: symmmetric: Iwork (nrow), unsym: Iwork (nrow+ncol) */
    ok = ok && etreeC (A->stype ? S:F, Parent, status) ;


    /* postorder the etree (required by cholmod_rowcolcounts) */
    /* workspace: Iwork (2*nrow) */
    ok = ok && postOrderC(Parent, n, NULL, Post, status) == n ;


    if (do_rowcolcounts)
    {
        ok = ok && rowcolcounts(A->stype ? F:S, fset, fsize, Parent,
                                Post, NULL, ColCount, First, Level,
                                fl, aatfl, lnz, status) ;
    }

    allocateAC(A1,0,0,0,FALSE);
    delete A1;
    allocateAC(A2,0,0,0,FALSE);
    delete A2;
    return (ok) ;
}


/* ========================================================================== */
/* === cholmod_analyze_p2 =================================================== */
/* ========================================================================== */

/* Ordering and analysis for sparse Cholesky or sparse QR.  */

BCSC* buildingInspectionSets
        (
                /* ---- input ---- */
                int for_whom,       /* FOR_SPQR     (0): for SPQR but not GPU-accelerated
                           FOR_CHOLESKY (1): for Cholesky (GPU or not)
                           FOR_SPQRGPU  (2): for SPQR with GPU acceleration */
                CSC *A,    /* matrix to order and analyze */
                int *UserPerm,    /* user-provided permutation, size A->nrow */
                int *fset,        /* subset of 0:(A->ncol)-1 */
                int *nrelax,
                double *zrelax,
                size_t fsize,    /* size of fset */
                //Outputs
                int *&prunePtr,
                int *&pruneSet,
                /* --------------- */
                int status
        )
{
    int *First, *Level, *Work4n, *Cmember, *CParent, *ColCount, *Lperm, *Parent,
            *Post, *Perm, *Lparent, *Lcolcount ;
    BCSC *L ;
    int n, ncol, k, ordering, method, nmethods,  default_strategy, uncol,
            skip_analysis, skip_best ;
    int amd_backup ;
    size_t s ;
    int ok = TRUE ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */


    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    n = A->nrow ;
    ncol = A->ncol ;
    uncol = (A->stype == 0) ? (A->ncol) : 0 ;

    /* ---------------------------------------------------------------------- */
    /* set the default strategy */
    /* ---------------------------------------------------------------------- */


    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */
    s = mult_size_t(n, 4, &ok);
    if(!ok)
        return FALSE;
    Work4n = new int[s]() ;
   // Work4n += 2*((size_t) n) + uncol ;
    Parent = Work4n ;
    First  = Work4n + n ;
    Level  = Work4n + 2*((size_t) n) ;
    Post   = Work4n + 3*((size_t) n) ;


    /* ---------------------------------------------------------------------- */
    /* allocate more workspace, and an empty simplicial symbolic factor */
    /* ---------------------------------------------------------------------- */

    L = new BCSC;
    L->Perm     = new int[n]();
    L->ColCount = new int[n]();
    L->ordering = CHOLMOD_AMD;
    Lparent  = new int[n]();
    Lperm = L->Perm ;
    Lcolcount = L->ColCount ;
/*
 * Ordering
 */
    double info[20]={0};
    double Control[2];
    Control [0] = 10; //TODO check later //AMD_Dense
    Control [1] = TRUE; //AMD_AGGRESSIVE
    amd_order(ncol,A->p,A->i,Lperm,NULL,info);

    if (!buildingSymbolicStructures(A, L->ordering, Lperm, fset, fsize,
                                    Lparent, Post, Lcolcount, First,
                                    Level, status))
    {
        /* out of memory, or method failed */
        //FREE_WORKSPACE_AND_RETURN ;
        return FALSE;
    }

    /* postorder the etree, weighted by the column counts */

    if (postOrderC(Lparent, n, Lcolcount, Post,status) == n)
    {
        /* use First and Level as workspace [ */
        int *Wi = First, *InvPost = Level ;
        int newchild, oldchild, newparent, oldparent ;

        for (k = 0 ; k < n ; k++)
        {
            Wi [k] = Lperm [Post [k]] ;
        }
        for (k = 0 ; k < n ; k++)
        {
            Lperm [k] = Wi [k] ;
        }

        for (k = 0 ; k < n ; k++)
        {
            Wi [k] = Lcolcount [Post [k]] ;
        }
        for (k = 0 ; k < n ; k++)
        {
            Lcolcount [k] = Wi [k] ;
        }
        for (k = 0 ; k < n ; k++)
        {
            InvPost [Post [k]] = k ;
        }

        /* updated Lparent needed only for supernodal case */
        for (newchild = 0 ; newchild < n ; newchild++)
        {
            oldchild = Post [newchild] ;
            oldparent = Lparent [oldchild] ;
            newparent = (oldparent == EMPTY) ? EMPTY : InvPost [oldparent] ;
            Wi [newchild] = newparent ;
        }
        for (k = 0 ; k < n ; k++)
        {
            Lparent [k] = Wi [k] ;
        }
        /* done using Iwork as workspace ] */

        /* L is now postordered, no longer in natural ordering */
        if (L->ordering == CHOLMOD_NATURAL)
        {
            L->ordering = CHOLMOD_POSTORDERED ;
        }
    }

    /* ---------------------------------------------------------------------- */
    /* supernodal analysis, if requested or if selected automatically */
    /* ---------------------------------------------------------------------- */


    int *Sparent = new int[n]();
    CSC *S, *F, *A2, *A1 ;

    permute_matrices (A, L->ordering, Lperm, fset, fsize, TRUE,
                      &A1, &A2, &S, &F, status) ;
    super_symbolic2(for_whom, S, F, Lparent, L,nrelax,zrelax,Sparent, status) ;

    prunePtr = new int[L->nsuper+1]();
    pruneSet = new int[L->ssize];
    int *col2Sup = new int[F->nrow];
    // Computing col2sup
    for (int i = 0; i < L->nsuper; ++i) {
        int k1 = L->super[i];
        int k2 = L->super[i+1];
        ASSERT(k1 <= F->nrow && k2 <= F->nrow);
        for (int j = k1; j < k2; ++j) {
            col2Sup[j] = i;
        }
    }
    //PRINT1 (("status %d\n", Common->status)) ;
    getBlockedPruneSet(L->nsuper,S->p,S->i,col2Sup,Sparent,L->super,
                       prunePtr,pruneSet);

    delete []col2Sup;
    allocateAC(A1,0,0,0,FALSE);
    allocateAC(A2,0,0,0,FALSE);

    /* ---------------------------------------------------------------------- */
    /* free temporary matrices and workspace, and return result L */
    /* ---------------------------------------------------------------------- */

    delete []Sparent;
    delete []Work4n;
    return L;
}
//#endif
#endif //CHOLOPENMP_LSPARSITY_H
