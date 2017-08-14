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
        else if (A->stype > 0)
        {
            /* symmetric upper case: F = pattern of triu (A)', S = A */
            /* workspace: Iwork (nrow) */
            if (do_rowcolcounts)
            {
                /* F not needed for symmetric case if do_rowcolcounts FALSE */
                A1 = ptranspose(A, 0, NULL, fset, fsize, status) ;
            }
            F = A1 ;
            S = A ;
        }
        else
        {
            /* unsymmetric case: F = pattern of A (:,f)',  S = A */
            /* workspace: Iwork (nrow if no fset, MAX(nrow,ncol) if fset) */
            A1 = ptranspose(A, 0, NULL, fset, fsize, status) ;
            F = A1 ;
            S = A ;
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
        else if (A->stype > 0)
        {
            /* symmetric upper case: F = triu (A (p,p))' and S = F' */
            /* workspace: Iwork (2*nrow) */
            A1 = ptranspose(A, 0, Perm, NULL, 0, status) ;
            F = A1 ;
            /* workspace: Iwork (nrow) */
            A2 = ptranspose(A1, 0, NULL, NULL, 0, status) ;
            S = A2 ;
        }
        else
        {
            /* unsymmetric case:     F = A (p,f)'         and S = F' */
            /* workspace: Iwork (nrow if no fset, MAX(nrow,ncol) if fset) */
            A1 = ptranspose(A, 0, Perm, fset, fsize, status) ;
            F = A1 ;
            if (do_rowcolcounts)
            {
                /* S not needed for unsymmetric case if do_rowcolcounts FALSE */
                /* workspace: Iwork (nrow) */
                A2 = ptranspose(A1, 0, NULL, NULL, 0, status) ;
            }
            S = A2 ;
        }
    }

    /* If any cholmod_*transpose fails, one or more matrices will be NULL */
    *A1_handle = A1 ;
    *A2_handle = A2 ;
    *S_handle = S ;
    *F_handle = F ;
    return (TRUE) ;
}

/* ========================================================================== */
/* === cholmod_analyze_ordering ============================================= */
/* ========================================================================== */

/* Given a matrix A and its fill-reducing permutation, compute the elimination
 * tree, its (non-weighted) postordering, and the number of nonzeros in each
 * column of L.  Also computes the flop count, the total nonzeros in L, and
 * the nonzeros in A (Common->fl, Common->lnz, and Common->anz).
 *
 * The column counts of L, flop count, and other statistics from
 * cholmod_rowcolcounts are not computed if ColCount is NULL.
 *
 * workspace: Iwork (2*nrow if symmetric, 2*nrow+ncol if unsymmetric),
 *	Flag (nrow), Head (nrow+1)
 */

int analyze_ordering
        (
                /* ---- input ---- */
                CSC *A,	/* matrix to analyze */
                int ordering,	/* ordering method used */
                int *Perm,		/* size n, fill-reducing permutation to analyze */
                int *fset,		/* subset of 0:(A->ncol)-1 */
                size_t fsize,	/* size of fset */
                /* ---- output --- */
                int *Parent,	/* size n, elimination tree */
                int *Post,		/* size n, postordering of elimination tree */
                int *ColCount,	/* size n, nnz in each column of L */

                /* ---- workspace  */
                int *First,		/* size n workspace for cholmod_postorder */
                int *Level,		/* size n workspace for cholmod_postorder */
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
#if 0
    for (int i = 0; i < n; ++i) {
        printf("%d :", Parent[i]);
    }
    printf("\n");
#endif

    /* postorder the etree (required by cholmod_rowcolcounts) */
    /* workspace: Iwork (2*nrow) */
    ok = ok && postOrderC(Parent, n, NULL, Post, status) == n ;

    /* cholmod_postorder doesn't set Common->status if it returns < n */
    //status = (!ok && status ) ?
                    // CHOLMOD_INVALID : Common->status ;

    /* analyze LL'=S or SS' or S(:,f)*S(:,f)' */
    /* workspace:
     *	if symmetric:   Flag (nrow), Iwork (2*nrow)
     *	if unsymmetric: Flag (nrow), Iwork (2*nrow+ncol), Head (nrow+1)
     */
    if (do_rowcolcounts)
    {
        ok = ok && rowcolcounts(A->stype ? F:S, fset, fsize, Parent,
                                          Post, NULL, ColCount, First, Level,
                                fl,
                                aatfl,
                                lnz, status) ;
    }

    /* free temporary matrices and return result */
   // CHOLMOD(free_sparse) (&A1, Common) ;
   // CHOLMOD(free_sparse) (&A2, Common) ;
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

BCSC* analyze_p2
        (
                /* ---- input ---- */
                int for_whom,       /* FOR_SPQR     (0): for SPQR but not GPU-accelerated
                           FOR_CHOLESKY (1): for Cholesky (GPU or not)
                           FOR_SPQRGPU  (2): for SPQR with GPU acceleration */
                CSC *A,	/* matrix to order and analyze */
                int *UserPerm,	/* user-provided permutation, size A->nrow */
                int *fset,		/* subset of 0:(A->ncol)-1 */
                int *nrelax,
                double *zrelax,
                size_t fsize,	/* size of fset */
                int* &prunePtr,
                int* &pruneSet,
                /* --------------- */
                int status
        )
{
    double lnz_best ;
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

/*    RETURN_IF_NULL_COMMON (NULL) ;
    RETURN_IF_NULL (A, NULL) ;
    RETURN_IF_XTYPE_INVALID (A, CHOLMOD_PATTERN, CHOLMOD_ZOMPLEX, NULL) ;
    Common->status = CHOLMOD_OK ;
    status = CHOLMOD_OK ;
    Common->selected = EMPTY ;
    Common->called_nd = FALSE ;*/

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    n = A->nrow ;
    ncol = A->ncol ;
    uncol = (A->stype == 0) ? (A->ncol) : 0 ;

    /* ---------------------------------------------------------------------- */
    /* set the default strategy */
    /* ---------------------------------------------------------------------- */

    lnz_best = (double) EMPTY ;
    skip_best = FALSE ;
  //  nmethods = MIN (Common->nmethods, CHOLMOD_MAXMETHODS) ;
   // nmethods = MAX (0, nmethods) ;



/*    default_strategy = (nmethods == 0) ;
    if (default_strategy)
    {
        *//* default strategy: try UserPerm, if given.  Try AMD for A, or AMD
         * to order A*A'.  Try METIS for the symmetric case only if AMD reports
             * a high degree of fill-in and flop count.  METIS is not tried if the
             * Partition Module isn't installed.   If Common->default_nesdis is
             * TRUE, then NESDIS is used as the 3rd ordering instead. *//*
        Common->method [0].ordering = CHOLMOD_GIVEN ;*//* skip if UserPerm NULL *//*
        Common->method [1].ordering = CHOLMOD_AMD ;
        Common->method [2].ordering =
                (Common->default_nesdis ? CHOLMOD_NESDIS : CHOLMOD_METIS) ;
        amd_backup = FALSE ;
#ifndef NPARTITION
        nmethods = 3 ;
#else
        nmethods = 2 ;
#endif
    }
    else
    {
        *//* If only METIS and NESDIS are selected, or if 2 or more methods are
         * being tried, then enable AMD backup *//*
        amd_backup = (nmethods > 1) || (nmethods == 1 &&
                                        (Common->method [0].ordering == CHOLMOD_METIS ||
                                         Common->method [0].ordering == CHOLMOD_NESDIS)) ;
    }*/

#ifdef NSUPERNODAL
    /* CHOLMOD Supernodal module not installed, just do simplicial analysis */
    Common->supernodal = CHOLMOD_SIMPLICIAL ;
#endif

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    /* Note: enough space needs to be allocated here so that routines called by
     * cholmod_analyze do not reallocate the space.
     */

    /* s = 6*n + uncol */
   /* s = CHOLMOD(mult_size_t) (n, 6, &ok) ;
    s = CHOLMOD(add_size_t) (s, uncol, &ok) ;
    if (!ok)
    {
//        ERROR (CHOLMOD_TOO_LARGE, "problem too large") ;
        return (NULL) ;
    }

    CHOLMOD(allocate_work) (n, s, 0, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
        return (NULL) ;	    *//* out of memory *//*
    }*/
 //   ASSERT (CHOLMOD(dump_work) (TRUE, TRUE, 0, Common)) ;

    /* ensure that subsequent routines, called by cholmod_analyze, do not
     * reallocate any workspace.  This is set back to FALSE in the
     * FREE_WORKSPACE_AND_RETURN macro, which is the only way this function
     * returns to its caller. */
   // Common->no_workspace_reallocate = TRUE ;

    /* Use the last 4*n int's in Iwork for Parent, First, Level, and Post, since
     * other CHOLMOD routines will use the first 2n+uncol space.  The ordering
     * routines (cholmod_amd, cholmod_colamd, cholmod_ccolamd, cholmod_metis)
     * are an exception.  They can use all 6n + ncol space, since the contents
     * of Parent, First, Level, and Post are not needed across calls to those
     * routines. */
    //s = 6*n + uncol;
    s=4*n;
    Work4n = new int[s]() ;
   // Work4n += 2*((size_t) n) + uncol ;
    Parent = Work4n ;
    First  = Work4n + n ;
    Level  = Work4n + 2*((size_t) n) ;
    Post   = Work4n + 3*((size_t) n) ;

    /* note that this assignment means that cholmod_nested_dissection,
     * cholmod_ccolamd, and cholmod_camd can use only the first 4n+uncol
     * space in Common->Iwork */
    Cmember = Post ;
    CParent = Level ;

    /* ---------------------------------------------------------------------- */
    /* allocate more workspace, and an empty simplicial symbolic factor */
    /* ---------------------------------------------------------------------- */

    //L = CHOLMOD(allocate_factor) (n, Common) ;//TODO
    L = new BCSC;
    L->Perm     = new int[n]();
    L->ColCount = new int[n]();
    L->ordering = CHOLMOD_AMD;
    Lparent  = new int[n]();//CHOLMOD(malloc) (n, sizeof (int), Common) ;
  //  Perm     = new int[n]();//CHOLMOD(malloc) (n, sizeof (int), Common) ;
  //  ColCount = new int[n]();//CHOLMOD(malloc) (n, sizeof (int), Common) ;
    /*if (Common->status < CHOLMOD_OK)
    {
        *//* out of memory *//*
        FREE_WORKSPACE_AND_RETURN ;
    }*/
    Lperm = L->Perm ;
    Lcolcount = L->ColCount ;
    //Common->anz = EMPTY ;

    /* ---------------------------------------------------------------------- */
    /* try all the requested ordering options and backup to AMD if needed */
    /* ---------------------------------------------------------------------- */
    double info[20]={0};
    double Control[2];
    Control [0] = 10; //TODO check later //AMD_Dense
    Control [1] = TRUE; //AMD_AGGRESSIVE
    /*Common->method [i].prune_dense2 = -1 ;	*//* COLAMD dense row control *//*
    320  	Common->method [i].aggressive = TRUE ;	*//* aggressive absorption *//*
    321  	Common->method [i].order_for_lu = FALSE ;*//* order for Cholesky not LU */
    amd_order(ncol,A->p,A->i,Lperm,NULL,info);
#if 0
    for (int i = 0; i < n; ++i) {
        printf("%d, ", Lperm[i]);
    }
    printf("\n");
#endif
//////// Can be added later
    skip_best=TRUE;

    /* ---------------------------------------------------------------------- */
    /* do the analysis for AMD, if skipped */
    /* ---------------------------------------------------------------------- */

/*    Common->fl  = Common->method [Common->selected].fl  ;
    Common->lnz = Common->method [Common->selected].lnz ;
    ASSERT (Common->lnz >= 0) ;*/

    if (skip_best)
    {
        if (!analyze_ordering (A, L->ordering, Lperm, fset, fsize,
                               Lparent, Post, Lcolcount, First,
                               Level, status))
        {
            /* out of memory, or method failed */
            //FREE_WORKSPACE_AND_RETURN ;
            return FALSE;
        }
    }

    /* ---------------------------------------------------------------------- */
    /* postorder the etree, weighted by the column counts */
    /* ---------------------------------------------------------------------- */

    //if (Common->postorder)
    if(true)//TODO
    {
        /* combine the fill-reducing ordering with the weighted postorder */
        /* workspace: Iwork (2*nrow) */
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
    }

    /* ---------------------------------------------------------------------- */
    /* supernodal analysis, if requested or if selected automatically */
    /* ---------------------------------------------------------------------- */

#ifndef NSUPERNODAL
  /*  if (Common->supernodal > CHOLMOD_AUTO
        || (Common->supernodal == CHOLMOD_AUTO &&
            Common->lnz > 0 &&
            (Common->fl / Common->lnz) >= Common->supernodal_switch))*/
    int supernodal = TRUE;
    int *Sparent = new int[n]();
    if(supernodal)
    {
        CSC *S, *F, *A2, *A1 ;

        permute_matrices (A, L->ordering, Lperm, fset, fsize, TRUE,
                          &A1, &A2, &S, &F, status) ;

        /* workspace: Flag (nrow), Head (nrow), Iwork (5*nrow) */
        /*int nrelax[3] = {4,16,48};//TODO
        double zrelax[3] = {0.8,0.1,0.05};*/
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

        //CHOLMOD(free_sparse) (&A1, Common) ;
        //CHOLMOD(free_sparse) (&A2, Common) ;
        delete []col2Sup;
        allocateAC(A1,0,0,0,FALSE);
        allocateAC(A2,0,0,0,FALSE);
    }
#endif

    /* ---------------------------------------------------------------------- */
    /* free temporary matrices and workspace, and return result L */
    /* ---------------------------------------------------------------------- */

    //FREE_WORKSPACE_AND_RETURN ;
    delete []Sparent;
    delete []Work4n;
    return L;
}
//#endif
#endif //CHOLOPENMP_LSPARSITY_H
