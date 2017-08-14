//
// Created by kazem on 7/25/17.
//

#ifndef CHOLOPENMP_ETREE_H
#define CHOLOPENMP_ETREE_H


#include "def.h"
#include "SparseUtils.h"

/* ========================================================================== */
/* === update_etree ========================================================= */
/* ========================================================================== */

void update_etree
        (
                /* inputs, not modified */
                int k,		/* process the edge (k,i) in the input graph */
                int i,
                /* inputs, modified on output */
                int Parent [ ],	/* Parent [t] = p if p is the parent of t */
                int Ancestor [ ]	/* Ancestor [t] is the ancestor of node t in the
			   partially-constructed etree */
        )
{
    int a ;
    for ( ; ; )		/* traverse the path from k to the root of the tree */
    {
        a = Ancestor [k] ;
        if (a == i)
        {
            /* final ancestor reached; no change to tree */
            return ;
        }
        /* perform path compression */
        Ancestor [k] = i ;
        if (a == EMPTY)
        {
            /* final ancestor undefined; this is a new edge in the tree */
            Parent [k] = i ;
            return ;
        }
        /* traverse up to the ancestor of k */
        k = a ;
    }
}


/* ========================================================================== */
/* === cholmod_etree ======================================================== */
/* ========================================================================== */

/* Find the elimination tree of A or A'*A */

int etreeC
        (
                /* ---- input ---- */
                CSC *A,
                /* ---- output --- */
                int *Parent,	/* size ncol.  Parent [j] = p if p is the parent of j */
                /* --------------- */
                //cholmod_common *Common
                int &status
        )
{
    int *Ap, *Ai, *Anz, *Ancestor, *Prev, *Iwork ;
    int i, j, jprev, p, pend, nrow, ncol, packed, stype ;
    size_t s ;
    int ok = TRUE ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    /*RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (A, FALSE) ;
    RETURN_IF_NULL (Parent, FALSE) ;
    RETURN_IF_XTYPE_INVALID (A, CHOLMOD_PATTERN, CHOLMOD_ZOMPLEX, FALSE) ;*/
    //Common->status = CHOLMOD_OK ;
    status = TRUE;

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    stype = A->stype ;

    /* s = A->nrow + (stype ? 0 : A->ncol) */
    s = add_size_t (A->nrow, (stype ? 0 : A->ncol), &ok) ;
    if (!ok)
    {
        //ERROR (CHOLMOD_TOO_LARGE, "problem too large") ;
        return (FALSE) ;
    }

    //CHOLMOD(allocate_work) (0, s, 0, Common) ;
    /*if (Common->status < CHOLMOD_OK)
    {
        return (FALSE) ;	*//* out of memory *//*
    }*/

    //ASSERT (CHOLMOD(dump_sparse) (A, "etree", Common) >= 0) ;
    //Iwork = Common->Iwork ;
    Iwork = new int[s]();

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    ncol = A->ncol ;	/* the number of columns of A */
    nrow = A->nrow ;	/* the number of rows of A */
    Ap = A->p ;		/* size ncol+1, column pointers for A */
    Ai = A->i ;		/* the row indices of A */
    Anz = A->nz ;	/* number of nonzeros in each column of A */
    packed = A->packed ;
    Ancestor = Iwork ;	/* size ncol (i/i/l) */

    for (j = 0 ; j < ncol ; j++)
    {
        Parent [j] = EMPTY ;
        Ancestor [j] = EMPTY ;
    }

    /* ---------------------------------------------------------------------- */
    /* compute the etree */
    /* ---------------------------------------------------------------------- */

    if (stype > 0)
    {

        /* ------------------------------------------------------------------ */
        /* symmetric (upper) case: compute etree (A) */
        /* ------------------------------------------------------------------ */

        for (j = 0 ; j < ncol ; j++)
        {
            /* for each row i in column j of triu(A), excluding the diagonal */
            p = Ap [j] ;
            pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
            for ( ; p < pend ; p++)
            {
                i = Ai [p] ;
                if (i < j)
                {
                    update_etree (i, j, Parent, Ancestor) ;
                }
            }
        }

    }
    else if (stype == 0)
    {

        /* ------------------------------------------------------------------ */
        /* unsymmetric case: compute etree (A'*A) */
        /* ------------------------------------------------------------------ */

        Prev = Iwork + ncol ;	/* size nrow (i/i/l) */
        for (i = 0 ; i < nrow ; i++)
        {
            Prev [i] = EMPTY ;
        }
        for (j = 0 ; j < ncol ; j++)
        {
            /* for each row i in column j of A */
            p = Ap [j] ;
            pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
            for ( ; p < pend ; p++)
            {
                /* a graph is constructed dynamically with one path per row
                 * of A.  If the ith row of A contains column indices
                 * (j1,j2,j3,j4) then the new graph has edges (j1,j2), (j2,j3),
                 * and (j3,j4).  When at node i of this path-graph, all edges
                 * (jprev,j) are considered, where jprev<j */
                i = Ai [p] ;
                jprev = Prev [i] ;
                if (jprev != EMPTY)
                {
                    update_etree (jprev, j, Parent, Ancestor) ;
                }
                Prev [i] = j ;
            }
        }

    }
    else
    {

        /* ------------------------------------------------------------------ */
        /* symmetric case with lower triangular part not supported */
        /* ------------------------------------------------------------------ */

       // ERROR (CHOLMOD_INVALID, "symmetric lower not supported") ;
        return (FALSE) ;
    }

//    ASSERT (CHOLMOD(dump_parent) (Parent, ncol, "Parent", Common)) ;
    delete []Iwork;
    return (TRUE) ;
}

/* ========================================================================== */
/* === cholmod_etree ======================================================== */
/* ========================================================================== */

/* Find the elimination tree of A or A'*A */

int etree1
        (
                /* ---- input ---- */
               // cholmod_sparse *A,
                int ncol,
                int nrow,
                int *Ap,
                int *Ai,
                int *Anz,
                int stype,
                int packed,
                /* ---- output --- */
                int *Parent	/* size ncol.  Parent [j] = p if p is the parent of j */
                /* --------------- */
                //cholmod_common *Common
        )
{
    int  *Ancestor, *Prev, *Iwork ;
    int i, j, jprev, p, pend;
    size_t s ;
    int ok = TRUE ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

/*    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (A, FALSE) ;
    RETURN_IF_NULL (Parent, FALSE) ;
    RETURN_IF_XTYPE_INVALID (A, CHOLMOD_PATTERN, CHOLMOD_ZOMPLEX, FALSE) ;
    Common->status = CHOLMOD_OK ;*/

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    //stype = A->stype ;

    /* s = A->nrow + (stype ? 0 : A->ncol) */
/*    s = CHOLMOD(add_size_t) (A->nrow, (stype ? 0 : A->ncol), &ok) ;
    if (!ok)
    {
        ERROR (CHOLMOD_TOO_LARGE, "problem too large") ;
        return (FALSE) ;
    }

    CHOLMOD(allocate_work) (0, s, 0, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
        return (FALSE) ;	*//* out of memory *//*
    }

    ASSERT (CHOLMOD(dump_sparse) (A, "etree", Common) >= 0) ;
    Iwork = Common->Iwork ;*/

    s = nrow + (stype ? 0 : ncol);
    Iwork = new int[s]();
    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    //ncol = A->ncol ;	/* the number of columns of A */
    //nrow = A->nrow ;	/* the number of rows of A */
    //Ap = A->p ;		/* size ncol+1, column pointers for A */
    //Ai = A->i ;		/* the row indices of A */
    //Anz = A->nz ;	/* number of nonzeros in each column of A */
    //packed = A->packed ;
    Ancestor = Iwork ;	/* size ncol (i/i/l) */

    for (j = 0 ; j < ncol ; j++)
    {
        Parent [j] = EMPTY ;
        Ancestor [j] = EMPTY ;
    }

    /* ---------------------------------------------------------------------- */
    /* compute the etree */
    /* ---------------------------------------------------------------------- */

    if (stype > 0)
    {

        /* ------------------------------------------------------------------ */
        /* symmetric (upper) case: compute etree (A) */
        /* ------------------------------------------------------------------ */

        for (j = 0 ; j < ncol ; j++)
        {
            /* for each row i in column j of triu(A), excluding the diagonal */
            p = Ap [j] ;
            pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
            for ( ; p < pend ; p++)
            {
                i = Ai [p] ;
                if (i < j)
                {
                    update_etree (i, j, Parent, Ancestor) ;
                }
            }
        }

    }
    else if (stype == 0)
    {

        /* ------------------------------------------------------------------ */
        /* unsymmetric case: compute etree (A'*A) */
        /* ------------------------------------------------------------------ */

        Prev = Iwork + ncol ;	/* size nrow (i/i/l) */
        for (i = 0 ; i < nrow ; i++)
        {
            Prev [i] = EMPTY ;
        }
        for (j = 0 ; j < ncol ; j++)
        {
            /* for each row i in column j of A */
            p = Ap [j] ;
            pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
            for ( ; p < pend ; p++)
            {
                /* a graph is constructed dynamically with one path per row
                 * of A.  If the ith row of A contains column indices
                 * (j1,j2,j3,j4) then the new graph has edges (j1,j2), (j2,j3),
                 * and (j3,j4).  When at node i of this path-graph, all edges
                 * (jprev,j) are considered, where jprev<j */
                i = Ai [p] ;
                jprev = Prev [i] ;
                if (jprev != EMPTY)
                {
                    update_etree (jprev, j, Parent, Ancestor) ;
                }
                Prev [i] = j ;
            }
        }

    }
    else
    {

        /* ------------------------------------------------------------------ */
        /* symmetric case with lower triangular part not supported */
        /* ------------------------------------------------------------------ */

        printf( "symmetric lower not supported") ;
        return (FALSE) ;
    }

    //ASSERT (CHOLMOD(dump_parent) (Parent, ncol, "Parent", Common)) ;
    return (TRUE) ;
}






#endif //CHOLOPENMP_ETREE_H
