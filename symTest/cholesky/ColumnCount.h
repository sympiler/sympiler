//
// Created by kazem on 7/25/17.
//

#ifndef CHOLOPENMP_COLUMNCOUNT_H
#define CHOLOPENMP_COLUMNCOUNT_H

#include "def.h"
#include "Transpose.h"

int initialize_node  /* initial work for kth node in postordered etree */
        (
                int k,		/* at the kth step of the algorithm (and kth node) */
                int Post [ ],	/* Post [k] = i, the kth node in postordered etree */
                int Parent [ ],	/* Parent [i] is the parent of i in the etree */
                int ColCount [ ],	/* ColCount [c] is the current weight of node c */
                int PrevNbr [ ]	/* PrevNbr [u] = k if u was last considered at step k */
        )
{
    int p, parent ;
    /* determine p, the kth node in the postordered etree */
    p = Post [k] ;
    /* adjust the weight if p is not a root of the etree */
    parent = Parent [p] ;
    if (parent != EMPTY)
    {
        ColCount [parent]-- ;
    }
    /* flag node p to exclude self edges (p,p) */
    PrevNbr [p] = k ;
    return (p) ;
}



/* edge (p,u) is being processed.  p < u is a descendant of its ancestor u in
 * the etree.  node p is the kth node in the postordered etree.  */

 void process_edge
        (
                int p,		/* process edge (p,u) of the matrix */
                int u,
                int k,		/* we are at the kth node in the postordered etree */
                int First [ ],	/* First [i] = k if the postordering of first
			 * descendent of node i is k */
                int PrevNbr [ ],	/* u was last considered at step k = PrevNbr [u] */
                int ColCount [ ],	/* ColCount [c] is the current weight of node c */
                int PrevLeaf [ ],	/* s = PrevLeaf [u] means that s was the last leaf
			 * seen in the subtree rooted at u.  */
                int RowCount [ ],	/* RowCount [i] is # of nonzeros in row i of L,
			 * including the diagonal.  Not computed if NULL. */
                int SetParent [ ],
                int Level [ ]	 /* Level [i] = length of path from node i to root */
        )
{
    int prevleaf, q, s, sparent ;
    if (First [p] > PrevNbr [u])
    {
        /* p is a leaf of the subtree of u */
        ColCount [p]++ ;
        prevleaf = PrevLeaf [u] ;
        if (prevleaf == EMPTY)
        {
            /* p is the first leaf of subtree of u; RowCount will be incremented
             * by the length of the path in the etree from p up to u. */
            q = u ;
        }
        else
        {
            /* q = FIND (prevleaf): find the root q of the
             * SetParent tree containing prevleaf */
            for (q = prevleaf ; q != SetParent [q] ; q = SetParent [q])
            {
                ;
            }
            /* the root q has been found; re-traverse the path and
             * perform path compression */
            s = prevleaf ;
            for (s = prevleaf ; s != q ; s = sparent)
            {
                sparent = SetParent [s] ;
                SetParent [s] = q ;
            }
            /* adjust the RowCount and ColCount; RowCount will be incremented by
             * the length of the path from p to the SetParent root q, and
             * decrement the ColCount of q by one. */
            ColCount [q]-- ;
        }
        if (RowCount != NULL)
        {
            /* if RowCount is being computed, increment it by the length of
             * the path from p to q */
            RowCount [u] += (Level [p] - Level [q]) ;
        }
        /* p is a leaf of the subtree of u, so mark PrevLeaf [u] to be p */
        PrevLeaf [u] = p ;
    }
    /* flag u has having been processed at step k */
    PrevNbr [u] = k ;
}


void finalize_node    /* compute UNION (p, Parent [p]) */
        (
                int p,
                int Parent [ ],	/* Parent [p] is the parent of p in the etree */
                int SetParent [ ]	/* see process_edge, above */
        )
{
    /* all nodes in the SetParent tree rooted at p now have as their final
     * root the node Parent [p].  This computes UNION (p, Parent [p]) */
    if (Parent [p] != EMPTY)
    {
        SetParent [p] = Parent [p] ;
    }
}

int rowcolcounts
        (
                /* ---- input ---- */
                CSC *A,	/* matrix to analyze */
                int *fset,		/* subset of 0:(A->ncol)-1 */
                size_t fsize,	/* size of fset */
                int *Parent,	/* size nrow.  Parent [i] = p if p is the parent of i */
                int *Post,		/* size nrow.  Post [k] = i if i is the kth node in
			 * the postordered etree. */
                /* ---- output --- */
                int *RowCount,	/* size nrow. RowCount [i] = # entries in the ith row of
			 * L, including the diagonal. */
                int *ColCount,	/* size nrow. ColCount [i] = # entries in the ith
			 * column of L, including the diagonal. */
                int *First,		/* size nrow.  First [i] = k is the least postordering
			 * of any descendant of i. */
                int *Level,		/* size nrow.  Level [i] is the length of the path from
			 * i to the root, with Level [root] = 0. */
                /* --------------- */
             int &fl,
                int &aatfl,
                int &lnz,
                //   cholmod_common *Common
                int status
        )
{
    double ff ;
    int *Ap, *Ai, *Anz, *PrevNbr, *SetParent, *Head, *PrevLeaf, *Anext, *Ipost,
            *Iwork ;
    int i, j, r, k, len, s, p, pend, inew, stype, nf, anz, parent,
            nrow, ncol, packed, use_fset, jj ;
    size_t w ;
    int ok = TRUE ;

    stype = A->stype ;
    if (stype > 0)
    {

        return (FALSE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    nrow = A->nrow ;	/* the number of rows of A */
    ncol = A->ncol ;	/* the number of columns of A */

    /* w = 2*nrow + (stype ? 0 : ncol) */
    w = mult_size_t (nrow, 2, &ok) ;
    w = add_size_t (w, (stype ? 0 : ncol), &ok) ;
    if (!ok)
    {
        return (FALSE) ;
    }

    Ap = A->p ;	/* size ncol+1, column pointers for A */
    Ai = A->i ;	/* the row indices of A, of size nz=Ap[ncol+1] */
    Anz = A->nz ;
    packed = A->packed ;

    /* ---------------------------------------------------------------------- */
    /* get workspace */
    /* ---------------------------------------------------------------------- */
    Iwork = new int[w]();
    //Iwork = Common->Iwork ;
    SetParent = Iwork ;		    /* size nrow (i/i/l) */
    PrevNbr   = Iwork + nrow ;	    /* size nrow (i/i/l) */
    Anext     = Iwork + 2*((size_t) nrow) ;    /* size ncol (i/i/l) (unsym only) */
    //PrevLeaf  = Common->Flag ;	    /* size nrow */
    //Head      = Common->Head ;	    /* size nrow+1 (unsym only)*/
    PrevLeaf = new int[nrow];
    Head = new int[nrow+1];

    /* ---------------------------------------------------------------------- */
    /* find the first descendant and level of each node in the tree */
    /* ---------------------------------------------------------------------- */

    /* First [i] = k if the postordering of first descendent of node i is k */
    /* Level [i] = length of path from node i to the root (Level [root] = 0) */

    for (i = 0 ; i < nrow ; i++)
    {
        First [i] = EMPTY ;
    }

    /* postorder traversal of the etree */
    for (k = 0 ; k < nrow ; k++)
    {
        /* node i of the etree is the kth node in the postordered etree */
        i = Post [k] ;

        /* i is a leaf if First [i] is still EMPTY */
        /* ColCount [i] starts at 1 if i is a leaf, zero otherwise */
        ColCount [i] = (First [i] == EMPTY) ? 1 : 0 ;

        /* traverse the path from node i to the root, stopping if we find a
         * node r whose First [r] is already defined. */
        len = 0 ;
        for (r = i ; (r != EMPTY) && (First [r] == EMPTY) ; r = Parent [r])
        {
            First [r] = k ;
            len++ ;
        }
        if (r == EMPTY)
        {
            /* we hit a root node, the level of which is zero */
            len-- ;
        }
        else
        {
            /* we stopped at node r, where Level [r] is already defined */
            len += Level [r] ;
        }
        /* re-traverse the path from node i to r; set the level of each node */
        for (s = i ; s != r ; s = Parent [s])
        {
            Level [s] = len-- ;
        }
    }


    fl = 0.0 ;
    if (stype == 0)
    {
        /* [ use PrevNbr [0..nrow-1] as workspace for Ipost */
        Ipost = PrevNbr ;
        /* Ipost [i] = k if i is the kth node in the postordered etree. */
        for (k = 0 ; k < nrow ; k++)
        {
            Ipost [Post [k]] = k ;
        }
        use_fset = (fset != NULL) ;
        if (use_fset)
        {
            nf = fsize ;
            /* clear Anext to check fset */
            for (j = 0 ; j < ncol ; j++)
            {
                Anext [j] = -2 ;
            }
            /* find the first postordered row in each column of A (post,f)
             * and place the column in the corresponding link list */
            for (jj = 0 ; jj < nf ; jj++)
            {
                j = fset [jj] ;
                if (j < 0 || j > ncol || Anext [j] != -2)
                {
                    /* out-of-range or duplicate entry in fset */
                 //   ERROR (CHOLMOD_INVALID, "fset invalid") ;
                    return (FALSE) ;
                }
                /* flag column j as having been seen */
                Anext [j] = EMPTY ;
            }
            /* fset is now valid */
           // ASSERT (CHOLMOD(dump_perm) (fset, nf, ncol, "fset", Common)) ;
        }
        else
        {
            nf = ncol ;
        }
        for (jj = 0 ; jj < nf ; jj++)
        {
            j = (use_fset) ? (fset [jj]) : jj ;
            /* column j is in the fset; find the smallest row (if any) */
            p = Ap [j] ;
            pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
            ff = (double) MAX (0, pend - p) ;
            fl += ff*ff + ff ;
            if (pend > p)
            {
                k = Ipost [Ai [p]] ;
                for ( ; p < pend ; p++)
                {
                    inew = Ipost [Ai [p]] ;
                    k = MIN (k, inew) ;
                }
                /* place column j in link list k */
                ASSERT (k >= 0 && k < nrow) ;
                Anext [j] = Head [k] ;
                Head [k] = j ;
            }
        }

    }

    /* ---------------------------------------------------------------------- */
    /* compute the row counts and node weights */
    /* ---------------------------------------------------------------------- */

    if (RowCount != NULL)
    {
        for (i = 0 ; i < nrow ; i++)
        {
            RowCount [i] = 1 ;
        }
    }
    for (i = 0 ; i < nrow ; i++)
    {
        PrevLeaf [i] = EMPTY ;
        PrevNbr [i] = EMPTY ;
        SetParent [i] = i ;	/* every node is in its own set, by itself */
    }



    /* ------------------------------------------------------------------ */
    /* symmetric case: LL' = A */
    /* ------------------------------------------------------------------ */

    /* also determine the number of entries in triu(A) */
    anz = nrow ;
    for (k = 0 ; k < nrow ; k++)
    {
        /* j is the kth node in the postordered etree */
        j = initialize_node (k, Post, Parent, ColCount, PrevNbr) ;

        /* for all nonzeros A(i,j) below the diagonal, in column j of A */
        p = Ap [j] ;
        pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
        for ( ; p < pend ; p++)
        {
            i = Ai [p] ;
            if (i > j)
            {
                /* j is a descendant of i in etree(A) */
                anz++ ;
                process_edge (j, i, k, First, PrevNbr, ColCount,
                              PrevLeaf, RowCount, SetParent, Level) ;
            }
        }
        /* update SetParent: UNION (j, Parent [j]) */
        finalize_node (j, Parent, SetParent) ;
    }
   // Common->anz = anz ;


    for (j = 0 ; j < nrow ; j++)
    {
        parent = Parent [j] ;
        if (parent != EMPTY)
        {
            /* add the ColCount of j to its parent */
            ColCount [parent] += ColCount [j] ;
        }
    }

    aatfl = fl ;
    lnz = 0. ;
    fl = 0 ;
    for (j = 0 ; j < nrow ; j++)
    {
        ff = (double) (ColCount [j]) ;
        lnz += ff ;
        fl += ff*ff ;
    }
    return (TRUE) ;
}

#endif //CHOLOPENMP_COLUMNCOUNT_H
