//
// The modified version of CHOLMOD code
//

#ifndef CHOLOPENMP_POSTORDER_H
#define CHOLOPENMP_POSTORDER_H

#include "def.h"
#include "SparseUtils.h"

int postOrderC	/* return # of nodes postordered */
        (
                /* ---- input ---- */
                int *Parent,	/* size n. Parent [j] = p if p is the parent of j */
                size_t n,
                int *Weight,	/* size n, optional. Weight [j] is weight of node j */
                /* ---- output --- */
                int *Post,		/* size n. Post [k] = j is kth in postordered tree */
                /* --------------- */
               // cholmod_common *Common
                int status
        )
{
    int *Head, *Next, *Pstack, *Iwork ;
    int j, p, k, w, nextj ;
    size_t s ;
    int ok = TRUE ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    if(Parent==NULL || n<0)
       return EMPTY;

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    /* s = 2*n */
    s = mult_size_t (n, 3, &ok) ;
    if (!ok)
    {
        return (EMPTY) ;
    }


    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */
    Iwork = new int[s+1]();
    Head = Iwork;
    Next  = Iwork + n + 1 ;		/* size n (i/i/l) */
    Pstack = Iwork + 2*n + 1 ;	/* size n (i/i/l) */

    for (int i = 0; i < n + 1; ++i) {
        Head[i] = EMPTY;
    }
    /* ---------------------------------------------------------------------- */
    /* construct a link list of children for each node */
    /* ---------------------------------------------------------------------- */

    if (Weight == NULL)
    {

        /* in reverse order so children are in ascending order in each list */
        for (j = n-1 ; j >= 0 ; j--)
        {
            p = Parent [j] ;
            if (p >= 0 && p < ((int) n))
            {
                /* add j to the list of children for node p */
                Next [j] = Head [p] ;
                Head [p] = j ;
            }
        }

        /* Head [p] = j if j is the youngest (least-numbered) child of p */
        /* Next [j1] = j2 if j2 is the next-oldest sibling of j1 */

    }
    else
    {

        /* First, construct a set of link lists according to Weight.
         *
         * Whead [w] = j if node j is the first node in bucket w.
         * Next [j1] = j2 if node j2 follows j1 in a link list.
         */

        int *Whead = Pstack ;	    /* use Pstack as workspace for Whead [ */

        for (w = 0 ; w < ((int) n) ; w++)
        {
            Whead [w] = EMPTY ;
        }
        /* do in forward order, so nodes that ties are ordered by node index */
        for (j = 0 ; j < ((int) n) ; j++)
        {
            p = Parent [j] ;
            if (p >= 0 && p < ((int) n))
            {
                w = Weight [j] ;
                w = MAX (0, w) ;
                w = MIN (w, ((int) n) - 1) ;
                /* place node j at the head of link list for weight w */
                Next [j] = Whead [w] ;
                Whead [w] = j ;
            }
        }

        /* traverse weight buckets, placing each node in its parent's list */
        for (w = n-1 ; w >= 0 ; w--)
        {
            for (j = Whead [w] ; j != EMPTY ; j = nextj)
            {
                nextj = Next [j] ;
                /* put node j in the link list of its parent */
                p = Parent [j] ;
                ASSERT (p >= 0 && p < ((int) n)) ;
                Next [j] = Head [p] ;
                Head [p] = j ;
            }
        }

        /* Whead no longer needed ] */
        /* Head [p] = j if j is the lightest child of p */
        /* Next [j1] = j2 if j2 is the next-heaviest sibling of j1 */
    }

    /* ---------------------------------------------------------------------- */
    /* start a DFS at each root node of the etree */
    /* ---------------------------------------------------------------------- */

    k = 0 ;
    for (j = 0 ; j < ((int) n) ; j++)
    {
        if (Parent [j] == EMPTY)
        {
            /* j is the root of a tree; start a DFS here */
            k = dfsC (j, k, Post, Head, Next, Pstack) ;
        }
    }

    /* this would normally be EMPTY already, unless Parent is invalid */
    for (j = 0 ; j < ((int) n) ; j++)
    {
        Head [j] = EMPTY ;
    }

    delete []Iwork;
    return (k) ;
}



#endif //CHOLOPENMP_POSTORDER_H
