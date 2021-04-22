//
// The modified version of CHOLMOD code
//

#ifndef CHOLOPENMP_POSTORDER_H
#define CHOLOPENMP_POSTORDER_H

#include "def.h"
#include "SparseUtils.h"

int postOrderC /* return # of nodes postordered */
  (int *eTree, size_t n, int *Weight, int *Post)
{
    int *Head, *Next, *Pstack, *tempSpace ;
    int j, p, k, w, nextj ;
    size_t s ;
    int ok = TRUE ;
    /* s = 2*n */
    s = mult_size_t (n, 3, &ok) ;
    if (!ok)
    {
        return (EMPTY) ;
    }
    tempSpace = new int[s+1]();
    Head = tempSpace;
    Next  = tempSpace + n + 1 ;		/* size n (i/i/l) */
    Pstack = tempSpace + 2*n + 1 ;	/* size n (i/i/l) */

    for (unsigned i = 0; i < n + 1; ++i) {
        Head[i] = EMPTY;
    }
    if (Weight == NULL)
    {
        for (j = n-1 ; j >= 0 ; j--)
        {
            p = eTree [j] ;
            if (p >= 0 && p < ((int) n))
            {
                /* add j to the list of children for node p */
                Next [j] = Head [p] ;
                Head [p] = j ;
            }
        }
    }
    else
    {
        int *Whead = Pstack ;	    /* use Pstack as workspace for Whead [ */

        for (w = 0 ; w < ((int) n) ; w++)
        {
            Whead [w] = EMPTY ;
        }
        for (j = 0 ; j < ((int) n) ; j++)
        {
            p = eTree [j] ;
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
        for (w = n-1 ; w >= 0 ; w--)
        {
            for (j = Whead [w] ; j != EMPTY ; j = nextj)
            {
                nextj = Next [j] ;
                /* put node j in the link list of its parent */
                p = eTree [j] ;
                ASSERT (p >= 0 && p < ((int) n)) ;
                Next [j] = Head [p] ;
                Head [p] = j ;
            }
        }
    }
    k = 0 ;
    for (j = 0 ; j < ((int) n) ; j++)
    {
        if (eTree [j] == EMPTY)
        {
            k = dfsC (j, k, Post, Head, Next, Pstack) ;
        }
    }
    for (j = 0 ; j < ((int) n) ; j++)
    {
        Head [j] = EMPTY ;
    }
    delete []tempSpace;
    return (k) ;
}


int *postOrder (const int *parent, int n)
{
/* post order a forest
 * Obtained from CSparse library
 * */
    int j, k = 0, *post, *w, *head, *next, *stack ;
    if (!parent) return (NULL) ;                        /* check inputs */
    //post = cs_malloc (n, sizeof (csi)) ;                /* allocate result */
    //w = cs_malloc (3*n, sizeof (csi)) ;                 /* get workspace */
    post = new int[n];
    w = new int[3*n];
    //if (!w || !post) return (cs_idone (post, NULL, w, 0)) ;
    if (!w || !post)
        return NULL;
    head = w ; next = w + n ; stack = w + 2*n ;
    for (j = 0 ; j < n ; j++) head [j] = -1 ;           /* empty linked lists */
    for (j = n-1 ; j >= 0 ; j--)            /* traverse nodes in reverse order*/
    {
        if (parent [j] == -1) continue ;    /* j is a root */
        next [j] = head [parent [j]] ;      /* add j to list of its parent */
        head [parent [j]] = j ;
    }
    for (j = 0 ; j < n ; j++)
    {
        if (parent [j] != -1) continue ;    /* skip j if it is not a root */
        k = tdfs(j, k, head, next, post, stack) ;
    }
    // return (cs_idone (post, NULL, w, 1)) ;  /* success; free w, return post */
    delete []w;
    return post;
}


#endif //CHOLOPENMP_POSTORDER_H
