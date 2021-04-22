//
// Created by kazem on 7/25/17.
//

#ifndef CHOLOPENMP_DFS_H
#define CHOLOPENMP_DFS_H

#include "def.h"

/* non-recursive version for actual use */

int dfsC		/* return the new value of k */
        (
                int p,		/* start the DFS at a root node p */
                int k,		/* start the node numbering at k */
                int *Post,	/* Post ordering, modified on output */
                int *Head ,	/* Head [p] = youngest child of p; EMPTY on output */
                int *Next ,	/* Next [j] = sibling of j; unmodified */
                int *Pstack 	/* workspace of size n, undefined on input or output */
        )
{
    int j, phead ;

    /* put the root node on the stack */
    Pstack [0] = p ;
    phead = 0 ;
    while (phead >= 0)
    {
        p = Pstack [phead] ;
        j = Head [p] ;
        if (j == EMPTY)
        {
            phead-- ;
            Post [k++] = p ;	/* order node p as the kth node */
        }
        else
        {
            Head [p] = Next [j] ;
            Pstack [++phead] = j ;
        }
    }
    return (k) ;	/* the next node will be numbered k */
}

//From CSPARSE
/* depth-first-search of the graph of a matrix, starting at node j */
int dfs (int j, int *Gp, int *Gi, int top, int *xi, int *pstack, const int *pinv)
{
    int i, p, p2, done, jnew, head = 0 ;
    //if (!CS_CSC (G) || !xi || !pstack) return (-1) ;    /* check inputs */
    //Gp = G->p ; Gi = G->i ;
    xi [0] = j ;                /* initialize the recursion stack */

    while (head >= 0)
    {
        j = xi [head] ;         /* get j from the top of the recursion stack */
        jnew = pinv ? (pinv [j]) : j ;
        if (!CS_MARKED (Gp, j))
        {
            CS_MARK (Gp, j) ;       /* mark node j as visited */
            pstack [head] = (jnew < 0) ? 0 : CS_UNFLIP (Gp [jnew]) ;
        }
        done = 1 ;                  /* node j done if no unvisited neighbors */
        p2 = (jnew < 0) ? 0 : CS_UNFLIP (Gp [jnew+1]) ;
        for (p = pstack [head] ; p < p2 ; p++)  /* examine all neighbors of j */
        {
            i = Gi [p] ;            /* consider neighbor node i */
            if (CS_MARKED (Gp, i)) continue ;   /* skip visited node i */
            pstack [head] = p ;     /* pause depth-first search of node j */
            xi [++head] = i ;       /* start dfs at node i */
            done = 0 ;              /* node j is not done */
            break ;                 /* break, to start dfs (i) */
        }
        if (done)               /* depth-first search at node j is done */
        {
            head-- ;            /* remove j from the recursion stack */
            xi [--top] = j ;    /* and place in the output stack */
        }
    }
    return (top) ;
}


/* depth-first search and postorder of a tree rooted at node j */
int tdfs(int j, int k, int *head, const int *next, int *post, int *stack)
{
    int i, p, top = 0 ;
    if (!head || !next || !post || !stack) return (-1) ;    /* check inputs */
    stack [0] = j ;                 /* place j on the stack */
    while (top >= 0)                /* while (stack is not empty) */
    {
        p = stack [top] ;           /* p = top of stack */
        i = head [p] ;              /* i = youngest child of p */
        if (i == -1)
        {
            top-- ;                 /* p has no unordered children left */
            post [k++] = p ;        /* node p is the kth postordered node */
        }
        else
        {
            head [p] = next [i] ;   /* remove i from children of p */
            stack [++top] = i ;     /* start dfs on child node i */
        }
    }
    return (k) ;
}


#endif //CHOLOPENMP_DFS_H
