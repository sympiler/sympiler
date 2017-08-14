//
// Created by kazem on 7/28/17.
//

#ifndef CHOLOPENMP_COLCOUNTS_H
#define CHOLOPENMP_COLCOUNTS_H

#include "Transpose.h"

/* consider A(i,j), node j in ith row subtree and return lca(jprev,j) */
int leaf (int i, int j, const int *first, int *maxfirst, int *prevleaf,
          int *ancestor, int *jleaf)
{
    int q, s, sparent, jprev ;
    if (!first || !maxfirst || !prevleaf || !ancestor || !jleaf) return (-1) ;
    *jleaf = 0 ;
    if (i <= j || first [j] <= maxfirst [i]) return (-1) ;  /* j not a leaf */
    maxfirst [i] = first [j] ;      /* update max first[j] seen so far */
    jprev = prevleaf [i] ;          /* jprev = previous leaf of ith subtree */
    prevleaf [i] = j ;
    *jleaf = (jprev == -1) ? 1: 2 ; /* j is first or subsequent leaf */
    if (*jleaf == 1) return (i) ;   /* if 1st leaf, q = root of ith subtree */
    for (q = jprev ; q != ancestor [q] ; q = ancestor [q]) ;
    for (s = jprev ; s != q ; s = sparent)
    {
        sparent = ancestor [s] ;    /* path compression */
        ancestor [s] = q ;
    }
    return (q) ;                    /* q = least common ancester (jprev,j) */
}

#define HEAD(k,j) (ata ? head [k] : j)
#define NEXT(J)   (ata ? next [J] : -1)
#define  CS_MIN(a,b) (((a) < (b)) ? (a) : (b))
void init_ata (int n, int m, int *ATp, int *ATi, const int *post, int *w, int **head, int **next)
{
    int i, k, p ;
    *head = w+4*n, *next = w+5*n+1 ;
    for (k = 0 ; k < n ; k++) w [post [k]] = k ;    /* invert post */
    for (i = 0 ; i < m ; i++)
    {
        for (k = n, p = ATp[i] ; p < ATp[i+1] ; p++) k = CS_MIN (k, w [ATi[p]]);
        (*next) [i] = (*head) [k] ;     /* place row i in linked list k */
        (*head) [k] = i ;
    }
}
int *counts (int n, int m, int *Ap, int *Ai, double *Ax, const int *parent, const int *post, int ata)
{
    int i, j, k, J, s, p, q, jleaf, *ATp, *ATi, *maxfirst, *prevleaf,
            *ancestor, *head = NULL, *next = NULL, *colcount, *w, *first, *delta ;
    //cs *AT ;
    if (!Ap || !parent || !post) return (NULL) ;    /* check inputs */
    //m = A->m ; n = A->n ;
    s = 4*n + (ata ? (n+m+1) : 0) ;
    //delta = colcount = cs_malloc (n, sizeof (int)) ;    /* allocate result */
    delta = colcount = new int[n];
    w = new int[s];                   /* get workspace */
    ATp = new int[n+1]; ATi = new int[Ap[n]];
    transpose (n,m,Ap,Ai,Ax,0,ATp,ATi,NULL);                          /* AT = A' */
    if (!ATp || !ATi || !colcount || !w) return (colcount) ;
    ancestor = w ; maxfirst = w+n ; prevleaf = w+2*n ; first = w+3*n ;
    for (k = 0 ; k < s ; k++) w [k] = -1 ;      /* clear workspace w [0..s-1] */
    for (k = 0 ; k < n ; k++)                   /* find first [j] */
    {
        j = post [k] ;
        delta [j] = (first [j] == -1) ? 1 : 0 ;  /* delta[j]=1 if j is a leaf */
        for ( ; j != -1 && first [j] == -1 ; j = parent [j]) first [j] = k ;
    }
//    ATp = AT->p ; ATi = AT->i ;
    if (ata) init_ata (n, m, ATp, ATi, post, w, &head, &next) ;
    for (i = 0 ; i < n ; i++) ancestor [i] = i ; /* each node in its own set */
    for (k = 0 ; k < n ; k++)
    {
        j = post [k] ;          /* j is the kth node in postordered etree */
        if (parent [j] != -1) delta [parent [j]]-- ;    /* j is not a root */
        for (J = HEAD (k,j) ; J != -1 ; J = NEXT (J))   /* J=j for LL'=A case */
        {
            for (p = ATp [J] ; p < ATp [J+1] ; p++)
            {
                i = ATi [p] ;
                q = leaf (i, j, first, maxfirst, prevleaf, ancestor, &jleaf);
                if (jleaf >= 1) delta [j]++ ;   /* A(i,j) is in skeleton */
                if (jleaf == 2) delta [q]-- ;   /* account for overlap in q */
            }
        }
        if (parent [j] != -1) ancestor [j] = parent [j] ;
    }
    for (j = 0 ; j < n ; j++)           /* sum up delta's of each child */
    {
        if (parent [j] != -1) colcount [parent [j]] += colcount [j] ;
    }
    delete [] w;
    delete [] ATi;
    delete [] ATp;
    return (colcount) ;    /* success: free workspace */
}
#endif //CHOLOPENMP_COLCOUNTS_H
