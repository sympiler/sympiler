//
// Created by kazem on 7/25/17.
//

#ifndef CHOLOPENMP_ETREE_H
#define CHOLOPENMP_ETREE_H


#include "def.h"
#include "SparseUtils.h"

void update_etree
        (

                int k,
                int i,

                int *curTree ,
                int *Ancestor
        )
{
    int a ;
    for ( ; ; )		/* traverse the path from k to the root of the tree */
    {
        a = Ancestor [k] ;
        if (a == i)
        {
            return ;
        }
        Ancestor [k] = i ;
        if (a == EMPTY)
        {
            curTree [k] = i ;
            return ;
        }
        k = a ;
    }
}



/* Find the elimination tree of A
 * */

int etree_sym
  (CSC *A, int *Parent)
{
    int *Ap, *Ai, *Anz, *Ancestor, *Iwork ;
    int i, j, p, pend,ncol, packed, stype ;
    size_t s ;
    int ok = TRUE ;

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    stype = A->stype ;
    s = add_size_t (A->nrow, (stype ? 0 : A->ncol), &ok) ;
    if (!ok)
    {
        return (FALSE) ;
    }

    Iwork = new int[s]();

    ncol = A->ncol ;
    Ap = A->p ;
    Ai = A->i ;
    Anz = A->nz ;
    packed = A->packed ;
    Ancestor = Iwork ;

    for (j = 0 ; j < ncol ; j++)
    {
        Parent [j] = EMPTY ;
        Ancestor [j] = EMPTY ;
    }

    /* ---------------------------------------------------------------------- */
    /* compute the etree */
    /* ---------------------------------------------------------------------- */

    if(stype<=0)
        return FALSE;
    for (j = 0 ; j < ncol ; j++){
        /* for each row i in column j of triu(A), excluding the diagonal */
        p = Ap [j] ;
        pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
        for ( ; p < pend ; p++){
            i = Ai [p] ;
            if (i < j){
                update_etree (i, j, Parent, Ancestor) ;
            }
        }
    }

    delete []Iwork;
    return (TRUE) ;
}




int *etree (int n, int* Ap, int* Ai, int ata){
    int i, k, p, m, inext, *w, *parent, *ancestor, *prev ;
    if (n<0 || Ap == NULL || Ai == NULL) //check inputs
        return 0;
    m = n;
    parent = new int[n]; //result allocation
    w = new int[n + (ata ? m : 0)]; // get workspace
    if (w == NULL || parent == NULL)
        return 0;

    ancestor = w ; prev = w + n ;
    if (ata) for (i = 0 ; i < m ; i++) prev [i] = -1 ;
    for (k = 0 ; k < n ; k++){
        parent [k] = -1 ;                   /* node k has no parent yet */
        ancestor [k] = -1 ;                 /* nor does k have an ancestor */
        for (p = Ap [k] ; p < Ap [k+1] ; p++){
            i = ata ? (prev [Ai [p]]) : (Ai [p]) ;
            for ( ; i != -1 && i < k ; i = inext){
                inext = ancestor [i] ;              /* inext = ancestor of i */
                ancestor [i] = k ;                  /* path compression */
                if (inext == -1) parent [i] = k ;   /* no anc., parent is k */
            }
            if (ata) prev [Ai [p]] = k ;
        }
    }
    delete []w;
    return parent;
}



#endif //CHOLOPENMP_ETREE_H
