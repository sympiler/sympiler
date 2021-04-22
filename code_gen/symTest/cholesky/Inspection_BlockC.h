//
// Created by kazem on 7/28/17.
//

#ifndef CHOLOPENMP_INSPECTION_BLOCKC_H
#define CHOLOPENMP_INSPECTION_BLOCKC_H
/* ========================================================================== */
/* === subtree ============================================================== */
/* ========================================================================== */

#include "def.h"
#include "SparseUtils.h"


void subtree
        (
                /* inputs, not modified: */
                int j,		/* j = k for symmetric case */
                int k,
                int Ap [ ],
                int Ai [ ],
                int Anz [ ],
                int SuperMap [ ],
                int Sparent [ ],
                int mark,
                int sorted,         /* true if the columns of A are sorted */
                int k1,             /* only consider A (0:k1-1,k) */

                /* input/output: */
                int Flag [ ],
                int Ls [ ],
                int Lpi2 [ ]
        )
{
    int p, pend, i, si ;
    p = Ap [j] ;
    //pend = (Anz == NULL) ? (Ap [j+1]) : (p + Anz [j]) ;
    pend = (Ap [j+1]) ;

    for ( ; p < pend ; p++)
    {
        i = Ai [p] ;
        if (i < k1)
        {

            for (si = SuperMap [i] ; Flag [si] < mark ; si = Sparent [si])
            {
                ASSERT (si <= SuperMap [k]) ;
                Ls [Lpi2 [si]++] = k ;
                Flag [si] = mark ;
            }
        }
        else if (sorted)
        {
            break ;
        }
    }
}




int supernodalAnalysis
  (CSC *A, int *Parent, BCSC *L, int *nrelax, double *zrelax,
   int *Sparent)
{
    double zrelax0, zrelax1, zrelax2, xxsize ;
    int *Wi, *Wj, *Super, *Snz, *Ap, *Ai, *Flag, *Head, *Ls, *Lpi, *Lpx,
             *Anz, *SuperMap, *Merged, *Nscol, *Zeros,
            *ColCount, *Lpi2, *Lsuper, *Iwork ;
    int nsuper, n, j, k, s, mark, parent, p, k1, k2, nscol,
            nsrow, stype, ssize, xsize, sparent, ss, nscol0, nscol1, ns, nfsuper, newzeros, totzeros,
            merge, snext, nrelax0, nrelax1, nrelax2, Asorted ;
    size_t w ;
    int ok = TRUE ;


    stype = A->stype ;
    if (stype < 0)
    {

        return (FALSE) ;
    }
    if (stype == 0)
    {
        return (FALSE);
    }
    n = A->nrow ;
    /* w = 5*n */
    w = mult_size_t(n, 4, &ok) ;
    if (!ok) {
        return (FALSE) ;
    }
    Ap = A->p ;
    Ai = A->i ;
    Anz = A->nz ;

    ColCount = L->ColCount ;
    nrelax0 = nrelax [0] ;
    nrelax1 = nrelax [1] ;
    nrelax2 = nrelax [2] ;

    zrelax0 = zrelax [0] ;
    zrelax1 = zrelax [1] ;
    zrelax2 = zrelax [2] ;
    Iwork = new int[w]();
    Wi      = Iwork ;	    /* size n (i/l/l).  Lpi2 is i/l/l */
    Wj      = Iwork + n ;   /* size n (i/l/l).  Zeros is i/l/l */
    Snz     = Iwork + 2*((size_t) n) ; /* size nfsuper <= n [ */
    Merged  = Iwork + 3*((size_t) n) ; /* size nfsuper <= n [ */


    Flag = new int[n]();//Common->Flag ;   /* size n */
    Head = new int[n+1]();//Common->Head ;   /* size n+1 */
    for (j = 0 ; j < n ; j++)
    {
        Wi [j] = 0 ;
    }
    for (j = 0 ; j < n ; j++)
    {
        parent = Parent [j] ;
        if (parent != EMPTY)
        {
            Wi [parent]++ ;
        }
    }
    Super = Head ;  /* use Head [0..nfsuper] as workspace for Super list ( */
    nfsuper = (n == 0) ? 0 : 1 ;	/* number of fundamental supernodes */
    Super [0] = 0 ;

    for (j = 1 ; j < n ; j++)
    {
        if (Parent [j-1] != j	    /* parent of j-1 is not j */
            || (ColCount [j-1] != ColCount [j] + 1) /* j-1 not subset of j*/
            || Wi [j] > 1	    /* j has more than one child */

                )
        {
            /* j is the leading node of a supernode */
            Super [nfsuper++] = j ;
        }
    }
    Super [nfsuper] = n ;
    Nscol = Wi ; /* use Wi as size-nfsuper workspace for Nscol [ */
    SuperMap = Wj ;	/* use Wj as workspace for SuperMap [ */

    for (s = 0 ; s < nfsuper ; s++)
    {
        for (k = Super [s] ; k < Super [s+1] ; k++)
        {
            SuperMap [k] = s ;
        }
    }

    for (s = 0 ; s < nfsuper ; s++)
    {
        j = Super [s+1] - 1 ;	/* last node in supernode s */
        parent = Parent [j] ;	/* parent of last node */
        Sparent [s] = (parent == EMPTY) ? EMPTY : SuperMap [parent] ;
    }

    Zeros = Wj ;   /* use Wj for Zeros, workspace of size nfsuper [ */

    /* relaxed amalgamation */

    for (s = 0 ; s < nfsuper ; s++)
    {
        Merged [s] = EMPTY ;			/* s not merged into another */
        Nscol [s] = Super [s+1] - Super [s] ;	/* # of columns in s */
        Zeros [s] = 0 ;				/* # of zero entries in s */
        ASSERT (s <= Super [s]) ;
        Snz [s] = ColCount [Super [s]] ;  /* # of entries in leading col of s */
    }

    for (s = nfsuper-2 ; s >= 0 ; s--)
    {
        double lnz1 ;
        ss = Sparent [s] ;
        if (ss == EMPTY)
        {
            continue ;
        }

        /* find the current parent of s (perform path compression as needed) */
        for (ss = Sparent [s] ; Merged [ss] != EMPTY ; ss = Merged [ss]) ;
        sparent = ss ;
        /* ss is the current parent of s */
        for (ss = Sparent [s] ; Merged [ss] != EMPTY ; ss = snext)
        {
            snext = Merged [ss] ;
            Merged [ss] = sparent ;
        }

        /* if s+1 is not the current parent of s, do not merge */
        if (sparent != s+1)
        {
            continue ;
        }

        nscol0 = Nscol [s] ;	/* # of columns in s */
        nscol1 = Nscol [s+1] ;	/* # of columns in s+1 */
        ns = nscol0 + nscol1 ;
        totzeros = Zeros [s+1] ;	/* current # of zeros in s+1 */
        lnz1 = (double) (Snz [s+1]) ;	/* # entries in leading column of s+1 */

        /* determine if supernodes s and s+1 should merge */
        if (ns <= nrelax0)
        {
            merge = TRUE ;
        }
        else
        {
            /* use double to avoid integer overflow */
            double lnz0 = Snz [s] ;	/* # entries in leading column of s */
            double xnewzeros = nscol0 * (lnz1 + nscol0 - lnz0) ;

            /* use int for the final update of Zeros [s] below */
            newzeros = nscol0 * (Snz [s+1] + nscol0 - Snz [s]) ;
            ASSERT (newzeros == xnewzeros) ;

            if (xnewzeros == 0)
            {
                merge = TRUE ;
            }
            else
            {
                double xtotzeros = ((double) totzeros) + xnewzeros ;
                double xns = (double) ns ;
                double xtotsize  = (xns * (xns+1) / 2) + xns * (lnz1 - nscol1) ;
                double z = xtotzeros / xtotsize ;
                totzeros += newzeros ;
                merge = ((ns <= nrelax1 && z < zrelax0) ||
                         (ns <= nrelax2 && z < zrelax1) ||
                         (z < zrelax2)) &&
                        (xtotsize < Int_max / sizeof (double)) ;
            }
        }
        if (merge)
        {
            Zeros [s] = totzeros ;
            Merged [s+1] = s ;
            Snz [s] = nscol0 + Snz [s+1] ;
            Nscol [s] += Nscol [s+1] ;
        }
    }

    nsuper = 0 ;
    for (s = 0 ; s < nfsuper ; s++)
    {
        if (Merged [s] == EMPTY)
        {
            Super [nsuper] = Super [s] ;
            Snz [nsuper] = Snz [s] ;
            nsuper++ ;
        }
    }
    Super [nsuper] = n ;

    for (s = 0 ; s < nsuper ; s++)
    {
        for (k = Super [s] ; k < Super [s+1] ; k++)
        {
            SuperMap [k] = s ;
        }
    }
#if 0
    printf("number of SN is %d: ", nsuper);
    for (int i = 0; i < n; ++i) {
        printf("%d ", SuperMap[i]);
    }
    printf("-----\n");
#endif
    /* construct the relaxed supernodal etree */

    for (s = 0 ; s < nsuper ; s++)
    {
        j = Super [s+1] - 1 ;	/* last node in supernode s */
        parent = Parent [j] ;	/* parent of last node */
        Sparent [s] = (parent == EMPTY) ? EMPTY : SuperMap [parent] ;
    }

    ssize = 0 ;
    xsize = 0 ;
    xxsize = 0 ;
    for (s = 0 ; s < nsuper ; s++)
    {
        nscol = Super [s+1] - Super [s] ;
        nsrow = Snz [s] ;
        ASSERT (nscol > 0) ;
        ssize += nsrow ;
        xsize += nscol * nsrow ;
        xxsize += ((double) nscol) * ((double) nsrow) ;
        if (ssize < 0 ||xxsize > Int_max){
            return (FALSE) ;
        }
        ASSERT (ssize > 0) ;
    }
    xsize = MAX (1, xsize) ;
    ssize = MAX (1, ssize) ;

    L->ssize = ssize ;
    L->xsize = xsize ;
    L->nsuper = nsuper ;

    allocateLC(L, true);
    Lpi = L->pi ;
    Lpx = L->px ;
    Ls = L->s ;
    Ls [0] = 0 ;    /* flag for cholmod_check_factor; supernodes are defined */
    Lsuper = L->super ;

    for (s = 0 ; s <= nsuper ; s++)
    {
        Lsuper [s] = Super [s] ;
    }
    Super = Lsuper ;	    /* Super is now the list of relaxed supernodes */

    p = 0 ;
    for (s = 0 ; s < nsuper ; s++)
    {
        Lpi [s] = p ;
        p += Snz [s] ;
    }
    Lpi [nsuper] = p ;
    ASSERT ((int) (L->ssize) == MAX (1,p)) ;

    /* ---------------------------------------------------------------------- */
    /* construct pointers for supernodal values (L->px) */
    /* ---------------------------------------------------------------------- */
    Lpx [0] = 0 ;
    p = 0 ;
    for (s = 0 ; s < nsuper ; s++)
    {
        nscol = Super [s+1] - Super [s] ;   /* number of columns in s */
        nsrow = Snz [s] ;           /* # of rows, incl triangular part*/
        Lpx [s] = p ;               /* pointer to numerical part of s */
        p += nscol * nsrow ;
    }
    Lpx [s] = p ;
    ASSERT ((int) (L->xsize) == MAX (1,p)) ;
    Lpi2 = new int[n+1]();
    for (s = 0 ; s < nsuper ; s++)
    {
        Lpi2 [s] = Lpi [s] ;
    }

    Asorted = A->sorted ;
    mark = EMPTY;
    for (s = 0 ; s < nsuper ; s++)
    {
        k1 = Super [s] ;
        k2 = Super [s+1] ;
        for (k = k1 ; k < k2 ; k++)
        {
            Ls [Lpi2 [s]++] = k ;
        }

        /* compute nonzero pattern each row k1 to k2-1 */
        for (k = k1 ; k < k2 ; k++)
        {
            mark++;
            if(mark <= 0){
                for (int i = 0 ; i < n ; i++)
                {
                    Flag [i] = EMPTY ;
                }
                mark = 0;
            }
            Flag [s] = mark ;
            ASSERT (s == SuperMap [k]) ;

            /* traverse the row subtree for each nonzero in A or AA' */
            subtree (k, k, Ap, Ai, Anz, SuperMap, Sparent, mark,
                         Asorted, k1, Flag, Ls, Lpi2) ;

        }
    }
#if 1
    for (s = 0 ; s < nsuper ; s++)
    {
        ASSERT (Lpi2 [s] == Lpi [s+1]) ;
    }
#endif

    delete []Iwork;
    delete []Flag;
    delete []Head;

    return (TRUE) ;
}
#endif //CHOLOPENMP_INSPECTION_BLOCKC_H
