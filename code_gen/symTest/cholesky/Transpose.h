//
// Created by kazem on 7/25/17.
//

#ifndef CHOLOPENMP_TRANSPOSE_H
#define CHOLOPENMP_TRANSPOSE_H

#include "def.h"
#include "SparseUtils.h"

int transpose_sym_real
  (CSC *A, int *Perm, CSC *F, int *Wi, int *Pinv)
{
    double *Ax,*Fx ;
    int *Ap, *Anz, *Ai, *Fj ;
    int p, pend, packed, fp, upper, permute, jold, n, i, j, iold ;

    permute = (Perm != NULL) ;
    n = A->nrow ;
    Ap = A->p ;		/* size A->ncol+1, column pointers of A */
    Ai = A->i ;		/* size nz = Ap [A->ncol], row indices of A */
    Ax = A->x ;		/* size nz, real values of A */
    Anz = A->nz ;
    packed = A->packed ;
    upper = (A->stype > 0) ;

    Fj = F->i ;		/* size nz, column indices of F */
    Fx = F->x ;		/* size nz, real values of F */
    if (permute)
    {
        if (upper)
        {
            /* permuted, upper */
            for (j = 0 ; j < n ; j++)
            {
                jold = Perm [j] ;
                p = Ap [jold] ;
                pend = (packed) ? Ap [jold+1] : p + Anz [jold] ;
                for ( ; p < pend ; p++)
                {
                    iold = Ai [p] ;
                    if (iold <= jold)
                    {
                        i = Pinv [iold] ;
                        if (i < j)
                        {
                            fp = Wi [i]++ ;
                            Fj [fp] = j ;
                            Fx[fp] = Ax[p];
                        }
                        else
                        {
                            fp = Wi [j]++ ;
                            Fj [fp] = i ;
                            Fx[fp] = Ax[p];
                        }
                    }
                }
            }
        }
        else
        {
            /* permuted, lower */
            for (j = 0 ; j < n ; j++)
            {
                jold = Perm [j] ;
                p = Ap [jold] ;
                pend = (packed) ? Ap [jold+1] : p + Anz [jold] ;
                for ( ; p < pend ; p++)
                {
                    iold = Ai [p] ;
                    if (iold >= jold)
                    {
                        i = Pinv [iold] ;
                        if (i > j)
                        {
                            fp = Wi [i]++ ;
                            Fj [fp] = j ;
                            Fx[fp] = Ax[p];
                        }
                        else
                        {
                            fp = Wi [j]++ ;
                            Fj [fp] = i ;
                            Fx[fp] = Ax[p];
                        }
                    }
                }
            }
        }
    }
    else
    {
        if (upper)
        {
            /* unpermuted, upper */
            for (j = 0 ; j < n ; j++)
            {
                p = Ap [j] ;
                pend = (packed) ? Ap [j+1] : p + Anz [j] ;
                for ( ; p < pend ; p++)
                {
                    i = Ai [p] ;
                    if (i <= j)
                    {
                        fp = Wi [i]++ ;
                        Fj [fp] = j ;
                        Fx[fp] = Ax[p];
                    }
                }
            }
        }
        else
        {
            /* unpermuted, lower */
            for (j = 0 ; j < n ; j++)
            {
                p = Ap [j] ;
                pend = (packed) ? Ap [j+1] : p + Anz [j] ;
                for ( ; p < pend ; p++)
                {
                    i = Ai [p] ;
                    if (i >= j)
                    {
                        fp = Wi [i]++ ;
                        Fj [fp] = j ;
                        Fx[fp] = Ax[p];
                    }
                }
            }
        }
    }

    return (TRUE) ;
}

int transpose_sym
  (CSC *A, int values, int *Perm, CSC *F)
{
    int *Ap, *Ai, *Fp, *Wi, *Pinv, *Iwork ;
    int p, pend, upper, permute, jold, n, i, j, k, iold ;
    size_t s ;
    int ok = TRUE ;
    if (A->nrow != A->ncol || A->stype == 0)
    {
        return (FALSE) ;
    }
    if (A->nrow != F->ncol || A->ncol != F->nrow)
    {
        return (FALSE) ;
    }

    permute = (Perm != NULL) ;
    n = A->nrow ;
    Ap = A->p ;		/* size A->ncol+1, column pointers of A */
    Ai = A->i ;		/* size nz = Ap [A->ncol], row indices of A */
    upper = (A->stype > 0) ;

    Fp = F->p ;		/* size A->nrow+1, row pointers of F */

    s = add_size_t (n, ((Perm != NULL) ? n : 0), &ok) ;
    if (!ok)
    {
        return (FALSE) ;
    }

    Iwork = new int[s]();
    Wi   = Iwork ;	    /* size n (i/l/l) */
    Pinv = Iwork + n ;	    /* size n (i/i/l) , unused if Perm NULL */
    if (permute)
    {
        for (i = 0 ; i < n ; i++)
        {
            Pinv [i] = EMPTY ;
        }
        for (k = 0 ; k < n ; k++)
        {
            i = Perm [k] ;
            if (i < 0 || i > n || Pinv [i] != EMPTY)
            {
                return (FALSE) ;
            }
            Pinv [i] = k ;
        }
    }

    for (i = 0 ; i < n ; i++)
    {
        Wi [i] = 0 ;
    }


    if (permute)
    {
        if (upper)
        {
            /* packed, permuted, upper */
            for (j = 0 ; j < n ; j++)
            {
                jold = Perm [j] ;
                pend = Ap [jold+1] ;
                for (p = Ap [jold] ; p < pend ; p++)
                {
                    iold = Ai [p] ;
                    if (iold <= jold)
                    {
                        i = Pinv [iold] ;
                        Wi [MIN (i, j)]++ ;
                    }
                }
            }
        }
        else
        {
            /* packed, permuted, lower */
            for (j = 0 ; j < n ; j++)
            {
                jold = Perm [j] ;
                pend = Ap [jold+1] ;
                for (p = Ap [jold] ; p < pend ; p++)
                {
                    iold = Ai [p] ;
                    if (iold >= jold)
                    {
                        i = Pinv [iold] ;
                        Wi [MAX (i, j)]++ ;
                    }
                }
            }
        }
    }
    else
    {
        if (upper)
        {
            /* packed, unpermuted, upper */
            for (j = 0 ; j < n ; j++)
            {
                pend = Ap [j+1] ;
                for (p = Ap [j] ; p < pend ; p++)
                {
                    i = Ai [p] ;
                    if (i <= j)
                    {
                        Wi [i]++ ;
                    }
                }
            }
        }
        else
        {
            /* packed, unpermuted, lower */
            for (j = 0 ; j < n ; j++)
            {
                pend = Ap [j+1] ;
                for (p = Ap [j] ; p < pend ; p++)
                {
                    i = Ai [p] ;
                    if (i >= j)
                    {
                        Wi [i]++ ;
                    }
                }
            }
        }
    }

    p = 0 ;
    for (i = 0 ; i < n ; i++)
    {
        Fp [i] = p ;
        p += Wi [i] ;
    }
    Fp [n] = p ;
    for (i = 0 ; i < n ; i++)
    {
        Wi [i] = Fp [i] ;
    }

    if (p > (int) (F->nzmax))
    {
        return (FALSE) ;
    }

    ok = FALSE ;
    if (values == 0 || F->xtype == MATRIX_PATTERN)
    {
        ok = transpose_sym_real(A, Perm, F, Wi, Pinv);
    }
    else if (F->xtype == MATRIX_REAL)
    {
        ok = transpose_sym_real(A, Perm, F, Wi, Pinv);
    }
    if (ok)
    {
        F->sorted = !permute ;
        F->packed = TRUE ;
        F->stype = - SIGN (A->stype) ;	/* flip the stype */
    }
    delete []Iwork;
    return (ok) ;
}

CSC *transposeMatrix
  (CSC *A, int values, int *Perm, int *fset, size_t fsize)
{
    int *Ap, *Anz ;
    CSC *F ;
    int ncol, use_fset, j, jj, fnz, packed, stype, nf ;
    int ok = TRUE ;

    nf = fsize ;
    stype = A->stype ;
    ncol = A->ncol ;

    if (stype != 0)
    {
        use_fset = FALSE ;
        if (Perm != NULL)
        {
            mult_size_t (A->nrow, 2, &ok) ;
        }
    }
    else
    {
        use_fset = (fset != NULL) ;
        if (use_fset)
        {
            MAX (A->nrow, A->ncol) ;
        }
    }

    if (!ok)
    {
        return (NULL) ;
    }

    Ap = A->p ;
    Anz = A->nz ;
    packed = A->packed ;

    if (stype != 0)
    {
        fnz = getNNZ(A->ncol, A->p);
    }
    else
    {
        nf = (use_fset) ? nf : ncol ;
        if (use_fset)
        {
            fnz = 0 ;
            for (jj = 0 ; jj < nf ; jj++)
            {
                j = fset [jj] ;
                if (j >= 0 && j < ncol)
                {
                    fnz += packed ? (Ap [j+1] - Ap [j]) : MAX (0, Anz [j]) ;
                }
            }
        }
        else
        {
            fnz = getNNZ(A->ncol, A->p);
        }
    }

    F = new CSC;
    allocateAC(F,A->nrow,fnz,-SIGN(A->stype),TRUE);

    if (stype != 0)
    {
        ok = transpose_sym(A, values, Perm, F);
    }
    else
    {
        printf("unsym is not supported.");
    }
    return (F) ;
}





//***********CSPARSE rountines
double cumsum (int *p, int *c, int n)
{
    int i, nz = 0 ;
    double nz2 = 0 ;
    if (!p || !c) return (-1) ;     /* check inputs */
    for (i = 0 ; i < n ; i++)
    {
        p [i] = nz ;
        nz += c [i] ;
        nz2 += c [i] ;              /* also in double to avoid csi overflow */
        c [i] = p [i] ;             /* also copy p[0..n-1] back into c[0..n-1]*/
    }
    p [n] = nz ;
    return (nz2) ;                  /* return sum (c [0..n-1]) */
}

int transpose (int n, int m, int *Ap, int *Ai, double *Ax, int values,
               int *Cp, int *Ci, double *Cx)
{
    int p, q, j, *w ;
    if (!Ai || !Ap) return (FALSE) ;    /* check inputs */

    w = new int[m]();
    if (!Cp || !w) return -1 ;       /* out of memory */
    //Cp = C->p ; Ci = C->i ; Cx = C->x ;
    for (p = 0 ; p < Ap [n] ; p++) w [Ai [p]]++ ;          /* row counts */
    cumsum (Cp, w, m) ;                                 /* row pointers */
    for (j = 0 ; j < n ; j++)
    {
        for (p = Ap [j] ; p < Ap [j+1] ; p++)
        {
            Ci [q = w [Ai [p]]++] = j ; /* place A(i,j) as entry C(j,i) */
            if (Cx) Cx [q] = Ax [p] ;
        }
    }
    delete []w;
    return 1 ;  /* success;  */
}



#endif //CHOLOPENMP_TRANSPOSE_H
