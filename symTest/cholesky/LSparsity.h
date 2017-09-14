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

#include "amd.h"
int permute_matrices
  (CSC *A, int *Perm, int do_rowcolcounts, CSC **A1_handle,
                       CSC **A2_handle, CSC **S_handle, CSC **F_handle,
                       int &status)
{
    CSC *A1, *A2, *S, *F ;

    *A1_handle = NULL ;
    *A2_handle = NULL ;
    *S_handle = NULL ;
    *F_handle = NULL ;
    A1 = NULL ;

    A2 = transposeMatrix(A, 0, Perm, NULL, 0);
    S = A2 ;
    if (do_rowcolcounts)
    {
        A1 = transposeMatrix(A2, 0, NULL, NULL, 0);
    }
    F = A1 ;

    *A1_handle = A1 ;
    *A2_handle = A2 ;
    *S_handle = S ;
    *F_handle = F ;
    return (TRUE) ;
}



int insGraph
  (CSC *A, int *Perm, int *fset, size_t fsize, int *Parent, int *Post,
               int *ColCount, int *First, int *Level, int &status)
{
    CSC *A1, *A2, *S, *F ;
    int n, ok, do_rowcolcounts ;
    int fl, aatfl, lnz;

    n = A->nrow ;

    do_rowcolcounts = (ColCount != NULL) ;
    ok = permute_matrices(A, Perm, do_rowcolcounts, &A1, &A2, &S, &F,
                          status);
    ok = ok && etree_sym(A->stype ? S : F, Parent);
#if 0
    for (int i = 0; i < n; ++i) {
        printf("%d :", Parent[i]);
    }
    printf("\n");
#endif

    ok = ok && postOrderC(Parent, n, NULL, Post) == n ;

    if (do_rowcolcounts)
    {
        ok = ok && rowcolcounts(A->stype ? F:S, fset, fsize, Parent,
                                          Post, NULL, ColCount, First, Level,
                                fl,
                                aatfl,
                                lnz, status) ;
    }


    allocateAC(A1,0,0,0,FALSE);
    delete A1;
    allocateAC(A2,0,0,0,FALSE);
    delete A2;
    return (ok) ;
}



BCSC *buildLFactor
  (CSC *A, int *fset, int *nrelax, double *zrelax, size_t fsize,
                     int *&prunePtr, int *&pruneSet, int status)
{
    int *First, *Level, *Work4n, *Lperm,
            *Post, *Lparent, *Lcolcount ;
    BCSC *L ;
    int n, ncol, k ;
    size_t s ;

    n = A->nrow ;
    ncol = A->ncol ;
    /*
     * Allocate Workspace
     */

    s=4*n;
    Work4n = new int[s]() ;
   // Work4n += 2*((size_t) n) + uncol ;
    First  = Work4n + n ;
    Level  = Work4n + 2*((size_t) n) ;
    Post   = Work4n + 3*((size_t) n) ;

    L = new BCSC;
    L->Perm     = new int[n]();
    L->ColCount = new int[n]();
    Lparent  = new int[n]();//CHOLMOD(malloc) (n, sizeof (int), Common) ;
    Lperm = L->Perm ;
    Lcolcount = L->ColCount ;

    double info[20]={0};

    amd_order(ncol,A->p,A->i,Lperm,NULL,info);

    if (!insGraph(A, Lperm, fset, fsize, Lparent, Post, Lcolcount, First,
                  Level, status))
    {
        return FALSE;
    }

    if (postOrderC(Lparent, n, Lcolcount, Post) == n)
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
    }

    int *Sparent = new int[n]();
    CSC *S, *F, *A2, *A1 ;

    permute_matrices(A, Lperm, TRUE, &A1, &A2, &S, &F,
                     status);

    supernodalAnalysis(S, Lparent, L, nrelax, zrelax, Sparent);

    prunePtr = new int[L->nsuper+1]();
    pruneSet = new int[L->ssize];
    int *col2Sup = new int[F->nrow];
    // Computing col2sup
    for (unsigned i = 0; i < L->nsuper; ++i) {
        int k1 = L->super[i];
        int k2 = L->super[i+1];
        for (int j = k1; j < k2; ++j) {
            col2Sup[j] = i;
        }
    }
    getBlockedPruneSet(L->nsuper,S->p,S->i,col2Sup,Sparent,L->super,
                       prunePtr,pruneSet);

    delete []col2Sup;
    allocateAC(A1,0,0,0,FALSE);
    allocateAC(A2,0,0,0,FALSE);


    delete []Sparent;
    delete []Work4n;
    return L;
}
//#endif
#endif //CHOLOPENMP_LSPARSITY_H
