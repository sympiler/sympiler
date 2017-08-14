//
// Created by kazem on 10/06/17.
//

#ifndef SYMPILER_PROJ_NUMERICALUTILS_H
#define SYMPILER_PROJ_NUMERICALUTILS_H

#include <chrono>

#define CS_FLIP(i) (-(i)-2)
#define CS_UNFLIP(i) (((i) < 0) ? CS_FLIP(i) : (i))
#define CS_MARKED(w,j) (w [j] < 0)
#define CS_MARK(w,j) { w [j] = CS_FLIP (w [j]) ; }
#define HEAD(k,j) (ata ? head [k] : j)
#define NEXT(J)   (ata ? next [J] : -1)
#define  CS_MIN(a,b) (((a) < (b)) ? (a) : (b))
#define NULL __null

#define EMPTY -1


/* depth-first-search of the graph of a matrix, starting at node j */
int dfs (int j, int *Gp, int *Gi, int top, int *xi, int *pstack, const int *pinv);

/* xi [top...n-1] = nodes reachable from graph of G*P' via nodes in B(:,k).
 * xi [n...2n-1] used as workspace */
int reach (int n, int* Gp, int *Gi, int* Bp, int* Bi, int k, int *pruneVal, const int *pinv);

int reach_sn (int n, int* Gp, int *Gi, int* Bp, int* Bi,
              int k, int *PBset, const int *pinv, int sn,
              int *col2sup );

int* superNodeDetection(int *Parent, int *ColCount, int n,
                    int *col2sup, int &bsSize, int &avgBSize);

//// Elimination Tree computation
int *etree (int n, int* Ap, int* Ai, int ata);

//// Post order
int tdfs(int j, int k, int *head, const int *next, int *post, int *stack);

int *postOrder (const int *parent, int n);

//// Column count computations

double cumsum (int *p, int *c, int n);

int transpose (int n, int m, int *Ap, int *Ai, double *Ax,
               int values, int *Cp, int *Ci, double *Cx);

int transposePattern (int n, int m, int *Ap, int *Ai,
               int values, int *Cp, int *Ci, double *Cx);

/* consider A(i,j), node j in ith row subtree and return lca(jprev,j) */
int leaf (int i, int j, const int *first, int *maxfirst, int *prevleaf,
          int *ancestor, int *jleaf);

void init_ata (int n, int m, int *ATp, int *ATi, const int *post,
               int *w, int **head, int **next);

int *counts (int n, int m, int *Ap, int *Ai, const int *parent,
             const int *post, int ata);

//// Prune-set calculation
/* find nonzero pattern of Cholesky L(k,1:k-1) using etree and triu(A(:,k)) */
int ereach (int n, int *Ap, int *Ai, int k, const int *parent, int *s, int *w);

/* find nonzero pattern of Cholesky L(k,1:k-1) using etree and triu(A(:,k)) */
int ereach_sn (int n, int *Ap, int *Ai, int col1,int col2, int *col2sup,
               const int *parent, int *s, int *w);

bool pruning(int n, int* c, int* r, int *ET, int* &prunePtr, int* &pruneSet);

#endif //SYMPILER_PROJ_NUMERICALUTILS_H
