//
// Created by kazem on 7/25/17.
//

#ifndef CHOLOPENMP_ORDERING_H
#define CHOLOPENMP_ORDERING_H

bool permute (int n, int m, int *Ap, int *Ai, double *Ax, const int *pinv, const int *q,
              int *Cp, int *Ci, double *Cx)
/* C = A(p,q) where p and q are permutations of 0..m-1 and 0..n-1.
 * pinv is the inverse permutation of P
 * */
{
    int t, j, k, nz = 0 ;
    if (!Ap || !Ai) return false;    /* check inputs */
    //C = cs_spalloc (m, n, Ap [n], values && Ax != NULL, 0) ;  /* alloc result */
    if (!Cp || !Ci) return false;   /* out of memory */
    for (k = 0 ; k < n ; k++)
    {
        Cp [k] = nz ;                   /* column k of C is column q[k] of A */
        j = q ? (q [k]) : k ;
        for (t = Ap [j] ; t < Ap [j+1] ; t++)
        {
            if (Cx) Cx [nz] = Ax [t] ;  /* row i of A is row pinv[i] of C */
            Ci [nz++] = pinv ? (pinv [Ai [t]]) : Ai [t] ;
        }
    }
    Cp [n] = nz ;                       /* finalize the last column of C */
    return true;
}

#endif //CHOLOPENMP_ORDERING_H
