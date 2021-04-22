//
// Created by kazem on 7/27/17.
//

#ifndef CHOLOPENMP_SPARSEUTILS_H
#define CHOLOPENMP_SPARSEUTILS_H

#include "def.h"

long int getNNZ
  (int ncol, int *Ap)
{
    size_t nz ;
    nz = Ap [ncol] ;
    return (nz) ;
}

size_t add_size_t (size_t a, size_t b, int *ok)
{
    size_t s = a + b ;
    (*ok) = (*ok) && (s >= a) ;
    return ((*ok) ? s : 0) ;
}


size_t mult_size_t(size_t a, size_t k, int *ok)
{
    size_t p = 0, s ;
    while (*ok)
    {
        if (k % 2)
        {
            p = p + a ;
            (*ok) = (*ok) && (p >= a) ;
        }
        k = k / 2 ;
        if (!k) return (p) ;
        s = a + a ;
        (*ok) = (*ok) && (s >= a) ;
        a = s ;
    }
    return (0) ;
}

#endif //CHOLOPENMP_SPARSEUTILS_H
