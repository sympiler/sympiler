/*******************************************************************************
* Copyright 2004-2017 Intel Corporation All Rights Reserved.
*
* The source code,  information  and material  ("Material") contained  herein is
* owned by Intel Corporation or its  suppliers or licensors,  and  title to such
* Material remains with Intel  Corporation or its  suppliers or  licensors.  The
* Material  contains  proprietary  information  of  Intel or  its suppliers  and
* licensors.  The Material is protected by  worldwide copyright  laws and treaty
* provisions.  No part  of  the  Material   may  be  used,  copied,  reproduced,
* modified, published,  uploaded, posted, transmitted,  distributed or disclosed
* in any way without Intel's prior express written permission.  No license under
* any patent,  copyright or other  intellectual property rights  in the Material
* is granted to  or  conferred  upon  you,  either   expressly,  by implication,
* inducement,  estoppel  or  otherwise.  Any  license   under such  intellectual
* property rights must be express and approved by Intel in writing.
*
* Unless otherwise agreed by Intel in writing,  you may not remove or alter this
* notice or  any  other  notice   embedded  in  Materials  by  Intel  or Intel's
* suppliers or licensors in any way.
*******************************************************************************/

/*
*   Content : Intel(R) MKL PARDISO C example
*
********************************************************************************
*/
/* -------------------------------------------------------------------- */
/* Example program to show the use of the "PARDISO" routine */
/* on symmetric linear systems */
/* -------------------------------------------------------------------- */
/* This program can be downloaded from the following site: */
/* www.pardiso-project.org */
/* */
/* (C) Olaf Schenk, Department of Computer Science, */
/* University of Basel, Switzerland. */
/* Email: olaf.schenk@unibas.ch */
/* -------------------------------------------------------------------- */
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <cstring>
#include <sstream>
#include <iostream>
#include <chrono>
#include <mkl.h>
#include "omp.h"
#include "../../util/Util.h"


MKL_INT main(int argc, char *argv[]) {
#ifdef COLDCACHE
    if(argc<3){
        std::cout<<"Two arguments are needed, input matrix path and a "
                "large matrix path\n";
        return -1;
    }
#else
    if(argc<2){
        printf("Please enter a path for the input matrix");
        return -1;
    }
#endif
    std::string inpath;
    inpath = argv[1];
    printf("%s \n",inpath.c_str());
    int tmp =1;
    if (argc >= 4)
        tmp = atoi(argv[3]);
    int num_thread = tmp; //Default is serial

    /* Matrix data. */
    double  *a ;
    MKL_INT    *ja ;
    MKL_INT    *ia ;
    MKL_INT n, nnz;

    if (!readMatrix(inpath,n,nnz,ia,ja,a)){
        printf("Input file is not loaded! \n");
        return -1;
    }
#ifdef COLDCACHE
    std::string waste=argv[2];
    std::ifstream wasteFile;
    wasteFile.open(waste);
    if(wasteFile.fail()){
        std::cout<<"The input file for simulating cold-cache conditions does "
                "not exist\n";
        return -1;
    }
#endif

    MKL_INT mtype = 2;       /* Real symmetric positive  matrix */
    /* RHS and solution vectors. */
    //double b[8], x[8];
    double *b = (double*)mkl_calloc(n,sizeof(double),64);
    double *x = (double*)mkl_calloc(n,sizeof(double),64);
    MKL_INT nrhs = 1;     /* Number of right hand sides. */
    /* Internal solver memory pointer pt, */
    /* 32-bit: int pt[64]; 64-bit: long int pt[64] */
    /* or void *pt[64] should be OK on both architectures */
    long int pt[64];
    /* Pardiso control parameters. */
    MKL_INT iparm[64];
    MKL_INT maxfct, mnum, phase, error, msglvl;
    /* Auxiliary variables. */
    MKL_INT i;
    double ddum;          /* Double dummy */
    MKL_INT idum;         /* Integer dummy. */
/* -------------------------------------------------------------------- */
/* .. Setup Pardiso control parameters. */
/* -------------------------------------------------------------------- */
    for ( i = 0; i < 64; i++ )
    {
        iparm[i] = 0;
    }
    iparm[0] = 1;         /* No solver default */
    iparm[1] = 0;         /* Fill-in reordering from METIS */
    iparm[3] = 0;         /* No iterative-direct algorithm */
    iparm[4] = 0;         /* No user fill-in reducing permutation */
    iparm[5] = 0;         /* Write solution into x */
    iparm[6] = 0;         /* Not in use */
    iparm[7] = 2;         /* Max numbers of iterative refinement steps */
    iparm[8] = 0;         /* Not in use */
    iparm[9] = 8;        /* Perturb the pivot elements with 1E-8 */
    iparm[10] = 0;        /* Use nonsymmetric permutation and scaling MPS */
    iparm[11] = 0;        /* A^TX=B */
    iparm[12] = 0;        /* Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
    iparm[13] = 0;        /* Output: Number of perturbed pivots */
    iparm[14] = 0;        /* Not in use */
    iparm[15] = 0;        /* Not in use */
    iparm[16] = 0;        /* Not in use */
    iparm[17] = -1;       /* Output: Number of nonzeros in the factor LU */
    iparm[18] = -1;       /* Output: Mflops for LU factorization */
    iparm[19] = 0;        /* Output: Numbers of CG Iterations */

    //
    iparm[26] = 1;
    iparm[34] = 1;
    maxfct = 1;           /* Maximum number of numerical factorizations. */
    mnum = 1;         /* Which factorization to use. */
    msglvl = 1;           /* Print statistical information in file */
    error = 0;            /* Initialize error flag */
/* -------------------------------------------------------------------- */
/* .. Initialize the internal solver memory pointer. This is only */
/* necessary for the FIRST call of the PARDISO solver. */
/* -------------------------------------------------------------------- */
    for ( i = 0; i < 64; i++ )
    {
        pt[i] = 0;
    }
/* -------------------------------------------------------------------- */
/* .. Reordering and Symbolic Factorization. This step also allocates */
/* all memory that is necessary for the factorization. */
/* -------------------------------------------------------------------- */
//    pardisoinit (pt,  &mtype, iparm);
/*    iparm[35] = 0;*/
    //iparm[1] = 0;
    omp_set_num_threads(num_thread);
    phase = 11;
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
    if ( error != 0 )
    {
        printf ("\nERROR during symbolic factorization: %d", error);
        exit (1);
    }
#ifdef VERBOSE
    printf ("\nReordering completed ... ");
    printf ("\nNumber of nonzeros in factors = %d", iparm[17]);
    printf ("\nNumber of factorization MFLOPS = %d", iparm[18]);
#endif
/* -------------------------------------------------------------------- */
/* .. Numerical factorization. */
/* -------------------------------------------------------------------- */
#ifdef COLDCACHE
    enableColdCache(1200,wasteFile);
#endif
    phase = 22;
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
    if ( error != 0 )
    {
        printf ("\nERROR during numerical factorization: %d", error);
        exit (2);
    }
#ifdef VERBOSE
    printf ("\nFactorization completed ... ");
#endif
/* -------------------------------------------------------------------- */
/* .. Back substitution and iterative refinement. */
/* -------------------------------------------------------------------- */
    phase = 33;
    iparm[7] = 2;         /* Max numbers of iterative refinement steps. */
    /* Set right hand side to one. */
    for ( i = 0; i < n; i++ )
    {
        b[i] = 1;
    }
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, b, x, &error);
    if ( error != 0 )
    {
        printf ("\nERROR during solution: %d", error);
        exit (3);
    }
#ifdef VERBOSE
    printf ("\nSolve completed ... \n");
#endif

/* -------------------------------------------------------------------- */
/* .. Termination and release of memory. */
/* -------------------------------------------------------------------- */
    phase = -1;           /* Release internal memory. */
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n, &ddum, ia, ja, &idum, &nrhs,
             iparm, &msglvl, &ddum, &ddum, &error);

    mkl_free(b);
    mkl_free(x);
    return 0;
}