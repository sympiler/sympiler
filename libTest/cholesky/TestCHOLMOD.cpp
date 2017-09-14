#include <iostream>
#include <fstream>
#include <chrono>
#include <cholmod.h>
#include <cholmod_function.h>
#ifdef COLDCACHE
#include "../../util/Util.h"
#endif

#define CPUTIME (SuiteSparse_time ( ))
#define TRUE 1
using namespace std;



int main(int argc, char *argv[]) {
#ifdef COLDCACHE
    if(argc<3){
        cout<<"Two arguments are needed, input matrix path and a "
                "large matrix path\n";
        return -1;
    }
#else
    if(argc<2){
        printf("Please enter a path for the input matrix");
        return -1;
    }
#endif
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> elapsed_seconds;
    double CHOLMODNumTime = 0 ;

    std::string f1 = argv[1];
#ifdef COLDCACHE
    std::string waste=argv[2];
    ifstream wasteFile;
    wasteFile.open(waste);
    if(wasteFile.fail()){
        cout<<"The input file for simulating cold-cache conditions does "
                "not exist\n";
        return -1;
    }
#endif


    FILE *f ;
    cholmod_sparse *Ac1 ;
    cholmod_dense *B;
//    double one [2], zero [2] ;
    cholmod_common Common, *cm ;
    cholmod_factor *L1 ;
    double *Bx ;
    int  xtype;
    int ver [3] ;

    /* ---------------------------------------------------------------------- */
    /* get the file containing the input matrix */
    /* ---------------------------------------------------------------------- */

    if ((f = fopen (argv [1], "r")) == NULL)
    {
        return -1;
    }


    /* ---------------------------------------------------------------------- */
    /* start CHOLMOD and set parameters */
    /* ---------------------------------------------------------------------- */

    cm = &Common ;
    cholmod_start (cm) ;
    CHOLMOD_FUNCTION_DEFAULTS ;     /* just for testing (not required) */

    /* use default parameter settings, except for the error handler.  This
     * demo program terminates if an error occurs (out of memory, not positive
     * definite, ...).  It makes the demo program simpler (no need to check
     * CHOLMOD error conditions).  This non-default parameter setting has no
     * effect on performance. */
   // cm->error_handler = my_handler ;

    /* Note that CHOLMOD will do a supernodal LL' or a simplicial LDL' by
     * default, automatically selecting the latter if flop/nnz(L) < 40. */

    cm->final_asis = TRUE; // Keep the output blocked
    cm->nmethods=1;
    cm->method[0].ordering = CHOLMOD_AMD ;
    cm->postorder = TRUE;
    cm->supernodal_switch = 0.0001; //Always Supernodal

    /*cm->nrelax [0] = 0 ;
    cm->nrelax [1] = 0 ;
    cm->nrelax [2] = 0 ;*/
    /* ---------------------------------------------------------------------- */
    /* create basic scalars */
    /* ---------------------------------------------------------------------- */


    /* ---------------------------------------------------------------------- */
    /* read in a matrix */
    /* ---------------------------------------------------------------------- */

    cholmod_version (ver) ;
    SuiteSparse_version (ver) ;
    Ac1 = cholmod_read_sparse (f, cm) ;

    xtype = Ac1->xtype ;
#ifndef NMATRIXOPS
     cholmod_norm_sparse (Ac1, 0, cm) ;
#endif

    if (Ac1->nrow > Ac1->ncol)
    {
        /* Transpose A so that A'A+beta*I will be factorized instead */
        cholmod_sparse *C = cholmod_transpose (Ac1, 2, cm) ;
        cholmod_free_sparse (&Ac1, cm) ;
        Ac1 = C ;
    }

    /* ---------------------------------------------------------------------- */
    /* create an arbitrary right-hand-side */
    /* ---------------------------------------------------------------------- */

    int n = Ac1->nrow ;
    B = cholmod_zeros (n, 1, xtype, cm) ;
    Bx = (double *)B->x ;

#if GHS
    {
	/* b = A*ones(n,1), used by Gould, Hu, and Scott in their experiments */
	cholmod_dense *X0 ;
	X0 = cholmod_ones (A->ncol, 1, xtype, cm) ;
	cholmod_sdmult (A, 0, one, zero, X0, B, cm) ;
	cholmod_free_dense (&X0, cm) ;
    }
#else
    if (xtype == MATRIX_REAL)
    {
        /* real case */
        for (int i = 0 ; i < n ; i++)
        {
            double x = n ;
            Bx [i] = 1 + i / x ;
        }
    }
    else
    {
        /* complex case */
        for (int i = 0 ; i < n ; i++)
        {
            double x = n ;
            Bx [2*i  ] = 1 + i / x ;		/* real part of B(i) */
            Bx [2*i+1] = (x/2 - i) / (3*x) ;	/* imag part of B(i) */
        }
    }
#endif

#ifndef NMATRIXOPS
     cholmod_norm_dense (B, 0, cm) ;	/* max norm */
#endif

    /* ---------------------------------------------------------------------- */
    /* analyze and factorize */
    /* ---------------------------------------------------------------------- */

    L1 = cholmod_analyze (Ac1, cm) ;
#ifdef COLDCACHE
    enableColdCache(1200,wasteFile);
#endif
    start = std::chrono::system_clock::now();
    cholmod_factorize (Ac1, L1, cm) ;
    end = std::chrono::system_clock::now();
    elapsed_seconds = end-start;
    CHOLMODNumTime=elapsed_seconds.count();
    long NNZL = L1->xsize;
    //cholmod_print_factor (L1, "L", cm) ;
    cholmod_free_factor (&L1, cm) ;
   // cholmod_free_dense (&X, cm) ;

    /* ---------------------------------------------------------------------- */
    /* free matrices and finish CHOLMOD */
    /* ---------------------------------------------------------------------- */
    cholmod_free_sparse (&Ac1, cm) ;
    //cholmod_free_dense (&B, cm) ;
    cholmod_finish (cm) ;
    cout<<f1<<","<<CHOLMODNumTime <<","<< NNZL <<",";

    return 0;
}