#include <iostream>
#include <fstream>
#include <chrono>
#include <chrono>
#include <cholmod.h>
#include <cholmod_function.h>
#include "Ordering.h"
#include "Inspection_Prune.h"
#include "Inspection_Block.h"
#include "Util.h"
#include "PB_Cholesky.h"
#include "LSparsity.h"
#include "chol_gen.h"
#define CPUTIME (SuiteSparse_time ( ))
#undef DEBUG
using namespace std;
int main(int argc, char *argv[]) {

    if(argc<2){
        printf("Please enter a path for the input matrix");
        return -1;
    }

    std::string f1 = argv[1];
    int *col, *row;
    int *colL, *rowL;
    double *valL;
    double  *y, *val, *x;
    int n, NNZ;
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> elapsed_seconds;
    double duration4 = 0 ,duration3 = 0, duration2=0, duration1=0;
    if (!readMatrix(f1,n,NNZ,col,row,val))
        return -1;
    std::string waste="/home/kazem/UFDB/rajat21.mtx";
    ifstream wasteFile;
    wasteFile.open(waste);
    if(wasteFile.fail())


    int factorSize=0;
    double timing[4];//for time measurement

    ifstream spFile1;
//    spFile1.open(fName1);
    int *prunePtr, *pruneSet;
    double timingChol[4]={.0,.0,.0,.0};//for time measurement
    int nrelax[3] = {4,16,48};//TODO
    double zrelax[3] = {0.8,0.1,0.05};
    int status=0;
    double *contribs;
    int super_max = 164; //tunig parameter for the max size of supernodes TODO: find it in analysis
    int col_max = n;
    int *col2sup=new int[n]();
    //int *blockSet;
    contribs = new double[super_max * col_max]();
    int *li_ptr = new int[n+1];
    int *map = new int[n]();
    colL = new int[n + 1]();
    CSC *Amat = new CSC;
    Amat->nzmax = NNZ; Amat->ncol=Amat->nrow=n;
    Amat->stype=-1;Amat->xtype=CHOLMOD_REAL;Amat->packed=TRUE;
    Amat->p = col; Amat->i = row; Amat->x=val; Amat->nz = NULL;
    Amat->sorted = TRUE;
    ///nrelax[0]=0;nrelax[1]=0;nrelax[2]=0;
    BCSC *L = analyze_p2(1,Amat,NULL,NULL,nrelax,zrelax,
                                    n,prunePtr,pruneSet, status);
#if 0
    cout<<"\n";
    for (int j = 0; j < L->nsuper; ++j) {
        int colLength = L->pi[j+1]-L->pi[j];
        int supWid = L->super[j+1]-L->super[j];
        cout<<"Supernode: "<<j<<" Len and Wid are: "
            <<colLength<<","<<supWid<<"\n";
        for (int i = L->pi[j]; i < L->pi[j+1]; ++i) {
            cout<<L->s[i]<<",";
        }
        cout<<"\n";
    }
#endif
    //Some conversion for sympiler

    int colLength=0;

    for (int j = 0; j < L->nsuper; ++j) {
        int curCol = L->super[j];
        int nxtCol = L->super[j+1];
        colLength = L->pi[j+1]-L->pi[j];
        //int supWid = L->super[j+1]-L->super[j];
        for (int i = curCol+1; i < nxtCol+1; ++i) {
            li_ptr[i-1] = L->pi[j];
            colL[i]= colL[curCol] + (i-curCol)*colLength;
        }
    }
    li_ptr[n] = L->pi[L->nsuper];
    colL[n] = colL[n-1] + colLength;
    valL = new double[L->xsize]();
#if 0

    for (int j = 0; j < 10; ++j) {
        int curCol = L->super[j];
        int nxtCol = L->super[j+1];
        for (int k = li_ptr[curCol]; k < li_ptr[nxtCol]; ++k) {
            cout<<L->s[k]<<",";

        }
        cout<<"\n";
    }
#endif

    CSC *A1 = ptranspose(Amat,2,L->Perm,NULL,0,status);
    CSC *A2 = ptranspose(A1,2,NULL,NULL,0,status);
   // enableColdCache(1200,wasteFile);
    /*double *tempVec = new double[n]();
    double *finger = new double[n*168]();*/
    start = std::chrono::system_clock::now();
    /*cholesky_left_sn_07(n,A2->p,A2->i,A2->x,colL,L->s,li_ptr,valL,
                        L->super,L->nsuper, timingChol,
                        prunePtr,pruneSet,map, contribs);*/
    Chol(n,A2->p,A2->i,A2->x,NULL,
         n,colL,L->s,valL,li_ptr,
         map,contribs,
         prunePtr,pruneSet,
         0, L->super,L->nsuper);

    end = std::chrono::system_clock::now();
    elapsed_seconds = end-start;
    duration2=elapsed_seconds.count();


    double resid [4], t, ta, tf, ts [3], tot, bnorm, xnorm, anorm, rnorm, fl,
            anz, axbnorm, rnorm2, resid2, rcond ;
    FILE *f ;
    cholmod_sparse *Ac1 ;
    cholmod_dense *X = NULL, *B, *W, *R = NULL ;
    double one [2], zero [2], minusone [2], beta [2], xlnz ;
    cholmod_common Common, *cm ;
    cholmod_factor *L1 ;
    double *Bx, *Rx, *Xx ;
    int i, isize, xsize, ordering, xtype, s, ss, lnz ;
    int trial, method, L_is_super ;
    int ver [3] ;

    ts[0] = 0.;
    ts[1] = 0.;
    ts[2] = 0.;

    /* ---------------------------------------------------------------------- */
    /* get the file containing the input matrix */
    /* ---------------------------------------------------------------------- */

    if ((f = fopen (argv [1], "r")) == NULL)
    {

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

    cm->final_asis = TRUE;
    cm->nmethods=1;
    cm->method[0].ordering = CHOLMOD_AMD ;
    cm->postorder = TRUE;
    cm->supernodal_switch = 0.0001;

    /*cm->nrelax [0] = 0 ;
    cm->nrelax [1] = 0 ;
    cm->nrelax [2] = 0 ;*/
    /* ---------------------------------------------------------------------- */
    /* create basic scalars */
    /* ---------------------------------------------------------------------- */

    zero [0] = 0 ;
    zero [1] = 0 ;
    one [0] = 1 ;
    one [1] = 0 ;
    minusone [0] = -1 ;
    minusone [1] = 0 ;
    beta [0] = 1e-6 ;
    beta [1] = 0 ;

    /* ---------------------------------------------------------------------- */
    /* read in a matrix */
    /* ---------------------------------------------------------------------- */

//    printf ("\n---------------------------------- cholmod_demo:\n") ;
    cholmod_version (ver) ;
//    printf ("cholmod version %d.%d.%d\n", ver [0], ver [1], ver [2]) ;
    SuiteSparse_version (ver) ;
//    printf ("SuiteSparse version %d.%d.%d\n", ver [0], ver [1], ver [2]) ;
    Ac1 = cholmod_read_sparse (f, cm) ;

    xtype = Ac1->xtype ;
    anorm = 1 ;
#ifndef NMATRIXOPS
    anorm = cholmod_norm_sparse (Ac1, 0, cm) ;
//    printf ("norm (A,inf) = %g\n", anorm) ;
//    printf ("norm (A,1)   = %g\n", cholmod_norm_sparse (Ac1, 1, cm)) ;
#endif
//    cholmod_print_sparse (Ac1, "A", cm) ;

    if (Ac1->nrow > Ac1->ncol)
    {
        /* Transpose A so that A'A+beta*I will be factorized instead */
        cholmod_sparse *C = cholmod_transpose (Ac1, 2, cm) ;
        cholmod_free_sparse (&Ac1, cm) ;
        Ac1 = C ;
//        printf ("transposing input matrix\n") ;
    }

    /* ---------------------------------------------------------------------- */
    /* create an arbitrary right-hand-side */
    /* ---------------------------------------------------------------------- */

    n = Ac1->nrow ;
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
    if (xtype == CHOLMOD_REAL)
    {
        /* real case */
        for (i = 0 ; i < n ; i++)
        {
            double x = n ;
            Bx [i] = 1 + i / x ;
        }
    }
    else
    {
        /* complex case */
        for (i = 0 ; i < n ; i++)
        {
            double x = n ;
            Bx [2*i  ] = 1 + i / x ;		/* real part of B(i) */
            Bx [2*i+1] = (x/2 - i) / (3*x) ;	/* imag part of B(i) */
        }
    }
#endif

//    cholmod_print_dense (B, "B", cm) ;
    bnorm = 1 ;
#ifndef NMATRIXOPS
    bnorm = cholmod_norm_dense (B, 0, cm) ;	/* max norm */
//    printf ("bnorm %g\n", bnorm) ;
#endif

    /* ---------------------------------------------------------------------- */
    /* analyze and factorize */
    /* ---------------------------------------------------------------------- */

    //t = CPUTIME ;
    L1 = cholmod_analyze (Ac1, cm) ;
    //ta = CPUTIME - t ;
    ta = MAX (ta, 0) ;

    //printf ("Analyze: flop %g lnz %g\n", cm->fl, cm->lnz) ;

    if (Ac1->stype == 0)
    {
    //    printf ("Factorizing A*A'+beta*I\n") ;
        //t = CPUTIME ;
        cholmod_factorize_p (Ac1, beta, NULL, 0, L1, cm) ;
        //tf = CPUTIME - t ;
        tf = MAX (tf, 0) ;
    }
    else
    {
        enableColdCache(1200,wasteFile);
        //start = std::chrono::system_clock::now();
        t = CPUTIME ;
        cholmod_factorize (Ac1, L1, cm) ;
        tf = CPUTIME - t ;
        tf = MAX (tf, 0) ;
        /*end = std::chrono::system_clock::now();
        elapsed_seconds = end-start;
        duration3=elapsed_seconds.count();*/
    }

    int *Lpx = static_cast<int*> (L1->px);
#ifdef VERIFY
    int *Lpi = static_cast<int*> (L1->pi);
    int *Lsuper = static_cast<int*> (L1->super);
    int *Ls = static_cast<int*> (L1->s);
    double *Lx = static_cast<double*> (L1->x);
    int *LPerm = static_cast<int*> (L1->Perm);
    int cnt=0;
   // cout<<"\n";
    ASSERT(L1->nsuper == L->nsuper);
    for (int j = 0; j < L1->nsuper; ++j) {
        int colLength = Lpi[j+1]-Lpi[j];
        int supWid = Lsuper[j+1]-Lsuper[j];
        //cout<<"Supernode: "<<j<<" Len and Wid are: "
        //    <<colLength<<","<<supWid<<"\n";
//        for (int k = Lpi[j]; k < Lpi[j+1]; ++k) {
//            cout<<Ls[k]<<",";
//        }
//        cout<<"\n";
        for (int i = Lpx[j]; i < Lpx[j+1]; ++i) {
            if(Lx[i] - valL[i] > 0.001){
                cnt++;
               // cout<<"SN: "<<j<<";";
            }
            //cout<<Lx[i]<<",";
        }
      //  cout<<"\n";
    }
   if(cnt>0)
       cout<<"#"<<cnt<<";"<<"\n";
#endif
    L_is_super = L->is_super ;

    //cholmod_print_factor (L1, "L", cm) ;
    int NNZL = Lpx[L1->nsuper];
    cholmod_free_factor (&L1, cm) ;
    cholmod_free_dense (&X, cm) ;

    /* ---------------------------------------------------------------------- */
    /* free matrices and finish CHOLMOD */
    /* ---------------------------------------------------------------------- */

    cholmod_free_sparse (&Ac1, cm) ;
    cholmod_free_dense (&B, cm) ;
    cholmod_finish (cm) ;

    allocateAC(Amat,0,0,0,FALSE);
    allocateLC(L,FALSE);

    cout<<f1<<","<<duration2<<","<<tf <<","<< NNZL <<",";


#if DEBUG > 0
    for (int i = n-10; i < n; ++i) {
            std::cout<<i<<":\n";
            for (int m = colL[i],cnt=0; m < colL[i+1]; ++m, ++cnt) {
                if(!std::isfinite(valL[m])) {
                    std::cout << "Error in col "<< i;
                    return -1;
                }
                if(rowL[li_ptr[i]+cnt] >= i )
                    std::cout<<valL[m]<<",";
            }
            std::cout<<"\n";
        }
        std::cout<<"\n";
#endif
    delete []col2sup;
    delete []prunePtr; delete []pruneSet;
    delete []contribs;
    delete []map;
    delete []valL;
    delete []colL;
    delete []li_ptr;
    return 0;
}