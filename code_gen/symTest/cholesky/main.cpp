#include <iostream>
#include <fstream>
#include <chrono>
#include <chrono>
#include <cholmod.h>
#include <cholmod_function.h>
#include "Inspection_Prune.h"
#include "PB_Cholesky.h"
#include "LSparsity.h"
#include "chol_gen.h"
#include "../../util/Util.h"

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

    std::string f1 = argv[1];
    int *col, *row;
    int *colL;
    double *valL;
    double  *val;
    int n, NNZ;
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> elapsed_seconds;
    double SympilerTime=0;
    if (!readMatrix(f1,n,NNZ,col,row,val))
        return -1;
#ifdef COLDCACHE
    std::string waste=argv[2];
    ifstream wasteFile;
    wasteFile.open(waste);
    if(wasteFile.fail()){
        std::cout<<"The input file for simulating cold-cache conditions does "
                "not exist\n";
    return -1;
    }
#endif

    int *prunePtr, *pruneSet;
    double timingChol[4]={.0,.0,.0,.0};//for time measurement
    int nrelax[3] = {4,16,48};//TODO
    double zrelax[3] = {0.8,0.1,0.05};
    int status=0;
    double *contribs;
    int super_max = 164; // TODO: find it in analysis phase
    int col_max = n;
    int *col2sup=new int[n]();
    //int *blockSet;
    contribs = new double[super_max * col_max]();
    int *li_ptr = new int[n+1];
    int *map = new int[n]();
    colL = new int[n + 1]();
    CSC *Amat = new CSC;
    Amat->nzmax = NNZ; Amat->ncol=Amat->nrow=n;
    Amat->stype=-1;Amat->xtype=MATRIX_REAL;Amat->packed=TRUE;
    Amat->p = col; Amat->i = row; Amat->x=val; Amat->nz = NULL;
    Amat->sorted = TRUE;
    BCSC *L = buildLFactor(Amat, NULL, nrelax, zrelax, n, prunePtr,
                           pruneSet, status);
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
    //Some conversion for src

    int colLength=0;

    for (unsigned j = 0; j < L->nsuper; ++j) {
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


    CSC *A1 = transposeMatrix(Amat, 2, L->Perm, NULL, 0);
    CSC *A2 = transposeMatrix(A1, 2, NULL, NULL, 0);
#ifdef COLDCACHE
    enableColdCache(1200,wasteFile);
#endif
    start = std::chrono::system_clock::now();
/*    cholesky_left_sn_07(n,A2->p,A2->i,A2->x,colL,L->s,li_ptr,valL,
                        L->super,L->nsuper, timingChol,
                        prunePtr,pruneSet,map, contribs);*/
    Chol(n,A2->p,A2->i,A2->x,NULL,
         n,colL,L->s,valL,li_ptr,
         map,contribs,
         prunePtr,pruneSet,
         0, L->super,L->nsuper);

    end = std::chrono::system_clock::now();
    elapsed_seconds = end-start;
    SympilerTime=elapsed_seconds.count();


#ifdef VERIFY
    FILE *f ;
    cholmod_sparse *Ac1 ;
    cholmod_dense *X = NULL, *B;
    cholmod_common Common, *cm ;
    cholmod_factor *L1 ;
    double *Bx ;
    int i,  xtype ;
    int ver [3] ;


    /* ---------------------------------------------------------------------- */
    /* get the file containing the input matrix */
    /* ---------------------------------------------------------------------- */

    if ((f = fopen (argv [1], "r")) == NULL)
    {
        printf("The input file is not valid.\n");
    }


    /* ---------------------------------------------------------------------- */
    /* start CHOLMOD and set parameters */
    /* ---------------------------------------------------------------------- */
    cm = &Common ;
    cholmod_start (cm) ;
    CHOLMOD_FUNCTION_DEFAULTS ;     /* just for testing (not required) */

    cm->final_asis = TRUE;
    cm->nmethods=1;
    cm->method[0].ordering = CHOLMOD_AMD ;
    cm->postorder = TRUE;
    cm->supernodal_switch = 0.0001;

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
    if (xtype == MATRIX_REAL)
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

#ifndef NMATRIXOPS
    cholmod_norm_dense (B, 0, cm) ;	/* max norm */
#endif

    L1 = cholmod_analyze (Ac1, cm) ;
    //start = std::chrono::system_clock::now();
    cholmod_factorize (Ac1, L1, cm) ;
    /*end = std::chrono::system_clock::now();
    elapsed_seconds = end-start;
    duration3=elapsed_seconds.count();*/

    int *Lpx = static_cast<int*> (L1->px);
    double *Lx = static_cast<double*> (L1->x);
    int cnt=0;
    ASSERT(L1->nsuper == L->nsuper);
    for (unsigned j = 0; j < L1->nsuper; ++j) {
        for (int i = Lpx[j]; i < Lpx[j+1]; ++i) {
            if(Lx[i] - valL[i] > 0.001){
                cnt++;
            }
        }
    }
   if(cnt>0)
       cout<<"#"<<cnt<<";"<<"\n";

    int NNZL = Lpx[L1->nsuper];
    cholmod_free_factor (&L1, cm) ;
    cholmod_free_dense (&X, cm) ;

    cholmod_free_sparse (&Ac1, cm) ;
    cholmod_free_dense (&B, cm) ;
    cholmod_finish (cm) ;

#endif
    allocateAC(Amat,0,0,0,FALSE);
    allocateLC(L,FALSE);
    cout<<f1<<","<<SympilerTime<<","<< NNZL <<",";

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

