#include <iostream>
#include <fstream>
#include <chrono>
#include <cholmod.h>
#include <cholmod_function.h>
#include <algorithm>
#include <iomanip>

#ifdef COLDCACHE
#include "../../util/Util.h"
#include "../../symTest/cholesky/def.h"

#endif

#define CPUTIME (SuiteSparse_time ( ))
#define TRUE 1
using namespace std;



int main(int argc, char *argv[]) {
#ifdef COLDCACHE
 if(argc<3){
  cout<<"Two arguments are needed, input matrix path and "
    "the output path for the generated triangular matrix \n";
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

 std::string f1 = argv[1];
 std::string f2 = argv[2];
 std::string f1Rev(f1);
 std::reverse(f1Rev.begin(),f1Rev.end());
 std::size_t pos = f1Rev.find("/");
 std::string inputName=f1Rev.substr(4,pos-4);
 std::reverse(inputName.begin(),inputName.end());
 inputName+="_trns.mtx";
 f2+=inputName;
 std::ofstream ofs;
 ofs.open (f2, std::ofstream::out );

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

 cm->final_asis = FALSE; // Keep the output blocked
 cm->final_super=FALSE;
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
 /* real case */
 for (int i = 0 ; i < n ; i++)
 {
  double x = n ;
  Bx [i] = 1 + i / x ;
 }

#endif

#ifndef NMATRIXOPS
 cholmod_norm_dense (B, 0, cm) ;	/* max norm */
#endif

 /* ---------------------------------------------------------------------- */
 /* analyze and factorize */
 /* ---------------------------------------------------------------------- */

 L1 = cholmod_analyze (Ac1, cm) ;
 cholmod_factorize (Ac1, L1, cm) ;

 // cholmod_free_dense (&X, cm) ;
 ofs<<"%%MatrixMarket matrix coordinate real symmetric\n"
   "%-------------------------------------------------------------------------------\n"
   "% Generated to test Sympiler \n"
   "% author: Kazem Cheshmi\n"
   "%-------------------------------------------------------------------------------\n";

 int *Lpi = static_cast<int*> (L1->i);
 int *Lp = static_cast<int*> (L1->p);
 double *Lx = static_cast<double*> (L1->x);
 long NNZ = Lp[L1->n];
 ofs<<L1->n<<" "<<L1->n << " "<<NNZ<<"\n";
 for (unsigned j = 0; j < L1->n; ++j) {
  for (int i = Lp[j]; i < Lp[j+1]; ++i) {
   ofs<< Lpi[i]+1 << " "<< j+1 << " "<< std::fixed << std::setprecision(12)<<Lx[i]<<"\n";
  }
 }
 cholmod_free_factor (&L1, cm) ;

 /* ---------------------------------------------------------------------- */
 /* free matrices and finish CHOLMOD */
 /* ---------------------------------------------------------------------- */
 cholmod_free_sparse (&Ac1, cm) ;
 //cholmod_free_dense (&B, cm) ;
 cholmod_finish (cm) ;

 return 0;
}