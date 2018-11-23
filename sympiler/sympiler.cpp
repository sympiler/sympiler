/**
* Original Sympiler code will be here
*/


#include <stdio.h>
#include "Triangular.h"
#include "Factorization.h"

//#undef TRN
#undef CHOL
#define TRN

using namespace Sympiler::Internal;
using namespace Sympiler;

/// entry to our framework
int main(int argc, char *argv[]) {
    if(argc < 2){
        printf("The input matrix path is missing!\n");
        return -1;
    }
    std::string matname=argv[1];

#ifdef TRN
    //Sparse l(10000,halide_type_t(halide_type_float,64),100,2,"");
    Sparse l(Float(64),matname);
    Sparse rhs(Float(64),"../../data/cbuckle_rhs.mtx");
    //Dense rhsDense(halide_type_t(halide_type_float,64),100,1,"");
    Triangular trns(l,rhs);
    trns.sympile_to_c("../../symGen/triang");
#endif
#ifdef CHOL
    Sparse A(Float(64),matname);
    Cholesky chol(A);
    chol.sympile_to_c("../../symGen/chol");
#endif
    printf("Sympiler code generated, but not yet run.\n");
    return 0;
}
