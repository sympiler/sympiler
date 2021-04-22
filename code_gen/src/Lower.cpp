//
// Created by kazem on 5/2/17.
//
#include "Lower.h"

namespace Sympiler {
namespace Internal {

Stmt lower(Kernel *ker, SymbolicObject *sym){
    Stmt s;
    //TODO in future, we may decide whether used a blocked format or not based on the symbolic information for demo:
    // if(sym->...) ker->block=true
    s = ker->baseCode();
    if(sym->IsVSBlock()){
        ker->isVSBlock= true;//baseline code will be a blocked code
        s = ker->VSBlockIG(sym);
    }
    if(sym->IsVIPrune()){
        ker->isVIPrune = true;
        s = ker->VIPruneIG(sym);
    }
    return s;
}

}
}