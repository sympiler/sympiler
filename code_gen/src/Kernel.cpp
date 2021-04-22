//
// Created by kazem on 4/18/17.
//

#include <fstream>
#include "Kernel.h"
#include "Output.h"
#include "Optimization.h"
#include "Lower.h"
#include "VIPrune.h"
#include "IROperator.h"

namespace Sympiler {
namespace Internal {

Kernel::Kernel():isSetArgInserted(false), isVIPrune(false),
                    isVSBlock(false){

}
Kernel::Kernel(std::string name):name(name),isSetArgInserted(false),
    isVSBlock(false), isVIPrune(false){
    //for now each kernel function returns int
    retVal = Argument(name,Argument::Kind::InputScalar,
                      halide_type_t(halide_type_int,32),0);
}

std::string Kernel::Name(){
    return name;
}


Stmt Kernel::baseCode(){
    return loweredKer;
}

Stmt Kernel::VIPruneIG(SymbolicObject *sym) {//TODO
    Stmt s;
    return s;
}

Stmt Kernel::VSBlockIG(SymbolicObject *sym) {
    Stmt s;
    return s;
}

Stmt Kernel::LoweredKer(){
    return loweredKer;
}

void Kernel::sympile_to_c(std::string fName, Target t) {

}


}
}