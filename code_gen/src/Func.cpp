//
// Created by kazem on 24/01/16.
//

#include "Func.h"
#include "Module.h"
#include "Output.h"
#include "Optimization.h"
namespace Sympiler {
namespace Internal {
Func::Func() { }

/*Func::Func(Param p1) {
    params.push_back(p1);
}*/

void Func::compile_to_c(std::string fName) {
    Target t;//with default values
    compile_to_c(fName, t);
}

void Func::compile_to_c(std::string fName, Target t) {
    //convert to module
    Target target(t);
    optimized=Optimization(*this);
    Module module(fName, target);
//    module.append(*this);
    //Output::compile_to_source_c(module, fName);
}

void Func::operator=(Expr e) {
    exprs.push_back(e);
}
/*template<class Var, class... Args>
void Func::operator() (Var t, Args... args){

}*/

}
}