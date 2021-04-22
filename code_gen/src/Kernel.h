//
// Created by kazem on 4/18/17.
//

#ifndef SYMPILER_PROJ_KERNEL_H
#define SYMPILER_PROJ_KERNEL_H

#include <vector>
#include "Inspector.h"
#include "Matrix.h"
#include "Target.h"
#include "Argument.h"



namespace Sympiler {
namespace Internal {

class Kernel {
protected:
    bool isSetArgInserted;
    std::string name;
    Stmt loweredKer; //The triangular code after possible IG transformations.
    bool isTile, isPeel, isUnroll;
    std::vector<Expr> args;
    Expr ret;

public:
    bool isVSBlock, isVIPrune;
    std::vector<Argument> argType;
    Argument retVal;
    Kernel();
    Kernel(std::string name);
    std::string Name();
    std::vector<Expr> exprs;
    Stmt LoweredKer();
    virtual Stmt baseCode(); //Generates the baseline code for the kernel
    virtual Stmt VSBlockIG(SymbolicObject *sym);
    virtual Stmt VIPruneIG(SymbolicObject *sym);
    virtual void sympile_to_c(std::string fName, Target t);
};


}
}


#endif //SYMPILER_PROJ_KERNEL_H
