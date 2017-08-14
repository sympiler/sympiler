/**
 * Describing Func as a basic block of the language.
 */

#ifndef DSLPROJECT_FUNC_H
#define DSLPROJECT_FUNC_H


#include "Var.h"
#include "Target.h"
#include "Expr.h"
namespace Sympiler {
namespace Internal {
//Declarations, not relevant to this file. TODO: This classed is not used
extern Expr operator+(Expr a, Expr b);

class Func { //TODO: More features needs to be added as we understand the domain.
   // std::vector<Param> params;
    std::string name;

public:
    Func();

    //Func(Param p1);

    std::string Name(){ return name;};
    void compile_to_c(std::string fName);

    void compile_to_c(std::string fName, Target t);

    void operator=(Expr e);

    std::vector<Expr> exprs;
    Stmt optimized;

    template<class Var, class... Args>
    void operator()(Var t, Args... args);
    //void comp(T var, Args... args);

};

}
}
#endif //DSLPROJECT_FUNC_H
