//
// Created by kazem on 28/01/16.
//

#ifndef DSLPROJECT_VAR_H
#define DSLPROJECT_VAR_H


#include <string>
#include "Expr.h"
#include "IR.h"
namespace Sympiler {
namespace Internal {
class Var {
    std::string _name;
    Type _t;

public:
    Var();

    Var(std::string);

    const std::string &name() const { return _name; }

    operator Expr() const {
        return Variable::make(_t,name());
    }

};
}
}

#endif //DSLPROJECT_VAR_H
