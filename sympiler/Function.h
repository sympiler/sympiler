//
// Created by kazem on 5/15/17.
//

#ifndef SYMPILER_PROJ_FUNCTION_H
#define SYMPILER_PROJ_FUNCTION_H

#include "Expr.h"

namespace Sympiler {
namespace Internal {

    struct FunctionContents;

class Function {

public:
/** Construct a new function with no definitions and no name. This
     * constructor only exists so that you can make vectors of
     * functions, etc.
     */
    Function();

    /** Reconstruct a Function from a FunctionContents pointer. */
    //Function(const IntrusivePtr<FunctionContents> &c) : contents(c) {}

    /** Construct a new function with the given name */
    Function(const std::string &n);

    /** Add a pure definition to this function. It may not already
     * have a definition. All the free variables in 'value' must
     * appear in the args list. 'value' must not depend on any
     * reduction domain */
    void define(const std::vector<std::string> &args, std::vector<Expr> values);

    /** Add an update definition to this function. It must already
     * have a pure definition but not an update definition, and the
     * length of args must match the length of args used in the pure
     * definition. 'value' must depend on some reduction domain, and
     * may contain variables from that domain as well as pure
     * variables. Any pure variables must also appear as Variables in
     * the args array, and they must have the same name as the pure
     * definition's argument in the same index. */
    void define_update(const std::vector<Expr> &args, std::vector<Expr> values);

};
}
}


#endif //SYMPILER_PROJ_FUNCTION_H
