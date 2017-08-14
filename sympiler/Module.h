//
// Created by kazem on 01/02/16.
//

#ifndef DSLPROJECT_MODULE_H
#define DSLPROJECT_MODULE_H


#include <vector>
#include "Target.h"
#include "Kernel.h"

namespace Sympiler {
namespace Internal {
class Module {
    Target target;
    std::string mName;

public:
    std::vector<Kernel> kers;

    Module();

    Module(std::string fName, Target target);

    void append(Kernel ker);

};
}
}

#endif //DSLPROJECT_MODULE_H
