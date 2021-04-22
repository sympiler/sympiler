//
// Created by kazem on 01/02/16.
//

#include "Module.h"
namespace Sympiler {
namespace Internal {
Module::Module() {

}

Module::Module(std::string fName, Target t) {
    mName = fName;
    target = t;
}

void Module::append(Kernel func) {
    kers.push_back(func);
}
}
}