//
// Created by george on 2020-02-16.
//

#ifndef FUSION_VIPEEL_H
#define FUSION_VIPEEL_H

#include "IRMutator.h"
#include "Substitute.h"

namespace Sympiler {
namespace Internal {

class VIPeel : public IRMutator {
    int num_peel;

    void visit(const For *for_loop);

public:
    VIPeel(int peel):
        num_peel(peel){}

};
}
}

#endif //FUSION_VIPEEL_H
