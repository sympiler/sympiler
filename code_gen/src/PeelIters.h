//
// Created by Shujian Qian on 2021-05-16.
//

#ifndef SYMPILER_PEELITERS_H
#define SYMPILER_PEELITERS_H

#include "IRMutator.h"

namespace Sympiler{
namespace Internal {

/**
 * \brief A visitor that peels iterations of a for loop.
 */
class PeelIters : public IRMutator {
private:
    const int *m_iters_to_peel; /**< The array of iterations to peel. */
    const int m_num_iters_to_peel; /**< Number of iterations to peel. */

    /**
     * \brief Mutates a for loop with type Sympiler::Internal::ForType::Peeled.
     * \param for_loop
     */
    void visit(const For *for_loop);

public:
    /**
     * \brief A constructor.
     * \param iters_to_peel The array of iterations to peel.
     * \param num_iters_to_peel Number of iterations to peel.
     */
    PeelIters(const int *iters_to_peel, int num_iters_to_peel):
        m_iters_to_peel(iters_to_peel), m_num_iters_to_peel(num_iters_to_peel) {}

};

}
}

#endif //SYMPILER_PEELITERS_H
