/**
 * Stand-alone functions managing output related tasks.
 */

#ifndef DSLPROJECT_OUTPUT_H
#define DSLPROJECT_OUTPUT_H


#include <string>

#include "Module.h"

namespace Sympiler {
namespace Internal {
class Output {
public:
    static void compile_to_source_c(Module& module, std::ostream &stream);
    static void compile_to_source_c(Module& module, std::string name, bool isBlock);

};
}
}


#endif //DSLPROJECT_OUTPUT_H
