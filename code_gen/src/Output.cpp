
#include <fstream>
#include "Output.h"
#include "CodeGen_C.h"
namespace Sympiler {
namespace Internal {

void Output::compile_to_source_c(Module& module, std::ostream &stream) {
    //std::ofstream file(fName + "_gen.cpp");
    CodeGen_C cg_c(stream,false);
    cg_c.compile(module);
    //file.close();
}

void Output::compile_to_source_c(Module& module, std::string name, bool isBlock){
    //Printing the header file, skipped this for now
    /*std::ofstream fileHeader(name + ".h");
    CodeGen_C cg_h(fileHeader,true);
    cg_h.compile(module);*/
    //Printing the body file
    std::ofstream fileBody(name + ".h");
    CodeGen_C cg_b(fileBody,false, isBlock);
    cg_b.compile(module);
}

}
}