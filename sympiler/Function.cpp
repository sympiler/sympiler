//
// Created by kazem on 5/15/17.
//

#include "Function.h"
namespace Sympiler{
namespace Internal{

struct FunctionContents {
    //mutable RefCount ref_count;
    std::string name;
    std::vector<std::string> args;
    std::vector<Expr> values;
    std::vector<Type> output_types;
    //Schedule schedule;

    //std::vector<UpdateDefinition> updates;

    std::string debug_file;

    //std::vector<Parameter> output_buffers;

    //std::vector<ExternFuncArgument> extern_arguments;
    std::string extern_function_name;

    bool trace_loads, trace_stores, trace_realizations;

    bool frozen;

    FunctionContents() : trace_loads(false), trace_stores(false), trace_realizations(false), frozen(false) { }
};
}
}
