cmake_minimum_required(VERSION 2.8)

project(sympiler_proj)

set(CMAKE_CXX_FLAGS "-fno-rtti -g -fpic -Wno-dev")
if(UNIX)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")
endif()

add_executable(sympiler_proj sympiler.cpp Func.cpp Func.h Var.cpp Var.h Module.cpp Module.h
        CodeGen_C.cpp CodeGen_C.h Target.cpp Target.h Output.cpp Output.h IR.cpp IR.h Expr.cpp Expr.h
        IROperator.cpp IROperator.h Schedule.cpp Schedule.h IRPrinter.h IRPrinter.cpp IRVisitor.h IRVisitor.cpp
        Type.h Type.cpp Optimization.h Optimization.cpp IntrusivePtr.h IRMutator.h IRMutator.cpp ScheduleCompute.cpp
        ScheduleCompute.h Kernel.cpp Kernel.h Inspector.cpp Inspector.h Matrix.cpp Matrix.h Util.cpp Util.h Scope.h
        Lower.cpp Lower.h Argument.h VIPrune.cpp VIPrune.h Substitute.cpp Substitute.h IREquality.cpp IREquality.h Function.cpp Function.h Factorization.cpp Factorization.h NumericalUtils.h FactorizationInspector.cpp FactorizationInspector.h NumericalUtils.cpp Triangular.cpp Triangular.h)
#add_executable(Substitute.cpp Substitute.h) #TODO will add later

target_link_libraries(sympiler_proj

        )
