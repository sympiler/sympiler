cmake_minimum_required(VERSION 3.0)

project(sympile)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-rtti -g -fpic -Wno-dev")
if(UNIX)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")
endif()

add_executable(sympile sympiler.cpp Func.cpp Func.h Var.cpp Var.h Module.cpp Module.h
        CodeGen_C.cpp CodeGen_C.h Target.cpp Target.h Output.cpp Output.h IR.cpp IR.h Expr.cpp Expr.h
        IROperator.cpp IROperator.h Schedule.cpp Schedule.h IRPrinter.h IRPrinter.cpp IRVisitor.h IRVisitor.cpp
        IRProfiler.cpp IRProfiler.h
        Type.h Type.cpp Optimization.h Optimization.cpp IntrusivePtr.h IRMutator.h IRMutator.cpp ScheduleCompute.cpp
        ScheduleCompute.h Kernel.cpp Kernel.h Inspector.cpp Inspector.h Matrix.cpp Matrix.h Util.cpp Util.h Scope.h
        Lower.cpp Lower.h Argument.h Substitute.cpp Substitute.h IREquality.cpp IREquality.h Function.cpp Function.h
        Factorization.cpp Factorization.h NumericalUtils.h FactorizationInspector.cpp FactorizationInspector.h
        NumericalUtils.cpp Triangular.cpp Triangular.h
        VIPrune.cpp VIPrune.h
        VIPeel.cpp VIPeel.h
        VILoopFuse.cpp VILoopFuse.h PeelIters.cpp PeelIters.h)


target_link_libraries(sympile

        )
