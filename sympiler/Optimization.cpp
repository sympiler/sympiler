//
// Created by kazem on 17/02/16.
//

#include "Optimization.h"
#include "ScheduleCompute.h"

namespace Sympiler {
namespace Internal {
Stmt Optimization(Func &output) { 
    //FIXME: assuming it is only a function, we need to change it for future.
    //create a for-loop for each func
    Stmt s,t;
    Expr i0 = Variable::make(halide_type_t(halide_type_int,16), output.Name()+"f0");
    Expr i1 = Variable::make(halide_type_t(halide_type_int,16),output.Name()+"f1");
    Expr imm1 = IntImm::make(halide_type_t(halide_type_int,16),1);
    Expr Lpi0 = Pointer::make(halide_type_t(halide_type_int,32),"Lp",i0); //Lp[i0]
    Expr Lip = Pointer::make(halide_type_t(halide_type_int,32),"Li",i1); //Li[i1]
    Expr Lxp = Pointer::make(halide_type_t(halide_type_float,64),"Lx",i1); //Lx[i1]
    Expr i2 = Add::make(i0 , imm1);
    Expr Lpi0_1 = Pointer::make(halide_type_t(halide_type_int,32),"Lp",i2);//Lp[i0+1]
    Expr rhs = Pointer::make(halide_type_t(halide_type_float,64),"x",i0); //x[i0]
    Expr xLiP = Pointer::make(halide_type_t(halide_type_int,32),"x",Lip); //x[Li[i1]]
    //Expr diag = Variable::make("diag");
    Expr diag=Div::make(Lpi0,rhs);
    //Load::make(halide_type_t(halide_type_float,32),);

    //Expr diag = Div::make(Lpi0,rhsDense);
    //xLiP =Sub::make(xLiP, Mul::make(Lxp,Variable::make("diag")));
    Expr imm0 = IntImm::make(halide_type_t(halide_type_int,16),0);
    Expr n = Variable::make(halide_type_t(halide_type_int,16),"n");

    //s=Evaluate::make(Let::make("x[Li[p]]", Sub::make(xLiP, Mul::make(Lxp,diag)),0));
    s=Store::make(xLiP, "st1", Sub::make(xLiP, Mul::make(Lxp,rhs)));
    s=For::make(output.Name()+"f1",Lpi0,Lpi0_1,ForType::Serial,DeviceAPI::Host, s);
    //s = LetStmt::make("diag",diag,s);
    t=Store::make(rhs,"st2" ,diag);
    s=Block::make(t,s);

    s=For::make(output.Name()+"f0",imm0,n,ForType::Serial,DeviceAPI::Host, s);
    //ScheduleCompute sc(s,f);
    //f.optimized =sc.mutate(s);

    return s;
}
}
}