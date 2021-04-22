//
// Created by kazem on 8/8/17.
//

#include "Triangular.h"
#include "IROperator.h"
#include "VIPrune.h"
#include "Lower.h"
#include "Module.h"
#include "Output.h"

namespace Sympiler {
namespace Internal {
Triangular::Triangular() : Kernel() {
    loweredKer = fwdSolve();
}

Triangular::Triangular(Matrix &Lmat, Matrix &rhsmat) : Kernel("trns") {
    L = &Lmat;
    rhs = &rhsmat;
    //argType.insert(argType.end(),Lmat.getDecl().begin(),Lmat.getDecl().end());
    std::vector<Expr> el;
    std::vector<Argument> al;
    Lmat.getDecl(el, al);
    args.insert(args.end(), el.begin(), el.end());
    argType.insert(argType.end(), al.begin(), al.end());

    el.clear();
    al.clear();
    rhsmat.getDecl(el, al);
    args.insert(args.end(), el.begin(), el.end());
    argType.insert(argType.end(), al.begin(), al.end());



}

Triangular::~Triangular() {
    delete rhsDense;
}

Stmt Triangular::copy2Dense() {//copy rhs to rhsDense
    Stmt s;
    //FIXME For now assuming zero, will be changed for multiple RHS
    rhsCol = rhs->accessCol(make_zero(Int(32)));
    Expr upperBound = rhs->accessCol(
            Add::make(rhsCol, make_one(Int(32))));
    Expr loopIdx = Variable::make(Int(32), "copy");
    s = Store::make(rhsDense->accessNNZ(rhs->accessRowIdx(loopIdx)),
                    "copyToDense", rhs->accessNNZ(loopIdx));
    s = For::make("copy", rhsCol, upperBound, ForType::Serial,
                  DeviceAPI::Host, s);
    return s;
}

Stmt Triangular::fwdSolve() {
    /*
     * builds the forward solve code using Sympiler IR.
     */
    Stmt s, t;

    Expr i0 = Variable::make(halide_type_t(halide_type_int, 16),
                             Kernel::name + "f0");
    Expr i1 = Variable::make(halide_type_t(halide_type_int, 16),
                             Kernel::name + "f1");
    Expr imm1 = IntImm::make(halide_type_t(halide_type_int, 16), 1);
    Expr Lxp = L->accessNNZ(i1); //Lx[i1]
    Expr i2 = Add::make(i0, imm1);
    Expr rhsDiag = rhsDense->diagonal(i0); //x[i0]
    Expr xLiP = rhsDense->accessNNZ(L->accessRowIdx(i1)); //x[Li[i1]]
    Expr diag = Div::make( rhsDense->diagonal(i0),L->diagonal(i0));
    Expr imm0 = IntImm::make(halide_type_t(halide_type_int, 16), 0);

   // s = Block::make(Free::make(rhsDense->Name()),s);
    s = Store::make(xLiP, "st1",
                    Sub::make(xLiP, Mul::make(Lxp, rhsDiag)));
    s = For::make(Kernel::name + "f1", Add::make(L->accessCol(i0),make_one(Int(16))),
                  L->accessCol(i2),
                  ForType::Serial, DeviceAPI::Host, s);
    t = Store::make(rhsDiag, "st2", diag);
    s = Block::make(t, s);
    s = For::make(Kernel::name + "f0", imm0, L->Order(),
                  ForType::Pruned,
                  DeviceAPI::Host, s);
    s = Block::make(copy2Dense(), s);
    //s = Allocate::make(rhsDense->Name(),Float(64),{rhsDense->Order()},const_true(),s);
    /*Stmt block;
    t = Store::make(rhsDiag,"st2",diag);
    block=Block::make(block,t);
    t =s;
    block = Block::make(block,t);
    s = For::make(Kernel::name+"f0",imm0,L->Order(),ForType::Pruned,
                  DeviceAPI::Host, t);*/
    loweredKer = s;
    return s;
}

Stmt Triangular::baseCode() {
    return fwdSolve();
}

Stmt Triangular::VIPruneIG(SymbolicObject *sym) {
    //VIPrune vi_prune(sym->setPtr,sym->setVal,sym->setSize);
    VIPrune vi_prune(make_zero(halide_type_t(halide_type_int, 32)),
                     sym->setVal,
                     Variable::make(halide_type_t(halide_type_int, 32),
                                    sym->setSize));
    //Adding the arguments for the inspection set
    if (!isSetArgInserted) {
        argType.push_back(
                Argument(sym->setPtr, Argument::Kind::InputBuffer,
                         halide_type_t(halide_type_int, 32), 1));
        argType.push_back(
                Argument(sym->setVal, Argument::Kind::InputBuffer,
                         halide_type_t(halide_type_int, 32), 1));
        argType.push_back(
                Argument(sym->setSize, Argument::Kind::InputScalar,
                         halide_type_t(halide_type_int, 32), 1));
        isSetArgInserted = true;//We need only one set, other
        //Let's keep the prototype fixed for easy use
        argType.push_back(
                Argument(sym->blk2Col, Argument::Kind::InputBuffer,
                         halide_type_t(halide_type_int, 32), 1));
        argType.push_back(Argument(sym->blkNo, Argument::Kind::InputScalar,
                                   halide_type_t(halide_type_int, 32), 1));
    }
    loweredKer = vi_prune.mutate(loweredKer);
    return loweredKer;
}

Stmt Triangular::VSBlockIG(SymbolicObject *sym) {
    //Adding the arguments for the inspection set
    if (!isSetArgInserted) {
        argType.push_back(
                Argument(sym->setPtr, Argument::Kind::InputBuffer,
                         halide_type_t(halide_type_int, 32), 1));
        argType.push_back(
                Argument(sym->setVal, Argument::Kind::InputBuffer,
                         halide_type_t(halide_type_int, 32), 1));
        argType.push_back(
                Argument(sym->setSize, Argument::Kind::InputScalar,
                         halide_type_t(halide_type_int, 32), 1));
        isSetArgInserted = true;//We need only one set, other
        //Adding arguments for block to column mapping
        argType.push_back(
                Argument(sym->blk2Col, Argument::Kind::InputBuffer,
                         halide_type_t(halide_type_int, 32), 1));
        argType.push_back(Argument(sym->blkNo, Argument::Kind::InputScalar,
                                   halide_type_t(halide_type_int, 32), 1));
    }
    //Generating code for the blocked triangular solver
    Stmt s;
    std::string tmpVec = "tempVec";
    Expr retVal = Variable::make(halide_type_t(halide_type_int, 32),
                                 "retval");
    Expr zro = make_zero(halide_type_t(halide_type_int, 32));
    Expr one = make_one(halide_type_t(halide_type_int, 32));
    Expr i0 = Variable::make(halide_type_t(halide_type_int, 16),
                             Kernel::name + "f0");
    Expr i1 = Variable::make(halide_type_t(halide_type_int, 16),
                             Kernel::name + "f1");

    Expr suptoCol = Pointer::make(halide_type_t(halide_type_int, 32),
                                  sym->blk2Col, i0);
    Expr nxtCol = Pointer::make(halide_type_t(halide_type_int, 32),
                                sym->blk2Col, i0);
    Expr curColtmp = Pointer::make(halide_type_t(halide_type_int, 32),
                                   sym->blk2Col,
                                   Sub::make(i0, one));
    Expr curCol = Select::make(NE::make(i0, zro), curColtmp, zro);
    Expr supWdt = Sub::make(nxtCol, curCol);
    Expr nSupR = Sub::make(L->accessRowIdxPntr(nxtCol),
                           L->accessRowIdxPntr(curCol));
    Expr LpCur = L->accessCol(curCol);

    Expr initFor = Add::make(L->accessRowIdxPntr(curCol), supWdt);
    Expr k = Sub::make(i1, initFor);

    Expr tempVec = Pointer::make(halide_type_t(halide_type_float, 64),
                                 tmpVec, k);

    Expr lii1 = L->accessRowIdx(i1);
    Expr subtract = Sub::make(rhsDense->accessNNZ(lii1), tempVec);
    Stmt cpyBack = Store::make(rhsDense->accessNNZ(lii1), "copyBack",
                               subtract);
    Stmt resetTemp = Store::make(tempVec,"resetTmp",make_zero(Int(16)));
    //Stmt loopBody = Block::make(resetTemp,loopBody);
    //loopBody = Block::make(cpyBack, loopBody);
    s = Block::make(resetTemp,s);
    s = Block::make(cpyBack,s);
    s = For::make(Kernel::name + "f1", initFor,
                  L->accessRowIdxPntr(nxtCol),
                  ForType::Serial, DeviceAPI::Host, s);

    Expr newIdx = Add::make(LpCur, supWdt);
    std::vector<Expr> matvecArgs;
    matvecArgs.push_back(nSupR);
    matvecArgs.push_back(Sub::make(nSupR, supWdt));
    matvecArgs.push_back(supWdt);
    matvecArgs.push_back(L->accessNNZ(newIdx));
    matvecArgs.push_back(rhsDense->accessNNZ(curCol));
    Expr tempVec1 = Pointer::make(halide_type_t(halide_type_float, 64),
                            tmpVec, zro);
    matvecArgs.push_back(tempVec1);
    std::vector<Argument> matvecArgums;
    matvecArgums.push_back(Argument(false));matvecArgums.push_back(Argument(false));
    matvecArgums.push_back(Argument(false));matvecArgums.push_back(Argument(true));
    matvecArgums.push_back(Argument(true));matvecArgums.push_back(Argument(true));

    Stmt matvec = CallX::make("matvec", matvecArgs, matvecArgums);
    s = Block::make(matvec, s);

    std::vector<Expr> trnsArgs;//adding input arguments of trsm
    trnsArgs.push_back(nSupR);
    trnsArgs.push_back(supWdt);
    trnsArgs.push_back(L->accessNNZ(LpCur));
    trnsArgs.push_back(rhsDense->accessNNZ(curCol));

    std::vector<Argument> trnsArgums;
    trnsArgums.emplace_back(Argument(false));trnsArgums.emplace_back(Argument(false));
    trnsArgums.emplace_back(Argument(true));trnsArgums.emplace_back(Argument(true));

    Stmt trsmCall = CallX::make("trsm_blas", trnsArgs, trnsArgums); //TODO use Let
    //s = LetStmtPtr::make(retVal,trsmCall,s);
    s = Block::make(trsmCall, s);

    s = For::make(Kernel::name + "f0", zro, sym->getBlockNo(),
                  ForType::Pruned, DeviceAPI::Host, s);
    s = Block::make(copy2Dense(), s);
    s = Allocate::make(tmpVec, halide_type_t(halide_type_float, 64),
                       {Mul::make(L->Order(),make_const(Int(16),164))}, const_true(), s);
    loweredKer = s;
    return loweredKer;
}

void Triangular::sympile_to_c(std::string fName, Target t) {
    //Tuning params;
    Target target(t);
    L->readPattern();
    rhs->readPattern();
    // we need put the final answer in a dense vector
    rhsDense = new Dense(*rhs);
    std::vector<Expr> el;
    std::vector<Argument> al;
    rhsDense->getDecl(el, al);
    args.insert(args.end(), el.begin(), el.end());
    argType.insert(argType.end(), al.begin(), al.end());

    TriangularInspector trnInspect(L, rhs);
    SymbolicObject *set = trnInspect.strategy(t.params);
    loweredKer = lower(this, set);
    Module module(fName, target);
    module.append(*this);
    //std::ofstream genCodeFile(fName + "_gen.cpp");
    Output::compile_to_source_c(module, fName + "_gen",set->IsVSBlock());
    //genCodeFile.close();
}

void Triangular::sympile_to_c(std::string fName) {
    Target t;
    sympile_to_c(fName, t);
}
}
}