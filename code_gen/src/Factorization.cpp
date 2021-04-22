//
// Created by kazem on 5/29/17.
//

#include <fstream>
#include "Factorization.h"
#include "IROperator.h"
#include "Module.h"
#include "Output.h"
#include "Lower.h"
#include "VIPrune.h"
#include "FactorizationInspector.h"

namespace Sympiler {
namespace Internal {

Cholesky::Cholesky(Sparse &a):Factorization(a),
                              tmpVec("tempVec"), finger("finger") {
    //L = new Sparse(a.NNZ(),a.getType(),a.Order_Num(),a.Dim(),"");
    L = new Sparse(a.getType(),a.Path());
    std::vector<Expr> el;
    std::vector<Argument> al;
    a.getDecl(el,al);
    args.insert(args.end(),el.begin(),el.end());
    argType.insert(argType.end(),al.begin(),al.end());
    el.clear();al.clear();
    L->getDecl(el,al);
    args.insert(args.end(),el.begin(),el.end());
    argType.insert(argType.end(),al.begin(),al.end());
    el.clear();al.clear();

    argType.emplace_back(Argument(tmpVec,Argument::Kind::InputBuffer,
                         Int(32),0));
//    Expr fingerExpr = Variable::make(Float(64),finger);
    argType.emplace_back(Argument(finger,Argument::Kind::InputBuffer,
                                  Float(64),0));

}

Cholesky::~Cholesky(){
    delete L;
}

Stmt Cholesky::baseCode() {
    Stmt s,t;
    Expr curCol = Variable::make(Int(32),Kernel::name+"f0");
    std::string tmpVec = "tempVec", finger = "finger";
    s = Free::make(tmpVec);
    s = Block::make(Free::make(finger),s);
    factCol = FactCol(Factorization::A,curCol,tmpVec);
    s = Block::make(factCol,s);
    update = Update(Factorization::A,make_zero(Int(32)),curCol,tmpVec,finger);
    s = Block::make(update,s);
    uncompressCol = UncompressCol(curCol,Factorization::A,tmpVec,ForType::Serial);
    s = Block::make(uncompressCol,s);
    s = For::make(Kernel::name+"0",make_zero(Int(32)),
                  A->Order(),ForType::Serial,DeviceAPI::Host, s);
    s = Allocate::make(tmpVec,Float(64),{A->Order()},const_true(),s);
    s = Allocate::make(finger,Float(64),{A->Order()},const_true(),s);
    Kernel::loweredKer=s;
    return Kernel::loweredKer;
}

Stmt Cholesky::UncompressCol(Expr col, Matrix *A, std::string tmp, ForType sched) {
    Stmt uncomp;
    Expr temp = Variable::make(A->getType(),tmp);
    Expr i0 = Variable::make(Int(32),"uncmpress");//TODO add an static member to make the loop variables unique
    Expr lhs = Pointer::make(A->getType(),tmp,A->accessRowIdx(i0));
    uncomp = Store::make(lhs,"stOp",A->accessNNZ(i0));
    Expr colp1 = Add::make(A->accessCol(col), make_one(Int(32)));
    uncomp = For::make("c0",A->accessCol(col),A->accessCol(colp1),sched,
                       DeviceAPI::Host, uncomp);
    return uncomp;
}

Stmt Cholesky::Update(Matrix *A, Expr lbCol, Expr ubCol,
                      std::string tmp, std::string extra){
    Stmt upd;
     //f[lR[l]] -= lValues[l] * lValues[finger[spCol]];
    Expr i0 = Variable::make(Int(32),"u0");
    Expr i1 = Variable::make(Int(32),"u1");
    Expr lhs = Pointer::make(A->getType(),tmp,A->accessRowIdx(i1));
    Expr finger = Pointer::make(A->getType(),extra,i0);//finger[spcol]
    Expr contrib = Mul::make(A->accessNNZ(i1),A->accessNNZ(finger));//lValues[l] * lValues[finger[spCol]]
    upd = Store::make(lhs,"stUp",Sub::make(lhs,contrib));
    Expr colp1 = Add::make(L->accessCol(i0), make_one(Int(32)));//i0+1
    upd = For::make("u1",L->accessCol(i0),L->accessCol(colp1),ForType::Serial,
                       DeviceAPI::Host, upd);
    upd = For::make("u0", lbCol, ubCol, ForType::Pruned,
                    DeviceAPI::Host, upd);
    return upd;
}

Stmt Cholesky::FactCol(Matrix *A, Expr col, std::string tmp) {
        Expr i0 = Variable::make(Int(32),"fa0");
        Stmt fact;
        std::vector<Expr> sqrtArgs;
        Expr diag = Pointer::make(L->getType(),tmp,col);
        sqrtArgs.push_back(diag);
        Expr sq = Call::make(Float(64),"sqrt",sqrtArgs);//FIXME
        Expr curTemp = Pointer::make(L->getType(),tmp,L->accessRowIdx(i0));
        Expr div = Div::make(curTemp,sq);//div = f[lR[i0]] / tmpSqrt;
        Stmt stFa = Store::make(L->accessNNZ(i0),"stFa",div);//lValues[i0]=div
        fact = Store::make(curTemp,"makeZero",
                           make_zero(Int(32)));
        fact = Block::make(stFa,fact);
        Expr colp1 = Add::make(L->accessCol(col), make_one(Int(32)));//i0+1
        fact = For::make("fa0",L->accessCol(col),L->accessCol(colp1),ForType::Serial,
                         DeviceAPI::Host, fact);
        return fact;
}

Stmt Cholesky::blockedUncompressCol(Expr curCol, Expr nxtCol, Matrix *A,
                                    std::string tmp, ForType sched) {
    Stmt uncomp;
    //int nSupR = Li_ptr[nxtCol]-Li_ptr[curCol];
    //nSupR = Sub::make(A->accessRowIdx(nxtCol),A->accessRowIdx(curCol));
    Expr temp = Variable::make(A->getType(),tmp);
    Expr i0 = Variable::make(Int(32),"c0");//TODO add an static member to make the loop variables unique
    Expr i1 = Variable::make(Int(32),"c1");
    Expr i2 = Variable::make(Int(32),"c2");
    Expr cnt = Variable::make(Int(32),"uncmpress");
    Expr lhs = Pointer::make(L->getType(),tmp,L->accessRowIdx(i0));
    Expr cntP1 = Sub::make(i0,L->accessRowIdxPntr(curCol));
    //Copy to L
    //lValues[lC[i]+map[r[j]]] = values[j];
    Expr idx2 = Add::make(L->accessCol(i1),Pointer::make(A->getType(),tmp,A->accessRowIdx(i2)));
    Stmt loopNest2 = Store::make(L->accessNNZ(idx2),"copyToL",A->accessNNZ(i2));
    loopNest2 = For::make("c2",A->accessCol(i1),
                          A->accessCol(Add::make(i1,make_one(Int(32)))),sched,
                       DeviceAPI::Host, loopNest2);
    loopNest2 = For::make("c1",curCol,nxtCol,sched, DeviceAPI::Host, loopNest2);

    //Map
   // uncomp = Store::make(cnt,"CntSt",cntP1);
    uncomp = Block::make(Store::make(lhs,"stOp",cntP1), uncomp);
    uncomp = For::make("c0",L->accessRowIdxPntr(curCol),L->accessRowIdxPntr(nxtCol),sched,
                       DeviceAPI::Host, uncomp);
    //Stmt initCnt = Store::make(cnt,"initCnt",make_zero(Int(16)));
    //uncomp = Block::make(initCnt,uncomp);
    uncomp = Block::make(uncomp,loopNest2);//joining two loops
    return uncomp;
}

Stmt Cholesky::blockedUpdate(SymbolicObject *sym,Matrix *A, Expr curCol, Expr ubCol,
                  std::string tmp, std::string extra){
    Stmt upd;
    //f[lR[l]] -= lValues[l] * lValues[contribs[spCol]];
    Expr curSNode = Variable::make(Int(32),Kernel::name+"f0");
    Expr i0 = Variable::make(Int(32),"u0");
    Expr i1 = Variable::make(Int(32),"u1");
    Expr lhs = Pointer::make(L->getType(),tmp,L->accessRowIdx(i1));
    Expr contribs = Variable::make(Float(64),extra);
    Expr one = Variable::make(Float(64),"one"); //TODO make it global
    Expr zero = Variable::make(Float(64),"zero");

    //int cSN = blockSet[lSN];//first col of current SN
    Expr cSN = Pointer::make(Int(32),sym->blk2Col,i0);
    Expr i0P1 = Add::make(i0,make_one(Int(32)));//i0+1
    //int cNSN = blockSet[lSN+1];//first col of Next SN
    Expr cNSN = Pointer::make(Int(32),sym->blk2Col,i0P1);
    //int Li_ptr_cNSN = Li_ptr[cNSN];
    Expr Li_ptr_cNSN = L->accessRowIdxPntr(cNSN);
    //int Li_ptr_cSN = Li_ptr[cSN];
    Expr Li_ptr_cSN = L->accessRowIdxPntr(cSN);
    //int nSNRCur=Li_ptr_cNSN-Li_ptr_cSN;
    Expr nSNRCur = Sub::make(Li_ptr_cNSN,Li_ptr_cSN);
    //int  supWdt=cNSN-cSN;//The width of current src SN
    //supWdt = Sub::make(cNSN,cSN);
    Expr sw = Variable::make(Bool(),"sw");

    //cur=&lValues[lC[curCol]];
    Expr cur = L->accessNNZ( L->accessCol(i0) );


    Expr cond1 = And::make(GE::make(L->accessRowIdx(i1),curCol) , sw );
    Expr cond2 = And::make(LT::make(L->accessRowIdx(i1),Add::make(curCol,supWdt)) , Not::make(sw) );
    Expr bnd = Sub::make(i1,Li_ptr_cSN);
    Expr lb = Variable::make(Int(32),"lb1");
    Expr ub = Variable::make(Int(32),"ub1");


    //nSupRs=Li_ptr_cNSN-Li_ptr_cSN-lb;
    Expr nSupRs = Sub::make(nSNRCur,lb);
    //int ndrow1=ub-lb+1;
    Expr ndrow1 = Add::make( Sub::make(ub,lb),make_one(Int(32)) );
    //int ndrow3 = nSupRs-ndrow1;
    Expr ndrow3 = Sub::make(nSupRs,ndrow1);
    //src=&lValues[lC[cSN]+lb];//first element of src supernode starting from row lb
    Expr src = L->accessNNZ( Add::make(L->accessCol(cSN),lb) );
    //double *srcL = &lValues[lC[cSN]+ub+1];
    Expr ubP1 = Add::make(ub,make_one(Int(32)));
    Expr srcL = L->accessNNZ( Add::make(L->accessCol(cSN),ubP1) );


    //upd = Store::make(lhs,"stUp",Sub::make(lhs,contrib));
    Expr i2 = Variable::make(Int(32),"uc0");
    Expr i3 = Variable::make(Int(32),"uc1");
    //int col=map[lR[Li_ptr_cSN+i+lb]];
    Expr idxTmp = Add::make(Li_ptr_cSN,lb);
    Expr colTmp = Pointer::make( Int(32),tmp,L->accessRowIdx(Add::make(idxTmp,i2)) );
    //int cRow= lR[Li_ptr_cSN+j+lb];
    Expr cRow = L->accessRowIdx(Add::make(idxTmp,i3));
    //cur[col*nSupR+map[cRow]]
    Expr cpyBack = Add::make( Mul::make(colTmp,nSupR), Pointer::make(Int(32),tmp,cRow) );
    //cpyBack = Add::make(cpyBack,L->accessCol(i0));
    //Expr lhsCpy = Pointer::make(Float(64),"cur",cpyBack);//
    Expr lhsCpy =  L->accessNNZ(Add::make(cpyBack,L->accessCol(curCol)));//
    //contribs[i*nSupRs+j]
    Expr rhsCpy = Pointer::make(Float(64),extra, Add::make( Mul::make(i2,nSupRs), i3) );
    Expr supWdts = Sub::make(cNSN,cSN);

    //cur[col*nSupR+map[cRow]] -= contribs[i*nSupRs+j];
    upd = Store::make(lhsCpy,"cpBack",Sub::make(lhsCpy,rhsCpy));
    upd = For::make("uc1",i2,nSupRs,ForType::Serial,
                    DeviceAPI::Host, upd);
    upd = For::make("uc0",make_zero(Int(32)),ndrow1,ForType::Serial,
                    DeviceAPI::Host, upd);


    /*if(ndrow3>0){
        dgemm_("N","C",&ndrow3,&ndrow1,&supWdt,one,srcL,&nSNRCur,
               src,&nSNRCur,zero,contribs+ndrow1,&nSupRs );
    }*/
    std::vector<Expr> gemmArgs;
    gemmArgs.push_back(ndrow3); gemmArgs.push_back(ndrow1);
    gemmArgs.push_back(supWdts); gemmArgs.push_back(one);
    gemmArgs.push_back(srcL); gemmArgs.push_back(nSNRCur);
    gemmArgs.push_back(src); gemmArgs.push_back(nSNRCur);
    gemmArgs.push_back(zero); gemmArgs.push_back( Pointer::make(Float(64), "finger",ndrow1) );
    gemmArgs.push_back(nSupRs);

    std::vector<Argument> gemmArgums;
    gemmArgums.emplace_back(Argument(true));gemmArgums.emplace_back(Argument(true));
    gemmArgums.emplace_back(Argument(true));gemmArgums.emplace_back(Argument(false));
    gemmArgums.emplace_back(Argument(true));gemmArgums.emplace_back(Argument(true));
    gemmArgums.emplace_back(Argument(true));gemmArgums.emplace_back(Argument(true));
    gemmArgums.emplace_back(Argument(false));gemmArgums.emplace_back(Argument(true));
    gemmArgums.emplace_back(Argument(true));
    Stmt gemm = CallX::make("gemm",gemmArgs,gemmArgums);
    upd = Block::make(IfThenElse::make( GT::make(ndrow1,make_zero(Int(32))),gemm ), upd);

    //dsyrk_("L","N",&ndrow1,&supWdt,one,src,&nSNRCur,zero,contribs,&nSupRs);
    std::vector<Expr> syrkArgs;
    syrkArgs.push_back(ndrow1); syrkArgs.push_back(supWdts);
    syrkArgs.push_back(one); syrkArgs.push_back(src);
    syrkArgs.push_back(nSNRCur);syrkArgs.push_back(zero);
    syrkArgs.push_back(contribs); syrkArgs.push_back(nSupRs);

    std::vector<Argument> syrkArgums;
    syrkArgums.emplace_back(Argument(true));syrkArgums.emplace_back(Argument(true));
    syrkArgums.emplace_back(Argument(false));syrkArgums.emplace_back(Argument(true));
    syrkArgums.emplace_back(Argument(true));syrkArgums.emplace_back(Argument(false));
    syrkArgums.emplace_back(Argument(false));syrkArgums.emplace_back(Argument(true));
    upd = Block::make( CallX::make("syrk",syrkArgs,syrkArgums), upd );

    Stmt lbStmt=Block::make(Store::make(lb,"lbSet",bnd),lbStmt);
    lbStmt = Block::make(Store::make(sw,"swreset",make_zero(Int(32))), lbStmt) ;
    Stmt ubStmt=Store::make(ub,"ubSet",bnd);
    //Finding the
    Stmt updBodyu1 = IfThenElse::make(cond2,ubStmt);
    updBodyu1 = Block::make( IfThenElse::make(cond1,lbStmt) ,updBodyu1);
    updBodyu1 = For::make("u1",Li_ptr_cSN,Li_ptr_cNSN,ForType::Serial,
                    DeviceAPI::Host, updBodyu1);
    upd = Block::make(updBodyu1,upd);
    //upd = Block::make( Store::make(sw,"swInit",make_one(Int(32))), upd );//sw=true
    upd = Block::make( Store::make(sw,"swInit",make_one(Int(32))), upd );//sw=true
    upd = Block::make( Store::make(lb,"lbInit",make_zero(Int(32))), upd );//sw=true
    upd = Block::make( Store::make(ub,"ubInit",make_zero(Int(32))), upd );//sw=true

    upd = For::make("u0", Sub::make(curSNode,make_one(Int(32))),
                    curSNode, ForType::Pruned, DeviceAPI::Host, upd);
    //upd = Block::make( Store::make(,"cur",cur),upd);
    return upd;
}

Stmt Cholesky::blockedFactCol(Matrix *A, Expr col, std::string tmp) {
    Stmt s;
    //Cholesky_col(nSupR,supWdt,cur);
    Expr cur = L->accessNNZ( L->accessCol(col) );
    Expr curPWdth = L->accessNNZ(Add::make(L->accessCol(col),supWdt));
    Expr info = Variable::make(Int(32),"info");
    std::vector<Expr> choleskyArgs;
     choleskyArgs.push_back(supWdt);
    choleskyArgs.push_back(cur);choleskyArgs.push_back(nSupR); choleskyArgs.push_back(info);
    //int rowNo=nSupR-supWdt;
    //dtrsm_("R", "L", "C", "N", &rowNo, &supWdt,one,
    //       cur,&nSupR,&cur[supWdt],&nSupR);
    Expr one = Variable::make(Float(64),"one");
    Expr zero = Variable::make(Float(64),"zero");

    std::vector<Expr> trsmArgs;
    Expr rowNo = Sub::make(nSupR,supWdt);
    trsmArgs.push_back(rowNo); trsmArgs.push_back(supWdt);
    trsmArgs.push_back(one); trsmArgs.push_back(cur);
    trsmArgs.push_back(nSupR); trsmArgs.push_back(curPWdth);
    trsmArgs.push_back(nSupR);

    std::vector<Argument> trsmArgums;
    trsmArgums.emplace_back(Argument(true));trsmArgums.emplace_back(Argument(true));
    trsmArgums.emplace_back(Argument(false));trsmArgums.emplace_back(Argument(true));
    trsmArgums.emplace_back(Argument(true));trsmArgums.emplace_back(Argument(true));
    trsmArgums.emplace_back(Argument(true));
    s = CallX::make("trsm",trsmArgs,trsmArgums);

    std::vector<Argument> potrfArgums;
    potrfArgums.emplace_back(Argument(true));potrfArgums.emplace_back(Argument(true));
    potrfArgums.emplace_back(Argument(true));potrfArgums.emplace_back(Argument(true));
    s = Block::make(CallX::make("potrf",choleskyArgs, potrfArgums),s);
    return s;
}

Stmt Cholesky::VSBlockIG(SymbolicObject *sym){
//Adding the arguments for the inspection set
    if(!isSetArgInserted){
        argType.push_back(Argument(sym->setPtr,Argument::Kind::InputBuffer,
                                   Int(32),1));
        argType.push_back(Argument(sym->setVal,Argument::Kind::InputBuffer,
                                   Int(32),1));
        argType.push_back(Argument(sym->setSize,Argument::Kind::InputScalar,
                                   Int(32),1));
        isSetArgInserted=true;//We need only one set, other
    }
    //Adding arguments for block to column mapping
    argType.push_back(Argument(sym->blk2Col,Argument::Kind::InputBuffer,
                               Int(32),1));
    argType.push_back(Argument(sym->blkNo,Argument::Kind::InputScalar,
                               Int(32),1));
    //Generating code for the blocked triangular solver
    Stmt s,t;
    Expr retVal = Variable::make(Int(32),"retval");
    Expr zro = make_zero(Int(32));
    Expr one = make_one(Int(32));
    Expr i0 = Variable::make(Int(32),Kernel::name+"f0");
    Expr suptoCol = Pointer::make(Int(32),sym->blk2Col,i0);
    Expr nxtCol = Pointer::make(Int(32),sym->blk2Col,i0);
    Expr curColtmp = Pointer::make(Int(32),sym->blk2Col,
                                   Sub::make(i0,one) );
    Expr curCol = Select::make( NE::make(i0,zro), curColtmp, zro);
    supWdt = Sub::make(nxtCol,curCol);
    nSupR = Sub::make(L->accessRowIdxPntr(nxtCol), L->accessRowIdxPntr(curCol));
    //curCol = Variable::make(Int(16),Kernel::name+"f0");
    //std::string tmpVec = "tempVec", finger = "finger";
   // s = Free::make(tmpVec);
   // s = Block::make(Free::make(finger),s);
    factCol = blockedFactCol(Factorization::A,curCol,tmpVec);
    s = Block::make(factCol,s);
    update = blockedUpdate(sym,Factorization::A,curCol,curCol,tmpVec,finger);
    s = Block::make(update,s);
    uncompressCol = blockedUncompressCol(curCol,nxtCol,Factorization::A,tmpVec,ForType::Serial);
    s = Block::make(uncompressCol,s);
    s = For::make(Kernel::name+"f0",make_one(Int(32)),
                  Add::make(Variable::make(Int(32),sym->blkNo),make_one(Int(32))),
                  ForType::Serial,DeviceAPI::Host, s);
  //  s = Allocate::make(tmpVec,Float(64),{A->Order()},const_true(),s);
  //  s = Allocate::make(finger,Float(64),{A->Order()},const_true(),s);
    //s = Allocate::make("sw",Bool(),{make_one(Int(16))},const_true(),s);

    Kernel::loweredKer=s;
    return Kernel::loweredKer;
}

Stmt Cholesky::VIPruneIG(SymbolicObject *sym){
    Expr zro = make_zero(Int(32));
    Expr one = make_one(Int(32));
    Expr i0 = Variable::make(Int(32),Kernel::name+"f0");
    Expr suptoCol = Pointer::make(Int(32),sym->blk2Col,i0);
    Expr nxtCol = Pointer::make(Int(32),sym->blk2Col,i0);
    Expr curColtmp = Pointer::make(Int(32),sym->blk2Col,
                                   Sub::make(i0,one) );
    Expr curCol = Select::make( NE::make(i0,zro), curColtmp, zro);

    VIPrune vi_prune(Pointer::make(Int(32),sym->setPtr,Sub::make(i0,make_one(Int(32))))
            ,sym->setVal, Pointer::make(Int(32),sym->setPtr,i0) );
    //Adding the arguments for the inspection set
    if(!isSetArgInserted){
        argType.push_back(Argument(sym->setPtr,Argument::Kind::InputBuffer,
                                   Int(32),1));
        argType.push_back(Argument(sym->setVal,Argument::Kind::InputBuffer,
                                   Int(32),1));
        argType.push_back(Argument(sym->setSize,Argument::Kind::InputScalar,
                                   Int(32),1));
        isSetArgInserted=true;//We need only one set, other
    }
    loweredKer = vi_prune.mutate(loweredKer);
    return loweredKer;
    //return vi_prune.mutate(loweredKer);
}

void Cholesky::sympile_to_c(std::string fName, Target t) {
    Target target(t);
    A->readPattern();
    L->mPattern->n = A->mPattern->n;
    L->mPattern->nnz = A->mPattern->nnz;
    CholeskyInspector cholInspect(A);
    SymbolicObject *set = cholInspect.strategy(t.params);
    set->setIsVIPrune(true);//VIPrune always happens for Cholesky
    loweredKer=lower(this,set);
    Module module(fName, target);
    module.append(*this);
    //std::ofstream genCodeFile(fName + "_gen.cpp");
    Output::compile_to_source_c(module, fName + "_gen", set->IsVSBlock());
   // genCodeFile.close();
}

void Cholesky::sympile_to_c(std::string fName) {
    Target t;
    sympile_to_c(fName,t);
}

}
}