//
// Created by kazem on 5/29/17.
//

#ifndef SYMPILER_PROJ_FACTORIZATION_H
#define SYMPILER_PROJ_FACTORIZATION_H

#include "Kernel.h"

namespace Sympiler {
namespace Internal {
class Factorization : public Kernel {
protected:
    Matrix *A;


public:
    Factorization():Kernel("Chol"){};
    Factorization(Matrix &a):Kernel("Chol"), A(&a){};

    ~Factorization(){};

    virtual Stmt baseCode()=0;

    virtual Stmt VSBlockIG(SymbolicObject *sym){}; //Inspector-guided transformations
    virtual Stmt VIPruneIG(SymbolicObject *sym){};

    virtual void sympile_to_c(std::string fName, Target t){};

    virtual void sympile_to_c(std::string fName){};
};

class Cholesky : public Factorization{
    Matrix *L;
    Stmt uncompressCol, update, factCol;
    Expr nSupR, supWdt;
    std::string tmpVec, finger;
    virtual Stmt VSBlockIG(SymbolicObject *sym);
    virtual Stmt VIPruneIG(SymbolicObject *sym);
    Stmt UncompressCol(Expr col, Matrix *A, std::string tmp,
                       ForType schedule);
    Stmt blockedUncompressCol( Expr col, Expr nxtCol,
                              Matrix *A, std::string tmp, ForType sched);

    Stmt Update(Matrix *A, Expr lbCol, Expr ubCol,
                std::string tmp, std::string extra);
    Stmt blockedUpdate(SymbolicObject *sym, Matrix *A, Expr curCol, Expr ubCol,
                        std::string tmp, std::string extra);

    Stmt FactCol(Matrix *A, Expr col, std::string tmp);
    Stmt blockedFactCol(Matrix *A, Expr col, std::string tmp);

public:
    Cholesky(){};
    Cholesky(Sparse &A);
    ~Cholesky();
    virtual Stmt baseCode();
    virtual void sympile_to_c(std::string fName, Target t);
    virtual void sympile_to_c(std::string fName);

};

}
}

#endif //SYMPILER_PROJ_FACTORIZATION_H
