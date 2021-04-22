//
// Created by kazem on 4/20/17.
//

#ifndef SYMPILER_PROJ_MATRIX_H
#define SYMPILER_PROJ_MATRIX_H

#include <string>
#include "Type.h"
#include "Expr.h"
#include "IR.h"
#include "Argument.h"

namespace Sympiler {
namespace Internal {

class MatrixPattern {

public:
    MatrixPattern();
    ~MatrixPattern();
    bool read();

    int n;//row number
    int m; // column number
    int nnz;
    int *p; //column pointer
    int *i; //row index
    //std::string path;
};

class Matrix {
protected:
    Type t;
    std::string name;
    std::string order;//The matrix order, assumed square
    std::string path;

    //int order_num;
    int dim;
    static int mNo;
public:
    Matrix();
    Matrix(Type , int , int , std::string , std::string ="A");
    Matrix(const Matrix& );
    ~Matrix(){};
    std::string Name(){ return name;};
    MatrixPattern *mPattern;
    Expr Order();
    std::string Path(){ return path;};
    Type getType(){ return t;};
    int Dim(){return dim;};
    int Order_Num(){ return mPattern->n;};
    void setOrder_Nym(int );
    static void incNo(){mNo++;};
    //virtual bool readPattern();
    virtual Expr diagonal(Expr e);
    virtual Expr accessRowIdx(Expr e);
    virtual Expr accessRowIdxPntr(Expr e);
    virtual Expr accessCol(Expr e);
    virtual Expr accessNNZ(Expr e);
    virtual Stmt allocateCol();
    virtual Stmt allocateRow();
    virtual Stmt allocateNNZ();
    virtual MatrixPattern* readPattern();
    virtual void getDecl(std::vector<Expr>& , std::vector<Argument>&){};
};

class Sparse:public Matrix{//CSC, will change once added other storages
    int nnz;
    std::string row, col, nz;//Access pointers
    std::string rowP;
    bool isBlock;
public:
    Sparse(int , Type , int , int , std::string , std::string="L");
    Sparse(Type , std::string);
    ~Sparse();
    int NNZ(){return nnz;};
    //virtual bool readPattern();
    virtual Expr diagonal(Expr e);
    virtual Expr accessRowIdx(Expr e);
    virtual Expr accessRowIdxPntr(Expr e);
    virtual Expr accessCol(Expr e);
    virtual Expr accessNNZ(Expr e);
    virtual Stmt allocateCol();
    virtual Stmt allocateRow();
    virtual Stmt allocateNNZ();
    virtual MatrixPattern* readPattern();
    virtual void getDecl(std::vector<Expr>& , std::vector<Argument>&);
};

class Dense:public Matrix{//FIXME Not fully implemented yet 1D dense array for now
public:
    Dense();
    Dense(Type , int , int , std::string , std::string ="A");
    Dense(Type , std::string );
    Dense(Sparse );
    Dense(const Matrix& m);
    virtual Expr diagonal(Expr e);
    virtual Expr accessNNZ(Expr e);
    virtual void getDecl(std::vector<Expr>& , std::vector<Argument>&);
};


}
}


#endif //SYMPILER_PROJ_MATRIX_H
