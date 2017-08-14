//
// Created by kazem on 4/18/17.
//

#ifndef SYMPILER_PROJ_INSPECTOR_H
#define SYMPILER_PROJ_INSPECTOR_H


#include <vector>
#include "Matrix.h"
#include "NumericalUtils.h"
#include "Target.h"

namespace Sympiler {
namespace Internal {



class SymbolicObject {
public:
    int order;
    int nnz;
    int nnzRes;
    int *setPointer;
    int *setValue;
    int setSizeNum;
    int *block2Col;
    int *col2Block;
    int BlockNo;
    int averageBlockSize;
    bool isVSBlock, isVIPrune;

public:
    std::string setPtr, setVal, setSize; //These are variables for defining the resulting set
    std::string blk2Col, blkNo;
    SymbolicObject() ;
    SymbolicObject(int, int);
    ~SymbolicObject();
    bool IsVSBlock(){ return isVSBlock;};
    bool IsVIPrune(){ return isVIPrune;};
    void setIsVIPrune(bool val){isVIPrune=val;};
    Expr getBlockNo(){ return Variable::make(halide_type_t(halide_type_int,32),blkNo); };
    //Expr getblk2Col(){ return Pointer::make(); };
};

class Inspector {
protected:
    SymbolicObject *sym;

public:
    Inspector();
    ~Inspector();
    virtual SymbolicObject* strategy(Tuning );

};

class TriangularInspector: public Inspector{
    Matrix *DG, *rhs;
    void superNodeDetection();
    void pruneSetDetection();
public:
    TriangularInspector(Matrix *DG, Matrix *rhs);
    ~TriangularInspector();
    virtual SymbolicObject* strategy(Tuning );

};



}
}


#endif //SYMPILER_PROJ_INSPECTOR_H
