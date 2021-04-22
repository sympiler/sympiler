//
// Created by kazem on 6/13/17.
//


#include "FactorizationInspector.h"
//#include "NumericalUtils.h"

namespace Sympiler {
namespace Internal {

FactorizationInspector::FactorizationInspector(Matrix *A):Inspector(), A(A) {
    int n=A->mPattern->n;
    eTree=new int[n];
    colCount=new int[n];
    postETree = new int[n];
}

FactorizationInspector::~FactorizationInspector() {
    delete[](eTree);
    delete[](colCount);
}

void FactorizationInspector::inspectionGraph(){

}

SymbolicObject* FactorizationInspector::strategy() {//TODO
    MatrixPattern *DGPattern = A->mPattern;

    return Inspector::sym;
}


/*
 * Cholesky Inspector Implementation
 */

CholeskyInspector::CholeskyInspector(Matrix *A):FactorizationInspector(A) {
}

CholeskyInspector::~CholeskyInspector() {
}

void CholeskyInspector::inspectionGraph(){
    MatrixPattern *APattern = A->mPattern;
    eTree = etree(APattern->n,APattern->p,APattern->i,0);//for now it is for symmetric case
    postETree = postOrder(eTree,APattern->n);
    colCount = counts(APattern->n,APattern->n,APattern->p,
                      APattern->i,eTree,postETree,0);
    sym=new SymbolicObject(APattern->n,APattern->nnz);
    sym->nnzRes=0;
    for (int i = 0; i < sym->order; ++i) {//couting nnz of L
        sym->nnzRes += colCount[i];
    }
}

void CholeskyInspector::blockDetection() {
    sym->col2Block = new int[A->mPattern->n];
    sym->block2Col = superNodeDetection(eTree,colCount,A->mPattern->n,
                                        sym->col2Block,sym->BlockNo,
                                        sym->averageBlockSize);
}

void CholeskyInspector::pruneCalculation() {
    sym->setValue = new int[sym->nnzRes];//We do it here since its size is different in different kernels
    int top,n=A->mPattern->n;
    int *xi = new int[2*A->mPattern->n];
    //For efficiency, We use two different versions of prune set detection
    if(!sym->isVSBlock){
        for (int colNo = 1; colNo < n; ++colNo) {
            top = ereach(n,A->mPattern->p,A->mPattern->i,
                         colNo,eTree,xi,xi+n);
            sym->setPointer[colNo+1] = sym->setPointer[colNo]+(n-top);
            for (int i = top, cnt=0; i < n; ++i, ++cnt) {
                sym->setValue[sym->setPtr[colNo]+cnt]=xi[i];
            }
        }
    }else{
        for (int colNo = 1; colNo < sym->BlockNo; ++colNo) {
            int curCol = colNo!=0 ?  sym->block2Col[colNo - 1] : 0;
            int nxtCol = sym->block2Col[colNo];
            top = ereach_sn(n,A->mPattern->p,A->mPattern->i,
                         curCol,nxtCol,sym->col2Block, eTree,xi,xi+n);
            sym->setPointer[colNo+1] = sym->setPointer[colNo]+(n-top);
            for (int i = top, cnt=0; i < n; ++i, ++cnt) {
                sym->setValue[sym->setPointer[colNo]+cnt] = xi[i];
            }
        }
    }
}

SymbolicObject* CholeskyInspector::strategy(Tuning params ) {
    inspectionGraph();
    blockDetection();
    sym->isVSBlock=true;//For now since we use block for Chokesky
    if(sym->averageBlockSize>params.BlockBound){
        sym->isVSBlock=true;
    }
    pruneCalculation();
    sym->isVIPrune = true;
    return Inspector::sym;
}

}
}
