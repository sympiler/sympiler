//
// Created by kazem on 4/18/17.
//

#include "Inspector.h"



namespace Sympiler {
namespace Internal {

SymbolicObject::SymbolicObject():
        isVIPrune(false),isVSBlock(false){
    setPtr = "setPtr";
    setVal = "setVal";
    setSize = "ub";
    blk2Col = "blk2Col";
    blkNo = "blkNo";
    BlockNo=0;
    averageBlockSize=0;
}

SymbolicObject::SymbolicObject(int ordr, int NNZ)
        :order(ordr), nnz(NNZ) {
    setPtr = "setPtr";
    setVal = "setVal";
    setSize = "ub";
    blk2Col = "blk2Col";
    blkNo = "blkNo";
    BlockNo=0;
    averageBlockSize=0;
    setPointer = new int[order+1]();
    block2Col = new int[order]();
    col2Block  = new int[order]();
    isVSBlock = false;
    isVIPrune = false;
}

SymbolicObject::~SymbolicObject() {
    delete[](setPointer);
    delete[] (setValue);
    delete[](col2Block);
    delete[] (block2Col);
}

SymbolicObject* Inspector::strategy(Tuning params) {

}

Inspector::Inspector() {
    sym=new SymbolicObject();

}

Inspector::~Inspector() {
    delete sym;
}

TriangularInspector::TriangularInspector(Matrix *DG, Matrix *rhs):
        Inspector(), DG(DG),rhs(rhs) {
    int n = DG->mPattern->n;
    int nnz = DG->mPattern->nnz;
    sym=new SymbolicObject(n,nnz);
    sym->nnzRes=nnz;//it is the same for triangular solver.

}

void TriangularInspector::superNodeDetection() {
    int prev, cur;
    int newRowSize=0, newNNZ=0, avgBSize=0;
    bool sim;
    int supNo=0;
    int n = DG->mPattern->n;
    int *col2sup = sym->col2Block, *col = DG->mPattern->p,
            *row = DG->mPattern->i, *sup2col = sym->block2Col;
    col2sup[0]=0;
    for (int i = 1; i < n; ++i) {
        sim=true;
        for (prev = col[i-1], cur = col[i]; 1; ) {
            if((row[prev] == i-1 && prev < col[i]) || (row[prev] == i && prev < col[i])){
                //skip diagonal block of prev col
                ++prev;
                continue;
            }
            if(row[cur]==i){//skip diagonal block of cur col
                ++cur;
                continue;
            }

            if(prev - col[i] != cur - col[i+1]){ // the off-diagonals length
                sim=false;
                break;
            }else if(prev - col[i]==0)
                break;
            if(row[prev] != row[cur]){//now off-diagonals
                sim=false;
                break;
            }else{
                ++prev;++cur;
                if (prev >= col[i] || cur >= col[i+1])
                    break;
            }
        }

        if(sim ){//col cur and nxt are similar
            col2sup[i]=supNo;
        }else{
            supNo++;
            col2sup[i]=supNo;
        }
    }
    supNo+= sim ;

    newNNZ=0;newRowSize=0;
    //newCol[0]=0;
    int curCol = 0, cnt=0, firstCol=0, tmpSize=0;
    for (int j = 0; j < supNo; ++j) {
        for (cnt=0; col2sup[curCol]==j && curCol<n; ++curCol, ++cnt);
        sup2col[j]=curCol;
        firstCol= j!=0 ? sup2col[j-1] : 0;
        avgBSize+=(curCol-firstCol);
        cnt--;//To find the last col offset of supernode
        tmpSize=col[firstCol+cnt+1]-col[firstCol+cnt]+cnt;
        newRowSize+=tmpSize;//diagonal dense part + off diagonal sparse one
        for (int i = firstCol; i < curCol; ++i) {
            newNNZ+=tmpSize;
            //newCol[i+1]=newNNZ;
        }
        if(curCol>=n)
            break;
    }
    sym->BlockNo=supNo;
    sym->averageBlockSize=n/supNo;
}

void  TriangularInspector::pruneSetDetection() {
    sym->setValue = new int[sym->order];//We do it here since its size is different in different kernels
    //For efficiency, We use two different versions of prune set detection
    //sym->setValue = new int[sym->order];
    if(!sym->isVSBlock){
        sym->setSizeNum = reach(DG->mPattern->n,DG->mPattern->p,DG->mPattern->i,
        rhs->mPattern->p,rhs->mPattern->i,0,sym->setValue, 0);
    }else{
        sym->setSizeNum = reach_sn(DG->mPattern->n,DG->mPattern->p,DG->mPattern->i,
              rhs->mPattern->p,rhs->mPattern->i,0,sym->setValue,
                 0,sym->BlockNo,sym->col2Block);
    }
    //TODO Future, for multiple right hand side
}

TriangularInspector::~TriangularInspector() {
}

SymbolicObject* TriangularInspector::strategy(Tuning params) {
    superNodeDetection();
    if(sym->averageBlockSize>params.BlockBound){
        sym->isVSBlock=true;
    }
    pruneSetDetection();
    if(sym->setSizeNum < sym->order-1) {//TODO: use a fraction of order?
        sym->isVIPrune = true;
    }
    return Inspector::sym;
}



}
}