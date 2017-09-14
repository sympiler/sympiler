//
// Created by kazem on 7/25/17.
//

#ifndef CHOLOPENMP_INSPECTION_PRUNE_H
#define CHOLOPENMP_INSPECTION_PRUNE_H

#include "Reach.h"

bool pruning(int n, int* c, int* r, int *ET, int* &prunePtr, int* &pruneSet){
    //Specifying the row patterns of L and pruneSet
    int top;
    int *xi = new int[2*n];
    prunePtr[0]=prunePtr[1]=0;
    for (int colNo = 1; colNo < n; ++colNo) {
        top = ereach(n,c,r,colNo,ET,xi,xi+n);
        prunePtr[colNo+1]=prunePtr[colNo]+(n-top);
        for (int i = top, cnt=0; i < n; ++i, ++cnt) {
            pruneSet[prunePtr[colNo]+cnt]=xi[i];
        }
    }
    return true;
}

bool getBlockedPruneSet(int n, int *Ap, int *Ai, int *col2sup,
                 const int *eTree, int *blockSet,
                 int *prunePtr, int *pruneSet){
    int top=0;
    int *xi = new int[2*n]();
    prunePtr[0]=0;
    for (int i = 0; i < n; ++i) {
        int curCol = blockSet[i];
        int nxtCol = blockSet[i+1];
       // ASSERT(curCol < 17361 && nxtCol < 17361);
        top=ereach_sn(n,Ap,Ai,curCol,nxtCol,col2sup, eTree,
                      xi,xi+n);
        //printf("%d \n",i);
        prunePtr[i+1]=prunePtr[i]+(n-top);
        for (int j = top,k=0; j < n; ++j, ++k) {
            pruneSet[prunePtr[i]+k] = xi[j];
           // ASSERT(prunePtr[i]+k < 100380);
            //   printf("%d, ", pruneSet[prunePtr[i]+k]);
        }
    }
#if 0
    std::ofstream spFile1;
    spFile1.open("/home/kazem/UFDB/sp2.txt");
    for (int l = 1; l < n; ++l) {
        spFile1<<"\n"<<l<<"\n";
        for (int i = prunePtr[l-1]; i <prunePtr[l] ; ++i) {
            spFile1<<pruneSet[i]<<",";
        }
    }
    spFile1.close();
#endif
    delete []xi;
    return true;
}

#endif //CHOLOPENMP_INSPECTION_PRUNE_H
