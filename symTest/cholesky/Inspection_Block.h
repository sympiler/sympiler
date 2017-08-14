//
// Created by kazem on 7/25/17.
//

#ifndef CHOLOPENMP_INSPECTION_BLOCK_H
#define CHOLOPENMP_INSPECTION_BLOCK_H

#include "Etree.h"
#include "PostOrder.h"
#include "colCounts.h"



int* supernodal(int* Etree, int* postETree, int* colCount, int n,
                int* col2sup, int& bsSize){
    int *blockSet = new int[n]();
    int blockSetSize=0;
    for (int i = 0; i < n; ++i) {
        int col1 = postETree[i];
        int col2 = postETree[i+1];
        if(Etree[col1] == col2 && colCount[col1]-1 == colCount[col2]){
            blockSet[blockSetSize+1]++;//increase the bound of SN blockSetSize
        }else{
            blockSet[blockSetSize+1] += blockSet[blockSetSize]+1;
            ++blockSetSize;//next SN
        }
    }
    bsSize=blockSetSize;
    for (int sNo = 0; sNo < blockSetSize; ++sNo) {
        for (int i = blockSet[sNo]; i < blockSet[sNo+1]; ++i) {
            col2sup[i]=sNo;
        }
    }
    return blockSet;
}

int* sNodeDetection(int* Parent, int* postETree, int* ColCount, int n,
                    int* col2sup, int& bsSize){
    int *Wi = new int[n]();
    int *blockSet = new int[n]();
    int blockSetSize=0, j=0, parent;
    /* ---------------------------------------------------------------------- */
    /* find the fundamental supernodes */
    /* ---------------------------------------------------------------------- */

    /* count the number of children of each node, using Wi [ */
    for (j = 0 ; j < n ; j++)
    {
        parent = Parent [j] ;
        if (parent != EMPTY)
        {
            Wi [parent]++ ;
        }
    }

    //blockSet = Head ;  /* use Head [0..bsSize] as workspace for blockSet list ( */

    /* column 0 always starts a new supernode */
    bsSize = (n == 0) ? 0 : 1 ;	/* number of fundamental supernodes */
    blockSet [0] = 0 ;

    for (j = 1 ; j < n ; j++)
    {
        /* check if j starts new supernode, or in the same supernode as j-1 */
        if (Parent [j-1] != j	    /* parent of j-1 is not j */
            || (ColCount [j-1] != ColCount [j] + 1) /* j-1 not subset of j*/
            || Wi [j] > 1	    /* j has more than one child */
                )
        {
            /* j is the leading node of a supernode */
            blockSet [bsSize++] = j ;
        }
    }
    blockSet [bsSize] = n ;

    /* contents of Wi no longer needed for child count ] */

    /* ---------------------------------------------------------------------- */
    /* find the mapping of fundamental nodes to supernodes */
    /* ---------------------------------------------------------------------- */

    //col2sup = Wj ;	/* use Wj as workspace for col2sup [ */
    //int *col2sup = new int[bsSize]();
    /* col2sup [k] = s if column k is contained in supernode s */
    for (int s = 0 ; s < bsSize ; s++)
    {
        for (int k = blockSet [s] ; k < blockSet [s+1] ; k++)
        {
            col2sup [k] = s ;
        }
    }
    return blockSet;
}

int* eTree2aTree(int* eTree, int n, int *blockSet, int* col2sup, int bsSize){
    int *aTree = new int[bsSize]();
    int parent,j;
    for (int s = 0 ; s < bsSize ; s++)
    {
        j = blockSet [s+1] - 1 ;	/* last node in supernode s */
        parent = eTree [j] ;	/* parent of last node */
        aTree [s] = (parent == EMPTY) ? EMPTY : col2sup [parent] ;
    }
    return aTree;
}

bool applyOrdering(int n, int *post, int *invPost, int *colCount, int *ETree){
    int k =0;
    int *Wi = new int[n]();
    int newchild, oldchild, newparent, oldparent ;

    for (k = 0 ; k < n ; k++){
        Wi [k] = colCount [post [k]] ;
    }
    for (k = 0 ; k < n ; k++){
        colCount [k] = Wi [k] ;
    }
    for (k = 0 ; k < n ; k++){
        invPost [post [k]] = k ;
    }
    /* updated ETree needed only for supernodal case */
    for (newchild = 0 ; newchild < n ; newchild++){
        oldchild = post [newchild] ;
        oldparent = ETree [oldchild] ;
        newparent = (oldparent == EMPTY) ? EMPTY : invPost [oldparent] ;
        Wi [newchild] = newparent ;
    }
    for (k = 0 ; k < n ; k++){
        ETree [k] = Wi [k] ;
    }
}
int *getBlockSet(int n, int* c, int* r, int* lC, int *colCnt,
                 int *AT, int *col2sup, int *postOrdr, int& blockSetSize){
    int *blockSet, *PET=new int[n]();
    int *postOrdrTmp, *invPosrOrdr=new int[n]();
    int *ET = etree(n,c,r,false);
#if DEBUG >0
    for (int i = 0; i < n; ++i) {
        std::cout<<ET[i]<<",";
    }
    std::cout<<"\n";
#endif
    postOrdrTmp = postOrder(ET,n);
    for (int k = 0; k < n; ++k) {
        PET[k] = ET[k];
        postOrdr[k] = postOrdrTmp[k];
    }
    applyOrdering(n,postOrdr,invPosrOrdr,colCnt,PET);
    blockSet = sNodeDetection(PET,invPosrOrdr,colCnt,n,col2sup,blockSetSize);
    int *ATtmp = eTree2aTree(PET,n,blockSet,col2sup,blockSetSize);
    for (int i = 0; i < blockSetSize; ++i) AT[i]=ATtmp[i];
    delete []ATtmp;
#if DEBUG >0
    printf("set\n");
    for (int j = 0; j < n; ++j) {
        printf("%d,",blockSet[j]);
    }
    printf("\netree\n");
    for (int j = 0; j < n; ++j) {
        printf("%d,",ET[j]);
    }
    printf("\npost\n");
    for (int j = 0; j < n; ++j) {
        printf("%d,",PET[j]);
    }
    printf("\nAssembly\n");
    for (int j = 0; j < blockSetSize; ++j) {
        printf("%d,",AT[j]);
    }
    printf("\n");
#endif
    return blockSet;
}
int *LFactorSparsity(int n, int* c, int* r, int* lC, int* lR,
                     int& factorSize){
    int *ET,*AT, *PET;
    int *postOrdr, *invPosrOrdr=new int[n]();
    int *colCnt;
    int *blockSet;
    int *col2sup = new int[n]();
    PET = new int[n]();
    int blockSetSize=0;
    ET = etree(n,c,r,false);
    for (int k = 0; k < n; ++k)
        PET[k]=ET[k];
    postOrdr = postOrder(ET,n);
    /*for (int l = 0; l < n; ++l) {
        invPosrOrdr[postOrdr[l]]=l;
    }
    for (int k = 0; k < n; ++k) {//Post order ETree
        PET[postOrdr[k]]=postOrdr[ET[k]];
    }*/
    colCnt = counts(n,n,c,r,NULL,ET,postOrdr,0);
/*    applyOrdering(n,postOrdr,invPosrOrdr,colCnt,PET);//creating post ordered ETree
    //blockSet = supernodal(ET,postOrdr,colCnt,n,col2sup,blockSetSize);
    blockSet = sNodeDetection(PET,invPosrOrdr,colCnt,n,col2sup,blockSetSize);
    AT = eTree2aTree(PET,n,blockSet,col2sup,blockSetSize);*/
#if DEBUG >0
    printf("\netree\n");
    for (int j = 0; j < n; ++j) {
        printf("%d,",ET[j]);
    }
    printf("\npost\n");
    for (int j = 0; j < n; ++j) {
        printf("%d,",PET[j]);
    }
    printf("\n");
#endif
    lC[0]=0;
    for (int i = 1; i < n+1; ++i) {
        lC[i] = lC[i-1]+colCnt[i-1];
        factorSize+=colCnt[i-1];
    }

    return(colCnt);
}


#endif //CHOLOPENMP_INSPECTION_BLOCK_H
