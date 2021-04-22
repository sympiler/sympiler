//
// Created by Kazem on 10/14/19.
//

#include <cstdint>
#include <cassert>
#include <climits>
#include <algorithm>
#include <omp.h>
#include <iostream>
#include "includes/lbc_utils.h"

namespace sym_lib{

 int make_unique(int *node2Par,
                 std::vector<int> &list,
                 int n, bool *ws){
  int min=INT32_MAX;
  for (int ii = 0; ii < list.size(); ) {
   int tmp = node2Par[list[ii]];
   if(!ws[tmp]){//if first time
    ws[tmp]= true;
    min=min>tmp?tmp:min;
    ii++;
   }else{//otherwise remove it
    list.erase(list.begin()+ii);
   }
  }
//Reset it for future use
  for (int ii = 0; ii < list.size(); ++ii) {
   int tmp = node2Par[list[ii]];
   ws[tmp] = false;
  }
  return min;//returns cluster with min number.
 }


 int dfs_CSC_CC (size_t  n, int j, int *Gp, int *Gi, int *marked, int top,
                 int *xi, int *pstack, std::vector<int> &clashedNodes,
                 const int *pinv) {
  int i, p, p2, done, jnew, head = 0 ;
  xi [0] = j ;                /* initialize the recursion stack */

  while (head >= 0) {
   j = xi [head] ;         /* get j from the top of the recursion stack */
   jnew = pinv ? (pinv [j]) : j ;
   if (!marked[jnew])  {
    marked[jnew] = 1;       /* mark node j as visited */
    pstack [head] =  Gp [jnew] ;
   }
   if(marked[jnew]==-1){//visited in previous CCs
    clashedNodes.push_back(jnew);
   }
   done = 1 ;                  /* node j done if no unvisited neighbors */
   p2 = Gp [jnew+1] ;
   for (p = pstack [head] ; p < p2 ; p++)  /* examine all neighbors of j */
   {
    i = Gi [p] ;            /* consider neighbor node i */
    if (marked[i]==-1)
     clashedNodes.push_back(i);//Another node from prev CCs.
    if (marked[i]) continue ;   /* skip visited node i */
    pstack [head] = p ;     /* pause depth-first search of node j */
    xi [++head] = i ;       /* start dfs at node i */
    done = 0 ;              /* node j is not done */
    break ;                 /* break, to start dfs (i) */
   }
   if (done)               /* depth-first search at node j is done */
   {
    head-- ;            /* remove j from the recursion stack */
    xi [--top] = j ;    /* and place in the output stack */
   }
  }
  return (top) ;
 }


 int dfs_BCSC_CC (size_t  n, int j, size_t *Gp, size_t *Gi_ptr, int *Gi,
                  const int *blk2Col,const int *col2Blk, int *marked,
                  int top, int *xi, int *pstack,
                  std::vector<int> &clashedNodes, const int *pinv){
  int i, p, done, jnew, head = 0 ;
  //bool *marked = new bool[n]();
  xi [0] = blk2Col[j] ;                /* initialize the recursion stack */
  while (head >= 0){
   j = xi [head];         /* get j from the top of the recursion stack */
   jnew=col2Blk[j];
   if (!marked[jnew]){//not visited before at all
    marked[jnew] = 1;       /* mark node j as visited */
    pstack [head] = Gi_ptr[j] ;
   }
   if(marked[jnew]==-1){//visited in previous CCs
    clashedNodes.push_back(jnew);
   }
   done = 1 ;                  /* node j done if no unvisited neighbors */

   int supWdth = blk2Col[jnew+1] - blk2Col[jnew];
   int nxtCol = blk2Col[jnew+1];
   int st=pstack [head];
   int p2 = Gi_ptr[nxtCol] ;
   for (int p1 = st ; p1 < p2 ; p1++)  /* examine all neighbors of j */
   {
    i = Gi [p1];            /* consider neighbor node i */
    int tmp = col2Blk[i];
    supWdth = blk2Col[tmp+1] - blk2Col[tmp];
    if (marked[tmp]==-1)
     clashedNodes.push_back(tmp);//Another node from prev CCs.
    if (marked[tmp]) continue ;   /* skip visited node i */
    pstack [head] = p1;     /* pause depth-first search of node j */
    xi [++head] = i ;       /* start dfs at node i */
    done = 0 ;              /* node j is not done */
    break ;                 /* break, to start dfs (i) */
   }
   if (done)               /* depth-first search at node j is done */
   {
    head-- ;            /* remove j from the recursion stack */
    xi [--top] = col2Blk[j] ;    /* and place in the output stack */
   }
  }
  return (top) ;
 }


 int merge_inner_part(std::vector <std::vector<int>> newLeveledParList,
                      int *inCost,
                      std::vector <std::vector<int>> &mergedLeveledParList,
                      int *outCost, int costThreshold) {
  int lClusterCnt = 0;
  int partNo = newLeveledParList.size();
  for (int i = 0; i < partNo;) {
   int curCost = 0;
   while (curCost < (1 * costThreshold) && i < partNo) {
    curCost += inCost[i];
    mergedLeveledParList[lClusterCnt].insert(
      mergedLeveledParList[lClusterCnt].end(),
      newLeveledParList[i].begin(), newLeveledParList[i].end());
    i++;
   }
   lClusterCnt++;
   outCost[lClusterCnt] = curCost;
  }
  return lClusterCnt;
 }


 bool cmp_cost(subTree a, subTree b) {
  return a.cost > b.cost;
 }


 int find_min(double *cost, int size) {
  double min = INT_MAX;
  int minBin = 0;
  for (int i = 0; i < size; ++i) {
   if (cost[i] < min) {
    min = cost[i];
    minBin = i;
   }
  }
  return minBin;
 }


 int worst_fit_bin_pack(const std::vector <std::vector<int>> &newLeveledParList,
                        double *inCost,
                        std::vector <std::vector<int>> &mergedLeveledParList,
                        double *outCost, int costThreshold, int numOfBins) {

  int lClusterCnt = 0;
  int partNo = newLeveledParList.size();
  //Sorting the subtree list
  std::vector<subTree> partList(partNo);
  for (int i = 0; i < partNo; ++i) {
   partList[i].cost = inCost[i];
   partList[i].nodeList.insert(partList[i].nodeList.begin(),
                               newLeveledParList[i].begin(), newLeveledParList[i].end());
#if 0
   for (int j = 0; j < newLeveledParList[i].size(); ++j) {
    assert(newLeveledParList[i][j] < 19602 );
   }
#endif
  }
  std::sort(partList.begin(), partList.end(), cmp_cost);
  int minBin = 0;
  for (int i = 0; i < partNo; i++) {
   minBin = find_min(outCost, numOfBins);
#if 0
   for (int j = 0; j < partList[i].nodeList.size(); ++j) {
    assert(partList[i].nodeList[j] < 19602 );
   }
#endif
   outCost[minBin] += partList[i].cost;
   mergedLeveledParList[minBin].insert(
     mergedLeveledParList[minBin].end(),
     partList[i].nodeList.begin(), partList[i].nodeList.end());
  }
  return numOfBins;
 }


 int height_partitioning_DAG_trng(int levelNo, int *levelPtr, int *node2Level,
                                  int originalHeight, int innerParts,
                                  int minLevelDist, int divRate,
                                  std::vector<int> &innerPartsSize,
                                  std::vector <std::vector<int>> &slackGroups,
                                  double *subTreeCost, int *partition2Level,
                                  bool sw) {
  int last_level = 0;
  for(int i = levelNo-1; i >= 0; i--) {
   int size = levelPtr[i+1] - levelPtr[i];
   if(size == 1)
    continue;
   else {
    last_level = i+1;
    break;
   }
  }
  int lClusterCnt = 0;
  int *accuSlackGroups = new int[levelNo];
  if (levelNo <= minLevelDist) {
   partition2Level[0] = 0;
   partition2Level[1] = levelNo;
   innerPartsSize.push_back(1);
   lClusterCnt = 1;
   delete []accuSlackGroups;
   return lClusterCnt;
  }
  //assign the nodes_ in normal level set
  for (int i = 0; i < levelNo; ++i) {
   accuSlackGroups[i] = levelPtr[i + 1] - levelPtr[i];
   assert(accuSlackGroups[i] >= 0);
  }
  int nthreads = innerParts>1 ? innerParts : omp_get_max_threads();

  partition2Level[0] = 0;

  /*if(minLevelDist<=0)
   minLevelDist=2;//default parameter*/
  int tmp = minLevelDist;
  if (tmp > partition2Level[lClusterCnt] && tmp < levelNo) {
   //Due to tuning parameter we need this check
   partition2Level[++lClusterCnt] = tmp;
   int size;
   if(accuSlackGroups[tmp-1] / 2 >= nthreads)
    size = nthreads;
   else if(accuSlackGroups[tmp-1] / 2 > 1)
    size = accuSlackGroups[tmp-1] / 2;
   else
    size = 1;
   innerPartsSize.push_back(size);
  }
  tmp += divRate;
  while (tmp < last_level) {//originalHeight - 1) {
   //Ensures a certain number of level in each partition
   int size;
   if(accuSlackGroups[tmp-1] >= nthreads)
    size = nthreads;
   else if(accuSlackGroups[tmp-1] > 1)
    size = accuSlackGroups[tmp-1];
   else
    size = 1;

   innerPartsSize.push_back(size);
   partition2Level[++lClusterCnt] = tmp;
   tmp += divRate;
  }
  partition2Level[++lClusterCnt] = originalHeight + 1;
  innerPartsSize.push_back(1);//The last partition has one element TODO is this true??
  delete[]accuSlackGroups;
  return lClusterCnt;
 }


 int height_partitioning(int levelNo, int *levelPtr, int *node2Level,
                         int originalHeight, int innerParts, int minLevelDist,
                         int divRate, std::vector<int> &innerPartsSize,
                         std::vector <std::vector<int>> &slackGroups,
                         double *subTreeCost, int *partition2Level,
                         bool sw) {
  int levelCut = 0, preLevelCut = 0, finalSeqNodes = 3;
  int lClusterCnt = 0, innerPartsTmp = innerParts;
  if(divRate <= 1){
   divRate = 10; // should be min 2.
   std::cout<<"div rate is set to 10, input is wrong\n";
  }
  //auto partition2Level = new int[levelNo+1]();
  int *accuSlackGroups = new int[levelNo];
  if (levelNo <= 2) {
   partition2Level[0] = 0;
   partition2Level[1] = levelNo;
   innerPartsSize.push_back(1);
   lClusterCnt = 1;
   delete []accuSlackGroups;
   return lClusterCnt;
  }
  //assigne the nodes_ in normal level set
  for (int i = 0; i < levelNo; ++i) {
   /*accuSlackGroups[i] = levelPtr[i+1]-levelPtr[i] +
     slackGroups[i].size();*/
   accuSlackGroups[i] = levelPtr[i + 1] - levelPtr[i];
   assert(accuSlackGroups[i] >= 0);
  }
  partition2Level[0] = 0;
  if (sw) {
   for (int i = minLevelDist; i < levelNo; i += minLevelDist) {
    lClusterCnt++;
    partition2Level[lClusterCnt] = i;
    //for leaves
    if (accuSlackGroups[i] >= divRate * innerPartsTmp) {
     innerPartsSize.push_back(innerPartsTmp);
    } else { //otherwise a divisor of divRate
     int tmp = accuSlackGroups[i] / 2;
     if (tmp > 1)
      innerPartsSize.push_back(tmp);
     else
      break;
    }
   }
   partition2Level[++lClusterCnt] = originalHeight + 1;
   innerPartsSize.push_back(1);//The last partition has one element

  } else {
   /*if(minLevelDist<=0)
    minLevelDist=2;//default parameter*/
   while (innerPartsTmp > 1) {
    while (accuSlackGroups[originalHeight - levelCut - 1] <= innerPartsTmp
           && levelCut < levelNo) {
     levelCut++;
     if(originalHeight - levelCut - 1 < 0 || originalHeight-levelCut-1>= levelNo)
      break;
    }
    //Ensures a certain number of level in each partition
    innerPartsSize.push_back(innerPartsTmp);
    int tmp = originalHeight - levelCut - minLevelDist;
    if (tmp > partition2Level[lClusterCnt] && tmp < levelNo) {
     //Due to tuning parameter we need this check
     partition2Level[++lClusterCnt] = tmp;
    }
    innerPartsTmp /= divRate;
    preLevelCut = levelCut;
    levelCut = 0; // starting from the root
   }
   partition2Level[++lClusterCnt] = originalHeight + 1;
   innerPartsSize.push_back(1);//The last partition has one element
  }
  delete[]accuSlackGroups;
  return lClusterCnt;
 }



}